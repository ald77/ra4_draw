#include "function_parser.hpp"

#include <cstdlib>
#include <cctype>

#include "utilities.hpp"
#include "named_func.hpp"

using namespace std;

using ScalarType = NamedFunc::ScalarType;
using VectorType = NamedFunc::VectorType;
using ScalarFunc = NamedFunc::ScalarFunc;
using VectorFunc = NamedFunc::VectorFunc;

FunctionParser::FunctionParser(const string &function_string):
  input_string_(function_string),
  tokens_(),
  tokenized_(false),
  solved_(false){
  ReplaceAll(input_string_, " ", "");
}

const string & FunctionParser::FunctionString() const{
  return input_string_;
}

FunctionParser & FunctionParser::FunctionString(const string &function_string){
  tokenized_ = false;
  solved_ = false;
  input_string_ = function_string;
  ReplaceAll(input_string_, " ", "");
  return *this;
}

const vector<Token> & FunctionParser::Tokens() const{
  return tokens_;
}

Token FunctionParser::ResolveAsToken() const{
  Solve();
  return tokens_.size() ? tokens_.at(0) : Token();
}

NamedFunc FunctionParser::ResolveAsNamedFunc() const{
  Solve();
  return tokens_.size() ? tokens_.at(0).function_
    : NamedFunc(input_string_,
                [](const Baby &){
                  return 0.;
                });
}

FunctionParser::FunctionParser(const vector<Token> &tokens):
  input_string_(""),
  tokens_(tokens),
  tokenized_(true),
  solved_(false){
  input_string_ = ConcatenateTokenStrings(0, tokens_.size());
}

void FunctionParser::Tokenize() const{
  if(tokenized_) return;
  size_t start = 0;
  while(start < input_string_.size()){
    char start_char = input_string_[start];
    if(Token::GetType(input_string_.substr(start, 1)) != Token::Type::unknown){
      tokens_.push_back(Token(input_string_.substr(start, 1)));
      ++start;
    }else if(Token::GetType(input_string_.substr(start, 2)) != Token::Type::unknown){
      //This branch will catch digraphs "<=", ">=", and "!="
      tokens_.push_back(Token(input_string_.substr(start, 2)));
      start+=2;
    }else if(start_char=='<'){
      //"<=" digraph already caught
      tokens_.push_back(Token(input_string_.substr(start, 1), Token::Type::less));
      ++start;
    }else if(start_char=='>'){
      //">=" digraph already caught
      tokens_.push_back(Token(input_string_.substr(start, 1), Token::Type::greater));
      ++start;
    }else if(start_char=='!'){
      //"!=" digraph already caught
      tokens_.push_back(Token(input_string_.substr(start, 1), Token::Type::logical_not));
      ++start;
    }else if(isalpha(start_char) || start_char == '_'){
      size_t count = 1;
      while(start+count < input_string_.size()
            && (isalnum(input_string_[start+count]) || input_string_[start+count] == '_')){
        ++count;
      }
      tokens_.push_back(Token(input_string_.substr(start, count), Token::Type::variable_name));
      start+=count;
    }else if(isdigit(start_char) || start_char == '.'){
      string from_start = input_string_.substr(start);
      char *cp = nullptr;
      double val = strtod(&from_start[0], &cp);
      string remaining = cp;
      if(val != 0. || remaining != from_start){
        size_t length = from_start.size() - remaining.size();
        if(length < 1) length = 1;
        tokens_.push_back(Token(from_start.substr(0, length), Token::Type::number));
        start+=length;
      }else{
        tokens_.push_back(Token(input_string_.substr(start, 1), Token::Type::unknown));
        ++start;
      }
    }else{
      tokens_.push_back(Token(input_string_.substr(start, 1), Token::Type::unknown));
      ++start;
    }
  }
  tokenized_ = true;
}

void FunctionParser::CheckForUnknowns() const{
  for(const auto &token: tokens_){
    if(token.type_ == Token::Type::unknown){
      ERROR("Function string \""+input_string_+"\" contains unknown token \""+token.string_rep_+"\".");
    }
  }
}

void FunctionParser::ResolveVariables() const{
  for(auto &token: tokens_){
    if(token.type_ == Token::Type::variable_name){
      token.function_ = Baby::GetFunction(token.string_rep_);
      token.type_ = token.function_.IsScalar() ? Token::Type::resolved_scalar : Token::Type::resolved_vector;
    }else if(token.type_ == Token::Type::number){
      char *cp = nullptr;
      NamedFunc::ScalarType val = strtod(&token.string_rep_[0], &cp);
      token.function_ = NamedFunc(token.string_rep_,
                                  [val](const Baby &){
                                    return val;
                                  });
      token.type_ = Token::Type::resolved_scalar;
    }
  }
}

void FunctionParser::EvaluateGroupings() const{
  for(size_t i_open = 0; i_open < tokens_.size(); ++i_open){
    size_t i_close = FindClose(i_open);
    if(i_close <= i_open || i_close >= tokens_.size()) continue;

    FunctionParser fp(vector<Token>(tokens_.cbegin()+i_open+1, tokens_.cbegin()+i_close));
    Token merged = fp.ResolveAsToken();

    CondenseTokens(i_open+1, i_close, merged);
  }
}

void FunctionParser::MergeParentheses() const{
  for(size_t i = 0; i+2 < tokens_.size(); ++i){
    Token &open = tokens_.at(i+0);
    Token &inner = tokens_.at(i+1);
    Token &close = tokens_.at(i+2);

    if(open.type_ != Token::Type::open_paren
       || (inner.type_ != Token::Type::resolved_scalar
           && inner.type_ != Token::Type::resolved_vector)
       || close.type_ != Token::Type::close_paren){
      continue;
    }

    string name = ConcatenateTokenStrings(i, i+3);
    NamedFunc merged_func = inner.function_;
    merged_func.Name(name);
    Token merged(merged_func);

    CondenseTokens(i, i+3, merged);
  }
}

void FunctionParser::ApplySubscripts() const{
  for(size_t i = 0; i+3 < tokens_.size(); ++i){
    Token &vec = tokens_.at(i);
    Token &open = tokens_.at(i+1);
    Token &sub = tokens_.at(i+2);
    Token &close = tokens_.at(i+3);

    if(vec.type_ != Token::Type::resolved_vector
       || open.type_ != Token::Type::open_square
       || sub.type_ != Token::Type::resolved_scalar
       || close.type_ != Token::Type::close_square){
      continue;
    }

    function<VectorFunc> vec_func = vec.function_.VectorFunction();
    function<ScalarFunc> sub_func = sub.function_.ScalarFunction();
    function<ScalarFunc> function = [vec_func,sub_func](const Baby &b){
      return vec_func(b).at(sub_func(b));
    };
    string name = ConcatenateTokenStrings(i, i+4);
    Token merged(NamedFunc(name, function));

    CondenseTokens(i, i+4, merged);
  }
}

void FunctionParser::DisambiguatePlusMinus() const{
  for(size_t i = 0; i < tokens_.size(); ++i){
    Token prev;
    if(i > 0){
      prev = tokens_.at(i-1);
    }else{
      prev.type_ = Token::Type::open_paren;
    }
    Token &cur = tokens_.at(i);
    if(cur.type_ != Token::Type::ambiguous_plus && cur.type_ != Token::Type::ambiguous_minus) continue;

    Token::Type binary_type = Token::Type::binary_plus;
    Token::Type unary_type = Token::Type::unary_plus;
    Token::Type ambiguous_type = Token::Type::ambiguous_plus;
    if(cur.type_ == Token::Type::ambiguous_minus){
      binary_type = Token::Type::binary_minus;
      unary_type = Token::Type::unary_minus;
      ambiguous_type = Token::Type::ambiguous_minus;
    }

    switch(prev.type_){
    case Token::Type::resolved_scalar:
    case Token::Type::resolved_vector:
    case Token::Type::number:
    case Token::Type::variable_name:
    case Token::Type::close_paren:
    case Token::Type::close_square:
      cur.type_ = binary_type;
      break;
    case Token::Type::binary_plus:
    case Token::Type::unary_plus:
    case Token::Type::ambiguous_plus:
    case Token::Type::binary_minus:
    case Token::Type::unary_minus:
    case Token::Type::ambiguous_minus:
    case Token::Type::multiply:
    case Token::Type::divide:
    case Token::Type::modulus:
    case Token::Type::equal:
    case Token::Type::not_equal:
    case Token::Type::greater:
    case Token::Type::less:
    case Token::Type::greater_equal:
    case Token::Type::less_equal:
    case Token::Type::logical_and:
    case Token::Type::logical_or:
    case Token::Type::logical_not:
    case Token::Type::open_paren:
    case Token::Type::open_square:
      cur.type_ = unary_type;
      break;
    case Token::Type::unknown:
    default:
      cur.type_ = ambiguous_type;
      break;
    }
  }
}

void FunctionParser::ApplyUnary() const{
  for(size_t i = 0; i+1 < tokens_.size(); ++i){
    Token &op = tokens_.at(i);
    Token &x = tokens_.at(i+1);

    if(x.type_ != Token::Type::resolved_scalar && x.type_ != Token::Type::resolved_vector) continue;

    Token merged;
    if(op.type_ == Token::Type::unary_plus){
      merged = Token(NamedFunc(+x.function_));
    }else if(op.type_ == Token::Type::unary_minus){
      merged = Token(NamedFunc(-x.function_));
    }else if(op.type_ == Token::Type::logical_not){
      merged = Token(NamedFunc(!x.function_));
    }else{
      continue;
    }

    CondenseTokens(i, i+2, merged);
    if(i!=0) i-=2;//Need to check previous token in case of multiple unary operators
  }
}

void FunctionParser::MultiplyAndDivide() const{
  for(size_t i = 0; i+2 < tokens_.size(); ++i){
    Token &a = tokens_.at(i);
    Token &op = tokens_.at(i+1);
    Token &b = tokens_.at(i+2);

    if(a.type_ != Token::Type::resolved_scalar && a.type_ != Token::Type::resolved_vector) continue;
    if(b.type_ != Token::Type::resolved_scalar && b.type_ != Token::Type::resolved_vector) continue;

    Token merged;
    if(op.type_ == Token::Type::multiply){
      merged = Token(NamedFunc(a.function_ * b.function_));
    }else if(op.type_ == Token::Type::divide){
      merged = Token(NamedFunc(a.function_ / b.function_));
    }else if(op.type_ == Token::Type::modulus){
      merged = Token(NamedFunc(a.function_ % b.function_));
    }else{
      continue;
    }

    CondenseTokens(i, i+3, merged);
    --i;//Need to recheck token in case of successive multiplications
  }
}

void FunctionParser::AddAndSubtract() const{
  for(size_t i = 0; i+2 < tokens_.size(); ++i){
    Token &a = tokens_.at(i);
    Token &op = tokens_.at(i+1);
    Token &b = tokens_.at(i+2);

    if(a.type_ != Token::Type::resolved_scalar && a.type_ != Token::Type::resolved_vector) continue;
    if(b.type_ != Token::Type::resolved_scalar && b.type_ != Token::Type::resolved_vector) continue;

    Token merged;
    if(op.type_ == Token::Type::binary_plus){
      merged = Token(NamedFunc(a.function_ + b.function_));
    }else if(op.type_ == Token::Type::binary_minus){
      merged = Token(NamedFunc(a.function_ - b.function_));
    }else{
      continue;
    }

    CondenseTokens(i, i+3, merged);
    --i;//Need to recheck token in case of successive additions
  }
}

void FunctionParser::LessGreater() const{
  for(size_t i = 0; i+2 < tokens_.size(); ++i){
    Token &a = tokens_.at(i);
    Token &op = tokens_.at(i+1);
    Token &b = tokens_.at(i+2);

    if(a.type_ != Token::Type::resolved_scalar && a.type_ != Token::Type::resolved_vector) continue;
    if(b.type_ != Token::Type::resolved_scalar && b.type_ != Token::Type::resolved_vector) continue;

    Token merged;
    if(op.type_ == Token::Type::greater){
      merged = Token(NamedFunc(a.function_ > b.function_));
    }else if(op.type_ == Token::Type::less){
      merged = Token(NamedFunc(a.function_ < b.function_));
    }else if(op.type_ == Token::Type::greater_equal){
      merged = Token(NamedFunc(a.function_ >= b.function_));
    }else if(op.type_ == Token::Type::less_equal){
      merged = Token(NamedFunc(a.function_ <= b.function_));
    }else{
      continue;
    }

    CondenseTokens(i, i+3, merged);
    --i;//Need to recheck token in case of successive comparisons
  }
}

void FunctionParser::EqualOrNot() const{
  for(size_t i = 0; i+2 < tokens_.size(); ++i){
    Token &a = tokens_.at(i);
    Token &op = tokens_.at(i+1);
    Token &b = tokens_.at(i+2);

    if(a.type_ != Token::Type::resolved_scalar && a.type_ != Token::Type::resolved_vector) continue;
    if(b.type_ != Token::Type::resolved_scalar && b.type_ != Token::Type::resolved_vector) continue;

    Token merged;
    if(op.type_ == Token::Type::equal){
      merged = Token(NamedFunc(a.function_ == b.function_));
    }else if(op.type_ == Token::Type::not_equal){
      merged = Token(NamedFunc(a.function_ != b.function_));
    }else{
      continue;
    }

    CondenseTokens(i, i+3, merged);
    --i;//Need to recheck token in case of successive comparisons
  }
}

void FunctionParser::And() const{
  for(size_t i = 0; i+2 < tokens_.size(); ++i){
    Token &a = tokens_.at(i);
    Token &op = tokens_.at(i+1);
    Token &b = tokens_.at(i+2);

    if(a.type_ != Token::Type::resolved_scalar && a.type_ != Token::Type::resolved_vector) continue;
    if(b.type_ != Token::Type::resolved_scalar && b.type_ != Token::Type::resolved_vector) continue;

    Token merged;
    if(op.type_ == Token::Type::logical_and){
      merged = Token(NamedFunc(a.function_ && b.function_));
    }else{
      continue;
    }

    CondenseTokens(i, i+3, merged);
    --i;//Need to recheck token in case of successive ANDs
  }
}

void FunctionParser::Or() const{
  for(size_t i = 0; i+2 < tokens_.size(); ++i){
    Token &a = tokens_.at(i);
    Token &op = tokens_.at(i+1);
    Token &b = tokens_.at(i+2);

    if(a.type_ != Token::Type::resolved_scalar && a.type_ != Token::Type::resolved_vector) continue;
    if(b.type_ != Token::Type::resolved_scalar && b.type_ != Token::Type::resolved_vector) continue;

    Token merged;
    if(op.type_ == Token::Type::logical_or){
      merged = Token(NamedFunc(a.function_ || b.function_));
    }else{
      continue;
    }

    CondenseTokens(i, i+3, merged);
    --i;//Need to recheck token in case of successive ORs
  }
}

void FunctionParser::CheckSolved() const{
  if(tokens_.size() == 0){
    DBG("No tokens found for " << (*this) << ".");
  }else if(tokens_.size() > 1){
    ostringstream oss;
    oss << "Could not condense " << (*this) << "." << flush;
    ERROR(oss.str());
  }else if(tokens_.at(0).type_ != Token::Type::resolved_vector
           && tokens_.at(0).type_ != Token::Type::resolved_scalar){
    ostringstream oss;
    oss << "Unknown token in " << (*this) << "." << flush;
    ERROR(oss.str());
  }
}

void FunctionParser::CleanupName() const{
  if(tokens_.size() == 1){
    tokens_.at(0).string_rep_ = input_string_;
    tokens_.at(0).function_.Name(input_string_);
  }
}

void FunctionParser::Solve() const{
  if(solved_) return;
  Tokenize();
  CheckForUnknowns();
  ResolveVariables();
  EvaluateGroupings();
  MergeParentheses();
  ApplySubscripts();
  DisambiguatePlusMinus();
  ApplyUnary();
  MultiplyAndDivide();
  AddAndSubtract();
  LessGreater();
  EqualOrNot();
  And();
  Or();
  CheckSolved();
  CleanupName();
  solved_ = true;
}

size_t FunctionParser::FindClose(size_t i_open_token) const{
  if(i_open_token > tokens_.size()) return tokens_.size();

  const auto &open_token = tokens_.at(i_open_token);
  if(open_token.type_ != Token::Type::open_paren
     && open_token.type_ != Token::Type::open_square){
    return i_open_token;
  }

  Token::Type open_type, close_type;
  if(open_token.type_ == Token::Type::open_square){
    open_type = Token::Type::open_square;
    close_type = Token::Type::close_square;
  }else{
    open_type = Token::Type::open_paren;
    close_type = Token::Type::close_paren;
  }

  int open_count = 1;
  size_t i_close_token;
  for(i_close_token = i_open_token+1; i_close_token < tokens_.size() && open_count > 0; ++i_close_token){
    const auto &token = tokens_.at(i_close_token);
    if(token.type_ == open_type){
      ++open_count;
    }else if(token.type_ == close_type){
      --open_count;
    }
  }
  if(open_count <= 0){
    return --i_close_token;
  }else{
    return tokens_.size();
  }
}

string FunctionParser::ConcatenateTokenStrings(size_t i_start, size_t i_end) const{
  string result = "";
  for(;i_start < i_end; ++i_start){
    result += tokens_.at(i_start).string_rep_;
  }
  return result;
}

void FunctionParser::CondenseTokens(size_t i_start, size_t i_end, const Token &replacement) const{
  if(i_end < i_start) return;
  size_t num_replaced = i_end-i_start;
  vector<Token> new_tokens(tokens_.size()+1-num_replaced);
  for(size_t i = 0; i < i_start; ++i){
    new_tokens.at(i) = tokens_.at(i);
  }
  new_tokens.at(i_start) = replacement;
  for(size_t i = i_start + 1; i < new_tokens.size(); ++i){
    new_tokens.at(i) = tokens_.at(i+num_replaced-1);
  }
  tokens_ = new_tokens;
}

ostream & operator << (ostream &stream, const FunctionParser &fp){
  stream << "FunctionParser::" << fp.FunctionString() << "::{";
  const auto &tokens = fp.Tokens();
  for(auto token = tokens.cbegin();
      token != tokens.cend();
      ++token){
    if(token != tokens.cbegin()) stream << ", ";
    stream << (*token);
  }
  stream << '}';
  return stream;
}
