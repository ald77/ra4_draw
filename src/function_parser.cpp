/*! \class FunctionParser

  \brief Converts a string into a NamedFunc

  A FunctionParser is initialized with a string representing a number, variable,
  function, cut, etc. and converts it to a NamedFunc. It first decomposes the
  string into components representing single numbers, variables, operators,
  etc. The components are stored as \link Token Tokens\endlink. The Tokens
  representing constants and variables in Baby are processed to obtain valid
  \link NamedFunc NamedFuncs\endlink. Operators are successively applied
  following standard order of operations to merge the \link Token Tokens\endlink
  into a single Token instance containing a NamedFunc which can return the value
  represented by the initial string.

  Parentheses and brackets are parsed recursively and can be arbitrarily nested.

  Currently has support for the basic arithmetic, logical, and comparison
  operators. Future versions may support ROOT's function syntax,
  e.g. Sum\$(jets_pt).
*/
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

/*!\brief Standard constructor from string representing a function

  \param[in] function_string String representing a number, variable, function,
  cut, etc.
*/
FunctionParser::FunctionParser(const string &function_string):
  input_string_(function_string),
  tokens_(),
  tokenized_(false),
  solved_(false){
  ReplaceAll(input_string_, " ", "");
}

/*!\brief Get string being parsed

  \return String being parsed
*/
const string & FunctionParser::FunctionString() const{
  return input_string_;
}

/*!\brief Set string to be parsed

  \param[in] function_string String to be parsed \return Reference to *this
*/
FunctionParser & FunctionParser::FunctionString(const string &function_string){
  tokenized_ = false;
  solved_ = false;
  input_string_ = function_string;
  ReplaceAll(input_string_, " ", "");
  return *this;
}

/*!\brief Get list of tokens as is at current stage of parsing
 */
const vector<Token> & FunctionParser::Tokens() const{
  return tokens_;
}

/*!\brief Parses provided string into a single Token
 */
Token FunctionParser::ResolveAsToken() const{
  Solve();
  return tokens_.size() ? tokens_.at(0) : Token();
}

/*!\brief Parses provided string into a single NamedFunc
 */
NamedFunc FunctionParser::ResolveAsNamedFunc() const{
  Solve();
  return tokens_.size() ? tokens_.at(0).function_
    : NamedFunc(input_string_,
                [](const Baby &){
                  return 0.;
                });
}

/*!\brief Constructs FunctionParser from list of \link Token Tokens\endlink

  Used by FunctionParser to recursively process lists of \link Token
  Tokens\endlink reprenting contents of parentheses or brackets

  \param[in] tokens List of \link Token Tokens\endlink. Should represent a
  single valid expression
*/
FunctionParser::FunctionParser(const vector<Token> &tokens):
  input_string_(""),
  tokens_(tokens),
  tokenized_(true),
  solved_(false){
  input_string_ = ConcatenateTokenStrings(0, tokens_.size());
}

/*!\brief Parses the string into a list of \link Token Tokens\endlink
 */
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

/*!\brief Checks that all \link Token Tokens\endlink were understood
 */
void FunctionParser::CheckForUnknowns() const{
  for(const auto &token: tokens_){
    if(token.type_ == Token::Type::unknown){
      ERROR("Function string \""+input_string_+"\" contains unknown token \""+token.string_rep_+"\".");
    }
  }
}

/*!\brief Generates NamedFunc for each Token representing a constant or Baby
  variable
*/
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

/*!\brief Recursively evaluates contents of parenthesis and brackets
 */
void FunctionParser::EvaluateGroupings() const{
  for(size_t i_open = 0; i_open < tokens_.size(); ++i_open){
    size_t i_close = FindClose(i_open);
    if(i_close <= i_open || i_close >= tokens_.size()) continue;

    FunctionParser fp(vector<Token>(tokens_.cbegin()+i_open+1, tokens_.cbegin()+i_close));
    Token merged = fp.ResolveAsToken();

    CondenseTokens(i_open+1, i_close, merged);
  }
}

/*!\brief Merges parenthesis \link Token Tokens\endlink with the contents

  Searches for patten {open paren}{value}{close paren} and replaces with single
  Token
*/
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

/*!\brief Merges subscript (bracket) \link Token Tokens\endlink with the
  contents

  Searches for patten {open bracket}{value}{close bracket} and replaces with
  single Token
*/
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

/*!\brief Determines whether "+" and "-" \link Token Tokens\endlink correspond
  to unary or binary operators
*/
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

/*!\brief Merges unary operator \link Token Tokens\endlink with their operands

  Searches for patten {unary operator}{value} and replaces with single Token
*/
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

/*!\brief Merges "*" and "/" \link Token Tokens\endlink with their operands

  Searches for patten {value}{"*" or "/"}{value} and replaces with single Token
*/
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

/*!\brief Merges binary "+" and "-" \link Token Tokens\endlink with their
  operands

  Searches for patten {value}{"+" or "-"}{value} and replaces with single Token
*/
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

/*!\brief Merges "<", "<=", ">", and ">=" \link Token Tokens\endlink with their
  operands

  Searches for patten {value}{comparison}{value} and replaces with single Token
*/
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

/*!\brief Merges "==" and "!=" \link Token Tokens\endlink with their operands

  Searches for patten {value}{comparison}{value} and replaces with single Token
*/
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

/*!\brief Merges "&&" \link Token Tokens\endlink with their operands

  Searches for patten {value}{&&}{value} and replaces with single Token
*/
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

/*!\brief Merges "||" \link Token Tokens\endlink with their operands

  Searches for patten {value}{||}{value} and replaces with single Token
*/
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

/*!\brief Check that we have a single Token with a valid NamedFun
 */
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

/*!\brief Restores string representation to initially provided string

  Token merging process often generates unnecessary parentheses in string, so
  just use the original string
*/
void FunctionParser::CleanupName() const{
  if(tokens_.size() == 1){
    tokens_.at(0).string_rep_ = input_string_;
    tokens_.at(0).function_.Name(input_string_);
  }
}

/*!\brief Runs full parse from start to finish, caching result
 */
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

/*!\brief Find position of closing parenthesis/bracket corresponding to given opening partner

  \param[in] i_open_token Position in current list of \link Token Tokens\endlink
  of opening parenthesis/bracket

  \return Position in current list of \link Token Tokens\endlink of closing
  parenthesis/bracket
*/
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

/*!\brief Concatenates string representation of \link Token Tokens\endlink in
  range [i_start, i_end)

  \param[in] i_start Position of starting Token (inclusive)

  \param[in] i_end Position of ending Token (exclusive)

  \return Concatentation of string representations of \link Token Tokens\endlink
  in range [i_start, i_end)
*/
string FunctionParser::ConcatenateTokenStrings(size_t i_start, size_t i_end) const{
  string result = "";
  for(;i_start < i_end; ++i_start){
    result += tokens_.at(i_start).string_rep_;
  }
  return result;
}

/*!\brief Replace \link Token Tokens\endlink in range [i_start, i_end) with
  replacement

  \param[in] i_start Position of starting Token (inclusive)

  \param[in] i_end Position of ending Token (exclusive)

  \param[in] replacement Token which replaces the range
*/
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

/*!\brief Print FunctionParser to output stream

  \param[in,out] stream Output stream to print to

  \param[in] fp FunctionParser to print

  \return Reference to stream
*/
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
