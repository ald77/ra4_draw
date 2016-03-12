#include "token.hpp"

using namespace std;

Token::Token(const string &function_string, Type type):
  function_(0.),
  string_rep_(function_string),
  type_(type){
  if(type_ == Type::unknown){
    type_ = GetType(function_string);
  }
}

Token::Token(const NamedFunc &function):
  function_(function),
  string_rep_(function.Name()),
  type_(function.IsScalar() ? Type::resolved_scalar : Type::resolved_vector){
}

Token::Type Token::GetType(const string &x){
  //Note that '<', '>', and '!' give "unknown" since they can be part of digraph with '='
  switch(x.size()){
  case 0: return Type::unknown;
  case 1:
    switch(x[0]){
    case '+': return Type::ambiguous_plus;
    case '-': return Type::ambiguous_minus;
    case '*': return Type::multiply;
    case '/': return Type::divide;
    case '%': return Type::modulus;
    case '(': return Type::open_paren;
    case ')': return Type::close_paren;
    case '[': return Type::open_square;
    case ']': return Type::close_square;
    default: return Type::unknown;
    }
  case 2:
    if(x == "=="){
      return Type::equal;
    }else if(x == "!="){
      return Type::not_equal;
    }else if(x == "<="){
      return Type::less_equal;
    }else if(x == ">="){
      return Type::greater_equal;
    }else if(x == "&&"){
      return Type::logical_and;
    }else if(x == "||"){
      return Type::logical_or;
    }else{
      return Type::unknown;
    }
  default:
    return Type::unknown;
  }
  return Type::unknown;
}

ostream & operator << (ostream &stream, const Token &token){
  stream << "Token{" << token.string_rep_ << "}{" << static_cast<unsigned>(token.type_) << '}';
  return stream;
}
