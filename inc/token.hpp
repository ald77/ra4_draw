#include <string>
#include <ostream>

#include "named_func.hpp"

struct Token{
  enum class Type{resolved_scalar, resolved_vector, //0-1
      number, variable_name, //2-3
      binary_plus, unary_plus, ambiguous_plus, //4-6
      binary_minus, unary_minus, ambiguous_minus, //7-9
      multiply, divide, modulus, //10-12
      equal, not_equal, greater, less, greater_equal, less_equal, //13-18
      logical_and, logical_or, logical_not, //19-21
      open_paren, close_paren, //22-23
      open_square, close_square, //24-25
      unknown};//26

  Token(const std::string &function_string="", Type type = Type::unknown);
  Token(const NamedFunc &function);
  Token(const Token &) = default;
  Token & operator=(const Token &) = default;
  Token(Token &&) = default;
  Token & operator=(Token &&) = default;
  ~Token() = default;

  static Type GetType(char x);
  static Type GetType(const std::string &x);

  NamedFunc function_;
  std::string string_rep_;
  Type type_;
};

std::ostream & operator << (std::ostream &stream, const Token &token);
