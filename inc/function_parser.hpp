#include <string>
#include <vector>
#include <ostream>

#include "token.hpp"

class NamedFunc;

class FunctionParser{
public:
  FunctionParser() = default;
  FunctionParser(const std::string &function_string);
  FunctionParser(const FunctionParser &) = default;
  FunctionParser & operator=(const FunctionParser &) = default;
  FunctionParser(FunctionParser &&) = default;
  FunctionParser & operator=(FunctionParser &&) = default;
  ~FunctionParser() = default;

  const std::string & FunctionString() const;
  FunctionParser & FunctionString(const std::string &function_string);

  const std::vector<Token> & Tokens() const;
  
  Token ResolveAsToken() const;
  NamedFunc ResolveAsNamedFunc() const;

private:
  std::string input_string_;
  mutable std::vector<Token> tokens_;
  mutable bool tokenized_, solved_;

  FunctionParser(const std::vector<Token> &tokens);

  void Tokenize() const;
  void CheckForUnknowns() const;
  void ResolveVariables() const;
  void EvaluateGroupings() const;
  void MergeParentheses() const;
  void ApplySubscripts() const;
  void DisambiguatePlusMinus() const;
  void ApplyUnary() const;
  void MultiplyAndDivide() const;
  void AddAndSubtract() const;
  void LessGreater() const;
  void EqualOrNot() const;
  void And() const;
  void Or() const;
  void CheckSolved() const;
  void CleanupName() const;
  
  void Solve() const;
  
  std::size_t FindClose(std::size_t i_open_token) const;
  std::string ConcatenateTokenStrings(std::size_t i_start, std::size_t i_end) const;
  void CondenseTokens(std::size_t i_start, std::size_t i_end, const Token &replacement) const;
};

std::ostream & operator << (std::ostream &stream, const FunctionParser &fp);
