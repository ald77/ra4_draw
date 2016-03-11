#ifndef H_NAMED_FUNC
#define H_NAMED_FUNC

#include <string>
#include <functional>
#include <ostream>
#include <vector>

#include "baby.hpp"

class NamedFunc{
public:
  using FuncType = std::vector<double>(const Baby &);

  NamedFunc(const std::string &name, const std::function<FuncType> &function, bool is_vector=false);
  NamedFunc(const std::string &function);
  NamedFunc(const char *function);
  NamedFunc(double x);
  NamedFunc(const NamedFunc &) = default;
  NamedFunc & operator=(const NamedFunc &) = default;
  NamedFunc(NamedFunc &&) = default;
  NamedFunc & operator=(NamedFunc &&) = default;
  ~NamedFunc() = default;  

  const std::string & Name() const;
  NamedFunc & Name(const std::string &name);

  const std::function<FuncType> & Function() const;
  NamedFunc & Function(const std::function<FuncType> & function, bool is_vector=false);
  bool IsVector() const;
  bool IsScalar() const;

  std::function<FuncType>::result_type operator()(const Baby &b);
  std::function<FuncType>::result_type::value_type operator()(const Baby &b, size_t);

  NamedFunc & operator += (const NamedFunc &func);
  NamedFunc & operator -= (const NamedFunc &func);
  NamedFunc & operator *= (const NamedFunc &func);
  NamedFunc & operator /= (const NamedFunc &func);
  NamedFunc & operator %= (const NamedFunc &func);

private:
  NamedFunc() = delete;
  std::string name_;
  std::function<FuncType> function_;
  bool is_vector_;

  void CleanName();
};

NamedFunc operator + (NamedFunc f, NamedFunc g);
NamedFunc operator - (NamedFunc f, NamedFunc g);
NamedFunc operator * (NamedFunc f, NamedFunc g);
NamedFunc operator / (NamedFunc f, NamedFunc g);
NamedFunc operator % (NamedFunc f, NamedFunc g);

NamedFunc operator + (NamedFunc f);
NamedFunc operator - (NamedFunc f);
  
NamedFunc operator == (NamedFunc f, NamedFunc g);
NamedFunc operator != (NamedFunc f, NamedFunc g);
NamedFunc operator > (NamedFunc f, NamedFunc g);
NamedFunc operator < (NamedFunc f, NamedFunc g);
NamedFunc operator >= (NamedFunc f, NamedFunc g);
NamedFunc operator <= (NamedFunc f, NamedFunc g);
  
NamedFunc operator && (NamedFunc f, NamedFunc g);
NamedFunc operator || (NamedFunc f, NamedFunc g);

NamedFunc operator ! (NamedFunc f);

std::ostream & operator<<(std::ostream &stream, NamedFunc function);

#define FUNC(x) NamedFunc(#x, [](const Baby &b){return std::vector<double>(1,x);})

#endif
