#ifndef H_NAMED_FUNC
#define H_NAMED_FUNC

#include <string>
#include <functional>
#include <ostream>

#include "baby.hpp"

class NamedFunc{
public:
  using FuncType = double(const Baby &);

  NamedFunc(const std::string &name, const std::function<FuncType> &function);
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
  NamedFunc & Function(const std::function<FuncType> & function);

  double operator()(const Baby &b);

  NamedFunc & operator += (const NamedFunc &func);
  NamedFunc & operator -= (const NamedFunc &func);
  NamedFunc & operator *= (const NamedFunc &func);
  NamedFunc & operator /= (const NamedFunc &func);

  NamedFunc operator + () const;
  NamedFunc operator - () const;
  NamedFunc operator ! () const;

private:
  NamedFunc() = delete;
  std::string name_;
  std::function<FuncType> function_;

  void CleanName();
};

NamedFunc operator+(NamedFunc a, NamedFunc b);
NamedFunc operator-(NamedFunc a, NamedFunc b);
NamedFunc operator*(NamedFunc a, NamedFunc b);
NamedFunc operator/(NamedFunc a, NamedFunc b);

NamedFunc operator&&(NamedFunc a, NamedFunc b);
NamedFunc operator||(NamedFunc a, NamedFunc b);

std::ostream & operator<<(std::ostream &stream, const NamedFunc &function);

#define FUNC(x) NamedFunc(#x, [](const Baby &b){return x;})

#endif
