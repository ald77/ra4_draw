#ifndef H_NAMED_FUNC
#define H_NAMED_FUNC

#include <string>
#include <functional>
#include <ostream>
#include <vector>

#include "TString.h"

#include "core/baby.hpp"

class NamedFunc{
public:
  using ScalarType = double;
  using VectorType = std::vector<ScalarType>;
  using ScalarFunc = ScalarType(const Baby &);
  using VectorFunc = VectorType(const Baby &);

  NamedFunc(const std::string &name,
            const std::function<ScalarFunc> &function);
  NamedFunc(const std::string &name,
            const std::function<VectorFunc> &function);
  NamedFunc(const std::string &function);
  NamedFunc(const char *function);
  NamedFunc(const TString &function);
  NamedFunc(ScalarType x);
  NamedFunc(const NamedFunc &) = default;
  NamedFunc & operator=(const NamedFunc &) = default;
  NamedFunc(NamedFunc &&) = default;
  NamedFunc & operator=(NamedFunc &&) = default;
  ~NamedFunc() = default;

  const std::string & Name() const;
  NamedFunc & Name(const std::string &name);
  std::string PlainName() const;
  std::string PrettyName() const;

  NamedFunc & Function(const std::function<ScalarFunc> &function);
  NamedFunc & Function(const std::function<VectorFunc> &function);
  const std::function<ScalarFunc> & ScalarFunction() const;
  const std::function<VectorFunc> & VectorFunction() const;

  bool IsScalar() const;
  bool IsVector() const;

  ScalarType GetScalar(const Baby &b) const;
  VectorType GetVector(const Baby &b) const;

  NamedFunc & operator += (const NamedFunc &func);
  NamedFunc & operator -= (const NamedFunc &func);
  NamedFunc & operator *= (const NamedFunc &func);
  NamedFunc & operator /= (const NamedFunc &func);
  NamedFunc & operator %= (const NamedFunc &func);

  NamedFunc operator [] (const NamedFunc &func) const;

private:
  NamedFunc() = delete;
  std::string name_;//!<String representation of the function
  std::function<ScalarFunc> scalar_func_;//<!Scalar function. Cannot be valid at same time as NamedFunc::vector_func_.
  std::function<VectorFunc> vector_func_;//<!Vector function. Cannot be valid at same time as NamedFunc::scalar_func_.

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

std::ostream & operator<<(std::ostream &stream, const NamedFunc &function);

bool HavePass(const NamedFunc::VectorType &v);
bool HavePass(const std::vector<NamedFunc::VectorType> &vv);

#endif
