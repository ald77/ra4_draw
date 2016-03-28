#ifndef H_NAMED_FUNC
#define H_NAMED_FUNC

#include <string>
#include <functional>
#include <ostream>
#include <vector>

#include "baby.hpp"

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
  NamedFunc(ScalarType x);
  NamedFunc(const NamedFunc &) = default;
  NamedFunc & operator=(const NamedFunc &) = default;
  NamedFunc(NamedFunc &&) = default;
  NamedFunc & operator=(NamedFunc &&) = default;
  ~NamedFunc() = default;  

  const std::string & Name() const;
  NamedFunc & Name(const std::string &name);
  std::string PlainName() const;
  
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

private:
  NamedFunc() = delete;
  std::string name_;
  std::function<ScalarFunc> scalar_func_;
  std::function<VectorFunc> vector_func_;
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

#endif
