#include "named_func.hpp"

#include <iostream>
#include <utility>

#include "utilities.hpp"
#include "function_parser.hpp"

using namespace std;

using ScalarType = NamedFunc::ScalarType;
using VectorType = NamedFunc::VectorType;
using ScalarFunc = NamedFunc::ScalarFunc;
using VectorFunc = NamedFunc::VectorFunc;

namespace{
  ScalarType MyModulus(ScalarType x, ScalarType y){
    return fmod(x,y);
  }
  
  template<typename Operator>
    static function<ScalarFunc> ApplyOp(const function<ScalarFunc> &f,
                                        const Operator &op){
    if(!static_cast<bool>(f)) return f;
    function<ScalarType(ScalarType)> op_c(op);
    return [f,op_c](const Baby &b){
      return op_c(f(b));
    };
  }

  template<typename Operator>
    static function<VectorFunc> ApplyOp(const function<VectorFunc> &f,
                                        const Operator &op){
    if(!static_cast<bool>(f)) return f;
    function<ScalarType(ScalarType)> op_c(op);
    return [f,op_c](const Baby &b){
      VectorType v = f(b);
      for(auto &x: v){
        x = op_c(x);
      }
      return v;
    };
  }

  template<typename Operator>
    static pair<function<ScalarFunc>, function<VectorFunc> > ApplyOp(const function<ScalarFunc> &sfa,
                                                                     const function<VectorFunc> &vfa,
                                                                     const function<ScalarFunc> &sfb,
                                                                     const function<VectorFunc> &vfb,
                                                                     const Operator &op){
    function<ScalarType(ScalarType,ScalarType)> op_c(op);
    function<ScalarFunc> sfo;
    function<VectorFunc> vfo;
    if(static_cast<bool>(sfa) && static_cast<bool>(sfb)){
      sfo = [sfa,sfb,op_c](const Baby &b){
        return op_c(sfa(b), sfb(b));
      };
    }else if(static_cast<bool>(sfa) && static_cast<bool>(vfb)){
      vfo = [sfa,vfb,op_c](const Baby &b){
        ScalarType sa = sfa(b);
        VectorType vb = vfb(b);
        VectorType vo(vb.size());
        for(size_t i = 0; i < vo.size(); ++i){
          vo.at(i) = op_c(sa, vb.at(i));
        }
        return vo;
      };
    }else if(static_cast<bool>(vfa) && static_cast<bool>(sfb)){
      vfo = [vfa,sfb,op_c](const Baby &b){
        VectorType va = vfa(b);
        ScalarType sb = sfb(b);
        VectorType vo(va.size());
        for(size_t i = 0; i < vo.size(); ++i){
          vo.at(i) = op_c(va.at(i), sb);
        }
        return vo;
      };
    }else if(static_cast<bool>(vfa) && static_cast<bool>(vfb)){
      vfo = [vfa,vfb,op_c](const Baby &b){
        VectorType va = vfa(b);
        VectorType vb = vfb(b);
        VectorType vo(va.size() > vb.size() ? vb.size() : va.size());
        for(size_t i = 0; i < vo.size(); ++i){
          vo.at(i) = op_c(va.at(i), vb.at(i));
        }
        return vo;
      };
    }
    return make_pair(sfo, vfo);
  }
}

NamedFunc::NamedFunc(const string &name,
                     const function<ScalarFunc> &function):
  name_(name),
  scalar_func_(function),
  vector_func_(){
  CleanName();
}

NamedFunc::NamedFunc(const string &name,
                     const function<VectorFunc> &function):
  name_(name),
  scalar_func_(),
  vector_func_(function){
  CleanName();
}

NamedFunc::NamedFunc(const string &function):
  NamedFunc(FunctionParser(function).ResolveAsNamedFunc()){
}

NamedFunc::NamedFunc(const char *function):
  NamedFunc(string(function)){
}

NamedFunc::NamedFunc(ScalarType x):
  name_(ToString(x)),
  scalar_func_([x](const Baby&){return x;}),
  vector_func_(){
}

const string & NamedFunc::Name() const{
  return name_;
}

NamedFunc & NamedFunc::Name(const string &name){
  name_ = name;
  CleanName();
  return *this;
}

string NamedFunc::PlainName() const{
  string plain = name_;
  ReplaceAll(plain, "b.", "");
  ReplaceAll(plain, "()", "");
  ReplaceAll(plain, ".", "p");
  ReplaceAll(plain, "(", "OP");
  ReplaceAll(plain, ")", "CP");
  ReplaceAll(plain, "[", "OB");
  ReplaceAll(plain, "]", "CB");
  ReplaceAll(plain, "{", "OC");
  ReplaceAll(plain, "}", "CC");
  ReplaceAll(plain, "+", "PLS");
  ReplaceAll(plain, "-", "MNS");
  ReplaceAll(plain, "*", "TMS");
  ReplaceAll(plain, "/", "DIV");
  ReplaceAll(plain, "%", "MOD");
  ReplaceAll(plain, "!", "NOT");
  ReplaceAll(plain, "&&", "AND");
  ReplaceAll(plain, "||", "OR");
  ReplaceAll(plain, "==", "EQL");
  ReplaceAll(plain, "<=", "GEQ");
  ReplaceAll(plain, ">=", "LEQ");
  ReplaceAll(plain, ">", "GTR");
  ReplaceAll(plain, "<", "LES");
  ReplaceAll(plain, "=", "EQL");
  ReplaceAll(plain, "&", "BITAND");
  ReplaceAll(plain, "|", "BITOR");
  ReplaceAll(plain, "^", "BITXOR");
  ReplaceAll(plain, "~", "BITNOT");
  ReplaceAll(plain, "__", "_");
  for(size_t i = 0; i < plain.size(); ++i){
    if(isalnum(plain.at(i)) || plain.at(i) == '.' || plain.at(i) == '_') continue;
    plain = plain.substr(0,i)+plain.substr(i+1);
  }
  return plain;
}

NamedFunc & NamedFunc::Function(const function<ScalarFunc> &f){
  if(!static_cast<bool>(f)) return *this;
  scalar_func_ = f;
  vector_func_ = function<VectorFunc>();
  return *this;
}

NamedFunc & NamedFunc::Function(const function<VectorFunc> &f){
  if(!static_cast<bool>(f)) return *this;
  scalar_func_ = function<ScalarFunc>();
  vector_func_ = f;
  return *this;
}

const function<ScalarFunc> & NamedFunc::ScalarFunction() const{
  return scalar_func_;
}

const function<VectorFunc> & NamedFunc::VectorFunction() const{
  return vector_func_;
}

bool NamedFunc::IsScalar() const{
  return static_cast<bool>(scalar_func_);
}

bool NamedFunc::IsVector() const{
  return static_cast<bool>(vector_func_);
}

ScalarType NamedFunc::GetScalar(const Baby &b) const{
  return scalar_func_(b);
}

VectorType NamedFunc::GetVector(const Baby &b) const{
  return vector_func_(b);
}

NamedFunc & NamedFunc::operator += (const NamedFunc &func){
  name_ = "("+name_ + ")+(" + func.name_ + ")";
  auto fp = ApplyOp(scalar_func_, vector_func_,
                    func.scalar_func_, func.vector_func_,
                    plus<ScalarType>());
  scalar_func_ = fp.first;
  vector_func_ = fp.second;
  return *this;
}

NamedFunc & NamedFunc::operator -= (const NamedFunc &func){
  name_ = "("+name_ + ")-(" + func.name_ + ")";
  auto fp = ApplyOp(scalar_func_, vector_func_,
                    func.scalar_func_, func.vector_func_,
                    minus<ScalarType>());
  scalar_func_ = fp.first;
  vector_func_ = fp.second;
  return *this;
}

NamedFunc & NamedFunc::operator *= (const NamedFunc &func){
  name_ = "("+name_ + ")*(" + func.name_ + ")";
  auto fp = ApplyOp(scalar_func_, vector_func_,
                    func.scalar_func_, func.vector_func_,
                    multiplies<ScalarType>());
  scalar_func_ = fp.first;
  vector_func_ = fp.second;
  return *this;
}

NamedFunc & NamedFunc::operator /= (const NamedFunc &func){
  name_ = "("+name_ + ")/(" + func.name_ + ")";
  auto fp = ApplyOp(scalar_func_, vector_func_,
                    func.scalar_func_, func.vector_func_,
                    divides<ScalarType>());
  scalar_func_ = fp.first;
  vector_func_ = fp.second;
  return *this;
}

NamedFunc & NamedFunc::operator %= (const NamedFunc &func){
  name_ = "("+name_ + ")%(" + func.name_ + ")";
  auto fp = ApplyOp(scalar_func_, vector_func_,
                    func.scalar_func_, func.vector_func_,
                    MyModulus);
  scalar_func_ = fp.first;
  vector_func_ = fp.second;
  return *this;
}

void NamedFunc::CleanName(){
  ReplaceAll(name_, " ", "");
}

NamedFunc operator+(NamedFunc f, NamedFunc g){
  return f+=g;
}

NamedFunc operator-(NamedFunc f, NamedFunc g){
  return f-=g;
}

NamedFunc operator*(NamedFunc f, NamedFunc g){
  return f*=g;
}

NamedFunc operator/(NamedFunc f, NamedFunc g){
  return f/=g;
}

NamedFunc operator%(NamedFunc f, NamedFunc g){
  return f%=g;
}

NamedFunc operator + (NamedFunc f){
  f.Name("+(" + f.Name() + ")");
  return f;
}

NamedFunc operator - (NamedFunc f){
  f.Name("-(" + f.Name() + ")");
  f.Function(ApplyOp(f.ScalarFunction(), negate<ScalarType>()));
  f.Function(ApplyOp(f.VectorFunction(), negate<ScalarType>()));
  return f;
}

NamedFunc operator == (NamedFunc f, NamedFunc g){
  f.Name("(" + f.Name() + ")==(" + g.Name() + ")");
  auto fp = ApplyOp(f.ScalarFunction(), f.VectorFunction(),
                    g.ScalarFunction(), g.VectorFunction(),
                    equal_to<ScalarType>());
  f.Function(fp.first);
  f.Function(fp.second);
  return f;
}

NamedFunc operator != (NamedFunc f, NamedFunc g){
  f.Name("(" + f.Name() + ")!=(" + g.Name() + ")");
  auto fp = ApplyOp(f.ScalarFunction(), f.VectorFunction(),
                    g.ScalarFunction(), g.VectorFunction(),
                    not_equal_to<ScalarType>());
  f.Function(fp.first);
  f.Function(fp.second);
  return f;
}

NamedFunc operator > (NamedFunc f, NamedFunc g){
  f.Name("(" + f.Name() + ")>(" + g.Name() + ")");
  auto fp = ApplyOp(f.ScalarFunction(), f.VectorFunction(),
                    g.ScalarFunction(), g.VectorFunction(),
                    greater<ScalarType>());
  f.Function(fp.first);
  f.Function(fp.second);
  return f;
}

NamedFunc operator < (NamedFunc f, NamedFunc g){
  f.Name("(" + f.Name() + ")<(" + g.Name() + ")");
  auto fp = ApplyOp(f.ScalarFunction(), f.VectorFunction(),
                    g.ScalarFunction(), g.VectorFunction(),
                    less<ScalarType>());
  f.Function(fp.first);
  f.Function(fp.second);
  return f;
}

NamedFunc operator >= (NamedFunc f, NamedFunc g){
  f.Name("(" + f.Name() + ")>=(" + g.Name() + ")");
  auto fp = ApplyOp(f.ScalarFunction(), f.VectorFunction(),
                    g.ScalarFunction(), g.VectorFunction(),
                    greater_equal<ScalarType>());
  f.Function(fp.first);
  f.Function(fp.second);
  return f;
}

NamedFunc operator <= (NamedFunc f, NamedFunc g){
  f.Name("(" + f.Name() + ")<=(" + g.Name() + ")");
  auto fp = ApplyOp(f.ScalarFunction(), f.VectorFunction(),
                    g.ScalarFunction(), g.VectorFunction(),
                    less_equal<ScalarType>());
  f.Function(fp.first);
  f.Function(fp.second);
  return f;
}

NamedFunc operator && (NamedFunc f, NamedFunc g){
  f.Name("(" + f.Name() + ")&&(" + g.Name() + ")");
  auto fp = ApplyOp(f.ScalarFunction(), f.VectorFunction(),
                    g.ScalarFunction(), g.VectorFunction(),
                    logical_and<ScalarType>());
  f.Function(fp.first);
  f.Function(fp.second);
  return f;
}

NamedFunc operator || (NamedFunc f, NamedFunc g){
  f.Name("(" + f.Name() + ")||(" + g.Name() + ")");
  auto fp = ApplyOp(f.ScalarFunction(), f.VectorFunction(),
                    g.ScalarFunction(), g.VectorFunction(),
                    logical_or<ScalarType>());
  f.Function(fp.first);
  f.Function(fp.second);
  return f;
}

NamedFunc operator ! (NamedFunc f){
  f.Name("!(" + f.Name() + ")");
  f.Function(ApplyOp(f.ScalarFunction(), logical_not<ScalarType>()));
  f.Function(ApplyOp(f.VectorFunction(), logical_not<ScalarType>()));
  return f;
}

ostream & operator<<(ostream &stream, const NamedFunc &function){
  stream << function.Name();
  return stream;
}
