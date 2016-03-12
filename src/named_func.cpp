#include "named_func.hpp"

#include <iostream>

#include "utilities.hpp"
#include "function_parser.hpp"

using namespace std;

namespace{
  double DoubleModulus(double x, double y){
    return fmod(x,y);
  }
  
  template<typename Operator>
    static function<NamedFunc::FuncType> ApplyOp(const function<NamedFunc::FuncType> &f, bool /*is_vector*/,
                                                 const Operator &op){
    return [f,op](const Baby &b){
      auto result = f(b);
      for(auto &x: result){
        op(x);
      }
      return result;
    };
  }
  
  template<typename Operator>
    static function<NamedFunc::FuncType> ApplyOp(const function<NamedFunc::FuncType> &f, bool f_is_vector,
                                                 const function<NamedFunc::FuncType> &g, bool g_is_vector,
                                                 const Operator &op){
    if(f_is_vector && g_is_vector){
      return [f,g,op](const Baby &b){
        auto result_f = f(b);
        auto result_g = g(b);
        size_t num_elements = min(result_f.size(), result_g.size());
        vector<function<NamedFunc::FuncType>::result_type::value_type> result(num_elements);
        for(size_t i = 0; i < num_elements; ++i){
          result[i] = op(result_f[i], result_g[i]);
        }
        return result;
      };
    }else if(f_is_vector && !g_is_vector){
      return [f,g,op](const Baby &b){
        auto result_f = f(b);
        auto vector_g = g(b);
        auto result_g = vector_g.size() ? vector_g[0] : 0.;
        size_t num_elements = result_f.size();
        vector<function<NamedFunc::FuncType>::result_type::value_type> result(num_elements);
        for(size_t i = 0; i < num_elements; ++i){
          result[i] = op(result_f[i], result_g);
        }
        return result;
      };
    }else if(!f_is_vector && g_is_vector){
      return [f,g,op](const Baby &b){
        auto vector_f = f(b);
        auto result_f = vector_f.size() ? vector_f[0] : 0.;
        auto result_g = g(b);
        size_t num_elements = result_g.size();
        vector<function<NamedFunc::FuncType>::result_type::value_type> result(num_elements);
        for(size_t i = 0; i < num_elements; ++i){
          result[i] = op(result_f, result_g[i]);
        }
        return result;
      };
    }else{
      return [f,g,op](const Baby &b){
        auto vector_f = f(b);
        auto result_f = vector_f.size() ? vector_f[0] : 0.;
        auto vector_g = g(b);
        auto result_g = vector_g.size() ? vector_g[0] : 0.;
        size_t num_elements = 1;
        vector<function<NamedFunc::FuncType>::result_type::value_type> result(num_elements);
        for(size_t i = 0; i < num_elements; ++i){
          result[i] = op(result_f, result_g);
        }
        return result;
      };
    }
  }
}

NamedFunc::NamedFunc(const string &name,
                     const function<FuncType> &function,
                     bool is_vector):
  name_(name),
  function_(function),
  is_vector_(is_vector){
  CleanName();
  }

NamedFunc::NamedFunc(const string &function):
  NamedFunc(FunctionParser(function).ResolveAsNamedFunc()){
}

NamedFunc::NamedFunc(const char *function):
  NamedFunc(string(function)){
}

NamedFunc::NamedFunc(double x):
  name_(ToString(x)),
  function_([x](const Baby&){return vector<double>(1,x);}),
  is_vector_(false){
}

const string & NamedFunc::Name() const{
  return name_;
}

NamedFunc & NamedFunc::Name(const string &name){
  name_ = name;
  CleanName();
  return *this;
}

const std::function<NamedFunc::FuncType> & NamedFunc::Function() const{
  return function_;
}

NamedFunc & NamedFunc::Function(const std::function<FuncType> & function, bool is_vector){
  function_ = function;
  is_vector_ = is_vector;
  return *this;
}

bool NamedFunc::IsVector() const{
  return is_vector_;
}

bool NamedFunc::IsScalar() const{
  return !is_vector_;
}

function<NamedFunc::FuncType>::result_type NamedFunc::operator()(const Baby &b){
  return function_(b);
}

function<NamedFunc::FuncType>::result_type::value_type NamedFunc::operator()(const Baby &b, size_t i){
  auto result = function_(b);
  return is_vector_ ? result[i] : (result.size() ? result[0] : 0.);
}

NamedFunc & NamedFunc::operator += (const NamedFunc &func){
  name_ = "("+name_ + ")+(" + func.name_ + ")";
  function_ = ApplyOp(function_, is_vector_,
                      func.function_, func.is_vector_,
                      std::plus<double>());
  is_vector_ = is_vector_ || func.is_vector_;
  return *this;
}

NamedFunc & NamedFunc::operator -= (const NamedFunc &func){
  name_ = "("+name_ + ")-(" + func.name_ + ")";
  function_ = ApplyOp(function_, is_vector_,
                      func.function_, func.is_vector_,
                      std::minus<double>());
  is_vector_ = is_vector_ || func.is_vector_;
  return *this;
}

NamedFunc & NamedFunc::operator *= (const NamedFunc &func){
  name_ = "("+name_ + ")*(" + func.name_ + ")";
  function_ = ApplyOp(function_, is_vector_,
                      func.function_, func.is_vector_,
                      std::multiplies<double>());
  is_vector_ = is_vector_ || func.is_vector_;
  return *this;
}

NamedFunc & NamedFunc::operator /= (const NamedFunc &func){
  name_ = "("+name_ + ")/(" + func.name_ + ")";
  function_ = ApplyOp(function_, is_vector_,
                      func.function_, func.is_vector_,
                      std::divides<double>());
  is_vector_ = is_vector_ || func.is_vector_;
  return *this;
}

NamedFunc & NamedFunc::operator %= (const NamedFunc &func){
  name_ = "("+name_ + ")%(" + func.name_ + ")";
  function_ = ApplyOp(function_, is_vector_,
                      func.function_, func.is_vector_,
                      DoubleModulus);
  is_vector_ = is_vector_ || func.is_vector_;
  return *this;
}

void NamedFunc::CleanName(){
  ReplaceAll(name_, " ", "");
  //ReplaceAll(name_, "b.", "");
  //ReplaceAll(name_, "()", "");
  //ReplaceAll(name_, ".", "p");
  //ReplaceAll(name_, "(", "OP");
  //ReplaceAll(name_, ")", "CP");
  //ReplaceAll(name_, "[", "OB");
  //ReplaceAll(name_, "]", "CB");
  //ReplaceAll(name_, "{", "OC");
  //ReplaceAll(name_, "}", "CC");
  //ReplaceAll(name_, "+", "PLS");
  //ReplaceAll(name_, "-", "MNS");
  //ReplaceAll(name_, "*", "TMS");
  //ReplaceAll(name_, "/", "DIV");
  //ReplaceAll(name_, "%", "MOD");
  //ReplaceAll(name_, "!", "NOT");
  //ReplaceAll(name_, "&&", "AND");
  //ReplaceAll(name_, "||", "OR");
  //ReplaceAll(name_, "==", "EQL");
  //ReplaceAll(name_, "<=", "GEQ");
  //ReplaceAll(name_, ">=", "LEQ");
  //ReplaceAll(name_, ">", "GTR");
  //ReplaceAll(name_, "<", "LES");
  //ReplaceAll(name_, "=", "EQL");
  //ReplaceAll(name_, "&", "BITAND");
  //ReplaceAll(name_, "|", "BITOR");
  //ReplaceAll(name_, "^", "BITXOR");
  //ReplaceAll(name_, "~", "BITNOT");
  //ReplaceAll(name_, "__", "_");
  //for(size_t i = 0; i < name_.size(); ++i){
  //  if(isalnum(name_.at(i)) || name_.at(i) == '.' || name_.at(i) == '_') continue;
  //  name_ = name_.substr(0,i)+name_.substr(i+1);
  //}
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
  f.Function(ApplyOp(f.Function(), f.IsVector(), std::negate<double>()),
             f.IsVector());
  return f;
}

NamedFunc operator == (NamedFunc f, NamedFunc g){
  f.Name("(" + f.Name() + ")==(" + g.Name() + ")");
  f.Function(ApplyOp(f.Function(), f.IsVector(),
                    g.Function(), g.IsVector(),
                    std::equal_to<double>()),
             f.IsVector() || g.IsVector());
  return f;
}

NamedFunc operator != (NamedFunc f, NamedFunc g){
  f.Name("(" + f.Name() + ")!=(" + g.Name() + ")");
  f.Function(ApplyOp(f.Function(), f.IsVector(),
                    g.Function(), g.IsVector(),
                    std::not_equal_to<double>()),
             f.IsVector() || g.IsVector());
  return f;
}

NamedFunc operator > (NamedFunc f, NamedFunc g){
  f.Name("(" + f.Name() + ")>(" + g.Name() + ")");
  f.Function(ApplyOp(f.Function(), f.IsVector(),
                    g.Function(), g.IsVector(),
                    std::greater<double>()),
             f.IsVector() || g.IsVector());
  return f;
}

NamedFunc operator < (NamedFunc f, NamedFunc g){
  f.Name("(" + f.Name() + ")<(" + g.Name() + ")");
  f.Function(ApplyOp(f.Function(), f.IsVector(),
                    g.Function(), g.IsVector(),
                    std::less<double>()),
             f.IsVector() || g.IsVector());
  return f;
}

NamedFunc operator >= (NamedFunc f, NamedFunc g){
  f.Name("(" + f.Name() + ")>=(" + g.Name() + ")");
  f.Function(ApplyOp(f.Function(), f.IsVector(),
                    g.Function(), g.IsVector(),
                    std::greater_equal<double>()),
             f.IsVector() || g.IsVector());
  return f;
}

NamedFunc operator <= (NamedFunc f, NamedFunc g){
  f.Name("(" + f.Name() + ")<=(" + g.Name() + ")");
  f.Function(ApplyOp(f.Function(), f.IsVector(),
                    g.Function(), g.IsVector(),
                    std::less_equal<double>()),
             f.IsVector() || g.IsVector());
  return f;
}

NamedFunc operator && (NamedFunc f, NamedFunc g){
  f.Name("(" + f.Name() + ")&&(" + g.Name() + ")");
  f.Function(ApplyOp(f.Function(), f.IsVector(),
                    g.Function(), g.IsVector(),
                    std::logical_and<double>()),
             f.IsVector() || g.IsVector());
  return f;
}

NamedFunc operator || (NamedFunc f, NamedFunc g){
  f.Name("(" + f.Name() + ")||(" + g.Name() + ")");
  f.Function(ApplyOp(f.Function(), f.IsVector(),
                    g.Function(), g.IsVector(),
                    std::logical_or<double>()),
             f.IsVector() || g.IsVector());
  return f;
}

NamedFunc operator ! (NamedFunc f){
  f.Name("!(" + f.Name() + ")");
  f.Function(ApplyOp(f.Function(), f.IsVector(), std::logical_not<double>()),
             f.IsVector());
  return f;
}

ostream & operator<<(ostream &stream, const NamedFunc &function){
  stream << function.Name();
  return stream;
}
