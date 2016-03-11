#include "named_func.hpp"

#include <iostream>

#include "TTreeFormula.h"

#include "utilities.hpp"

using namespace std;

NamedFunc::NamedFunc(const string &name,
                     const function<FuncType> &function):
  name_(name),
  function_(function){
  CleanName();
  }

NamedFunc::NamedFunc(const string &function):
  NamedFunc(function.c_str()){
}

NamedFunc::NamedFunc(const char *function):
  name_(function),
  function_(){
  map<TChain *, pair<shared_ptr<TTreeFormula>, int> > fmap;
  function_ = [fmap,function](const Baby &b) mutable {
    TChain *tree_ptr = b.GetTree().get();
    auto loc = fmap.find(tree_ptr);
    if(loc == fmap.end()){
      fmap[tree_ptr] = make_pair(make_shared<TTreeFormula>("",function,tree_ptr), -1);
      fmap.at(tree_ptr).first->SetQuickLoad(true);
    }
    pair<shared_ptr<TTreeFormula>, int> &entry = fmap.at(tree_ptr);
    shared_ptr<TTreeFormula> &ttf = entry.first;
    int &last_num = entry.second;
    int tree_num = tree_ptr->GetTreeNumber();
    if(tree_num != last_num){
      ttf->UpdateFormulaLeaves();
      last_num = tree_num;
    }
    return static_cast<double>(ttf->EvalInstance());
  };
}

NamedFunc::NamedFunc(double x):
  name_(ToString(x)),
  function_([x](const Baby&){return x;}){
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

NamedFunc & NamedFunc::Function(const std::function<FuncType> & function){
  function_ = function;
  return *this;
}

double NamedFunc::operator()(const Baby &b){
  return function_(b);
}

NamedFunc & NamedFunc::operator += (const NamedFunc &func){
  name_ = name_ + "PLS" + func.name_;
  auto func_a = function_;
  auto func_b = func.function_;
  function_ = [func_a, func_b](const Baby &b){return func_a(b)+func_b(b);};
  return *this;
}

NamedFunc & NamedFunc::operator -= (const NamedFunc &func){
  name_ = name_ + "MNS" + func.name_;
  auto func_a = function_;
  auto func_b = func.function_;
  function_ = [func_a, func_b](const Baby &b){return func_a(b)-func_b(b);};
  return *this;
}

NamedFunc & NamedFunc::operator *= (const NamedFunc &func){
  name_ = name_ + "TMS" + func.name_;
  auto func_a = function_;
  auto func_b = func.function_;
  function_ = [func_a, func_b](const Baby &b){return func_a(b)*func_b(b);};
  return *this;
}

NamedFunc & NamedFunc::operator /= (const NamedFunc &func){
  name_ = name_ + "DIV" + func.name_;
  auto func_a = function_;
  auto func_b = func.function_;
  function_ = [func_a, func_b](const Baby &b){return func_a(b)/func_b(b);};
  return *this;
}

NamedFunc NamedFunc::operator+() const{
  NamedFunc copy = *this;
  copy.name_ = "POS" + name_;
  auto orig = function_;
  copy.function_ = [orig](const Baby &b){return +orig(b);};
  return copy;
}

NamedFunc NamedFunc::operator-() const{
  NamedFunc copy = *this;
  copy.name_ = "NEG" + name_;
  auto orig = function_;
  copy.function_ = [orig](const Baby &b){return -orig(b);};
  return copy;
}

NamedFunc NamedFunc::operator!() const{
  NamedFunc copy = *this;
  copy.name_ = "NOT" + name_;
  auto orig = function_;
  copy.function_ = [orig](const Baby &b){return !orig(b);};
  return copy;
}

void NamedFunc::CleanName(){
  ReplaceAll(name_, " ", "");
  ReplaceAll(name_, "b.", "");
  ReplaceAll(name_, "()", "");
  ReplaceAll(name_, ".", "p");
  ReplaceAll(name_, "(", "OP");
  ReplaceAll(name_, ")", "CP");
  ReplaceAll(name_, "[", "OB");
  ReplaceAll(name_, "]", "CB");
  ReplaceAll(name_, "{", "OC");
  ReplaceAll(name_, "}", "CC");
  ReplaceAll(name_, "+", "PLS");
  ReplaceAll(name_, "-", "MNS");
  ReplaceAll(name_, "*", "TMS");
  ReplaceAll(name_, "/", "DIV");
  ReplaceAll(name_, "%", "MOD");
  ReplaceAll(name_, "!", "NOT");
  ReplaceAll(name_, "&&", "AND");
  ReplaceAll(name_, "||", "OR");
  ReplaceAll(name_, "==", "EQL");
  ReplaceAll(name_, "<=", "GEQ");
  ReplaceAll(name_, ">=", "LEQ");
  ReplaceAll(name_, ">", "GTR");
  ReplaceAll(name_, "<", "LES");
  ReplaceAll(name_, "=", "EQL");
  ReplaceAll(name_, "&", "BITAND");
  ReplaceAll(name_, "|", "BITOR");
  ReplaceAll(name_, "^", "BITXOR");
  ReplaceAll(name_, "~", "BITNOT");
  ReplaceAll(name_, "__", "_");
  for(size_t i = 0; i < name_.size(); ++i){
    if(isalnum(name_.at(i)) || name_.at(i) == '.' || name_.at(i) == '_') continue;
    name_ = name_.substr(0,i)+name_.substr(i+1);
  }
}

NamedFunc operator+(NamedFunc a, NamedFunc b){
  return a+=b;
}

NamedFunc operator-(NamedFunc a, NamedFunc b){
  return a-=b;
}

NamedFunc operator*(NamedFunc a, NamedFunc b){
  return a*=b;
}

NamedFunc operator/(NamedFunc a, NamedFunc b){
  return a/=b;
}

NamedFunc operator&&(NamedFunc a, NamedFunc b){
  string name = a.Name() + "AND" + b.Name();
  auto fa = a.Function();
  auto fb = b.Function();
  function<NamedFunc::FuncType> func([fa,fb](const Baby &baby){return fa(baby)&&fb(baby);});
  return NamedFunc(name,func);
}

NamedFunc operator||(NamedFunc a, NamedFunc b){
  string name = a.Name() + "OR" + b.Name();
  auto fa = a.Function();
  auto fb = b.Function();
  function<NamedFunc::FuncType> func([fa,fb](const Baby &baby){return fa(baby)||fb(baby);});
  return NamedFunc(name,func);
}

ostream & operator<<(ostream &stream, const NamedFunc &function){
  stream << function.Name();
  return stream;
}
