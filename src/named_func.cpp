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

NamedFunc::NamedFunc(const std::string &name,
                     const std::function<ScalarFunc> &function):
  name_(name),
  scalar_func_(function),
  vector_func_(){
  CleanName();
}

NamedFunc::NamedFunc(const std::string &name,
                     const std::function<VectorFunc> &function):
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

string NamedFunc::PrettyName() const{
  string pretty = name_;
  if(pretty == "1") return "";
  ReplaceAll(pretty, "1==1", "Full Sample");
  ReplaceAll(pretty, "el_tks_chg*lep_charge<0", "OS");
  ReplaceAll(pretty, "mu_tks_chg*lep_charge<0", "OS");
  ReplaceAll(pretty, "had_tks_chg*lep_charge<0", "OS");
  ReplaceAll(pretty, "Sum$(abs(mc_id)==11)","n^{true}_{e}");
  ReplaceAll(pretty, "Sum$(abs(mc_id)==13)","n^{true}_{#mu}");
  ReplaceAll(pretty, "Sum$(genels_pt>0)", "n^{true}_{e}");
  ReplaceAll(pretty, "Sum$(genmus_pt>0)", "n^{true}_{#mu}");
  ReplaceAll(pretty, "Sum$(mus_sigid&&mus_miniso<0.2)","n_{#mu}^{10}");
  ReplaceAll(pretty, "Sum$(els_sigid&&els_miniso<0.1)","n_{e}^{10}");
  ReplaceAll(pretty, "nvmus==1&&nmus==1&&nvels==0","1 #mu");
  ReplaceAll(pretty, "nvmus10==0&&nvels10==0", "0 leptons");
  ReplaceAll(pretty, "(nmus+nels)", "n_{lep}");
  ReplaceAll(pretty, "(nels+nmus)", "n_{lep}");
  ReplaceAll(pretty, "(nvmus+nvels)", "n^{veto}_{lep}");
  ReplaceAll(pretty, "nvmus==2&&nmus>=1","n_{#mu}#geq1, n^{veto}_{#mu}=2");
  ReplaceAll(pretty, "nvels==2&&nels>=1","n_{e}#geq1, n^{veto}_{e}=2");
  ReplaceAll(pretty, "(nvmus>=2||nvels>=2)","n^{veto}_{lep} #geq 2");
  ReplaceAll(pretty, "(mumu_m*(mumu_m>0)+elel_m*(elel_m>0))>80&&(mumu_m*(mumu_m>0)+elel_m*(elel_m>0))<100",
             "80<m_{ll}<100");
  ReplaceAll(pretty, "mumuv_m>80&&mumuv_m<100",
             "80<m_{ll}<100");
  ReplaceAll(pretty, "elelv_m>80&&elelv_m<100",
             "80<m_{ll}<100");
  ReplaceAll(pretty, "onht>350&&onmet>100&&","");
  ReplaceAll(pretty, "jets_islep[0]==0","");
  ReplaceAll(pretty, "(nels==0&&nmus==1)","n_{#mu}=1");
  ReplaceAll(pretty, "(nels==1&&nmus==0)","n_{#font[12]{e}}=1");
  ReplaceAll(pretty, "Max$(abs(els_eta)*(els_sigid&&els_miniso<0.1&&els_pt>20))<1.479","barrel #font[12]{e}");
  ReplaceAll(pretty, "Max$(abs(els_eta)*(els_sigid&&els_miniso<0.1&&els_pt>20))>1.479","endcap #font[12]{e}");

  ReplaceAll(pretty, "nleps", "n_{l}");
  ReplaceAll(pretty, "nvleps", "n_{l}");
  ReplaceAll(pretty, "nmus", "n_{#mu}");
  ReplaceAll(pretty, "nels", "n_{e}");
  ReplaceAll(pretty, "nvmus", "n^{veto}_{#mu}");
  ReplaceAll(pretty, "nvels", "n^{veto}_{e}");
  ReplaceAll(pretty, "ntruleps", "n^{true}_{l}");
  ReplaceAll(pretty, "_ra2b", "^{ra2b}");
  ReplaceAll(pretty, "npv", "n_{PV}");
  ReplaceAll(pretty, "mumu_pt1", "p_{T}^{#mu}");
  ReplaceAll(pretty, "elel_pt1", "p_{T}^{e}");

  ReplaceAll(pretty, "abs(mc_id)==1000006", "stop");
  ReplaceAll(pretty, "abs(mc_id)==1000022", "LSP");

  ReplaceAll(pretty, "onmet", "MET^{on}");
  ReplaceAll(pretty, "onht", "H_{T}^{on}");
  ReplaceAll(pretty, "njets30","n_{jets}^{30}");
  ReplaceAll(pretty, "els_pt","p^{e}_{T}");
  ReplaceAll(pretty, "mus_pt","p^{#mu}_{T}");
  ReplaceAll(pretty, "fjets_nconst","n_{const}^{fat jet}");
  ReplaceAll(pretty, "fjets_30_m[0]","m(J_{1})");
  ReplaceAll(pretty, "fjets_m[0]","m(J_{1})");
  ReplaceAll(pretty, "(fjets_pt*cosh(fjets_eta))","p_{fatjet}");
  ReplaceAll(pretty, "fjets_pt","p^{fatjet}_{T}");
  ReplaceAll(pretty, "jets_pt","p^{jet}_{T}");
  ReplaceAll(pretty, "mus_reliso","RelIso");
  ReplaceAll(pretty, "els_reliso","RelIso");
  ReplaceAll(pretty, "mus_miniso_tr15","MiniIso");
  ReplaceAll(pretty, "els_miniso_tr15","MiniIso");
  ReplaceAll(pretty, "njets","n_{jets}");
  ReplaceAll(pretty, "abs(lep_id)==13&&","");
  ReplaceAll(pretty, ">=", " #geq ");
  ReplaceAll(pretty, ">", " > ");
  ReplaceAll(pretty, "<=", " #leq ");
  ReplaceAll(pretty, "<", " < ");
  ReplaceAll(pretty, "&&", ", ");
  ReplaceAll(pretty, "==", " = ");
  ReplaceAll(pretty, "met", "MET");
  ReplaceAll(pretty, "ht_hlt", "H_{T}^{HLT}");
  ReplaceAll(pretty, "mht", "MHT");
  ReplaceAll(pretty, "ht", "H_{T}");
  ReplaceAll(pretty, "mt", "m_{T}");
  ReplaceAll(pretty, "ntks_chg==0", " ITV");
  ReplaceAll(pretty, "nbm","n_{b}");
  ReplaceAll(pretty, "nbl","n_{b,l}");
  ReplaceAll(pretty, "mj", " M_{J}");

  ReplaceAll(pretty, "el_tks_mt", "Track m_{T}");
  ReplaceAll(pretty, "mu_tks_mt", "Track m_{T}");
  ReplaceAll(pretty, "had_tks_mt", "Track m_{T}");
  ReplaceAll(pretty, "el_tks_pt", "Track p_{T}");
  ReplaceAll(pretty, "mu_tks_pt", "Track p_{T}");
  ReplaceAll(pretty, "had_tks_pt", "Track p_{T}");
  ReplaceAll(pretty, "el_tks_miniso", "Track miniso");
  ReplaceAll(pretty, "mu_tks_miniso", "Track miniso");
  ReplaceAll(pretty, "had_tks_miniso", "Track miniso");
  ReplaceAll(pretty, "el_tks_chg_miniso", "Track charged miniso");
  ReplaceAll(pretty, "mu_tks_chg_miniso", "Track charged miniso");
  ReplaceAll(pretty, "had_tks_chg_miniso", "Track charged miniso");
  return pretty;
}

NamedFunc & NamedFunc::Function(const std::function<ScalarFunc> &f){
  if(!static_cast<bool>(f)) return *this;
  scalar_func_ = f;
  vector_func_ = function<VectorFunc>();
  return *this;
}

NamedFunc & NamedFunc::Function(const std::function<VectorFunc> &f){
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
