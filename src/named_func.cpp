/*! \class NamedFunc

  \brief Combines a callabale function taking a Baby and returning a scalar or
  vector with its string representation.

  NamedFunc contains both a callabale function which takes a Baby as a parameter
  and returns either a scalar or a vector result. It also contains a string
  representation of that function. Typically, this is a C++/"TTree::Draw"-like
  expression, but can be manually set to any desired string independently of the
  callable function. Given only a C++ expression, it is able to dynamically
  generate the appropriate callable and fully compiled function. The string
  parsing is done just once by FunctionParser, and the resulting callable
  function stored for fast calling without reinterpreting the string each time.

  \link NamedFunc NamedFuncs\endlink can be manipulated in much the same way as
  C++ arithmetic types. For example, given \link NamedFunc NamedFuncs\endlink x
  and y, x+y will return a NamedFunc z whose function evaluates the functions of
  x and y and returns the sum of the results. It is important to note that z's
  function does not simply return the result it obtains by evaluating x and y on
  construction. Rather, it remembers the functions from x and y, and reevaluates
  the component addends and resulting sum each time z is called. This allows
  construction of arbitrarily complicated functions by applying standard C++
  operators to simple \link NamedFunc NamedFuncs\endlink. This ability is used
  heavily by FunctionParser to build a single NamedFunc from complex
  expressions. Currently, the operators "+" (unary and binary), "-" (unary and
  binary), "*", "/", "%", "+=", "-=", "*=", "/=", "%=", "==", "!=", ">", "<",
  ">=", "<=", "&&", "||", and "!" are supported. Bit-level operators "<<", ">>",
  "~", "^", "^=", "&=", and "|=" and not supported. The "<<" is used for
  printing to an output stream.

  The current implementation keeps both a scalar and vector function internally,
  only one of which is valid at any time. To the scalar function is evaluated
  with NamedFunc::GetScalar(), while the vector function is evaluated with
  NamedFunc::GetVector(). A possible alternative is to always return a vector
  (of length 1 in the case of a scalar result), and use a bool to determine if
  the result should be considered a scalar syntactically. This would allow
  NamedFunc to act as a true functor with a cleaner interface, but results in
  extra vectors being constructed (and often copied if care is not taken with
  results) even when evaluating a simple scalar value.

  \see FunctionParser for allowed expression syntax for constructing a
  NamedFunc.
*/
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
  /*!\brief Get a functor applying unary operator op to f

    \param[in] f Function which takes a Baby and returns a single value

    \param[in] op Unary operator to apply to f

    \return Functor which takes a Baby and returns the result of applying op to
    the result of f
  */
  template<typename Operator>
    static function<ScalarFunc> ApplyOp(const function<ScalarFunc> &f,
                                        const Operator &op){
    if(!static_cast<bool>(f)) return f;
    function<ScalarType(ScalarType)> op_c(op);
    return [f,op_c](const Baby &b){
      return op_c(f(b));
    };
  }

  /*!\brief Get a functor applying unary operator op to f

    \param[in] f Function which takes a Baby and returns a vector of values

    \param[in] op Unary operator to apply to f

    \return Functor which takes a Baby and returns the result of applying op to
    each element of the result of f
  */
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

  /*!\brief Get a functor applying binary operator op to operands (sfa or vfa)
    and (sfb or vfb)

    sfa and vfa are associated to the same NamedFunc, and sfb and vfb are
    associated to the same NamedFunc. Exactly one from each pair should be a
    valid function. Determines the valid function from each pair and applies
    binary operator op between the "a" function result on the left and the "b"
    function result on the right.

    \param[in] sfa Scalar function from the same NamedFunc as vfa

    \param[in] vfa Vector function from the same NamedFunc as sfa

    \param[in] sfb Scalar function from the same NamedFunc as vfb

    \param[in] vfb Vector function from the same NamedFunc as sfb

    \param[in] op Unary operator to apply to (sfa or vfa) and (sfb or vfb)

    \return Functor which takes a Baby and returns the result of applying op to
    (sfa or vfa) and (sfb or vfb)
  */
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

/*!\brief Constructor of a scalar NamedFunc

  \param[in] name Text representation of function

  \param[in] function Functor taking a Baby and returning a scalar
*/
NamedFunc::NamedFunc(const std::string &name,
                     const std::function<ScalarFunc> &function):
  name_(name),
  scalar_func_(function),
  vector_func_(){
  CleanName();
}

/*!\brief Constructor of a vector NamedFunc

  \param[in] name Text representation of function

  \param[in] function Functor taking a Baby and returning a vector
*/
NamedFunc::NamedFunc(const std::string &name,
                     const std::function<VectorFunc> &function):
  name_(name),
  scalar_func_(),
  vector_func_(function){
  CleanName();
  }

/*!\brief Constructor using FunctionParser to produce a real function from a
  string

  \param[in] function C++/"TTree::Draw"-like expression containing constants,
  Baby variables, operators, parenthesis, brackets, etc.
*/
NamedFunc::NamedFunc(const string &function):
  NamedFunc(FunctionParser(function).ResolveAsNamedFunc()){
}

/*!\brief Constructor using FunctionParser to produce a real function from a
  string

  \param[in] function C++/"TTree::Draw"-like expression containing constants,
  Baby variables, operators, parenthesis, brackets, etc.
*/
NamedFunc::NamedFunc(const char *function):
  NamedFunc(string(function)){
}

/*!\brief Constructor for NamedFunc returning a constant

  \param[in] x The constant to be returned
*/
NamedFunc::NamedFunc(ScalarType x):
  name_(ToString(x)),
  scalar_func_([x](const Baby&){return x;}),
  vector_func_(){
}

/*!\brief Get the string representation of this function

  \return The standard string representation of this function
*/
const string & NamedFunc::Name() const{
  return name_;
}

/*!\brief Set the string representation of this function

  \param[in] name String representation of the function
*/
NamedFunc & NamedFunc::Name(const string &name){
  name_ = name;
  CleanName();
  return *this;
}

/*!\brief Get name stripped of special symbols for use in file names

  \return String representation of function with characters unsuitable for file
  names removed
*/
string NamedFunc::PlainName() const{
  string plain = name_;
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

/*!\brief Get name with common symbols pretty-printed for use in titles

  \return Pretty-printed version of name for use in title
*/
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

/*!\brief Set function to given scalar function

  This function overwrites the scalar function and invalidates the vector
  function if set.

  \param[in] f Valid function taking a Baby and returning a scalar

  \return Reference to *this
*/
NamedFunc & NamedFunc::Function(const std::function<ScalarFunc> &f){
  if(!static_cast<bool>(f)) return *this;
  scalar_func_ = f;
  vector_func_ = function<VectorFunc>();
  return *this;
}

/*!\brief Set function to given vector function

  This function overwrites the vector function and invalidates the scalar
  function if set.

  \param[in] f Valid function taking a Baby and returning a vector

  \return Reference to *this
*/
NamedFunc & NamedFunc::Function(const std::function<VectorFunc> &f){
  if(!static_cast<bool>(f)) return *this;
  scalar_func_ = function<ScalarFunc>();
  vector_func_ = f;
  return *this;
}

/*!\brief Return the (possibly invalid) scalar function

  \return The (possibly invalid) scalar function associated to *this
*/
const function<ScalarFunc> & NamedFunc::ScalarFunction() const{
  return scalar_func_;
}

/*!\brief Return the (possibly invalid) vector function

  \return The (possibly invalid) vector function associated to *this
*/
const function<VectorFunc> & NamedFunc::VectorFunction() const{
  return vector_func_;
}

/*!\brief Check if scalar function is valid

  \return True if scalar function is valid; false otherwise.
*/
bool NamedFunc::IsScalar() const{
  return static_cast<bool>(scalar_func_);
}

/*!\brief Check if vectorr function is valid

  \return True if vector function is valid; false otherwise.
*/
bool NamedFunc::IsVector() const{
  return static_cast<bool>(vector_func_);
}

/*!\brief Evaluate scalar function with b as argument

  \param[in] b Baby to pass to scalar function

  \return Result of applying scalar function to b
*/
ScalarType NamedFunc::GetScalar(const Baby &b) const{
  return scalar_func_(b);
}

/*!\brief Evaluate vector function with b as argument

  \param[in] b Baby to pass to vector function

  \return Result of applying vector function to b
*/
VectorType NamedFunc::GetVector(const Baby &b) const{
  return vector_func_(b);
}

/*!\brief Add func to *this

  \param[in] func Function to be added to *this

  \return Reference to *this
*/
NamedFunc & NamedFunc::operator += (const NamedFunc &func){
  name_ = "("+name_ + ")+(" + func.name_ + ")";
  auto fp = ApplyOp(scalar_func_, vector_func_,
                    func.scalar_func_, func.vector_func_,
                    plus<ScalarType>());
  scalar_func_ = fp.first;
  vector_func_ = fp.second;
  return *this;
}

/*!\brief Subtract func from *this

  \param[in] func Function to be subtracted from *this

  \return Reference to *this
*/
NamedFunc & NamedFunc::operator -= (const NamedFunc &func){
  name_ = "("+name_ + ")-(" + func.name_ + ")";
  auto fp = ApplyOp(scalar_func_, vector_func_,
                    func.scalar_func_, func.vector_func_,
                    minus<ScalarType>());
  scalar_func_ = fp.first;
  vector_func_ = fp.second;
  return *this;
}

/*!\brief Multiply *this by func

  \param[in] func Function by which *this is multiplied

  \return Reference to *this
*/
NamedFunc & NamedFunc::operator *= (const NamedFunc &func){
  name_ = "("+name_ + ")*(" + func.name_ + ")";
  auto fp = ApplyOp(scalar_func_, vector_func_,
                    func.scalar_func_, func.vector_func_,
                    multiplies<ScalarType>());
  scalar_func_ = fp.first;
  vector_func_ = fp.second;
  return *this;
}

/*!\brief Divide *this by func

  \param[in] func Function by which to divide *this

  \return Reference to *this
*/
NamedFunc & NamedFunc::operator /= (const NamedFunc &func){
  name_ = "("+name_ + ")/(" + func.name_ + ")";
  auto fp = ApplyOp(scalar_func_, vector_func_,
                    func.scalar_func_, func.vector_func_,
                    divides<ScalarType>());
  scalar_func_ = fp.first;
  vector_func_ = fp.second;
  return *this;
}

/*!\brief Set *this to remainder of *this divided by func

  \param[in] func Function with respect to which to take remainder

  \return Reference to *this
*/
NamedFunc & NamedFunc::operator %= (const NamedFunc &func){
  name_ = "("+name_ + ")%(" + func.name_ + ")";
  auto fp = ApplyOp(scalar_func_, vector_func_,
                    func.scalar_func_, func.vector_func_,
                    static_cast<ScalarType (*)(ScalarType ,ScalarType)>(fmod));
  scalar_func_ = fp.first;
  vector_func_ = fp.second;
  return *this;
}

/*!\brief Strip spaces from name
 */
void NamedFunc::CleanName(){
  ReplaceAll(name_, " ", "");
}

/*!\brief Add two \link NamedFunc NamedFuncs\endlink

  \param[in] f Augend

  \param[in] g Addend

  \return NamedFunc which returns the sum of the results of f and g
*/
NamedFunc operator+(NamedFunc f, NamedFunc g){
  return f+=g;
}

/*!\brief Add a NamedFunc from another

  \param[in] f Minuend

  \param[in] g Subtrahend

  \return NamedFunc which returns the difference of the results of f and g
*/
NamedFunc operator-(NamedFunc f, NamedFunc g){
  return f-=g;
}

/*!\brief Multiply two \link NamedFunc NamedFuncs\endlink

  \param[in] f Multiplier

  \param[in] g Multiplicand

  \return NamedFunc which returns the product of the results of f and g
*/
NamedFunc operator*(NamedFunc f, NamedFunc g){
  return f*=g;
}

/*!\brief Divide two \link NamedFunc NamedFuncs\endlink

  \param[in] f Dividend

  \param[in] g Divisor

  \return NamedFunc which returns the quotient of the results of f and g
*/
NamedFunc operator/(NamedFunc f, NamedFunc g){
  return f/=g;
}

/*!\brief Get remainder form division of two \link NamedFunc NamedFuncs\endlink

  \param[in] f Dividend

  \param[in] g Divisor

  \return NamedFunc which returns the remainder from dividing the results of f
  and g
*/
NamedFunc operator%(NamedFunc f, NamedFunc g){
  return f%=g;
}

/*!\brief Applied unary plus operator. Acts as identity operation.

  \param[in] f NamedFunc to apply unary "+" to

  \return f
*/
NamedFunc operator + (NamedFunc f){
  f.Name("+(" + f.Name() + ")");
  return f;
}

/*!\brief Negates f

  \param[in] f NamedFunc to apply unary "-" to

  \return NamedFunc returing the negative of the result of f
*/
NamedFunc operator - (NamedFunc f){
  f.Name("-(" + f.Name() + ")");
  f.Function(ApplyOp(f.ScalarFunction(), negate<ScalarType>()));
  f.Function(ApplyOp(f.VectorFunction(), negate<ScalarType>()));
  return f;
}

/*!\brief Gets NamedFunc which tests for equality of results of f and g

  \param[in] f Left hand operand

  \param[in] g Right hand operand

  \return NamedFunc returning whether the results of f and g are equal
*/
NamedFunc operator == (NamedFunc f, NamedFunc g){
  f.Name("(" + f.Name() + ")==(" + g.Name() + ")");
  auto fp = ApplyOp(f.ScalarFunction(), f.VectorFunction(),
                    g.ScalarFunction(), g.VectorFunction(),
                    equal_to<ScalarType>());
  f.Function(fp.first);
  f.Function(fp.second);
  return f;
}

/*!\brief Gets NamedFunc which tests for inequality of results of f and g

  \param[in] f Left hand operand

  \param[in] g Right hand operand

  \return NamedFunc returning whether the results of f and g are not equal
*/
NamedFunc operator != (NamedFunc f, NamedFunc g){
  f.Name("(" + f.Name() + ")!=(" + g.Name() + ")");
  auto fp = ApplyOp(f.ScalarFunction(), f.VectorFunction(),
                    g.ScalarFunction(), g.VectorFunction(),
                    not_equal_to<ScalarType>());
  f.Function(fp.first);
  f.Function(fp.second);
  return f;
}

/*!\brief Gets NamedFunc which tests if result of f is greater than result of g

  \param[in] f Left hand operand

  \param[in] g Right hand operand

  \return NamedFunc returning whether the results of f is greater than result of
  g
*/
NamedFunc operator > (NamedFunc f, NamedFunc g){
  f.Name("(" + f.Name() + ")>(" + g.Name() + ")");
  auto fp = ApplyOp(f.ScalarFunction(), f.VectorFunction(),
                    g.ScalarFunction(), g.VectorFunction(),
                    greater<ScalarType>());
  f.Function(fp.first);
  f.Function(fp.second);
  return f;
}

/*!\brief Gets NamedFunc which tests if result of f is less than result of g

  \param[in] f Left hand operand

  \param[in] g Right hand operand

  \return NamedFunc returning whether the results of f is less than result of g
*/
NamedFunc operator < (NamedFunc f, NamedFunc g){
  f.Name("(" + f.Name() + ")<(" + g.Name() + ")");
  auto fp = ApplyOp(f.ScalarFunction(), f.VectorFunction(),
                    g.ScalarFunction(), g.VectorFunction(),
                    less<ScalarType>());
  f.Function(fp.first);
  f.Function(fp.second);
  return f;
}

/*!\brief Gets NamedFunc which tests if result of f is greater than or equal to
  result of g

  \param[in] f Left hand operand

  \param[in] g Right hand operand

  \return NamedFunc returning whether the results of f is greater than or equal
  to result of g
*/
NamedFunc operator >= (NamedFunc f, NamedFunc g){
  f.Name("(" + f.Name() + ")>=(" + g.Name() + ")");
  auto fp = ApplyOp(f.ScalarFunction(), f.VectorFunction(),
                    g.ScalarFunction(), g.VectorFunction(),
                    greater_equal<ScalarType>());
  f.Function(fp.first);
  f.Function(fp.second);
  return f;
}

/*!\brief Gets NamedFunc which tests if result of f is less than or equal to
  result of g

  \param[in] f Left hand operand

  \param[in] g Right hand operand

  \return NamedFunc returning whether the results of f is less than or equal to
  result of g
*/
NamedFunc operator <= (NamedFunc f, NamedFunc g){
  f.Name("(" + f.Name() + ")<=(" + g.Name() + ")");
  auto fp = ApplyOp(f.ScalarFunction(), f.VectorFunction(),
                    g.ScalarFunction(), g.VectorFunction(),
                    less_equal<ScalarType>());
  f.Function(fp.first);
  f.Function(fp.second);
  return f;
}

/*!\brief Gets NamedFunc which tests if results of both f and g are true

  \param[in] f Left hand operand

  \param[in] g Right hand operand

  \return NamedFunc returning whether the results of both f and g are true
*/
NamedFunc operator && (NamedFunc f, NamedFunc g){
  f.Name("(" + f.Name() + ")&&(" + g.Name() + ")");
  auto fp = ApplyOp(f.ScalarFunction(), f.VectorFunction(),
                    g.ScalarFunction(), g.VectorFunction(),
                    logical_and<ScalarType>());
  f.Function(fp.first);
  f.Function(fp.second);
  return f;
}

/*!\brief Gets NamedFunc which tests if result of f or g is true

  \param[in] f Left hand operand

  \param[in] g Right hand operand

  \return NamedFunc returning whether the results of f or g is true
*/
NamedFunc operator || (NamedFunc f, NamedFunc g){
  f.Name("(" + f.Name() + ")||(" + g.Name() + ")");
  auto fp = ApplyOp(f.ScalarFunction(), f.VectorFunction(),
                    g.ScalarFunction(), g.VectorFunction(),
                    logical_or<ScalarType>());
  f.Function(fp.first);
  f.Function(fp.second);
  return f;
}

/*!\brief Gets NamedFunct returning logical inverse of result of f

  \param[in] f Function whose result the logical not is applied to

  \return NamedFunc returning logical inverse of result of f
*/
NamedFunc operator ! (NamedFunc f){
  f.Name("!(" + f.Name() + ")");
  f.Function(ApplyOp(f.ScalarFunction(), logical_not<ScalarType>()));
  f.Function(ApplyOp(f.VectorFunction(), logical_not<ScalarType>()));
  return f;
}

/*!\brief Print NamedFunc to output stream

  \param[in,out] stream Output stream to print to

  \param[in] function NamedFunc to print
*/
ostream & operator<<(ostream &stream, const NamedFunc &function){
  stream << function.Name();
  return stream;
}
