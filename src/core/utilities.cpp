#include "core/utilities.hpp"

#include <cstdlib>
#include <cstring>
#include <cstdint>
#include <cmath>

#include <sstream>
#include <string>
#include <vector>
#include <memory>

#include <unistd.h>
#include <glob.h>
#include <libgen.h>
#include <sys/stat.h>

#include "TMath.h"
#include "TCanvas.h"
#include "TVector2.h"
#include "TH1D.h"
#include "RooStats/RooStatsUtils.h"

using namespace std;

mutex Multithreading::root_mutex;

set<string> Glob(const string &pattern){
  glob_t glob_result;
  glob(pattern.c_str(), GLOB_TILDE, nullptr, &glob_result);
  set<string> ret;
  for(size_t i=0; i<glob_result.gl_pathc; ++i){
    ret.emplace(realpath(glob_result.gl_pathv[i], nullptr));
  }
  globfree(&glob_result);
  return ret;
}

string Basename(const string &filename){
  vector<char> c(filename.cbegin(), filename.cend());
  c.push_back(0);
  return string(basename(&c.at(0)));
}

bool Contains(const string &str, const string &pat){
  return str.find(pat) != string::npos;
}

bool StartsWith(const string &str, const string &pat){
  return str.find(pat) == 0;
}

void ReplaceAll(string &str, const string &orig, const string &rep){
  size_t loc = 0;
  while ((loc = str.find(orig, loc)) != string::npos) {
    str.replace(loc, orig.length(), rep);
    loc += rep.length();
  }
}

string CopyReplaceAll(string str, const string &orig, const string &rep){
  ReplaceAll(str, orig, rep);
  return str;
}

string LeftStrip(string s){
  s.erase(s.begin(), find_if(s.begin(), s.end(), [](char c){return isspace(c)==0;}));
  return s;
}

string RightStrip(string s){
  s.erase(find_if(s.rbegin(), s.rend(), [](char c){return isspace(c)==0;}).base(), s.end());
  return s;
}

string Strip(string s){
  return LeftStrip(RightStrip(s));
}

bool FileExists(const string &path){
  struct stat buffer;
  return (stat (path.c_str(), &buffer) == 0);
}

string execute(const string &cmd){
  FILE *pipe = popen(cmd.c_str(), "r");
  if(!pipe) throw runtime_error("Could not open pipe.");
  const size_t buffer_size = 128;
  char buffer[buffer_size];
  string result = "";
  while(!feof(pipe)){
    if(fgets(buffer, buffer_size, pipe) != NULL) result += buffer;
  }

  pclose(pipe);
  return result;
}

string CodeToPlainText(string code){
  ReplaceAll(code, " ", "");
  ReplaceAll(code, "((leps_pt[0]>30&&leps_eta[0]<2.1&&leps_eta[0]>-2.1)||(leps_pt[1]>30&&leps_eta[1]<2.1&&leps_eta[1]>-2.1))", "SLtrig");
  ReplaceAll(code, "(mumu_m*(mumu_m>0&&mumu_pt1>30)+elel_m*(elel_m>0&&elel_pt1>30))>80&&(mumu_m*(mumu_m>0&&mumu_pt1>30)+elel_m*(elel_m>0&&elel_pt1>30))<100", "zmasswindow");
  ReplaceAll(code, "nbt==2&&nbm==2","2b");
  ReplaceAll(code, "nbt>=2&&nbm>=3","ge3b");
  ReplaceAll(code, "nbt>=2&&nbm==3&&nbl==3","3b");
  ReplaceAll(code, "nbt>=2&&nbm>=3&&nbl>=4","4b");
  ReplaceAll(code, "pass&&stitch&&nvleps==0&&ntks==0&&!low_dphi&&njets>=4&&njets<=5","preseln");
  ReplaceAll(code, "hig_am>100&&hig_am<=140&&hig_dm<=40","HIG");
  ReplaceAll(code, "(hig_am<=100||hig_am>140||hig_dm>40)","SBD");
  ReplaceAll(code, "mumu_pt*(mumu_pt>0)+elel_pt*(elel_pt>0)","zpt");
  ReplaceAll(code, ".", "p");
  ReplaceAll(code, "(", "");
  ReplaceAll(code, ")", "");
  ReplaceAll(code, "[", "");
  ReplaceAll(code, "]", "");
  ReplaceAll(code, "{", "");
  ReplaceAll(code, "}", "");
  ReplaceAll(code, "+", "p");
  ReplaceAll(code, "-", "m");
  ReplaceAll(code, "*", "x");
  ReplaceAll(code, "/", "d");
  ReplaceAll(code, "%", "_");
  ReplaceAll(code, "!", "n");
  ReplaceAll(code, "&&", "_");
  ReplaceAll(code, "||", "_");
  ReplaceAll(code, "==", "");
  ReplaceAll(code, "<=", "le");
  ReplaceAll(code, ">=", "ge");
  ReplaceAll(code, ">", "g");
  ReplaceAll(code, "<", "l");
  ReplaceAll(code, "=", "");
  ReplaceAll(code, "&", "_");
  ReplaceAll(code, "|", "_");
  ReplaceAll(code, "^", "_");
  ReplaceAll(code, "~", "_");
  ReplaceAll(code, "___", "__");
  for(size_t i = 0; i < code.size(); ++i){
    if(isalnum(code.at(i)) || code.at(i) == '.' || code.at(i) == '_') continue;
    code = code.substr(0,i)+code.substr(i+1);
  }

  return code;
}

string CodeToRootTex(string code){
  ReplaceAll(code, " ", "");
  ReplaceAll(code, "pass&&stitch", "ps");
  ReplaceAll(code, "&&1", "");
  ReplaceAll(code, "weight", "w");
  ReplaceAll(code, "nbt<=1&&nbm==2","TM+MM");
  ReplaceAll(code, "nbt==2&&nbm==2","2b");
  ReplaceAll(code, "nbt>=2&&nbm>=3","3+b");
  ReplaceAll(code, "nbt>=2&&nbm==3&&nbl==3","3b");
  ReplaceAll(code, "nbt>=2&&nbm>=3&&nbl>=4","4b");
  ReplaceAll(code, "hig_am>100&&hig_am<=140&&hig_dm<=40","HIG");
  ReplaceAll(code, "(hig_am<=100||hig_am>140||hig_dm>40)","SBD");
  ReplaceAll(code, "ht1l_stmin2l", "1l: H_{T}>500, 2l: H_{T} + p_{T}^{l,min}");
  ReplaceAll(code, "ht1l_stmax2l", "1l: H_{T}>500, 2l: H_{T} + p_{T}^{l,max}");
  ReplaceAll(code, "ht1l_stave2l", "1l: H_{T}>500, 2l: H_{T} + p_{T}^{l,ave}");
  ReplaceAll(code, "st", "S_{T}");

  ReplaceAll(code, "met>100&&met<=150", "100<met<=150");
  ReplaceAll(code, "met>100&&met<=200", "100<met<=200");
  ReplaceAll(code, "met>150&&met<=200", "150<met<=200");
  ReplaceAll(code, "met>200&&met<=300", "200<met<=300");
  ReplaceAll(code, "met>200&&met<=350", "200<met<=350");
  ReplaceAll(code, "met>300&&met<=500", "300<met<=500");
  ReplaceAll(code, "met>350&&met<=500", "350<met<=500");
  ReplaceAll(code, "met>200&&met<=500", "200<met<=500");
  ReplaceAll(code, "njets>=4&&njets<=5", "4<=njets<=5");
  ReplaceAll(code, "njets>=5&&njets<=7", "5<=njets<=7");
  ReplaceAll(code, "njets>=6&&njets<=8", "6<=njets<=8");
  ReplaceAll(code, "nbm>=1&&nbm<=2", "1<=nbm<=2");

  ReplaceAll(code, "1==1", "Full Sample");
  ReplaceAll(code, "el_tks_chg*lep_charge<0", "OS");
  ReplaceAll(code, "mu_tks_chg*lep_charge<0", "OS");
  ReplaceAll(code, "had_tks_chg*lep_charge<0", "OS");
  ReplaceAll(code, "Sum$(abs(mc_id)==11)","N^{true}_{e}");
  ReplaceAll(code, "Sum$(abs(mc_id)==13)","N^{true}_{#mu}");
  ReplaceAll(code, "Sum$(genels_pt>0)", "N^{true}_{e}");
  ReplaceAll(code, "Sum$(genmus_pt>0)", "N^{true}_{#mu}");
  ReplaceAll(code, "Sum$(mus_sigid&&mus_miniso<0.2)","N_{#mu}^{10}");
  ReplaceAll(code, "Sum$(els_sigid&&els_miniso<0.1)","N_{e}^{10}");
  ReplaceAll(code, "nvmus==1&&nmus==1&&nvels==0","1 #mu");
  ReplaceAll(code, "nvmus10==0&&nvels10==0", "0 leptons");
  ReplaceAll(code, "(nmus+nels)", "N_{lep}");
  ReplaceAll(code, "(nels+nmus)", "N_{lep}");
  ReplaceAll(code, "nveto", "N^{veto}_{tks}");
  ReplaceAll(code, "(nvmus+nvels)", "N^{veto}_{lep}");
  ReplaceAll(code, "nvmus==2&&nmus>=1","N_{#mu}#geq1, N^{veto}_{#mu}=2");
  ReplaceAll(code, "nvels==2&&nels>=1","N_{e}#geq1, N^{veto}_{e}=2");
  ReplaceAll(code, "(nvmus>=2||nvels>=2)","N^{veto}_{lep} #geq 2");
  ReplaceAll(code, "(mumu_pt*(mumu_pt>0)+elel_pt*(elel_pt>0))>150&&(mumu_pt*(mumu_pt>0)+elel_pt*(elel_pt>0))<=200",
                   "150<p_{T}^{Z}#leq200");
  ReplaceAll(code, "(mumu_pt*(mumu_pt>0)+elel_pt*(elel_pt>0))>200&&(mumu_pt*(mumu_pt>0)+elel_pt*(elel_pt>0))<=300",
                   "200<p_{T}^{Z}#leq300");
  ReplaceAll(code, "(mumu_pt*(mumu_pt>0)+elel_pt*(elel_pt>0))>200&&(mumu_pt*(mumu_pt>0)+elel_pt*(elel_pt>0))<=300",
                   "200<p_{T}^{Z}#leq300");
  ReplaceAll(code, "(mumu_pt*(mumu_pt>0)+elel_pt*(elel_pt>0))","p_{T}^{Z}");
  ReplaceAll(code, "(mumu_m*(mumu_m>0)+elel_m*(elel_m>0))>80&&(mumu_m*(mumu_m>0)+elel_m*(elel_m>0))<100",
             "80<m_{ll}<100");
  ReplaceAll(code, "mumuv_m>80&&mumuv_m<100",
             "80<m_{ll}<100");
  ReplaceAll(code, "elelv_m>80&&elelv_m<100",
             "80<m_{ll}<100");
  ReplaceAll(code, "onht>350&&onmet>100&&","");
  ReplaceAll(code, "jets_islep[0]==0","");
  ReplaceAll(code, "(nels==0&&nmus==1)","N_{#mu}=1");
  ReplaceAll(code, "(nels==1&&nmus==0)","N_{#font[12]{e}}=1");
  ReplaceAll(code, "Max$(abs(els_eta)*(els_sigid&&els_miniso<0.1&&els_pt>20))<1.479","barrel #font[12]{e}");
  ReplaceAll(code, "Max$(abs(els_eta)*(els_sigid&&els_miniso<0.1&&els_pt>20))>1.479","endcap #font[12]{e}");

  ReplaceAll(code, "!low_dphi", "high #Delta#phi");
  ReplaceAll(code, "hig_drmax", "#DeltaR^{max}_{bb}");
  ReplaceAll(code, "ntks", "N_{tks}");
  ReplaceAll(code, "nleps", "N_{lep}");
  ReplaceAll(code, "nvleps", "N_{lep}");
  ReplaceAll(code, "nmus", "N_{#mu}");
  ReplaceAll(code, "nels", "N_{e}");
  ReplaceAll(code, "nvmus", "N^{veto}_{#mu}");
  ReplaceAll(code, "nvels", "N^{veto}_{e}");
  ReplaceAll(code, "ntruleps", "N^{true}_{lep}");
  ReplaceAll(code, "_ra2b", "^{ra2b}");
  ReplaceAll(code, "npv", "N_{PV}");
  ReplaceAll(code, "mumu_pt1", "p_{T}^{#mu}");
  ReplaceAll(code, "elel_pt1", "p_{T}^{e}");

  ReplaceAll(code, "abs(mc_id)==1000006", "stop");
  ReplaceAll(code, "abs(mc_id)==1000022", "LSP");

  ReplaceAll(code, "onmet", "MET^{on}");
  ReplaceAll(code, "onht", "H_{T}^{on}");
  ReplaceAll(code, "njets30","N_{jets}^{30}");
  ReplaceAll(code, "els_pt","p^{e}_{T}");
  ReplaceAll(code, "mus_pt","p^{#mu}_{T}");
  ReplaceAll(code, "fjets_nconst","N_{const}^{fat jet}");
  ReplaceAll(code, "fjets_30_m[0]","m(J_{1})");
  ReplaceAll(code, "fjets_m[0]","m(J_{1})");
  ReplaceAll(code, "(fjets_pt*cosh(fjets_eta))","p_{fatjet}");
  ReplaceAll(code, "fjets_pt","p^{fatjet}_{T}");
  ReplaceAll(code, "jets_pt","p^{jet}_{T}");
  ReplaceAll(code, "mus_reliso","RelIso");
  ReplaceAll(code, "els_reliso","RelIso");
  ReplaceAll(code, "mus_miniso_tr15","MiniIso");
  ReplaceAll(code, "els_miniso_tr15","MiniIso");
  ReplaceAll(code, "njets","N_{jets}");
  ReplaceAll(code, "abs(lep_id)==13&&","");
  ReplaceAll(code, ">=", " #geq ");
  ReplaceAll(code, ">", " > ");
  ReplaceAll(code, "<=", " #leq ");
  ReplaceAll(code, "<", " < ");
  ReplaceAll(code, "&&", ", ");
  ReplaceAll(code, "==", " = ");
  ReplaceAll(code, "met", "E_{T}^{miss}");
  ReplaceAll(code, "ht_hlt", "H_{T}^{HLT}");
  ReplaceAll(code, "mht", "MHT");
  ReplaceAll(code, "ht", "H_{T}");
  ReplaceAll(code, "mt", "m_{T}");
  ReplaceAll(code, "ntks_chg==0", " ITV");
  ReplaceAll(code, "nbm","N_{b}");
  ReplaceAll(code, "nbl","N_{b,l}");
  ReplaceAll(code, "mj", " M_{J}");

  ReplaceAll(code, "el_tks_mt", "Track m_{T}");
  ReplaceAll(code, "mu_tks_mt", "Track m_{T}");
  ReplaceAll(code, "had_tks_mt", "Track m_{T}");
  ReplaceAll(code, "el_tks_pt", "Track p_{T}");
  ReplaceAll(code, "mu_tks_pt", "Track p_{T}");
  ReplaceAll(code, "had_tks_pt", "Track p_{T}");
  ReplaceAll(code, "el_tks_miniso", "Track miniso");
  ReplaceAll(code, "mu_tks_miniso", "Track miniso");
  ReplaceAll(code, "had_tks_miniso", "Track miniso");
  ReplaceAll(code, "el_tks_chg_miniso", "Track charged miniso");
  ReplaceAll(code, "mu_tks_chg_miniso", "Track charged miniso");
  ReplaceAll(code, "had_tks_chg_miniso", "Track charged miniso");
  ReplaceAll(code, "jetsys_nob_pt", "ISR p_{T}");
  ReplaceAll(code, "(", "");
  ReplaceAll(code, ")", "");

  return code;
}

string CodeToLatex(string code){
  code = CodeToRootTex(code);
  ReplaceAll(code, "#", "\\");
  ReplaceAll(code, "\\DeltaR", "\\Delta R");
  return code;
}

vector<string> Tokenize(const string& input,
                        const string& tokens){
  char* ipt(new char[input.size()+1]);
  memcpy(ipt, input.data(), input.size());
  ipt[input.size()]=static_cast<char>(0);
  char* ptr(strtok(ipt, tokens.c_str()));
  vector<string> output(0);
  while(ptr!=NULL){
    output.push_back(ptr);
    ptr=strtok(NULL, tokens.c_str());
  }
  return output;
}

string MakeDir(string prefix){
  prefix += "XXXXXX";
  char *dir_name = new char[prefix.size()];
  if(dir_name == nullptr) ERROR("Could not allocate directory name");
  strcpy(dir_name, prefix.c_str());
  mkdtemp(dir_name);
  prefix = dir_name;
  delete[] dir_name;
  return prefix;
}

string MakeTemp(string prefix){
  prefix += "XXXXXX";
  char *file_name = new char[prefix.size()];
  if(file_name == nullptr) ERROR("Could not allocate file name");
  strcpy(file_name, prefix.c_str());
  mkstemp(file_name);
  prefix = file_name;
  delete[] file_name;
  return prefix;
}

void AdjustDensityForBinWidth(TH1D &h){
  double entries = h.GetEntries();
  int nbins = h.GetNbinsX();
  double low = h.GetBinLowEdge(1);
  double high = h.GetBinLowEdge(nbins+1);
  double width = (high-low)/nbins;
  for(int bin = 1; bin <= nbins; ++bin){
    double content = h.GetBinContent(bin);
    double error = h.GetBinError(bin);
    double this_width = h.GetBinWidth(bin);
    double scale = width/this_width;
    h.SetBinContent(bin, content*scale);
    h.SetBinError(bin, error*scale);
  }
  h.SetEntries(entries);
}

void Normalize(TH1D &h, double normalization, bool norm_per_avg_width){
  int nbins = h.GetNbinsX();
  double low = h.GetBinLowEdge(1);
  double high = h.GetBinLowEdge(nbins+1);
  double width = (high-low)/nbins;
  if(norm_per_avg_width) normalization *= width;
  double integral = h.Integral("width");
  h.Scale(normalization/integral);
}

void MergeOverflow(TH1D &h, bool merge_underflow, bool merge_overflow){
  if(merge_underflow){
    h.SetBinContent(1, h.GetBinContent(0)+h.GetBinContent(1));
    h.SetBinContent(0, 0.);
    h.SetBinError(1, hypot(h.GetBinError(0), h.GetBinError(1)));
    h.SetBinError(0, 0.);
  }
  int nbins = h.GetNbinsX();
  if(merge_overflow){
    h.SetBinContent(nbins, h.GetBinContent(nbins)+h.GetBinContent(nbins+1));
    h.SetBinContent(nbins+1, 0.);
    h.SetBinError(nbins, hypot(h.GetBinError(nbins), h.GetBinError(nbins+1)));
    h.SetBinError(nbins+1, 0.);
  }
}

string FixedDigits(double x, int n_digits){
  int digits_left = max(floor(log10(x))+1., 0.);
  int digits_right = max(n_digits-digits_left, 0);

  double multiplier = pow(10., digits_right);

  ostringstream oss;
  oss << setprecision(numeric_limits<double>::digits10) << round(x*multiplier)/multiplier << flush;
  string out = oss.str();
  if(out.substr(0,2) == "0."){
    out = out.substr(1);
  }
  return out;
}

string FullTitle(const TH1 &h){
  return string(h.GetTitle())
    +";"+h.GetXaxis()->GetTitle()
    +";"+h.GetYaxis()->GetTitle();
}

TString HoursMinSec(float fseconds){
  long seconds = round(fseconds);
  int minutes((seconds/60)%60), hours(seconds/3600);
  TString hhmmss("");
  if(hours<10) hhmmss += "0";
  hhmmss += hours; hhmmss += ":";
  if(minutes<10) hhmmss += "0";
  hhmmss += minutes; hhmmss += ":";
  if((seconds%60)<10) hhmmss += "0";
  hhmmss += seconds%60; 

  return hhmmss;
}

TString AddCommas(double num){
  TString result(""); result += num;
  int posdot(result.First('.'));
  if(posdot==-1) posdot = result.Length();
  for(int ind(posdot-3); ind > 0; ind -= 3)
    result.Insert(ind, ",");
  return result;
}


TString RoundNumber(double num, int decimals, double denom){
  if(denom==0 || !isfinite(num) || !isfinite(denom)) return " - ";
  double neg = 1; if(num*denom<0) neg = -1;
  num /= neg*denom; num += 0.5*pow(10.,-decimals);
  if(abs(num) > 1e16) return "-";
  long num_int = static_cast<long>(num);
  long num_dec = static_cast<long>((1+num-num_int)*pow(10.,decimals));
  TString s_dec = ""; s_dec += num_dec; s_dec.Remove(0,1);
  TString result="";
  if(neg<0) result+="-";
  result+= num_int;
  if(decimals>0) {
    result+="."; result+=s_dec;
  }

  TString afterdot = result;
  afterdot.Remove(0,afterdot.First(".")+1);
  for(int i=0; i<decimals-afterdot.Length(); i++)
    result += "0";
  if(result.Length()>15) cout<<"num "<<num<<", denom "<<denom<<"  ---->  "<<result<<endl;
  return result;
}

// Code from http://www.hongliangjie.com/2012/12/19/how-to-generate-gamma-random-variables/
// Parameter b could be theta...
double gsl_ran_gamma(const double a, const double b, TRandom3 &rand){
  if (a < 1){
    double u = rand.Uniform(1);
    return gsl_ran_gamma(1.0 + a, b, rand) * pow (u, 1.0 / a);
  }

  double x, v, u;
  double d = a - 1.0 / 3.0;
  double c = (1.0 / 3.0) / sqrt (d);
  
  while (1) {
    do {
      x = rand.Gaus(0, 1.0);
      v = 1.0 + c * x;
    }
    while (v <= 0);
      
    v = v * v * v;
    u = rand.Uniform(1);

    if (u < 1 - 0.0331 * x * x * x * x) 
      break;

    if (log (u) < 0.5 * x * x + d * (1 - v + log (v)))
      break;
  }
    
  return b * d * v;
}

double intGaus(double mean, double sigma, double minX, double maxX){
  return (TMath::Erf((maxX-mean)/sigma/sqrt(2.))-TMath::Erf((minX-mean)/sigma/sqrt(2.)))/2.;
}

float deltaR(float eta1, float phi1, float eta2, float phi2){
  return hypot(TVector2::Phi_mpi_pi(phi2-phi1), eta2-eta1);
}

// Finds significance of observation Nbkg, for an expected background Nbkg+Eup_bkg-Edown_bkg
// The mean of the Poisson is Nbkg convolved with the asymmetric lognormal(Eup_bkg, Edown_bkg)
double Significance(double Nobs, double Nbkg, double Eup_bkg, double Edown_bkg){
  double precision = 0.03; // Desired precision on the p-value
  double Nmin = 1/pow(precision,2), Nmax = 5e7; // Nmax controls max sigmas achievable (5e7->5.5 sigma)
  double Nbelow=0, Nabove=0, Nequal=0;
  if(Edown_bkg<0) Edown_bkg = Eup_bkg; // If down uncertainty not specified, symmetric lognormal
  if(Edown_bkg>Nbkg) {
    cout<<"Down uncertainty ("<<Edown_bkg<<") has to be smaller than Nbkg ("<<Nbkg<<")"<<endl;
    return -999.;
  }
  //if(Nobs==Nbkg) return 0.;
  TRandom3 rand(1234);
  double mu, valG;
  while( (min(Nbelow,Nabove)+Nequal)<Nmin && (Nbelow+Nequal+Nabove)<Nmax){
    // Convolving expected bkg with log-normal
    mu = Nbkg;
    valG = rand.Gaus(0,1);
    if(mu==0) mu = fabs(valG)*Eup_bkg; // Apply 2-sided Gaussian uncertainty
    else if(valG>=0) mu *= exp(valG*log(1+Eup_bkg/Nbkg));
    else if(Edown_bkg<0.8*Nbkg) mu *= exp(-valG*log(1-Edown_bkg/Nbkg));
    else mu = max(0., mu + valG*Edown_bkg);
    // Finding if toy above the observed yield
    double valPois = rand.PoissonD(mu);
    if(valPois>Nobs) Nabove++;
    else if(valPois==Nobs) Nequal++;
    else Nbelow++;
  }

  if(Nabove+Nequal==0){
    cout<<"No toys above or at Nobs="<<Nobs<<" for Nbkg "<<Nbkg<<"+"<<Eup_bkg<<"-"<<Edown_bkg<<". Returning "
	<<RooStats::PValueToSignificance(1/Nmax)<<endl;
    return RooStats::PValueToSignificance(1/Nmax);
  }
  if(Nbelow+Nequal==0){
    cout<<"No toys below or at Nobs="<<Nobs<<" for Nbkg "<<Nbkg<<"+"<<Eup_bkg<<"-"<<Edown_bkg<<". Returning "
	<<RooStats::PValueToSignificance(1-1/Nmax)<<endl;
    return RooStats::PValueToSignificance(1-1/Nmax);
  }
  return RooStats::PValueToSignificance((Nabove+Nequal/2.)/(Nbelow+Nequal+Nabove));
}

// yields[Nobs][Nsam] has the entries for each sample for each observable going into kappa
// weights[Nobs][Nsam] has the average weight of each observable for each sample
// powers[Nobs] defines kappa = Product_obs{ Sum_sam{yields[sam][obs]*weights[sam][obs]}^powers[obs] }
double calcKappa(vector<vector<float> > &entries, vector<vector<float> > &weights,
                 vector<float> &powers, float &mSigma, float &pSigma, bool do_data,
                 bool verbose, double syst, bool do_plot, int nrep){
  TRandom3 rand(1234);
  int nbadk(0);
  vector<float> fKappas;
  double mean(0.), bignum(1e10);
  // Doing kappa variations
  for(int rep(0), irep(0); rep < nrep; rep++) {
    fKappas.push_back(1.);
    bool Denom_is0(false);
    for(unsigned obs(0); obs < powers.size(); obs++) {
      float observed(0.);
      for(unsigned sam(0); sam < entries[obs].size(); sam++) {
        // With a flat prior, the expected average of the Poisson with N observed is Gamma(N+1,1)
        // Rounding the expected yield for data stats
        if(do_data) observed += entries[obs][sam]*weights[obs][sam];
        else observed += gsl_ran_gamma(entries[obs][sam]+1,1,rand)*weights[obs][sam];
      } // Loop over samples
      //if(do_data) observed = gsl_ran_gamma(static_cast<int>(0.5+observed)+1,1,rand);
      if(do_data) observed = gsl_ran_gamma(observed+1,1,rand);
      if(observed <= 0 && powers[obs] < 0) Denom_is0 = true;
      else fKappas[irep] *= pow(observed, powers[obs]);
    } // Loop over number of observables going into kappa

    if(syst>=0){
      double factor = exp(rand.Gaus(0,log(1+syst)));
      fKappas[irep] *= factor;
    }
    if(Denom_is0 && fKappas[irep]==0) {
      fKappas.pop_back();
      nbadk++;
    }else {
      if(Denom_is0) fKappas[irep] = bignum;
      else mean += fKappas[irep];
      irep++;
    }
  } // Loop over fluctuations of kappa (repetitions)
  int ntot(nrep-nbadk);
  mean /= static_cast<double>(ntot);

  sort(fKappas.begin(), fKappas.end());
  double gSigma = intGaus(0,1,0,1);
  int iMedian((nrep-nbadk+1)/2-1);
  int imSigma(iMedian-static_cast<int>(gSigma*ntot)), ipSigma(iMedian+static_cast<int>(gSigma*ntot));
  float median(fKappas[iMedian]);
  mSigma = median-fKappas[imSigma]; pSigma = fKappas[ipSigma]-median;

  // Finding standard value
  float stdval(1.);
  bool infStd(false);
  for(unsigned obs(0); obs < powers.size(); obs++) {
    float stdyield(0.);
    if(verbose) cout<<obs<<": ";
    for(unsigned sam(0); sam < entries[obs].size(); sam++) {
      if(verbose) cout<<"Yield"<<sam<<" "<<entries[obs][sam]*weights[obs][sam]
                      <<", N"<<sam<<" "<<entries[obs][sam]
                      <<", avW"<<sam<<" "<<weights[obs][sam]<<". ";
      stdyield += entries[obs][sam]*weights[obs][sam];
    }
    if(verbose) cout<<"  ==> Total yield "<<stdyield<<endl;
    if(stdyield <= 0 && powers[obs] < 0) infStd = true;
    else stdval *= pow(stdyield, powers[obs]);
  } // Loop over number of observables going into kappa
  if(infStd) stdval = median;
  else {
    int istd(0);
    for(int rep(0); rep < ntot; rep++) 
      if(fKappas[rep]>stdval) {istd = rep; break;}
    imSigma = istd-static_cast<int>(gSigma*ntot);
    ipSigma = istd+static_cast<int>(gSigma*ntot);
    if(imSigma<0){ // Adjusting the length of the interval in case imSigma has less than 1sigma
      ipSigma += (-imSigma);
      imSigma = 0;
    }
    if(ipSigma>=ntot){ // Adjusting the length of the interval in case ipSigma has less than 1sigma
      imSigma -= (ipSigma-ntot+1);
      ipSigma = ntot-1;
    }
    mSigma = stdval-fKappas[imSigma]; pSigma = fKappas[ipSigma]-stdval;
  }

  TCanvas can;
  int nbins(100);
  double minH(stdval-3*mSigma), maxH(stdval+3*pSigma);
  if(minH < fKappas[0]) minH = fKappas[0];
  if(maxH > fKappas[ntot-1]) maxH = fKappas[ntot-1];
  TH1D histo("h","",nbins, minH, maxH);
  for(int rep(0); rep < ntot; rep++) 
    histo.Fill(fKappas[rep]);   
  //histo.SetBinContent(1, histo.GetBinContent(1)+nbadk);
  //histo.SetBinContent(nbins, histo.GetBinContent(nbins)+histo.GetBinContent(nbins+1));
  histo.Scale(1/histo.Integral());
  histo.SetMaximum(histo.GetMaximum()*1.2);
  histo.SetLineWidth(3);
  histo.Draw();
  histo.SetXTitle("Expected value");
  histo.SetYTitle("Probability");
  histo.Draw();
  if(do_plot) can.SaveAs("test.eps");

  double mode(histo.GetBinLowEdge(histo.GetMaximumBin()));
  if(verbose) cout<<"Std kappa = "<<stdval<<"+"<<pSigma<<"-"<<mSigma<<".   Mean = "<<mean
                  <<". Mode = "<<mode<<". Median = "<<median<<endl;

  return stdval;
}

set<string> attach_folder(string folder, set<string> &fileset) {
  set<string> fset = set<string>();
  for (auto &ifile: fileset) fset.insert(folder+ifile);
  return fset; 
}
