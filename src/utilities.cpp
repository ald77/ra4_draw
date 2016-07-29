#include "utilities.hpp"

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

#include "TMath.h"
#include "TCanvas.h"
#include "TVector2.h"
#include "TH1D.h"

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
  if(dir_name == nullptr) throw runtime_error("Could not allocate directory name");
  strcpy(dir_name, prefix.c_str());
  mkdtemp(dir_name);
  prefix = dir_name;
  delete[] dir_name;
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

//// Converting ROOT cuts to ROOT labels
TString cutsToLabel(TString cut){
  if(cut.Contains("mm_")){
    cut.ReplaceAll("mm_","");
    int ind;
    do{
      ind = cut.First('[');
      cut.Remove(ind, cut.First(']')-ind+1);
    }while(ind>=0);
  } // Cleaning up mismeasured cuts
  cut.ReplaceAll(" ","");
  cut.ReplaceAll("met>150&&met<=200", "150<met<=200");
  cut.ReplaceAll("met>200&&met<=350", "200<met<=350");
  cut.ReplaceAll("met>350&&met<=500", "350<met<=500");
  cut.ReplaceAll("njets>=5&&njets<=7", "5<=njets<=7");
  cut.ReplaceAll("njets>=6&&njets<=8", "6<=njets<=8");
  cut.ReplaceAll("nbm>=1&&nbm<=2", "1<=nbm<=2");

  cut.ReplaceAll("met","E_{T}^{miss}");
  cut.ReplaceAll("njets","N_{jets}");
  cut.ReplaceAll("nbm","N_{b}");
  cut.ReplaceAll("=="," = ");
  cut.ReplaceAll(">="," #geq ");
  cut.ReplaceAll("<="," #leq ");
  cut.ReplaceAll(">"," > ");
  cut.ReplaceAll("<"," < ");
  cut.ReplaceAll("&&",", ");

  return cut;
}



TString cuts2tex(TString cuts){
  cuts.ReplaceAll(" ", "");
  if(cuts.Contains("met>200")){
    if(cuts.Contains("met<=400")) {
      cuts.ReplaceAll("met<=400", "");
      cuts.ReplaceAll("met>200", "200<met<=400");
    }
    if(cuts.Contains("met>400")) {
      cuts.ReplaceAll("met>400", "");
      cuts.ReplaceAll("met>200", "met>400");
    }
  }
  //cuts.ReplaceAll("1&&", ""); 
  cuts.ReplaceAll("&&1", "");
  cuts.ReplaceAll("trig[0]", "\\text{HT350\\_MET100}");  cuts.ReplaceAll("trig[22]", "\\text{Ele27\\_eta2p1}");   
  cuts.ReplaceAll("trig[4]", "\\text{Mu15\\_VVVL}"); cuts.ReplaceAll("trig[8]", "\\text{Ele15\\_VVVL}"); 
  cuts.ReplaceAll("trig[14]", "\\text{MET170}");  cuts.ReplaceAll("trig[28]", "\\text{MET90}");  
  cuts.ReplaceAll("trig[12]", "\\text{HT800}");  

  cuts.ReplaceAll("&&&&", "&&");  cuts.ReplaceAll("&&", ", ");  
  cuts.ReplaceAll("elelv_m", "m_{ee}"); cuts.ReplaceAll("elel_m", "m_{ee}"); cuts.ReplaceAll("elelv_pt", "p^{ee}_T"); 
  cuts.ReplaceAll("mumuv_m", "m_{\\mu\\mu}"); cuts.ReplaceAll("mumu_m", "m_{\\mu\\mu}"); 
  cuts.ReplaceAll("mumuv_pt", "p^{\\mu\\mu}_T"); 
  cuts.ReplaceAll("ht_ra2", "H_T"); cuts.ReplaceAll("ht_clean", "H_T"); cuts.ReplaceAll("ht", "H_T"); 
  cuts.ReplaceAll("mj14", "M_J^{1.4}"); cuts.ReplaceAll("mj", "M_J"); cuts.ReplaceAll("met", "\\mathrm{MET}");  
  cuts.ReplaceAll("njets_ra2", "n_j");  cuts.ReplaceAll("njets_clean", "n_j");  cuts.ReplaceAll("njets", "N_j");  
  cuts.ReplaceAll("nbm", "N_b");  cuts.ReplaceAll("nleps", "n_{\\ell}"); 
  cuts.ReplaceAll("nvels", "n_e"); cuts.ReplaceAll("nels", "n_e");  
  cuts.ReplaceAll("nvmus", "n_\\mu"); cuts.ReplaceAll("nmus", "n_\\mu");  cuts.ReplaceAll("nveto", "N_{\\rm veto}");  
  cuts.ReplaceAll(">=", "\\geq ");  cuts.ReplaceAll("<=", " \\leq "); cuts.ReplaceAll("==", " = ");
  cuts.ReplaceAll("pass_jets","\\text{JetID}"); cuts.ReplaceAll("pass_ra2, ","");cuts.ReplaceAll("pass, ","");
  cuts.ReplaceAll("pass","\\text{all filters}");

  cuts = "$"+cuts+"$";
  return cuts;
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
