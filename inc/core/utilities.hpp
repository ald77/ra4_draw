#ifndef H_UTILITIES
#define H_UTILITIES

#include <string>
#include <vector>
#include <iostream>
#include <iomanip>
#include <sstream>
#include <limits>
#include <mutex>
#include <set>
#include <numeric>
#include <algorithm>

#include "TH1D.h"
#include "TRandom3.h"
#include "TLegend.h"
#include "TH1.h"
#include "TH2.h"
#include "TH3.h"

#define ERROR(x) do{throw std::runtime_error(std::string("Error in file ")+__FILE__+" at line "+std::to_string(__LINE__)+" (in "+__func__+"): "+x);}while(false)
#define DBG(x) do{std::cerr << "In " << __FILE__ << " at line " << __LINE__ << " (in function " << __func__ << "): " << x << std::endl;}while(false)

namespace Multithreading{
  extern std::mutex root_mutex;
}

std::set<std::string> Glob(const std::string &pattern);
std::string Basename(const std::string &filename);

bool Contains(const std::string &str, const std::string &pat);
bool StartsWith(const std::string &str, const std::string &pat);
void ReplaceAll(std::string &str, const std::string &orig, const std::string &rep);
std::string CopyReplaceAll(const std::string str, const std::string &orig, const std::string &rep);
std::string LeftStrip(std::string str);
std::string RightStrip(std::string str);
std::string Strip(std::string str);
std::string ChangeExtension(std::string path, const std::string &new_ext);

bool FileExists(const std::string &path);

std::string execute(const std::string &cmd);

std::string CodeToPlainText(std::string code);
std::string CodeToLatex(std::string code);
std::string CodeToRootTex(std::string code);

void getLegendBoxes(TLegend &leg, std::vector<std::vector<float> > &boxes);

template<typename T>
std::vector<std::size_t> SortPermutation(const std::vector<T>& vec){
  std::vector<std::size_t> p(vec.size());
  std::iota(p.begin(), p.end(), 0);
  std::sort(p.begin(), p.end(),
	    [&](std::size_t i, std::size_t j){ return vec[i] < vec[j]; });
  return p;
}
template<typename T>
std::vector<T> ApplyPermutation(const std::vector<T>& vec, const std::vector<std::size_t>& perm){
  if(vec.size() != perm.size()) ERROR("Bad permutation: vector and perm have to have the same size");
  std::vector<T> vsorted(vec.size());
  for(std::size_t ind=0; ind<vec.size(); ind++) vsorted[ind] = vec[perm[ind]];
  return vsorted;
}

std::vector<std::string> Tokenize(const std::string& input,
                                  const std::string& tokens=" ");

std::string MakeDir(std::string prefix);
std::string MakeTemp(std::string prefix);

void AdjustDensityForBinWidth(TH1D &h);
void Normalize(TH1D &h, double normalization, bool norm_per_avg_x);

void MergeOverflow(TH1D &h, bool merge_underflow, bool merge_overflow);

std::string FixedDigits(double x, int n_digits);

std::string FullTitle(const TH1 &h);

std::vector<double> ScaleBins(const TAxis &a, double scale);

TH1D ScaleAxes(const TH1 &h, double scale, const std::string &axes = "xyz");
TH2D ScaleAxes(const TH2 &h, double scale, const std::string &axes = "xyz");
TH3D ScaleAxes(const TH3 &h, double scale, const std::string &axes = "xyz");

void CopyStyle(const TH1 &hin, TH1 &hout);

template<typename T>
void Append(T &collection, const typename T::value_type &value){
  collection.insert(collection.end(), value);
}

template<typename T>
std::string ToString(const T& x){
  std::ostringstream oss;
  oss << std::setprecision(std::numeric_limits<double>::digits10) << x << std::flush;
  return oss.str();
}

template<typename T>
std::string ToLongString(const T& x){
  std::ostringstream oss;
  oss << std::setprecision(std::numeric_limits<double>::max_digits10) << x << std::flush;
  return oss.str();
}

TString HoursMinSec(float fseconds);
TString AddCommas(double num);
TString RoundNumber(double num, int decimals, double denom=1.);

double Significance(double Nobs, double Nbkg, double Eup_bkg, double Edown_bkg=-1.);
double gsl_ran_gamma (const double a, const double b, TRandom3 &rand);
double intGaus(double mean, double sigma, double minX, double maxX);
float deltaR(float eta1, float phi1, float eta2, float phi2);
double deltaPhi(double phi1, double phi2);
// yields[Nobs][Nsam] has the entries for each sample for each observable going into kappa
// weights[Nobs][Nsam] has the average weight of each observable for each sample
// powers[Nobs] defines kappa = Product_obs{ Sum_sam{yields[sam][obs]*weights[sam][obs]}^powers[obs] }
double calcKappa(std::vector<std::vector<float> > &entries, std::vector<std::vector<float> > &weights,
		 std::vector<float> &powers, float &mSigma, float &pSigma, bool do_data=false, 
		 bool verbose=false, double syst=-1., bool do_plot=false, int nrep=100000);

std::set<std::string> attach_folder(std::string folder, std::set<std::string> &fileset);

#endif
