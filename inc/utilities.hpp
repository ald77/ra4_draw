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

#include "TH1D.h"
#include "TRandom3.h"

#define ERROR(x) do{throw std::runtime_error(string("Error in file ")+__FILE__+" at line "+to_string(__LINE__)+" (in "+__func__+"): "+x);}while(false)
#define DBG(x) do{std::cerr << "In " << __FILE__ << " at line " << __LINE__ << " (in function " << __func__ << "): " << x << std::endl;}while(false)

namespace Multithreading{
  extern std::mutex root_mutex;
}

std::set<std::string> Glob(const std::string &pattern);
std::string Basename(const std::string &filename);

bool Contains(const std::string &str, const std::string &pat);
bool StartsWith(const std::string &str, const std::string &pat);
void ReplaceAll(std::string &str, const std::string &orig, const std::string &rep);

std::string execute(const std::string &cmd);

std::vector<std::string> Tokenize(const std::string& input,
                                  const std::string& tokens=" ");

std::string MakeDir(std::string prefix);

void AdjustDensityForBinWidth(TH1D &h);
void Normalize(TH1D &h, double normalization, bool norm_per_avg_x);

void MergeOverflow(TH1D &h, bool merge_underflow, bool merge_overflow);

std::string FixedDigits(double x, int n_digits);

std::string FullTitle(const TH1 &h);

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

TString RoundNumber(double num, int decimals, double denom=1.);
TString cuts2tex(TString cuts);
TString cutsToLabel(TString cut);

double gsl_ran_gamma (const double a, const double b, TRandom3 &rand);
double intGaus(double mean, double sigma, double minX, double maxX);
float deltaR(float eta1, float phi1, float eta2, float phi2);
// yields[Nobs][Nsam] has the entries for each sample for each observable going into kappa
// weights[Nobs][Nsam] has the average weight of each observable for each sample
// powers[Nobs] defines kappa = Product_obs{ Sum_sam{yields[sam][obs]*weights[sam][obs]}^powers[obs] }
double calcKappa(std::vector<std::vector<float> > &entries, std::vector<std::vector<float> > &weights,
		 std::vector<float> &powers, float &mSigma, float &pSigma, bool do_data=false, 
		 bool verbose=false, double syst=-1., bool do_plot=false, int nrep=100000);


#endif
