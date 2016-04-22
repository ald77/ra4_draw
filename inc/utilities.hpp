#ifndef H_UTILITIES
#define H_UTILITIES

#include <string>
#include <vector>
#include <iostream>
#include <iomanip>
#include <sstream>
#include <limits>
#include <mutex>

#include "TH1D.h"

#define ERROR(x) do{throw std::runtime_error(string("Error in file ")+__FILE__+" at line "+to_string(__LINE__)+" (in "+__func__+"): "+x);}while(false)
#define DBG(x) do{std::cerr << "In " << __FILE__ << " at line " << __LINE__ << " (in function " << __func__ << "): " << x << std::endl;}while(false)

namespace Multithreading{
  extern std::mutex root_mutex;
}

bool Contains(const std::string &str, const std::string &pat);
bool StartsWith(const std::string &str, const std::string &pat);
void ReplaceAll(std::string &str, const std::string &orig, const std::string &rep);

std::string execute(const std::string &cmd);

std::vector<std::string> Tokenize(const std::string& input,
                                  const std::string& tokens=" ");

std::string MakeDir(std::string prefix);

void Scale(TH1D &h, bool adjust_width = false, double normalization = -1.);
void MergeOverflow(TH1D &h, bool merge_underflow, bool merge_overflow);

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

#endif
