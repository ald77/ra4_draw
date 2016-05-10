#include "utilities.hpp"

#include <cstdlib>
#include <cstring>
#include <cstdint>

#include <sstream>
#include <string>
#include <vector>
#include <memory>

#include <unistd.h>

using namespace std;

mutex Multithreading::root_mutex;

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
