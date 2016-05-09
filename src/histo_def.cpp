#include "histo_def.hpp"

#include <algorithm>

using namespace std;

HistoDef::HistoDef(const vector<double> &bins,
                   const NamedFunc &var,
                   const string &x_title,
                   const NamedFunc &cut,
                   const NamedFunc &weight,
                   const set<double> &cut_vals):
  var_(var),
  cut_(cut),
  weight_(weight),
  x_title_(x_title),
  units_(""),
  cut_vals_(cut_vals),
  bins_(bins){
  sort(bins_.begin(), bins_.end());
  ParseUnits();
  }

HistoDef::HistoDef(size_t nbins,
                   double xmin,
                   double xmax,
                   const NamedFunc &var,
                   const string &x_title,
                   const NamedFunc &cut,
                   const NamedFunc &weight,
                   const set<double> &cut_vals):
  var_(var),
  cut_(cut),
  weight_(weight),
  x_title_(x_title),
  units_(""),
  cut_vals_(cut_vals),
  bins_(GetEdges(nbins, xmin, xmax)){
  ParseUnits();
}

size_t HistoDef::Nbins() const{
  return bins_.size()-1;
}

HistoDef & HistoDef::Bins(const std::vector<double> &bins){
  bins_ = bins;
  sort(bins_.begin(), bins_.end());
  return *this;
}

const vector<double> & HistoDef::Bins() const{
  return bins_;
}

string HistoDef::Name() const{
  return var_.PlainName() + "_CUT_" + cut_.PlainName() + "_WGT_" + weight_.PlainName();
}

string HistoDef::Title() const{
  bool cut = (cut_.Name() != "" && cut_.Name() != "1");
  bool weight = weight_.Name() != "weight";
  if(cut && weight){
    return cut_.Name()+" (weight="+weight_.Name()+")";
  }else if(cut){
    return cut_.Name();
  }else if(weight){
    return "weight="+weight_.Name();
  }else{
    return "";
  }
}

vector<double> HistoDef::GetEdges(size_t nbins, double xmin, double xmax){
  vector<double> edges(nbins+1);
  if(nbins != 0){
    double delta = (xmax-xmin)/nbins;
    for(size_t i = 0; i < nbins+1; ++i){
      edges.at(i) = i*delta + xmin;
    }
  }

  //Not necessary, but make sure that first and last edge are correct to available precision
  edges.front() = xmin;
  edges.back() = xmax;

  return edges;
}

void HistoDef::ParseUnits(){
  auto p1 = x_title_.rfind('[');
  auto p2 = x_title_.rfind(']');
  if(p1 >= p2 || p2 == string::npos) return;
  units_ = x_title_.substr(p1+1, p2-p1-1);
  x_title_ = x_title_.substr(0, p1);
  while(x_title_.back() == ' '){
    x_title_.pop_back();
  }
}
