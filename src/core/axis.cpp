#include "core/axis.hpp"

#include <algorithm>

using namespace std;

Axis::Axis(const vector<double> &bins,
           const NamedFunc &var,
           const string &title,
           const set<double> &cut_vals):
  var_(var),
  title_(title),
  units_(""),
  cut_vals_(cut_vals),
  bins_(bins){
  sort(bins_.begin(), bins_.end());
  ParseUnits();
  }

Axis::Axis(size_t nbins,
           double xmin, double xmax,
           const NamedFunc &var,
           const string &title,
           const set<double> &cut_vals):
  var_(var),
  title_(title),
  units_(""),
  cut_vals_(cut_vals),
  bins_(GetEdges(nbins, xmin, xmax)){
  ParseUnits();
  }

size_t Axis::Nbins() const{
  return bins_.size()-1;
}


Axis & Axis::Bins(const vector<double> &bins){
  bins_ = bins;
  sort(bins_.begin(), bins_.end());
  return *this;
}

Axis & Axis::Bins(size_t nbins, double xmin, double xmax){
  bins_ = GetEdges(nbins, xmin, xmax);
  return *this;
}

const vector<double> & Axis::Bins() const{
  return bins_;
}

double Axis::AvgBinWidth() const{
  if(bins_.size() < 2) return 0.;
  else return (bins_.back()-bins_.front())/(bins_.size()-1);
}

string Axis::Title() const{
  if(units_ == ""){
    return title_;
  }else{
    return title_+" ["+units_+"]";
  }
}

vector<double> Axis::GetEdges(size_t nbins, double xmin, double xmax){
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

void Axis::ParseUnits(){
  auto p1 = title_.rfind('[');
  auto p2 = title_.rfind(']');
  if(p1 >= p2 || p2 == string::npos) return;
  units_ = title_.substr(p1+1, p2-p1-1);
  title_ = title_.substr(0, p1);
  while(title_.back() == ' '){
    title_.pop_back();
  }
}
