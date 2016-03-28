#include "histo_def.hpp"

using namespace std;

HistoDef::HistoDef(const vector<double> &bins,
                   const NamedFunc &var,
                   const string &x_title,
                   const NamedFunc &cut,
                   const NamedFunc &weight):
  bins_(bins),
  var_(var),
  cut_(cut),
  weight_(weight),
  x_title_(x_title){
  }

HistoDef::HistoDef(size_t nbins,
                   double xmin,
                   double xmax,
                   const NamedFunc &var,
                   const std::string &x_title,
                   const NamedFunc &cut,
                   const NamedFunc &weight):
  bins_(GetEdges(nbins, xmin, xmax)),
  var_(var),
  cut_(cut),
  weight_(weight),
  x_title_(x_title){
  }

size_t HistoDef::GetNbins() const{
  return bins_.size()-1;
}

string HistoDef::GetName() const{
  return var_.PlainName() + "_CUT_" + cut_.PlainName() + "_WGT_" + weight_.PlainName();
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
