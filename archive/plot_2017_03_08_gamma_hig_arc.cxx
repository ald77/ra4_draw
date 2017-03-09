// Cheatsheet:
//    To obtain plots with no mismeasurement and kappa's with expected data stats:
//         ./run/hig/plot_systematics.exe --mm mc_as_data -l 36.2
//    There are 4 possibilities for the skim, requested with option -s. These are: search, zll, qcd, ttbar
//    Option -t plots the kappas with a tighter selection, see basecuts below, e.g.
//         ./run/hig/plot_systematics.exe --mm mc_as_data -t -s zll -l 36.2

#include <fstream>
#include <iostream>

#include "core/utilities.hpp"

using namespace std;

int main(){

  int nrep = 1000000;
  vector<vector<float> > entries, weights;
  vector<float> powers;
  float val, mSigma, pSigma;

  entries = vector<vector<float> >({{0}});
  powers = vector<float>({1});
  weights = vector<vector<float> >(entries.size(), vector<float>(entries[0].size(), 1));  
  val = calcKappa(entries, weights, powers, mSigma, pSigma, true, true, 0, true, nrep);

  entries = vector<vector<float> >({{5}});
  powers = vector<float>({1});
  weights = vector<vector<float> >(entries.size(), vector<float>(entries[0].size(), 1));  
  val = calcKappa(entries, weights, powers, mSigma, pSigma, true, true, 0, true, nrep);

  entries = vector<vector<float> >({{0},{5},{5}});
  powers = vector<float>({1,1,-1});
  weights = vector<vector<float> >(entries.size(), vector<float>(entries[0].size(), 1));  
  val = calcKappa(entries, weights, powers, mSigma, pSigma, true, true, 0, true, nrep);

  entries = vector<vector<float> >({{1},{5},{5}});
  powers = vector<float>({1,1,-1});
  weights = vector<vector<float> >(entries.size(), vector<float>(entries[0].size(), 1));  
  val = calcKappa(entries, weights, powers, mSigma, pSigma, true, true, 0, true, nrep);

  entries = vector<vector<float> >({{2},{39},{74}});
  powers = vector<float>({1,1,-1});
  weights = vector<vector<float> >(entries.size(), vector<float>(entries[0].size(), 1));  
  val = calcKappa(entries, weights, powers, mSigma, pSigma, true, true, 0, true, nrep);

}
