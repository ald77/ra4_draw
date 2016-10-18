// Class defining ABCD planes

#ifndef H_ABCD_METHOD
#define H_ABCD_METHOD

#include <iostream>
#include <vector>

#include "TString.h"


class abcd_method {
public:
  bool is2l, int_nbnj;
  TString method; 
  std::vector<TString> planecuts, abcdcuts, allcuts;
  std::vector<std::vector<TString> > bincuts;
  TString caption, basecuts, title, rd_letter; 

  size_t indexBin(size_t iplane, size_t ibin, size_t iabcd);
  void setLeptons();
  void setMj12();
  void printCuts();
  TString lowerNjets(TString &cut);
  void serializeCuts();
  void setIntNbNj(bool int_nbnj_b);


  abcd_method(TString imethod, std::vector<TString> iplanecuts, std::vector<TString> ibincuts, 
	      std::vector<TString> iabcdcuts, TString icaption="", TString ibasecuts="", TString ititle="");
  abcd_method(TString imethod, std::vector<TString> iplanecuts, std::vector<std::vector<TString> > ibincuts, 
	      std::vector<TString> iabcdcuts, TString icaption="", TString ibasecuts="", TString ititle="");
  ~abcd_method();

};

#endif
