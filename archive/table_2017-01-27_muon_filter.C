#include <fstream>
#include <iostream>
#include <string>
#include <vector>
#include <iomanip>  // setw

#include "TChain.h"
#include "TMath.h"
#include "TString.h"

using namespace std;


TString RoundNumber(double num, int decimals, double denom=1.){
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

void table_badmuons(){

  vector<TString> sfolders({"/cms2r0/babymaker/babies/2017_01_27/data/merged_database_standard/*root", "/cms2r0/babymaker/babies/2017_01_27/mc/unprocessed/fullbaby_SMS-T1tttt_mGluino-1200_mLSP-800*.root"});
  vector<TString> names({"Data", "T1tttt(1200,800)"});

  // TString baseline = "st>500&&nbm>=1&&njets>=4&&pass"; // Table we presented
  TString baseline = "st>500&&nbm>=1&&njets>=4&&pass&&nleps==1&&mt>140"; // Table we want
  TString filters = "pass_ra2_badmu&&met/met_calo<5";
  TString mbad = "(Sum$(mus_bad_dupl||mus_bad)>0)";
  //vector<TString> cuts({"nmus==1", "nels==1", filters+"&&nmus==1", filters+"&&nels==1"});
  vector<TString> cuts({"1", filters});
  vector<TString> mets({"met>100&&met<=200", "met>200&&met<=350", "met>350&&met<=500", "met>500"});
  vector<TString> metslbl({"100<met<=200", "200<met<=350", "350<met<=500", "met>500"});

  string file_name = "badmu_table.tex";
  std::ofstream file(file_name);
  file << "\\documentclass[10pt,oneside]{report}\n";
  file << "\\usepackage{graphicx,xspace,amssymb,amsmath,colordvi,colortbl,verbatim,multicol}\n";
  file << "\\usepackage{multirow, rotating}\n\n";
  file << "\\usepackage[active,tightpage]{preview}\n\n";

  file << "\\begin{document}\n";
  file << "\\begin{preview}\n";
  file << "  \\begin{tabular}{ l";
  for (size_t ind=0; ind<sfolders.size(); ind++) file << "rl";
  file << "}\n";
  file << "    \\hline\\hline\n";
  file<<"& \\multicolumn{2}{c}{Base filters} & \\multicolumn{2}{c}{Base + SUSY filters} \\\\ \n";
  file << "    \\hline\n";
  
  cout<<endl<<endl;
  for(size_t ind=0; ind<sfolders.size(); ind++){
    cout<<endl<<"Doing "<<sfolders[ind]<<endl;
    file<<" \\hline \n";
    file<<" \\multicolumn{5}{c}{"<<names[ind]<<"} \\\\ \n";
    file<<" \\hline \n";
    TChain ntu("tree"); ntu.Add(sfolders[ind]);
    TString fullbase = baseline;
    if (sfolders[ind].Contains("data")) fullbase += "&&trig_ra4";
    for(unsigned imet(0); imet<mets.size(); imet++){
      metslbl[imet].ReplaceAll("<met<=","$<E_{T}^{\\rm miss}\\leq$").ReplaceAll("met>","$E_{T}^{\\rm miss}>$");
      file<<metslbl[imet];
      for(auto cut:cuts){
       long all = ntu.GetEntries(fullbase+"&&"+mets[imet]+"&&"+cut);
       long bad = ntu.GetEntries(fullbase+"&&"+mets[imet]+"&&"+cut+"&&"+mbad);

       file<<" & "<<bad<<"/"<<all<<" &= "<<RoundNumber(bad*100,1,all)<<"\\%";
       cout<<setw(3)<<bad<<"/"<<setw(5)<<all<<" = "<<RoundNumber(bad*100,1,all)<<"%\t ";
      } // Loop over cuts
      file<<"\\\\ \n";
      cout<<endl;
    } // Loop over mets
  } // Loop over ntuples
  cout<<endl<<endl;
  file << "    \\hline\\hline\n";
  file << "  \\end{tabular}\n";
  file << "\\end{preview}\n";
  file << "\\end{document}\n";
  file << flush;
  file.close();
  cout << "pdflatex " << file_name << " &> /dev/null" << endl;
  exit(0);
}
