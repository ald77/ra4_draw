#include "hig/extract_yields.hpp"

#include <cstdlib>

#include <iostream>
#include <iomanip>
#include <sstream>
#include <fstream>
#include <algorithm>
#include <initializer_list>
#include <stdexcept>

#include <getopt.h>

#include "TFile.h"
#include "TLine.h"
#include "TCanvas.h"
#include "TPad.h"
#include "TLegend.h"
#include "TColor.h"
#include "TH2D.h"
#include "TStyle.h"
#include "TLatex.h"

#include "RooArgList.h"
#include "RooArgSet.h"
#include "RooRealVar.h"
#include "RooAbsData.h"

#include "core/utilities.hpp"

using namespace std;

namespace{
  string file_wspace("/cms2r0/babymaker/fits/2017_02_26/TChiHH/wspace_nor4_TChiHH_mGluino-1000_mLSP-1_xsecNom_nbTTML_lumi35p9_sig0_higgsCombine.root");
  string name_wspace("w");
  bool table_clean(false);
  bool r4_only(true);
}

int main(int argc, char *argv[]){
  GetOptionsExtract(argc, argv);

  string higgs_name(file_wspace);
  TFile higgs_file(higgs_name.c_str(),"read");
  if(!higgs_file.IsOpen()) ERROR("File "+higgs_name+" not produced");
  RooWorkspace *w = static_cast<RooWorkspace*>(higgs_file.Get(name_wspace.c_str()));
  if(w == nullptr) ERROR("Workspace "+name_wspace+" not found");

  string full_fit_name = higgs_name;
  ReplaceAll(full_fit_name, "higgsCombine", "mlfit");
  TFile fit_file(full_fit_name.c_str(),"read");
  if(!fit_file.IsOpen()) ERROR("Could not open "+full_fit_name);
  RooFitResult *fit_b = static_cast<RooFitResult*>(fit_file.Get("fit_b"));
  RooFitResult *fit_s = static_cast<RooFitResult*>(fit_file.Get("fit_s"));
  
  string prefix = (Contains(file_wspace,"nor4")?"nor4_":"");
  if(fit_b != nullptr){
    //PrintDebug(*w, *fit_b, prefix+"table_bkg_debug.tex");
    PrintTable(*w, *fit_b, prefix+"table_bkg_fit.tex");
    //MakeYieldPlot(*w, *fit_b, prefix+"bkg_plot.pdf");
    MakeCovarianceMatrix(*w, *fit_b, prefix+"matrix_covar.pdf");
  }
  if(fit_s != nullptr){
    //PrintDebug(*w, *fit_s, ChangeExtension(file_wspace, "_sig_debug.tex"));
    PrintTable(*w, *fit_s, ChangeExtension(file_wspace, "_sig_table.tex"));
    //MakeYieldPlot(*w, *fit_s, ChangeExtension(file_wspace, "_sig_plot.pdf"));
  }

  higgs_file.Close();
  fit_file.Close();

}

string GetSignalName(const RooWorkspace &w){
  TIter iter(w.allFunctions().createIterator());
  int size = w.allFunctions().getSize();
  RooAbsArg *arg = nullptr;
  int i = 0;
  while((arg = static_cast<RooAbsArg*>(iter())) && i < size){
    ++i;
    if(arg == nullptr) continue;
    string name = arg->GetName();
    if(name.substr(0,9) != "nsig_BLK_") continue;
    TIter iter2(arg->getVariables()->createIterator());
    int size2 = arg->getVariables()->getSize();
    RooAbsArg *arg2 = nullptr;
    int i2 = 0;
    while((arg2 = static_cast<RooAbsArg*>(iter2())) && i2 < size2){
      ++i2;
      if(arg2 == nullptr) continue;
      string name2 = arg2->GetName();
      auto pos2 = name2.find("_PRC_");
      if(pos2 != string::npos){
        iter2.Reset();
        iter.Reset();
        return name2.substr(pos2+5);
      }
    }
    iter2.Reset();
  }
  iter.Reset();
  return "signal";
}

string TexFriendly(const string &s){
  string out;
  for(size_t i = 0; i < s.size(); ++i){
    if(s.at(i) == '_'){
      out += "\\_";
    }else{
      out += s.at(i);
    }
  }
  return out;
}

void PrintDebug(RooWorkspace &w,
                const RooFitResult &f,
                const string &file_name){
  SetVariables(w, f);

  vector<string> var_names = GetVarNames(w);
  vector<string> func_names = GetFuncNames(w);

  ofstream out(file_name);
  out << "\\documentclass{article}\n";
  out << "\\usepackage{amsmath,graphicx,rotating,longtable}\n";
  out << "\\thispagestyle{empty}\n";
  out << "\\begin{document}\n";
  out << "\\begin{longtable}{rr}\n";
  out << "\\hline\\hline\n";
  out << "Variable & Fit Value\\\\\n";
  out << "\\hline\n";
  out << fixed << setprecision(2);
  for(const auto &var: var_names){
    RooRealVar *varo = w.var(var.c_str());
    if(varo == nullptr) continue;
    if(!varo->isConstant()){
      out << TexFriendly(var) << " & $" << varo->getVal() << "\\pm" << GetError(*varo, f) << "$\\\\\n";
    }else{
      out << TexFriendly(var) << " & $" << varo->getVal() << "$\\\\\n";
    }
  }
  for(const auto &func: func_names){
    RooAbsReal *funco = w.function(func.c_str());
    if(funco == nullptr) continue;
    if(!funco->isConstant()){
      out << TexFriendly(func) << " & $" << funco->getVal() << "\\pm" << GetError(*funco, f) << "$\\\\\n";
    }else{
      out << TexFriendly(func) << " & $" << funco->getVal() << "$\\\\\n";
    }
  }

  out << "\\hline\\hline\n";
  out << "\\end{longtable}\n";
  out << "\\end{document}\n";
  out << endl;
  out.close();
  cout<<"Saved "<<file_name.c_str()<<endl;
}

void PrintTable(RooWorkspace &w,
                const RooFitResult &f,
                const string &file_name){
  RooRealVar *r_var = SetVariables(w, f);
  if(r_var != nullptr && !r_var->isConstant()){
    cout<<"Signal strength: "<<r_var->getVal() << " + " << GetError(*r_var, f, 1) 
	<< " - " << GetError(*r_var, f, -1) << endl;
  }

  string sig_name = GetSignalName(w);
  vector<string> prc_names = GetProcessNames(w);
  vector<string> bin_names = GetPlainBinNames(w);

  bool dosig(Contains(file_name, "sig_table"));
  bool blind_all(Contains(file_name, "r4blinded"));
  bool blind_2b(Contains(file_name, "1bunblinded"));
  size_t digits(2), ncols(10);
  if(!dosig) ncols = 8;
  if(table_clean) {
    ncols--;
    digits = 1;
  }

  ofstream out(file_name);
  out << fixed << setprecision(digits);
  out << "\\documentclass{article}\n";
  out << "\\usepackage{amsmath,graphicx,rotating}\n";
  out << "\\usepackage[landscape]{geometry}\n";
  out << "\\thispagestyle{empty}\n";
  out << "\\begin{document}\n";
  out << "\\begin{table}\n";
  out << "\\centering\n";
  out << "\\resizebox{\\textwidth}{!}{\n";
  out << "\\begin{tabular}{l ";
  for(size_t i = 0; i < ncols-1; ++i) out << "r";
  out << "}\n";
  out << "\\hline\\hline\n\n";
  out << "Bin & ";
  for(const auto &prc_name: prc_names){
    out << prc_name << " & ";
  }
  out << "MC Bkg. "<<(dosig?"& Bkgnd. Pred. ":"")<<"& Signal "<<(dosig?"& Sig. Pred. ":"")
      <<"& Tot. Pred. & Obs.";
  if(!table_clean) out << " & $\\lambda$";
  out<<"\\\\\n";

  // out << "\\hline\n";
  for(const auto &bin_name: bin_names){
    if(Contains(bin_name, "r1")) {
      out << "\n\\hline\\hline"<<endl;
      if(Contains(bin_name, "lowmet")) out<<"\\multicolumn{"<<ncols<<"}{c}{$200<\\text{MET}\\leq 400$} \\\\ \\hline"<<endl;
      if(Contains(bin_name, "highmet")) out<<"\\multicolumn{"<<ncols<<"}{c}{$\\text{MET}>400$} \\\\ \\hline"<<endl;
    }
    string bin_tex(TexFriendly(bin_name));
    ReplaceAll(bin_tex, "lowmet\\_","");
    ReplaceAll(bin_tex, "highmet\\_","");
    ReplaceAll(bin_tex, "lownj\\_","$n_j\\leq8$, ");
    ReplaceAll(bin_tex, "highnj\\_","$n_j\\geq9$, ");
    ReplaceAll(bin_tex, "allnb","all $n_j,n_b$");
    ReplaceAll(bin_tex, "1b","$n_b=1$");
    ReplaceAll(bin_tex, "3b","$n_b\\geq3$");
    if(Contains(bin_name, "lowmet")) ReplaceAll(bin_tex, "2b","$n_b=2$");
    else ReplaceAll(bin_tex, "2b","$n_b\\geq2$");
    for(int ind(1); ind<=4; ind++){
      ReplaceAll(bin_tex, "r"+to_string(ind)+"\\_","R"+to_string(ind)+": ");
      ReplaceAll(bin_tex, "r"+to_string(ind)+"c\\_","R"+to_string(ind)+": ");
      ReplaceAll(bin_tex, "d"+to_string(ind)+"\\_","D"+to_string(ind)+": ");
    }
    out << bin_tex << " & ";
    for(const auto &prc_name: prc_names){
      out << GetMCYield(w, bin_name, prc_name) << " & ";
    }
    out << "$" << GetMCTotal(w, bin_name);
    if(!table_clean) out << "\\pm" << GetMCTotalErr(w, f, bin_name);
    out <<  "$ & ";

    if(dosig) out << "$" << GetBkgPred(w, bin_name) << "\\pm" << GetBkgPredErr(w, f, bin_name) <<  "$ & ";
    out << GetMCYield(w, bin_name, sig_name) << " & ";
    if(dosig) out << "$" << GetSigPred(w, bin_name) << "\\pm" << GetSigPredErr(w, f, bin_name) <<  "$ & ";
    if(!Contains(file_wspace, "nor4") || (Contains(bin_name, "hig_3b")||Contains(bin_name, "hig_4b")))
      out << "$" << GetTotPred(w, bin_name) << "^{+" << GetTotPredErr(w, f, bin_name,1) 
	  <<"}_{-"<< GetTotPredErr(w, f, bin_name,-1) <<  "}$";
    out<<" & ";
    if(Contains(bin_name,"4") && (blind_all || (!Contains(bin_name,"1b") && blind_2b))) out << "-- & ";
    else out << setprecision(0) << GetObserved(w, bin_name);
    out << setprecision(digits);
    if(!table_clean) out << "& $" << GetLambda(w, bin_name) << "\\pm" << GetLambdaErr(w, f, bin_name) <<  "$";
    out << "\\\\\n";
    if(Contains(bin_name, "r3") || Contains(bin_name, "d3")) out << "\\hline"<<endl;
  }
  out << "\n\\hline\\hline\n";
  out << "\\end{tabular}\n";
  out << "}\n";
  out << "\\end{table}\n";
  out << "\\end{document}\n";
  out << endl;
  out.close();
  cout<<"Saved "<<file_name.c_str()<<endl;
}

double GetMCYield(const RooWorkspace &w,
                  const string &bin_name,
                  const string &prc_name){
  TIter iter(w.allFunctions().createIterator());
  int size = w.allFunctions().getSize();
  RooAbsArg *arg = nullptr;
  int i = 0;
  while((arg = static_cast<RooAbsArg*>(iter())) && i < size){
    ++i;
    if(arg == nullptr) continue;
    string name = arg->GetName();
    if(name.substr(0,8) != "ymc_BLK_") continue;
    if(!(Contains(name, "_BIN_"+bin_name))) continue;
    if(!(Contains(name, "_PRC_"+prc_name))) continue;
    return static_cast<RooRealVar*>(arg)->getVal();
  }
  iter.Reset();
  return -1.;
}

double GetMCTotal(const RooWorkspace &w,
                  const string &bin_name){
  TIter iter(w.allFunctions().createIterator());
  int size = w.allFunctions().getSize();
  RooAbsArg *arg = nullptr;
  int i = 0;
  while((arg = static_cast<RooAbsArg*>(iter())) && i < size){
    ++i;
    if(arg == nullptr) continue;
    string name = arg->GetName();
    if(name.substr(0,8) != "ymc_BLK_") continue;
    if(!(Contains(name, "_BIN_"+bin_name))) continue;
    if(Contains(name, "_PRC_")) continue;
    return static_cast<RooRealVar*>(arg)->getVal();
  }
  iter.Reset();
  return -1.;
}

double GetMCTotalErr(RooWorkspace &w,
                     const RooFitResult &f,
                     const string &bin_name){
  TIter iter(w.allFunctions().createIterator());
  int size = w.allFunctions().getSize();
  RooAbsArg *arg = nullptr;
  int i = 0;
  while((arg = static_cast<RooAbsArg*>(iter())) && i < size){
    ++i;
    if(arg == nullptr) continue;
    string name = arg->GetName();
    if(name.substr(0,8) != "ymc_BLK_") continue;
    if(!(Contains(name, "_BIN_"+bin_name))) continue;
    if(Contains(name, "_PRC_")) continue;
    return GetError(*static_cast<RooAbsReal*>(arg), f);
  }
  iter.Reset();
  return -1.;
}

double GetBkgPred(const RooWorkspace &w,
                  const string &bin_name){
  TIter iter(w.allFunctions().createIterator());
  int size = w.allFunctions().getSize();
  RooAbsArg *arg = nullptr;
  int i = 0;
  while((arg = static_cast<RooAbsArg*>(iter())) && i < size){
    ++i;
    if(arg == nullptr) continue;
    string name = arg->GetName();
    if(name.substr(0,9) != "nbkg_BLK_") continue;
    if(!(Contains(name, "_BIN_"+bin_name))) continue;
    if(Contains(name, "_PRC_")) continue;
    return static_cast<RooRealVar*>(arg)->getVal();
  }
  iter.Reset();
  return -1.;
}

double GetBkgPredErr(RooWorkspace &w,
                     const RooFitResult &f,
                     const string &bin_name){
  TIter iter(w.allFunctions().createIterator());
  int size = w.allFunctions().getSize();
  RooAbsArg *arg = nullptr;
  int i = 0;
  while((arg = static_cast<RooAbsArg*>(iter())) && i < size){
    ++i;
    if(arg == nullptr) continue;
    string name = arg->GetName();
    if(name.substr(0,9) != "nbkg_BLK_") continue;
    if(!(Contains(name, "_BIN_"+bin_name))) continue;
    if(Contains(name, "_PRC_")) continue;
    return GetError(*static_cast<RooAbsReal*>(arg), f);
  }
  iter.Reset();
  return -1.;
}

double GetSigPred(const RooWorkspace &w,
                  const string &bin_name){
  TIter iter(w.allFunctions().createIterator());
  int size = w.allFunctions().getSize();
  RooAbsArg *arg = nullptr;
  int i = 0;
  while((arg = static_cast<RooAbsArg*>(iter())) && i < size){
    ++i;
    if(arg == nullptr) continue;
    string name = arg->GetName();
    if(name.substr(0,9) != "nsig_BLK_") continue;
    if(!(Contains(name, "_BIN_"+bin_name))) continue;
    if(Contains(name, "_PRC_")) continue;
    return static_cast<RooRealVar*>(arg)->getVal();
  }
  iter.Reset();
  return -1.;
}

double GetSigPredErr(RooWorkspace &w,
                     const RooFitResult &f,
                     const string &bin_name){
  TIter iter(w.allFunctions().createIterator());
  int size = w.allFunctions().getSize();
  RooAbsArg *arg = nullptr;
  int i = 0;
  while((arg = static_cast<RooAbsArg*>(iter())) && i < size){
    ++i;
    if(arg == nullptr) continue;
    string name = arg->GetName();
    if(name.substr(0,9) != "nsig_BLK_") continue;
    if(!(Contains(name, "_BIN_"+bin_name))) continue;
    if(Contains(name, "_PRC_")) continue;
    return GetError(*static_cast<RooAbsReal*>(arg), f);
  }
  iter.Reset();
  return -1.;
}

double GetTotPred(const RooWorkspace &w,
                  const string &bin_name){
  TIter iter(w.allFunctions().createIterator());
  int size = w.allFunctions().getSize();
  RooAbsArg *arg = nullptr;
  int i = 0;
  while((arg = static_cast<RooAbsArg*>(iter())) && i < size){
    ++i;
    if(arg == nullptr) continue;
    string name = arg->GetName();
    if(name.substr(0,9) != "nexp_BLK_") continue;
    if(!(Contains(name, "_BIN_"+bin_name))) continue;
    if(Contains(name, "_PRC_")) continue;
    return static_cast<RooRealVar*>(arg)->getVal();
  }
  iter.Reset();
  return -1.;
}

double GetTotPredErr(RooWorkspace &w,
                     const RooFitResult &f,
                     const string &bin_name, int errtype){
  TIter iter(w.allFunctions().createIterator());
  int size = w.allFunctions().getSize();
  RooAbsArg *arg = nullptr;
  int i = 0;
  while((arg = static_cast<RooAbsArg*>(iter())) && i < size){
    ++i;
    if(arg == nullptr) continue;
    string name = arg->GetName();
    if(name.substr(0,9) != "nexp_BLK_") continue;
    if(!(Contains(name, "_BIN_"+bin_name))) continue;
    if(Contains(name, "_PRC_")) continue;
    //cout<<name<<endl;
    return GetError(*static_cast<RooAbsReal*>(arg), f, errtype);
  }
  iter.Reset();
  return -1.;
}

double GetObserved(const RooWorkspace &w,
                   const string &bin_name){
  ostringstream oss;
  oss << "data_obs";
  oss << flush;
  RooAbsData *data = w.data(oss.str().c_str());
  if(data == nullptr) ERROR("Could not find dataset "+oss.str());
  const RooArgSet *args = data->get();
  if(args == nullptr) ERROR("Could not extract args");
  TIter iter(args->createIterator());
  int size = args->getSize();
  RooAbsArg *arg = nullptr;
  int i = 0;
  while((arg = static_cast<RooAbsArg*>(iter())) && i < size){
    ++i;
    if(arg == nullptr) continue;
    string name = arg->GetName();
    if(name.substr(0,9) != "nobs_BLK_") continue;
    if(!(Contains(name, "_BIN_"+bin_name))) continue;
    if(Contains(name, "_PRC_")) continue;
    return static_cast<RooRealVar*>(arg)->getVal();
  }
  iter.Reset();
  return -1.;
}

double GetLambda(const RooWorkspace &w,
                 const string &bin_name){
  TIter iter(w.allFunctions().createIterator());
  int size = w.allFunctions().getSize();
  RooAbsArg *arg = nullptr;
  int i = 0;
  while((arg = static_cast<RooAbsArg*>(iter())) && i < size){
    ++i;
    if(arg == nullptr) continue;
    string name = arg->GetName();
    if(name.substr(0,12) != "kappamc_BLK_") continue;
    if(!(Contains(name, "_BIN_"+bin_name))) continue;
    if(Contains(name, "_PRC_")) continue;
    return static_cast<RooRealVar*>(arg)->getVal();
  }
  iter.Reset();
  return -1.;
}

double GetLambdaErr(RooWorkspace &w,
                    const RooFitResult &f,
                    const string &bin_name){
  TIter iter(w.allFunctions().createIterator());
  int size = w.allFunctions().getSize();
  RooAbsArg *arg = nullptr;
  int i = 0;
  while((arg = static_cast<RooAbsArg*>(iter())) && i < size){
    ++i;
    if(arg == nullptr) continue;
    string name = arg->GetName();
    if(name.substr(0,12) != "kappamc_BLK_") continue;
    if(!(Contains(name, "_BIN_"+bin_name))) continue;
    if(Contains(name, "_PRC_")) continue;
    return GetError(*static_cast<RooAbsReal*>(arg), f);
  }
  iter.Reset();
  return -1.;
}

RooRealVar * SetVariables(RooWorkspace &w,
                          const RooFitResult &f){
  bool set_r = false;
  RooArgList pars = f.floatParsFinal();
  for(int ipar = 0; ipar < pars.getSize(); ++ipar){
    RooRealVar *fit_var = static_cast<RooRealVar*>(pars.at(ipar));
    if(fit_var == nullptr) continue;
    RooRealVar *w_var = w.var(fit_var->GetName());
    if(w_var == nullptr) continue;
    w_var->removeRange();
    w_var->setVal(fit_var->getVal());
    w_var->setError(fit_var->getError());
    if(fit_var->GetName() == string("r")) set_r = true;
  }
  vector<string> var_names = GetVarNames(w);
  vector<string> func_names = GetFuncNames(w);
  for(const auto &var: var_names){
    RooRealVar *varo = w.var(var.c_str());
    if(varo == nullptr) continue;
    if(!varo->isConstant()){
      varo->removeRange();
    }
  }

  RooRealVar *r_var = static_cast<RooRealVar*>(w.var("r"));
  if(r_var != nullptr){
    if(!set_r){
      r_var->setVal(0);
      r_var->setConstant(true);
    }else{
      r_var->setConstant(false);
    }
  }
  return r_var;
}

void MakeYieldPlot(RooWorkspace &w,
                   const RooFitResult &f,
                   const string &file_name){
  RooRealVar *r_var = SetVariables(w, f);

  vector<string> bin_names = GetBinNames(w);
  vector<string> prc_names = GetProcessNames(w);

  vector<vector<double> > component_yields = GetComponentYields(w, bin_names, prc_names);

  vector<TH1D> histos = MakeBackgroundHistos(component_yields, bin_names, prc_names);
  TH1D signal = MakeTotalHisto(w, f, bin_names);
  TGraphErrors band = MakeErrorBand(signal);
  TH1D obs = MakeObserved(w, bin_names);

  SetBounds(obs, signal, histos);

  TCanvas c;
  c.cd();
  TPad bot_pad("bot_pad", "bot_pad", 0., 0., 1., 0.4);
  bot_pad.SetFillColor(0); bot_pad.SetFillStyle(4000);
  bot_pad.SetMargin(0.1, 0., 0.5, 0.);
  bot_pad.Draw();
  c.cd();
  TPad mid_pad("mid_pad", "mid_pad", 0., 0.4, 1., 0.85);
  mid_pad.SetFillColor(0); mid_pad.SetFillStyle(4000);
  mid_pad.SetMargin(0.1, 0., 0.0, 0.);
  mid_pad.Draw();
  c.cd();
  TPad top_pad("top_pad", "top_pad", 0., 0.85, 1., 1.0);
  top_pad.SetFillColor(0); top_pad.SetFillStyle(4000);
  top_pad.SetMargin(0.1, 0., 0.0, 0.);
  top_pad.Draw();

  double font_size = 0.1;
  double offset = 0.5;

  mid_pad.cd();
  if(!Contains(file_wspace, "nor4")) mid_pad.SetLogy();
  signal.SetTitleSize(font_size, "Y");
  signal.SetTitleOffset(offset, "Y");
  signal.SetFillColor(2);
  signal.SetFillStyle(1001);
  signal.SetLineColor(2);
  signal.SetLineStyle(1);
  signal.SetLineWidth(0);
  signal.Draw("hist");
  for(auto h = histos.rbegin(); h!= histos.rend(); ++h){
    h->Draw("same");
  }

  double marker_size(1.4);
  obs.SetMarkerStyle(20); obs.SetMarkerSize(marker_size);
  band.Draw("02 same");
  obs.Draw("ex0 same");
  signal.Draw("same axis");

  top_pad.cd();
  TLegend l(0.1, 0., 1., 1.);
  l.SetNColumns(3);
  l.SetFillColor(0); l.SetFillStyle(4000);
  l.SetBorderSize(0);
  l.AddEntry(&obs, "Observed", "lep");
  l.AddEntry(&signal, "Signal", "f");
  ostringstream oss;
  oss << setprecision(2) << fixed;
  oss << "r=";
  if(r_var == nullptr){
    oss << "???";
  }else if(r_var->isConstant()){
    oss << r_var->getVal() << " (fixed)";
  }else{
    oss << r_var->getVal() << "#pm" << GetError(*r_var, f);
    cout<<"Signal strength: "<<r_var->getVal() << "#pm" << GetError(*r_var, f) << endl;
  }
  oss << flush;
  l.AddEntry(&obs, oss.str().c_str(), "");
  for(auto h = histos.crbegin(); h != histos.crend(); ++h){
    l.AddEntry(&(*h), h->GetName(), "f");
  }
  l.Draw("same");

  bot_pad.cd();
  TLine line; line.SetLineStyle(2);
  TGraphErrors obs_rat = MakeRatio(obs, signal);
  TGraphErrors pred_rat = MakeRatio(signal, signal);
  TH1D dumb = obs;
  obs_rat.SetMarkerStyle(20); obs_rat.SetMarkerSize(marker_size);
  obs_rat.SetMarkerColor(1);
  dumb.SetLineColor(0);
  dumb.SetLineWidth(0);
  dumb.SetFillColor(0);
  dumb.SetFillStyle(4000);
  dumb.SetMinimum(0.);
  dumb.SetMaximum(2.8);
  dumb.SetTitle(";;Obs/Pred ");
  dumb.GetXaxis()->LabelsOption("V");
  dumb.SetTitleSize(font_size, "Y");
  dumb.SetTitleOffset(offset, "Y");
  dumb.Draw();
  pred_rat.SetFillColor(kGray);
  pred_rat.SetFillStyle(3001);
  pred_rat.Draw("02 same");
  obs_rat.Draw("ep0 same");
  line.DrawLine(0.5, 1, 0.5+dumb.GetNbinsX(), 1);
  c.Print(file_name.c_str());
}

vector<string> GetVarNames(const RooWorkspace &w){
  vector<string> names;
  TIter iter(w.allVars().createIterator());
  int size = w.allVars().getSize();
  TObject *obj = nullptr;
  int i = 0;
  while((obj = iter()) && i < size){
    ++i;
    if(obj == nullptr) continue;
    string name = obj->GetName();
    Append(names, name);
  }
  iter.Reset();
  sort(names.begin(), names.end());
  return names;
}

vector<string> GetFuncNames(const RooWorkspace &w){
  vector<string> names;
  TIter iter(w.allFunctions().createIterator());
  int size = w.allFunctions().getSize();
  TObject *obj = nullptr;
  int i = 0;
  while((obj = iter()) && i < size){
    ++i;
    if(obj == nullptr) continue;
    string name = obj->GetName();
    Append(names, name);
  }
  iter.Reset();
  sort(names.begin(), names.end());
  return names;
}

vector<string> GetBinNames(const RooWorkspace &w){
  vector<string> names;
  TIter iter(w.allFunctions().createIterator());
  int size = w.allFunctions().getSize();
  RooAbsArg *arg = nullptr;
  int i = 0;
  while((arg = static_cast<RooAbsArg*>(iter())) && i < size){
    ++i;
    if(arg == nullptr) continue;
    string name = arg->GetName();
    if(name.substr(0,9) != "nexp_BLK_") continue;
    if((!Contains(name, "sig_3b") || !Contains(name, "sig_4b")) && Contains(file_wspace, "nor4")) continue;
    string bin_name = name.substr(5);
    Append(names, bin_name);
  }
  iter.Reset();
  reverse(names.begin(), names.end());
  return names;
}

vector<string> GetPlainBinNames(const RooWorkspace &w){
  vector<string> names;
  TIter iter(w.allFunctions().createIterator());
  int size = w.allFunctions().getSize();
  RooAbsArg *arg = nullptr;
  int i = 0;
  while((arg = static_cast<RooAbsArg*>(iter())) && i < size){
    ++i;
    if(arg == nullptr) continue;
    string name = arg->GetName();
    if(name.substr(0,9) != "nexp_BLK_") continue;
    auto bpos = name.find("_BIN_");
    auto ppos = name.find("_PRC_");
    if(bpos == string::npos) continue;
    string bin_name = name.substr(bpos+5, ppos-bpos-5);
    Append(names, bin_name);
  }
  iter.Reset();
  reverse(names.begin(), names.end());
  return names;
}

vector<string> GetProcessNames(const RooWorkspace &w){
  vector<string> names;
  TIter iter(w.allFunctions().createIterator());
  int size = w.allFunctions().getSize();
  RooAbsArg *arg = nullptr;
  int i = 0;
  while((arg = static_cast<RooAbsArg*>(iter())) && i < size){
    ++i;
    if(arg == nullptr) continue;
    string name = arg->GetName();
    if(name.substr(0,9) != "frac_BIN_") continue;
    auto prc_pos = name.find("_PRC_");
    if(prc_pos == string::npos) continue;
    string bin_name = name.substr(prc_pos+5);
    if(find(names.cbegin(), names.cend(), bin_name) != names.cend()) continue;
    Append(names, bin_name);
  }
  iter.Reset();
  return names;
}

vector<vector<double> > GetComponentYields(const RooWorkspace &w,
                                           const vector<string> &bin_names,
                                           const vector<string> &prc_names){
  vector<vector<double> > yields(bin_names.size());
  for(auto &bin: yields){
    bin = vector<double>(prc_names.size(), 0.);
  }
  for(size_t ibin = 0; ibin < yields.size(); ++ibin){
    const string &bin_name = bin_names.at(ibin);
    auto blk_pos = bin_name.find("_BIN_");
    if(blk_pos == string::npos) continue;
    string plain_name = bin_name.substr(blk_pos+5);
    for(size_t iprc = 0; iprc < yields.at(ibin).size(); ++iprc){
      const string &prc_name = prc_names.at(iprc);
      RooRealVar *nbkg_arg = static_cast<RooRealVar*>(w.function(("nbkg_"+bin_name).c_str()));
      if(nbkg_arg == nullptr) continue;
      RooRealVar *frac_arg = static_cast<RooRealVar*>(w.function(("frac_BIN_"+plain_name+"_PRC_"+prc_name).c_str()));
      if(frac_arg == nullptr) continue;
      yields.at(ibin).at(iprc) = nbkg_arg->getVal() * frac_arg->getVal();
    }
  }
  return yields;
}

vector<TH1D> MakeBackgroundHistos(const vector<vector<double> > &yields,
                                  const vector<string> &bin_names,
                                  const vector<string> &prc_names){
  if(yields.size() == 0){
    return vector<TH1D>();
  }
  vector<TH1D> histos(yields.at(0).size(),
                      TH1D("", ";;Yield ", yields.size(), 0.5, yields.size()+0.5));
  for(size_t ibin = 0; ibin < yields.size(); ++ibin){
    for(size_t iprc = 0; iprc < yields.at(ibin).size(); ++iprc){
      histos.at(iprc).SetBinContent(ibin+1, yields.at(ibin).at(iprc));
    }
  }

  for(size_t iprc = 0; iprc < histos.size(); ++iprc){
    TH1D &h = histos.at(iprc);
    h.SetName(prc_names.at(iprc).c_str());
    h.SetFillColor(iprc+3);
    h.SetLineColor(iprc+3);
    h.SetLineWidth(0);
    for(size_t ibin = 0; ibin < bin_names.size(); ++ibin){
      const string &name = bin_names.at(ibin);
      auto pos = name.find("_BIN_");
      if(pos == string::npos) continue;
      pos = name.find("4");
      if(pos != string::npos && Contains(file_wspace, "nor4")) continue;
      string label = name.substr(pos+5);
      h.GetXaxis()->SetBinLabel(ibin+1, label.c_str());
    }
  }
  sort(histos.begin(), histos.end(),
       [](const TH1D &a, const TH1D &b) -> bool{return a.Integral() < b.Integral();});

  for(size_t iprc = histos.size()-1; iprc < histos.size(); --iprc){
    TH1D &h = histos.at(iprc);
    for(size_t isum = iprc-1; isum < histos.size(); --isum){
      h.Add(&histos.at(isum));
    }
  }

  return histos;
}

TH1D MakeTotalHisto(RooWorkspace &w,
                    const RooFitResult &f,
                    const vector<string> &bin_names){
  TH1D h("signal", ";;Yield ", bin_names.size(), 0.5, bin_names.size()+0.5);
  h.SetFillColor(2);
  h.SetLineColor(2);
  h.SetLineWidth(0);

  for(size_t ibin = 0; ibin < bin_names.size(); ++ibin){
    const string &name = bin_names.at(ibin);
    auto pos = name.find("_BIN_");
    if(pos == string::npos) continue;
    string label = name.substr(pos+5);
    h.GetXaxis()->SetBinLabel(ibin+1, label.c_str());
    RooRealVar *var = static_cast<RooRealVar*>(w.function(("nexp_"+name).c_str()));
    if(var == nullptr) continue;
    h.SetBinContent(ibin+1, var->getVal());
    h.SetBinError(ibin+1, GetError(*var, f));
  }

  return h;
}

TH1D MakeObserved(const RooWorkspace &w,
                  const vector<string> &bin_names){
  TH1D h("observed", ";;Yield ", bin_names.size(), 0.5, bin_names.size()+0.5);
  h.SetLineColor(1);
  h.SetFillColor(0);
  h.SetFillStyle(4000);

  for(size_t ibin = 0; ibin < bin_names.size(); ++ibin){
    const string &name = bin_names.at(ibin);
    auto pos = name.find("_BIN_");
    if(pos == string::npos) continue;
    string label = name.substr(pos+5);
    h.GetXaxis()->SetBinLabel(ibin+1, label.c_str());
    RooRealVar *var = static_cast<RooRealVar*>(w.var(("nobs_"+name).c_str()));
    if(var == nullptr) continue;
    h.SetBinContent(ibin+1, var->getVal());
  }

  return h;
}

void SetBounds(TH1D &a,
               TH1D &b,
               std::vector<TH1D> &cs){
  double factor = 0.02;

  double hmax = GetMaximum(a, b, cs);
  double hmin = GetMinimum(a, b, cs);
  double lmax = log(hmax);
  double lmin = log(hmin);
  double log_diff = lmax-lmin;
  lmin -= factor*log_diff;
  lmax += factor*log_diff;
  hmin = exp(lmin);
  hmax = exp(lmax);
  if(!Contains(file_wspace, "nor4")){
    a.SetMinimum(hmin);
    a.SetMaximum(hmax);
    b.SetMinimum(hmin);
    b.SetMaximum(hmax);
    for(auto &c: cs){
      c.SetMinimum(hmin);
      c.SetMaximum(hmax);
    }
  } else {
    a.SetMaximum(hmax+1.1*sqrt(hmax));
    b.SetMaximum(hmax+1.1*sqrt(hmax));
    a.SetMinimum(0);
    b.SetMinimum(0);
  }
}

double GetMaximum(const TH1D &a,
                  const TH1D &b,
                  const vector<TH1D> &cs){
  double the_max = GetMaximum(a);
  double this_max = GetMaximum(b);
  if(this_max > the_max) the_max = this_max;
  for(const auto &c: cs){
    this_max = GetMaximum(c);
    if(this_max > the_max) the_max = this_max;
  }
  return the_max;
}

double GetMinimum(const TH1D &a,
                  const TH1D &b,
                  const vector<TH1D> &cs){
  double the_min = GetMinimum(a, 0.1);
  double this_min = GetMinimum(b, 0.1);
  if(this_min < the_min) the_min = this_min;
  for(const auto &c: cs){
    this_min = GetMinimum(c, 0.1);
    if(this_min < the_min) the_min = this_min;
  }
  return the_min;
}

double GetMaximum(const TH1D &h, double y){
  double the_max = -numeric_limits<double>::max();
  for(int bin = 1; bin <= h.GetNbinsX(); ++bin){
    double content = h.GetBinContent(bin);
    if(content > the_max){
      if(content < y){
        the_max = content;
      }else{
        the_max = y;
      }
    }
  }
  return the_max;
}

double GetMinimum(const TH1D &h, double y){
  double the_min = numeric_limits<double>::max();
  for(int bin = 1; bin <= h.GetNbinsX(); ++bin){
    double content = h.GetBinContent(bin);
    if(content < the_min){
      if(content > y){
        the_min = content;
      }else{
        the_min = y;
      }
    }
  }
  return the_min;
}

TGraphErrors MakeErrorBand(const TH1D &h){
  TGraphErrors g(h.GetNbinsX());
  for(int bin = 1; bin <= h.GetNbinsX(); ++bin){
    g.SetPoint(bin, h.GetBinCenter(bin), h.GetBinContent(bin));
    g.SetPointError(bin, 0.5, h.GetBinError(bin));
  }
  g.SetFillColor(kGray);
  g.SetFillStyle(3001);
  return g;
}

TGraphErrors MakeRatio(const TH1D &num, const TH1D &den){
  TGraphErrors g(num.GetNbinsX());
  double xerror(0.5);
  if(&num != &den) xerror = 0;
  for(int bin = 1; bin <= num.GetNbinsX(); ++bin){
    double x = num.GetBinCenter(bin);
    double nc = num.GetBinContent(bin);
    double dc = den.GetBinContent(bin);
    double ne = num.GetBinError(bin);
    double big_num = 0.5*numeric_limits<float>::max();
    if(dc != 0.){
      g.SetPoint(bin, x, nc/dc);
      g.SetPointError(bin, xerror, ne/dc);
    }else if(nc == 0.){
      g.SetPoint(bin, x, 1.);
      g.SetPointError(bin, xerror, big_num);
    }else{
      g.SetPoint(bin, x, nc > 0. ? big_num : -big_num);
      g.SetPointError(bin, xerror, big_num);
    }
  }
  return g;
}

void MakeCorrectionPlot(RooWorkspace &w,
                        const RooFitResult &f,
                        const string &file_name){
  SetVariables(w, f);

  vector<string> bin_names = GetBinNames(w);
  vector<string> prc_names = GetProcessNames(w);

  TCanvas c;
  c.cd();

  TH1D h("", ";;#lambda", bin_names.size(), 0.5, bin_names.size()+0.5);
  for(size_t ibin = 0; ibin < bin_names.size(); ++ibin){
    string bin = bin_names.at(ibin);
    auto pos = bin.find("_BIN_");
    string plain_bin = bin.substr(pos+5);
    h.GetXaxis()->SetBinLabel(ibin+1, plain_bin.c_str());
    h.SetBinContent(ibin+1, static_cast<RooRealVar*>(w.function(("kappamc_"+bin).c_str()))->getVal());
    h.SetBinError(ibin+1, GetError(*w.function(("kappamc_"+bin).c_str()), f));
  }
  h.GetXaxis()->LabelsOption("V");
  h.Draw();
  c.SetMargin(0.1, 0.05, 1./3., 0.05);
  c.Print(file_name.c_str());
}

void MakeCovarianceMatrix(RooWorkspace &w,
			  const RooFitResult &f,
			  string covar_file_name){
  SetVariables(w, f);
  const RooArgList &fpf = f.floatParsFinal();

  // Make list of parameter instances of cloneFunc in order of error matrix
  RooArgList paramList;
  vector<int> fpf_idx;

  vector<RooAbsReal*> yields;
  RooArgSet funcs = w.allFunctions();
  TIter iter(funcs.createIterator());
  int size = funcs.getSize();
  RooAbsReal *arg = nullptr;
  int i = 0;
  while((arg = static_cast<RooAbsReal*>(iter())) && i < size){
    ++i;
    if(arg == nullptr) continue;
    string name = arg->GetName();
    if(name.substr(0,9) != "nbkg_BLK_") continue;
    if(!Contains(name, "_BIN_")) continue;
    if(Contains(name, "_PRC_")) continue;
    if(r4_only && !(Contains(name, "hig_3b")||Contains(name, "hig_4b"))) continue;
    yields.push_back(arg);
  }

  vector<vector<double> > errors(fpf.getSize(), vector<double>(yields.size(), 0.));
  for(Int_t iparam = 0; iparam<fpf.getSize(); ++iparam){
    RooRealVar &rrv2 = static_cast<RooRealVar&>(*w.var(fpf.at(iparam)->GetName()));
    RooRealVar &rrv = static_cast<RooRealVar&>(fpf[iparam]);

    double lambda = 0.01;
    double cenVal = rrv.getVal();
    double minVal = rrv.getMin();
    double maxVal = rrv.getMax();
    double downVal = cenVal-lambda*fabs(rrv.getErrorLo());
    double upVal = cenVal+lambda*fabs(rrv.getErrorHi());

    string parname = rrv.GetName();
    // if(Contains(parname, "rx21_BLK_met3") || Contains(parname, "rx31_BLK_met3")){
    //   cout<<iparam<<" - "<<parname<<": cenVal "<<setw(12)<<cenVal<<", Hi "<<setw(12)<< rrv.getErrorHi()
    // 	  <<", Lo "<<setw(12)<<rrv.getErrorLo()
    // 	  <<", upVal "<<setw(12)<<upVal
    // 	  <<", downVal "<<setw(12)<<downVal<<", minVal "<<setw(12)<<minVal<<", maxVal "<<setw(12)<<maxVal<<endl;
    // }

    if(upVal-downVal >= maxVal-minVal){
      //Error bars bigger than variable range
      downVal = minVal;
      upVal = maxVal;
    }else if(downVal < minVal){
      upVal += minVal - downVal;
      downVal = minVal;
    }else if(upVal > maxVal){
      downVal -= upVal - maxVal;
      upVal = maxVal;
    }

    // /// Adam's original
    // rrv.setVal(upVal);
    // for(size_t iyield = 0; iyield<yields.size(); ++iyield){
    //   errors.at(iparam).at(iyield) = 0.5*yields.at(iyield)->getVal();
    // }
    // rrv.setVal(downVal);
    // for(size_t iyield = 0; iyield<yields.size(); ++iyield){
    //   errors.at(iparam).at(iyield) -= 0.5*yields.at(iyield)->getVal();
    // }
    // rrv.setVal(cenVal);

    /// Up variation
    rrv2.setVal(upVal);
    for(size_t iyield = 0; iyield<yields.size(); ++iyield){
      errors.at(iparam).at(iyield) = yields.at(iyield)->getVal()/lambda;
    }
    rrv2.setVal(cenVal);
    for(size_t iyield = 0; iyield<yields.size(); ++iyield){
      errors.at(iparam).at(iyield) -= yields.at(iyield)->getVal()/lambda;
      // if(Contains(parname, "BLK_met1") || Contains(parname, "BLK_met2") || Contains(parname, "BLK_met3")){
      // 	cout<<iparam<<" - "<<parname<<", Iyield "<<iyield<<": error "<<setw(12)<<errors.at(iparam).at(iyield)
      // 	    <<", val "<<setw(12)<<yields.at(iyield)->getVal()<<endl;
      // }
    }
    rrv2.setVal(cenVal);
  }

  vector<vector<double> > right(fpf.getSize(), vector<double>(yields.size(), 0.));
  for(Int_t iparam = 0; iparam<fpf.getSize(); ++iparam){
    for(size_t iyield = 0; iyield<yields.size(); ++iyield){
      right.at(iparam).at(iyield) = 0.;
      for(Int_t entry = 0; entry<fpf.getSize(); ++entry){
	right.at(iparam).at(iyield) += f.correlation(fpf.at(iparam)->GetName(),fpf.at(entry)->GetName())
	  * errors.at(entry).at(iyield);
      }
    }
  }

  vector<vector<double> > covar(yields.size(), vector<double>(yields.size(), 0.));
  for(size_t irow = 0; irow < yields.size(); ++irow){
    for(size_t icol = 0.; icol < yields.size(); ++icol){
      covar.at(irow).at(icol) = 0.;
      for(Int_t ientry = 0.; ientry < fpf.getSize(); ++ientry){
	// if(((irow==4&&icol==5) || (irow==6&&icol==7)) && 
	//    fabs(errors.at(ientry).at(irow) * right.at(ientry).at(icol))>0.0001)
	//   cout<<irow<<", "<<icol<<" - entry "<<ientry<<": error "<<setw(12)<<errors.at(ientry).at(irow)
	//       <<", right "<<setw(12)<<right.at(ientry).at(icol)
	//       <<" -> prod = "<<setw(12)<<errors.at(ientry).at(irow) * right.at(ientry).at(icol)<<endl;
	covar.at(irow).at(icol) += errors.at(ientry).at(irow) * right.at(ientry).at(icol);
      }
    }
  }

  TH2D h_covar("", "Covariance Matrix",
	       covar.size(), -0.5, covar.size()-0.5,
	       covar.size(), -0.5, covar.size()-0.5);
  TH2D h_corr("", "Correlation Matrix",
	      covar.size(), -0.5, covar.size()-0.5,
	      covar.size(), -0.5, covar.size()-0.5);
  float labelSizeX = 0.05, labelSizeY = 0.05, markerSize = 1.9;
  if(!r4_only) {
    labelSizeX = 0.03;
    labelSizeY = 0.035;
    markerSize = 0.7;
  }
  h_covar.SetLabelSize(labelSizeX, "x");
  h_covar.SetLabelSize(labelSizeY, "y");
  h_covar.SetMarkerSize(markerSize);
  h_covar.SetTickLength(0., "xy");
  for(size_t x = 0; x < yields.size(); ++x){
    string name = yields.at(x)->GetName();
    auto pos = name.find("_BIN_");
    name = name.substr(pos+5);
    name = PrettyBinName(name);
    h_covar.GetXaxis()->SetBinLabel(x+1, name.c_str());
    h_covar.GetYaxis()->SetBinLabel(x+1, name.c_str());
    h_corr.GetXaxis()->SetBinLabel(x+1, name.c_str());
    h_corr.GetYaxis()->SetBinLabel(x+1, name.c_str());
    for(size_t y = 0; y < yields.size(); ++y){
      h_covar.SetBinContent(x+1,y+1,covar.at(x).at(y));
      h_corr.SetBinContent(x+1,y+1,covar.at(x).at(y)/sqrt(covar.at(x).at(x)*covar.at(y).at(y)));
    }
  }
  h_covar.LabelsOption("vd","X");
  h_corr.LabelsOption("vd","X");

  const int bands = 255;
  const unsigned num = 2;
  int colors[bands];
  double stops[num] = {0., 1.};
  double red[num] =   {1., 207/255.};
  double green[num] = {1., 131/255.};
  double blue[num] =  {1., 132/255.};
  int fi = TColor::CreateGradientColorTable(num, stops, red, green, blue, bands);
  for(int ib = 0; ib < bands; ++ib){
    colors[ib] = fi+ib;
  }
  gStyle->SetNumberContours(bands);
  gStyle->SetPalette(bands, colors);

  float LeftMargin = 0.12, RightMargin = 0.15, BottomMargin = 0.15, TopMargin = 0.07;
  TString cmsPrel = "#font[62]{CMS} #scale[0.8]{#font[52]{Supplementary}}  #scale[0.73]{#font[82]{arXiv:xxxx.xxxxx}}";
  TString lumiEner = "#font[42]{35.9 fb^{-1} (13 TeV)}"; 
  TLatex cmslabel;  
  cmslabel.SetNDC(kTRUE);

  //////// Covariance matrix
  gStyle->SetOptStat(0);              // No Stats box
  gStyle->SetTextFont(42);
  TCanvas c("", "", 1024, 700);
  c.SetMargin(LeftMargin, RightMargin, BottomMargin, TopMargin);
  gStyle->SetPaintTextFormat("5.1f");
  h_covar.SetTitle(""); h_corr.SetTitle("");
  if(!r4_only) h_covar.GetXaxis()->LabelsOption("v");
  h_covar.GetZaxis()->SetTitle("Covariance");
  h_covar.GetZaxis()->SetTitleSize(0.045);
  h_covar.GetZaxis()->CenterTitle(true);
  if(r4_only) h_covar.GetXaxis()->SetLabelOffset(0.0085);
  h_covar.Draw("axis");
  h_covar.Draw("colz same");
  h_covar.Draw("text same");
  cmslabel.SetTextAlign(11); cmslabel.SetTextSize(0.045);
  cmslabel.DrawLatex(LeftMargin+0.005, 1-TopMargin+0.015, cmsPrel);
  cmslabel.SetTextAlign(31); cmslabel.SetTextSize(0.041);
  cmslabel.DrawLatex(1-RightMargin-0.005, 1-TopMargin+0.015, lumiEner);
  c.SaveAs(covar_file_name.c_str());

  const unsigned num2 = 3;
  int colors2[bands];
  double stops2[num2] = {0., 0.5, 1.};
  double red2[num2] =   {137/255., 1., 207/255.};
  double green2[num2] = {162/255., 1., 131/255.};
  double blue2[num2] =  {215/255., 1., 132/255.};
  int fi2 = TColor::CreateGradientColorTable(num2, stops2, red2, green2, blue2, bands);
  for(int ib = 0; ib < bands; ++ib){
    colors2[ib] = fi2+ib;
  }
  gStyle->SetNumberContours(bands);
  gStyle->SetPalette(bands, colors2);

  //////// Correlation matrix
  gStyle->SetPaintTextFormat("5.2f");
  c.SetLogz(false);
  h_corr.SetMinimum(-1.);
  h_corr.SetMaximum(1.);
  if(!r4_only) h_corr.GetXaxis()->LabelsOption("v");
  h_corr.SetLabelSize(labelSizeX, "x");
  h_corr.SetLabelSize(labelSizeY, "y");
  h_corr.SetMarkerSize(markerSize);
  h_corr.SetTickLength(0., "xy");
  h_corr.GetZaxis()->SetTitle("Correlation");
  h_corr.GetZaxis()->SetTitleSize(0.045);
  h_corr.GetZaxis()->CenterTitle(true);
  if(r4_only) h_corr.GetXaxis()->SetLabelOffset(0.0085);
  h_corr.Draw("colz");
  h_corr.Draw("text same");
  ReplaceAll(covar_file_name, "_covar.pdf", "_corr.pdf");
  cmslabel.SetTextAlign(11); cmslabel.SetTextSize(0.045);
  cmslabel.DrawLatex(LeftMargin+0.005, 1-TopMargin+0.015, cmsPrel);
  cmslabel.SetTextAlign(31); cmslabel.SetTextSize(0.041);
  cmslabel.DrawLatex(1-RightMargin-0.005, 1-TopMargin+0.015, lumiEner);
  c.SaveAs(covar_file_name.c_str());

  TString pname = "CMS-SUS-16-044_AuxFigure_2_CorrGlobalFit.root";
  TString fitname = "Global";
  if(Contains(covar_file_name, "nor4")) {
    pname = "CMS-SUS-16-044_AuxFigure_1_CorrPredFit.root";
    fitname = "Predictive";
  }
  TFile file(pname, "recreate");
  file.cd();
  h_corr.Write("CorrelationMatrix_"+fitname+"Fit");
  h_covar.Write("CovarianceMatrix_"+fitname+"Fit");
  file.Close();
  cout<<"Saved correlation matrix in "<<pname<<endl<<endl;
}

string PrettyBinName(string name){
  // ReplaceAll(name, "sbd_", "SBD, ");
  // ReplaceAll(name, "hig_", "HIG, ");
  ReplaceAll(name, "sbd_", "SBD, ");
  if(r4_only) ReplaceAll(name, "hig_", "");
  else ReplaceAll(name, "hig_", "HIG, ");
  ReplaceAll(name, "2b_", "2b, ");
  ReplaceAll(name, "3b_", "3b, ");
  ReplaceAll(name, "4b_", "4b, ");

  // ReplaceAll(name, "met0", "150<p_{T}^{miss}#leq 200");
  // ReplaceAll(name, "met1", "200<p_{T}^{miss}#leq 300");
  // ReplaceAll(name, "met2", "300<p_{T}^{miss}#leq 450");
  // ReplaceAll(name, "met3", "p_{T}^{miss}>450");

  ReplaceAll(name, "met0", "MET1");
  ReplaceAll(name, "met1", "MET2");
  ReplaceAll(name, "met2", "MET3");
  ReplaceAll(name, "met3", "MET4");
  return name;
}

double GetError(const RooAbsReal &var,
                const RooFitResult &f, int errtype){
  // Clone self for internal use
  RooAbsReal* cloneFunc = static_cast<RooAbsReal*>(var.cloneTree());
  RooArgSet* errorParams = cloneFunc->getObservables(f.floatParsFinal());
  RooArgSet* nset = cloneFunc->getParameters(*errorParams);

  // Make list of parameter instances of cloneFunc in order of error matrix
  RooArgList paramList;
  const RooArgList& fpf = f.floatParsFinal();
  vector<int> fpf_idx;
  for (int i=0; i<fpf.getSize(); i++) {
    RooAbsArg* par = errorParams->find(fpf[i].GetName());
    if (par) {
      paramList.add(*par);
      fpf_idx.push_back(i);
    }
  }

  string name = var.GetName();
  vector<double> errors(paramList.getSize());
  for (Int_t ivar=0; ivar<paramList.getSize(); ivar++) {
    RooRealVar& rrv = static_cast<RooRealVar&>(fpf[fpf_idx[ivar]]);

    double cenVal = rrv.getVal();
    double errVal = rrv.getError();
    if(errtype == 1) errVal = rrv.getErrorHi();
    if(errtype ==-1) errVal = rrv.getErrorLo();

    // Make Plus variation
    if(errtype==1) static_cast<RooRealVar*>(paramList.at(ivar))->setVal(cenVal+errVal);
    else if(errtype==-1) static_cast<RooRealVar*>(paramList.at(ivar))->setVal(cenVal);
    else static_cast<RooRealVar*>(paramList.at(ivar))->setVal(cenVal+0.5*errVal);
    double up = cloneFunc->getVal(nset);

    // Make Minus variation
    if(errtype==1) static_cast<RooRealVar*>(paramList.at(ivar))->setVal(cenVal);
    else if(errtype==-1) static_cast<RooRealVar*>(paramList.at(ivar))->setVal(cenVal+errVal);
    else static_cast<RooRealVar*>(paramList.at(ivar))->setVal(cenVal-0.5*errVal);
    double down = cloneFunc->getVal(nset);

    // string parname = rrv.GetName();
    // if(Contains(parname, "rx21_BLK_met3") || Contains(parname, "rx31_BLK_met3")){
    //   cout<<name<<" "<<ivar<<" - "<<parname<<": cenVal "<<setw(12)<<cenVal<<", errVal "<<setw(12)<<errVal
    // 	  <<", up "<<setw(12)<<up<<", down "<<setw(12)<<down<<", errtype "<<setw(12)<<errtype<<endl;
    // }

    errors.at(ivar) = (up-down);

    static_cast<RooRealVar*>(paramList.at(ivar))->setVal(cenVal);
  }
  // cout<<endl;

  // bool print_corr = false;
  // if(name.substr(0,9) == "nexp_BLK_" && Contains(name, "_BIN_hig_4b_met") && !Contains(name, "_PRC_")){
  //      cout<<"Doing hig_4b_met3 = "<<var.getVal() <<endl;
  //      print_corr = true;
  // }
  vector<double> right(errors.size());
  for(size_t i = 0; i < right.size(); ++i){
    right.at(i) = 0.;
    for(size_t j = 0; j < errors.size(); ++j){
      right.at(i) += f.correlation(paramList.at(i)->GetName(),paramList.at(j)->GetName())*errors.at(j);
      // string par1 = paramList.at(i)->GetName(), par2 = paramList.at(j)->GetName();
      // if(print_corr && (Contains(par1, "rx") || Contains(par1, "ry") || Contains(par1, "norm"))
      // 	 && (Contains(par2, "rx") || Contains(par2, "ry") || Contains(par2, "norm"))) 
      // 	cout<<setw(40)<<par1<<", "<<setw(40)<<par2<<": corr = "
      // 	    <<f.correlation(paramList.at(i)->GetName(),paramList.at(j)->GetName())<<endl;
    }
  }
  double sum = 0.;
  for(size_t i = 0; i < right.size(); ++i){
    sum += errors.at(i)*right.at(i);
  }

  if(cloneFunc != nullptr){
    delete cloneFunc;
    cloneFunc = nullptr;
  }
  if(errorParams != nullptr){
    delete errorParams;
    errorParams = nullptr;
  }
  if(nset != nullptr){
    delete nset;
    nset = nullptr;
  }

  return sqrt(sum);
}

void GetOptionsExtract(int argc, char *argv[]){
  while(true){
    static struct option long_options[] = {
      {"file_wspace", required_argument, 0, 'f'},
      {"name_wspace", required_argument, 0, 'w'},
      {"r4_only", no_argument, 0, '4'},
      {"table_clean", no_argument, 0, 'c'},
      {0, 0, 0, 0}
    };

    char opt = -1;
    int option_index;
    opt = getopt_long(argc, argv, "f:w:c4", long_options, &option_index);
    if( opt == -1) break;

    string optname;
    switch(opt){
    case 'w':
      name_wspace = optarg;
      break;
    case 'f':
      file_wspace = optarg;
      break;
    case '4':
      r4_only = false;
      break;
    case 'c':
      table_clean = true;
      break;
    default:
      printf("Bad option! getopt_long returned character code 0%o\n", opt);
      break;
    }
  }
}
