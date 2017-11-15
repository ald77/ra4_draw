#include <fstream>
#include <sstream>
#include <string>
#include <algorithm>

#include <unistd.h>
#include <getopt.h>

#include "TCanvas.h"
#include "TFile.h"
#include "TGraph.h"
#include "TGraphErrors.h"
#include "TH1D.h"
#include "TStyle.h"
#include "TLegend.h"
#include "TLatex.h"
#include "TLine.h"
#include "TCollection.h"
#include "TSystemDirectory.h"
#include "TSystemFile.h"
#include "TSystem.h"
#include "TRegexp.h"

#include "core/styles.hpp"
#include "core/utilities.hpp"
#include "core/plot_opt.hpp"

using namespace std;

namespace{
  TString filename = "txt/tdr/cfeb_dcfeb_loss.txt";
  PlotOpt opts("txt/plot_styles.txt", "OneCol1D");
}

void GetOptions(int argc, char *argv[]);

void multVector(vector<float> &vals, float factor){
  for(unsigned ind=0; ind<vals.size(); ind++) vals[ind] *= factor;
}

float toHours(TString time, TString ampm){
  time.ReplaceAll(".", " ");
  istringstream iss(time.Data());
  float hours, min, sec;
  TString decsec;
  iss >> hours >> min >> sec >> decsec;
  //cout<<hours<<":"<<min<<":"<<sec<<":"<<decsec<<" from "<<time<<endl;

  if(ampm=="AM") {
    if(hours>11.5) hours = 0;
  } else
    if (hours<11.5) hours += 12;
  
  decsec = "0."+decsec;
  return hours + min/60. + (sec+decsec.Atof())/(60*60);
}

// Returns list of directorites or files in folder
vector<TString> dirlist(const TString &folder,
                        const TString &inname,
                        const TString &tag){
  TRegexp regex_tag(tag,true), regex_inname(inname, true);
  TString pwd(gSystem->pwd());
  vector<TString> v_dirs;
  TSystemDirectory dir(folder, folder);
  TList *files = dir.GetListOfFiles();
  if (files) {
    TSystemFile *file;
    TString fname;
    TIter next(files);
    while ((file=static_cast<TSystemFile*>(next()))) {
      fname = file->GetName();
      if (inname=="dir") {
        if ((file->IsDirectory() && !fname.Contains(".") && fname.Contains(regex_tag))) v_dirs.push_back(fname);
      } else  if(fname.Contains(regex_inname)) v_dirs.push_back(fname);
    }
  } // if(files)
  gSystem->cd(pwd); // The TSystemDirectory object seems to change current folder
  return v_dirs;
}

TString readFile(TString fname, vector<vector<float> > &vals, vector<vector<float> > &times, vector<int> &channels,
                 TString &maratonID, TString &date){
  ;
  TString time_s, ampm, title;
  
  ifstream infile(fname);
  string line_s;
  int channel;
  float current;
  int nrows=0;
  getline(infile, line_s);
  while(getline(infile, line_s)){
    ReplaceAll(line_s, ",", " ");
    istringstream iss(line_s);
    iss >> date >> time_s >> ampm >> maratonID >> channel >> current;
    //cout << "  " << date << "  " << time_s << "  " << ampm << "  " << maratonID << "  " << channel << "  " << current<<endl;
    bool found=false;
    for(unsigned ind=0; ind<channels.size(); ind++){
      if(channel == channels[ind]) {
        vals[ind].push_back(current);
        times[ind].push_back(toHours(time_s, ampm));
        found = true;
      }
    }// Loop over channels
    if(!found){
      vals.push_back(vector<float>()); vals[vals.size()-1].push_back(current);
      times.push_back(vector<float>()); times[times.size()-1].push_back(toHours(time_s, ampm));
      channels.push_back(channel);
    }
    nrows++;
  } // Loop over rows
  infile.close();

  title = "Currents for Maraton #" + maratonID + " on " + date;
  return title;
}

void plotCurrent(TString file){

  TString title, maratonID, date;
  vector<int> channels;
  vector<vector<float> > currents, times;
  title = readFile(file, currents, times, channels, maratonID, date);
  float minY=1e11, maxY=-1e11;
  for(size_t ichan=0; ichan<currents.size(); ichan++){
    for(unsigned ind=0; ind<currents[ichan].size(); ind++){
      //cout<<"At "<<times[ichan][ind]<<" current was "<<currents[ichan][ind]<<endl;
      if(currents[ichan][ind] > maxY) maxY = currents[ichan][ind];
      if(currents[ichan][ind] < minY) minY = currents[ichan][ind];
    }
  }

  vector<int> colors({kOrange, kGreen+1, kBlack, kCyan+1, kBlue, kMagenta+1, kOrange+10, kRed+1});
  int linw = 3;
  TString pname;
  TCanvas can;
  //can.SetGrid(); 
  can.SetFillStyle(4000);
  // TLine line; line.SetLineWidth(2); line.SetLineStyle(2); line.SetLineColor(kOrange+4);
  double legSize = 0.04;
  double legX(opts.LeftMargin()+0.03), legY(1-opts.TopMargin()-0.02), legSingle = 0.06;
  double legW = 0.7, legH = legSingle*2;
  // TLatex label; label.SetTextSize(0.058); label.SetTextFont(42); label.SetTextAlign(13);
  // label.SetTextColor(kOrange+4);

  TH1D histoModel("histoModel", "", 24, 0, 24);
  histoModel.SetMinimum(minY);
  histoModel.SetMaximum(maxY*1.3);
  histoModel.GetXaxis()->CenterTitle(true);
  histoModel.GetYaxis()->CenterTitle(true);
  histoModel.GetXaxis()->SetLabelOffset(0.01);
  histoModel.SetXTitle("Time [24 hours]");
  histoModel.SetYTitle("Current [A]");
  histoModel.SetTitle(title);
  histoModel.Draw();

  TLegend legModel(legX, legY-legH, legX+legW, legY);
  legModel.SetTextSize(legSize); legModel.SetFillColor(0); 
  legModel.SetFillStyle(0); legModel.SetBorderSize(0);
  legModel.SetNColumns(3);
  vector<TGraph*> graphs;
  for(size_t ichan=0; ichan<currents.size(); ichan++){
    graphs.push_back(new TGraph(currents[ichan].size(), &(times[ichan][0]), &(currents[ichan][0])));
    graphs[ichan]->SetLineWidth(linw); graphs[ichan]->SetLineColor(colors[channels[ichan]]);
    graphs[ichan]->Draw("same");
    TString chan_s = "Channel "; chan_s += channels[ichan];
    legModel.AddEntry(graphs[ichan], chan_s, "l");
  } // Loop over Maraton channels

  legModel.Draw();

  histoModel.Draw("axis same");
  pname = "plots/maraton"+maratonID+"_"+date+".pdf";
  can.SaveAs(pname);
 

}

int main(int argc, char *argv[]){
  GetOptions(argc, argv);

  //// Setting plot style
  setPlotStyle(opts);
  gStyle->SetGridStyle(3);

  TString txtfolder = "txt/maratons/";
  vector<TString> files = dirlist(txtfolder, "*.csv", "");
  for(size_t file=0; file<files.size(); file++){
    plotCurrent(txtfolder+files[file]);
  }

}



void GetOptions(int argc, char *argv[]){
  while(true){
    static struct option long_options[] = {
      {"file", required_argument, 0, 'f'},
      {0, 0, 0, 0}
    };

    char opt = -1;
    int option_index;
    opt = getopt_long(argc, argv, "f:", long_options, &option_index);
    if( opt == -1) break;

    string optname;
    switch(opt){
    case 'f':
      filename = optarg;
      break;
    case 0:
      optname = long_options[option_index].name;
      if(optname == ""){
        printf("Bad option! Found option name %s\n", optname.c_str());
      }else{
        printf("Bad option! Found option name %s\n", optname.c_str());
      }
      break;
    default:
      printf("Bad option! getopt_long returned character code 0%o\n", opt);
      break;
    }
  }
}
