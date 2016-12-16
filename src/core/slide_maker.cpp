#include "core/slide_maker.hpp"

#include "core/utilities.hpp"

using namespace std;

SlideMaker::SlideMaker(string fname, string aspect_ratio):
  filename_(fname){
  outfile_.open("slides/"+filename_, ofstream::out);
  outfile_<<"\\documentclass[8pt,aspectratio="+aspect_ratio+"]{beamer}\n";
  outfile_<<"\\usepackage{graphicx}\n\n";
  outfile_<<"\\begin{document}\n";
}

void SlideMaker::AddSlide(vector<string> pnames, int ncols, string title){
  //set # of columns and rows
  unsigned nplots = pnames.size();
  if (ncols==-1) ncols = nplots/2;
  int nrows = nplots/ncols;
  if (nplots%ncols!=0) nrows +=1;

  string width = RoundNumber(1./ncols,2).Data();
  string height = RoundNumber(1./nrows,2).Data();
  if (title!="") height = RoundNumber(0.9/nrows,2).Data();
  outfile_<<"\\frame{\n";
  outfile_<<"  \\frametitle{"+title+"}\n";
  for (size_t i(0); i<nplots; i++){
    if (pnames[i]=="") continue;
    string line = "  \\includegraphics[height="+height+"\\textheight,width="
                  +width+"\\textwidth,keepaspectratio]{plots/"+pnames[i]+"}";
    if ((i+1)%unsigned(ncols)==0) line +="\\\\";
    outfile_<<line<<endl;
  }
  outfile_<<"}\n";
}

void SlideMaker::AddSlideWithReplace(string oldexp, string newexp, vector<string> pnames, int ncols, string title){
  for (auto &iname: pnames) ReplaceAll(iname, oldexp, newexp);
  AddSlide(pnames, ncols, title);
}

void SlideMaker::Close(){
  outfile_<<"\\end{document}\n";
  outfile_.close();

  execute("pdflatex -output-directory=slides slides/"+filename_+"  > /dev/null");
  ReplaceAll(filename_, ".tex",".pdf");
  cout<<endl<<"open slides/"+filename_<<endl;
}
