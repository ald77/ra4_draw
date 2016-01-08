//----------------------------------------------------------------------------
// styles - Class to set default plotting styles, read from a text file
//----------------------------------------------------------------------------

#ifndef STYLES_HH
#define STYLES_HH

#include "TH1.h"
#include "TPad.h"
#include "TString.h"

class styles {
public: 
  styles(TString group="Standard");  
  void readGroupStyle();
  void testGlobalStyle(bool fixY = true, float scale = 1000.); 
  void setGlobalStyle();
  void setDefaultStyle();
  void setHistoStyle(TH1 *h);
  void printValues();
  void moveYAxisLabel(TH1 *h, float maxi, bool isLog=false);
  void styleHist(TH1 *h, Int_t color = 1, Int_t fillstyle = 0,
                 Int_t symbol = 8,Double_t size = 0.7, Int_t width = 1);
  void setMarkers(TH1 *h, float Msize=0.6, int Mstyle=20) ;
  void setTitles(TH1 *h, TString xTitle="", TString yTitle="", TString Left="", TString Right="");
  void setTitleSizes(TH1 *h,  float size, float lsize, int font=62, 
                     float xoff=1., float yoff=1., int divisions=405);
  void parseStyleFile(TString group, TString fnames[], float *fvalues[], int nFloat, 
                      TString inames[], int *ivalues[], int nInt);
  void setGroup(TString group);

  TString confFile, Group;
  int nFont, nPads, nDivisions;
  int CanvasW, CanvasH;
  float LegendSize, TitleSize, LabelSize, xTitleOffset, yTitleOffset, zTitleOffset;
  float PadRightMargin, PadBottomMargin, PadLeftMargin, PadTopMargin;
};


#endif  // STYLES_HH
