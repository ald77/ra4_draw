// Class defining ABCD planes

#include "core/abcd_method.hpp"

using namespace std;

//// Constructor
abcd_method::abcd_method(TString imethod, vector<TString> iplanecuts, vector<vector<TString> > ibincuts, 
			 vector<TString> iabcdcuts, TString icaption, TString ibasecuts, TString ititle):
  method(imethod),
  planecuts(iplanecuts),
  abcdcuts(iabcdcuts),
  bincuts(ibincuts),
  caption(icaption),
  basecuts(ibasecuts),
  title(ititle){

  //// By default, all planes are signal
  signalplanes = vector<bool>(planecuts.size(), true);

  // Set up dilepton and region name
  is2l = (method.Contains("2l") || method.Contains("veto"));
  rd_letter = is2l?"D":"R";

  // Set up whether we integrate R1/R3 in Nb and Njets
  int_nbnj = true;

  serializeCuts();

  } // Constructor

abcd_method::abcd_method(TString imethod, vector<TString> iplanecuts, vector<TString> ibincuts, 
			 vector<TString> iabcdcuts, TString icaption, TString ibasecuts, TString ititle):
  method(imethod),
  planecuts(iplanecuts),
  abcdcuts(iabcdcuts),
  caption(icaption),
  basecuts(ibasecuts),
  title(ititle){

  //// By default, all planes are signal
  signalplanes = vector<bool>(planecuts.size(), true);

  // Pushing bincuts for each planecut. This will allow to have different bins in different planes
  for(size_t ind=0; ind<planecuts.size(); ind++) bincuts.push_back(ibincuts);

  // Set up dilepton and region name
  is2l = (method.Contains("2l") || method.Contains("veto"));
  rd_letter = is2l?"D":"R";

  // Set up whether we integrate R1/R3 in Nb and Njets
  int_nbnj = true;

  serializeCuts();

  } // Constructor

//// Setting the planes that are signal
void abcd_method::setFirstSignalBin(int firstSigBin){
  if(firstSigBin>=static_cast<int>(planecuts.size())) {
    cout<<"Tried to set firstSigBin to "<<firstSigBin<<", but there's only "<<planecuts.size()
	<<" MET bins. Leaving all bins as signal"<<endl;
    return;
  }
  int lastBkgBin = firstSigBin-1;
  if(lastBkgBin<0) lastBkgBin = static_cast<int>(planecuts.size())-1;
  for(int bin=0; bin<=lastBkgBin; bin++)
    signalplanes[bin] = false;
}


//// Setting the integration in Nb and Njets
void abcd_method::setIntNbNj(bool int_nbnj_b){
  int_nbnj = int_nbnj_b;
  serializeCuts();
}

// Setting up all the cuts serially
void abcd_method::serializeCuts(){
  allcuts.clear();
  for(size_t iplane=0; iplane < planecuts.size(); iplane++) {
    //// Finding the OR of all bin cuts to apply to R1/R3 if int_nbnj is true
    TString c_allnbnj = "(("+bincuts[iplane][0];
    for(size_t ibin=1; ibin < bincuts[iplane].size(); ibin++)
      c_allnbnj += ")||("+bincuts[iplane][ibin];
    c_allnbnj += "))";

    //// Serializing all the cuts
    for(size_t ibin=0; ibin < bincuts[iplane].size(); ibin++){
      for(size_t iabcd=0; iabcd < abcdcuts.size(); iabcd++){
	TString totcut = "";
	if(basecuts != "") totcut += basecuts +"  &&  ";
	totcut += planecuts[iplane] +"  &&  "+ abcdcuts[iabcd];

	// Setting up the nbm/njets cuts
	if(method.Contains("agg_")){
	  totcut.ReplaceAll("nj_all_1l", "nbm>=1&&njets>=6");
	}else if(int_nbnj && bincuts[iplane].size()>1){
	  totcut.ReplaceAll("nj_all_1l", c_allnbnj);
	  totcut.ReplaceAll("nj_all_2l", lowerNjets(c_allnbnj));
	} else {
	  totcut.ReplaceAll("nj_all_1l", bincuts[iplane][ibin]);
	  totcut.ReplaceAll("nj_all_2l", lowerNjets(bincuts[iplane][ibin]));
	}
	totcut.ReplaceAll("nj_1l", bincuts[iplane][ibin]);
	totcut.ReplaceAll("nj_2l", lowerNjets(bincuts[iplane][ibin]));

	allcuts.push_back(totcut);
      } // Loop over ABCD cuts
    } // Loop over bin cuts
  } // Loop over plane cuts
} //serializeCuts

//// Returns index in allcuts of a given bin in a plane and ABCD
size_t abcd_method::indexBin(size_t indplane, size_t indbin, size_t indabcd){
  if(indplane >= planecuts.size()) {
    cout<<"Plane index "<<indplane<<" greater than Nplanes = "<<planecuts.size()<<". Exiting"<<endl;
    exit(1);
  }
  size_t Nabcd = abcdcuts.size(), index = 0;
  for(size_t iplane=0; iplane < indplane; iplane++) index += Nabcd*bincuts[iplane].size();
  
  return index + Nabcd*indbin + indabcd;
}

//// Changes all cuts to use mj12
void abcd_method::setMj12(){
  for(size_t iabcd=0; iabcd < abcdcuts.size(); iabcd++) abcdcuts[iabcd].ReplaceAll("mj14", "mj");
  basecuts.ReplaceAll("mj14", "mj");
  serializeCuts();
}


//// Changes all cuts to use only electrons, muons, or emu
void abcd_method::setLeptons(){
  for(size_t iplane=0; iplane < planecuts.size(); iplane++) {
    if(method.Contains("_el"))  planecuts[iplane].ReplaceAll("nleps", "nels");
    if(method.Contains("_mu"))  planecuts[iplane].ReplaceAll("nleps", "nmus");
    if(method.Contains("_emu")) planecuts[iplane].ReplaceAll("nleps==2", "(nels==1&&nmus==1)");
    for(size_t ibin=0; ibin < bincuts[iplane].size(); ibin++){
      if(method.Contains("_el"))  bincuts[iplane][ibin].ReplaceAll("nleps", "nels");
      if(method.Contains("_mu"))  bincuts[iplane][ibin].ReplaceAll("nleps", "nmus");
      if(method.Contains("_emu")) bincuts[iplane][ibin].ReplaceAll("nleps==2", "(nels==1&&nmus==1)");
    } // Loop over bin cuts
  } // Loop over plane cuts
  for(size_t iabcd=0; iabcd < abcdcuts.size(); iabcd++){
    if(method.Contains("_el"))  abcdcuts[iabcd].ReplaceAll("nleps", "nels");
    if(method.Contains("_mu"))  abcdcuts[iabcd].ReplaceAll("nleps", "nmus");
    if(method.Contains("_emu")) abcdcuts[iabcd].ReplaceAll("nleps==2", "(nels==1&&nmus==1)");
  } // Loop over ABCD cuts
  if(method.Contains("_el"))  basecuts.ReplaceAll("nleps", "nels");
  if(method.Contains("_mu"))  basecuts.ReplaceAll("nleps", "nmus");
  if(method.Contains("_emu")) basecuts.ReplaceAll("nleps==2", "(nels==1&&nmus==1)");

  serializeCuts();
}


//// Prints the cut for all the ABCD planes
void abcd_method::printCuts(){
  cout<<endl<<endl<<"=================== Printing cuts for method "<<method<<" ==================="<<endl;
  
  for(size_t iplane=0; iplane < planecuts.size(); iplane++) {
    cout<<endl<<" **** Plane "<<planecuts[iplane]<<" ***"<<endl;
    for(size_t ibin=0; ibin < bincuts[iplane].size(); ibin++){
      for(size_t iabcd=0; iabcd < abcdcuts.size(); iabcd++)
	cout<<allcuts[indexBin(iplane, ibin, iabcd)]<<endl;
      cout<<endl;
    } // Loop over bin cuts
  } // Loop over plane cuts
}


//// Lowers the njets requirements for dilepton bins
TString abcd_method::lowerNjets(TString &cut){
  TString lowcut = cut;
  for(int nj=6; nj<=9; nj++){
    TString nj_s = ""; nj_s += nj;
    TString njlo_s = ""; njlo_s += nj-1;
    lowcut.ReplaceAll(nj_s, njlo_s);
  }
  return lowcut;
}


abcd_method::~abcd_method(){}
