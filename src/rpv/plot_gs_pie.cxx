#include <iostream>
#include <string>
#include <vector>
#include <memory>

#include "TError.h"
#include "TColor.h"

#include "core/baby.hpp"
#include "core/process.hpp"
#include "core/named_func.hpp"
#include "core/plot_maker.hpp"
#include "core/plot_opt.hpp"
#include "core/palette.hpp"
#include "core/hist1d.hpp"
#include "core/hist2d.hpp"
#include "core/table.hpp"


using namespace std;
using namespace PlotOptTypes;

NamedFunc::VectorType num_bjets_from_gluon(const Baby &b, int flavor=0, int ptrange=0, int status=0, int scalar=0);
NamedFunc::ScalarType num_gluon_split(const Baby &b, int flavor=0, int ptrange=0, int status=0);
NamedFunc::ScalarType num_nonGS_true_bjets(const Baby &b);
NamedFunc::ScalarType numHSb(const Baby &b); 

NamedFunc::VectorType gs_flavs(const Baby &b);
NamedFunc::VectorType gs_pts(const Baby &b);
NamedFunc::VectorType gs_cats(const Baby &b);
NamedFunc::VectorType gs_num_bjets(const Baby &b);
const float CSVM=0.8;

int main(){
  gErrorIgnoreLevel = 6000;
 
  double lumi = 35.9;

  //bool test=true;
  
  string mc = "/net/cms29/cms29r0/babymaker/babies/2017_01_27/mc/skim_st1200/";

  Palette colors("txt/colors.txt", "default");

  NamedFunc num_gs("num_gs",[&](const Baby &b){
      return num_gluon_split(b);
    });

  NamedFunc num_gs_bb("num_gs_bb",[&](const Baby &b){
      return num_gluon_split(b,5);
    });

  NamedFunc num_gs_bb30("num_gs_bb30",[&](const Baby &b){
      return num_gluon_split(b,5,6);
    });

  NamedFunc num_gs_bb60("num_gs_bb60",[&](const Baby &b){
      return num_gluon_split(b,5,7);
    });


  NamedFunc total_gs_numb_bb30("total_gs_numb_bb30",[&](const Baby &b){
      return num_bjets_from_gluon(b,5,6,0,1);
    });


  NamedFunc numNonGSTruebjets("numNonGSTruebjets", num_nonGS_true_bjets);
  NamedFunc num_HS_b("num_HS_b",numHSb);

  NamedFunc glu_cats("glu_cats",[&](const Baby &b){
      return gs_cats(b);
    });

  NamedFunc glu_pts("glu_pts",[&](const Baby &b){
      return gs_pts(b);
    });

  NamedFunc glu_flavs("glu_flavs",[&](const Baby &b){
      return gs_flavs(b);
    });

  NamedFunc glu_nb("glu_nb",[&](const Baby &b){
      return gs_num_bjets(b);
    });

  NamedFunc total_gs_numb("total_gs_numb",[&](const Baby &b){
      return num_bjets_from_gluon(b,0,0,0,1);
    });


     
  set<string> allfiles = {mc+"*_TTJets*Lept*.root",mc+"*_QCD_HT*.root"};
  set<string> ttonly = {mc+"*_TTJets*Lept*.root"};
  set<string> qcdonly = {mc+"*_QCD_HT*.root"};
  //{mc+"fullbaby_QCD_HT2000toInf_TuneCUETP8M1_13TeV-madgraphMLM-pythia8_RunIISpring16MiniAODv2-PUSpring16_80X_mcRun2_asymptotic_2016_miniAODv2_v0-v1_9_renorm_stG1200.root"};// //
    


  PlotMaker pm;

  vector<NamedFunc> sels = {"stitch_met&&pass&&st>1200&&njets>=4","stitch_met&&pass&&st>1200&&njets>=4&&nbm>=1","stitch_met&&pass&&st>1200&&njets>=4&&nbm>=1&&mj12>500"};
  vector<NamedFunc> leps = {"nleps==0","nleps==1"};
  vector<NamedFunc> gluon_pt_cuts = {"1",glu_pts<=30.,glu_pts>30.&&glu_pts<=100.,glu_pts>100.&&glu_pts<=250.,glu_pts>250.&&glu_pts<=500.,glu_pts>500.};
  vector<NamedFunc> gluon_flav_cuts = {"1",glu_flavs==4.,glu_flavs==5.};
 
  vector<TableRow> table_cuts;

  for(unsigned int icut=0; icut<sels.size(); icut++){
    for(unsigned int ilep=0; ilep<leps.size(); ilep++){
      for(unsigned int iflav=0; iflav<gluon_flav_cuts.size(); iflav++){
	for(unsigned int ipt=0; ipt<gluon_pt_cuts.size(); ipt++){
	  //table_cuts.push_back(TableRow(Form("$%i%i%i%i$",icut,ilep,iflav,ipt),sels[icut]&&leps[ilep]&&gluon_flav_cuts[iflav]&&gluon_pt_cuts[ipt],0,0,"weight"));  
	  table_cuts.push_back(TableRow("$"+CodeToLatex((sels[icut]&&leps[ilep]&&gluon_flav_cuts[iflav]&&gluon_pt_cuts[ipt]).Name())+"$",sels[icut]&&leps[ilep]&&gluon_flav_cuts[iflav]&&gluon_pt_cuts[ipt],0,0,"weight"));  
	}
      }
    }
  }
  
  map<string, vector<shared_ptr<Process> > > gluons;

  vector< set<string> > processes = {allfiles,qcdonly,ttonly};
  vector<string> procnames = {"cats","cats_QCD","cats_ttbar"};
  for(unsigned int iproc=0; iproc<processes.size();iproc++){
    
    gluons[procnames[iproc]] = vector<shared_ptr<Process> >();
    gluons[procnames[iproc]].push_back(Process::MakeShared<Baby_full>
				       ("One or more quark is lost", Process::Type::background, kRed-4,
					processes[iproc], glu_cats==1.));
    gluons[procnames[iproc]].push_back(Process::MakeShared<Baby_full>
				       ("Both quarks in same jet", Process::Type::background, kCyan-3, 
					processes[iproc], glu_cats==2.));
    gluons[procnames[iproc]].push_back(Process::MakeShared<Baby_full>
				       ("#splitline{Quarks form two jets,}{but not both tagged}", Process::Type::background, kOrange, 
					processes[iproc], glu_cats==3.));
    gluons[procnames[iproc]].push_back(Process::MakeShared<Baby_full>
				       ("Quark is lost to lepton cleaning", Process::Type::background, kMagenta+2, 
					processes[iproc], glu_cats==5.));
    gluons[procnames[iproc]].push_back(Process::MakeShared<Baby_full>
				       ("Quarks form two b-tags", Process::Type::background, kGreen-3, 
					processes[iproc], glu_cats==4.));
  }

  gluons["nb_bb"] = vector<shared_ptr<Process> >();
  gluons["nb_bb"].push_back(Process::MakeShared<Baby_full>
			    ("Zero b-tags", Process::Type::background, kRed-4,
			     allfiles, glu_nb==0.));
  gluons["nb_bb"].push_back(Process::MakeShared<Baby_full>
			    ("One b-tag", Process::Type::background, kCyan-3, 
			     allfiles, glu_nb==1.));
  gluons["nb_bb"].push_back(Process::MakeShared<Baby_full>
			    ("Two b-tags", Process::Type::background, kGreen-3, 
			     allfiles, glu_nb==2.));

  gluons["eff_double"] = vector<shared_ptr<Process> >();
  gluons["eff_double"].push_back(Process::MakeShared<Baby_full>
				 ("Not tagged", Process::Type::background, kRed-4,
				  allfiles, glu_cats==2. && glu_nb==0.));
  gluons["eff_double"].push_back(Process::MakeShared<Baby_full>	
				 ("Tagged", Process::Type::background, kGreen-3, 
				  allfiles, glu_cats==2. && glu_nb>0.));

  for(auto &ipr: gluons) pm.Push<Table>("chart_"+ipr.first,  table_cuts, ipr.second, true, true, true, false);

    



  vector<TableRow> reco_cats;
  vector<NamedFunc> recosels = {"ntruleps<=1&&stitch&&pass&&st>1200&&njets>=4","ntruleps<=1&&stitch&&pass&&st>1200&&njets>=6","ntruleps<=1&&stitch&&pass&&st>1200&&njets>=6&&mj12>500"};
  vector<NamedFunc> nb_sels = {"nbm==1","nbm==2","nbm==3","nbm>=4"};
 
    
  for(size_t icut=0; icut<recosels.size(); icut++){
    for(size_t ilep=0; ilep<leps.size(); ilep++){
      for(size_t inb=0; inb<nb_sels.size(); inb++){
	reco_cats.push_back(TableRow("$"+CodeToLatex((sels[icut]&&leps[ilep]&&nb_sels[inb]).Name())+"$",sels[icut]&&leps[ilep]&&nb_sels[inb],0,0,"weight"));  

      }
    }
  }
  map<string, vector<shared_ptr<Process> > > recos;
  recos["gs_frac"] = vector<shared_ptr<Process> >();
  recos["gs_frac"].push_back(Process::MakeShared<Baby_full>
			     ("Has no 30 GeV GS to bb", Process::Type::background, kOrange,
			      allfiles, num_gs_bb30==0.));
  recos["gs_frac"].push_back(Process::MakeShared<Baby_full>
			     ("Has one 30 GeV GS to bb", Process::Type::background, kGreen-3,
			      allfiles, num_gs_bb30==1.));
  recos["gs_frac"].push_back(Process::MakeShared<Baby_full>
			     ("Has two or more 30 GeV GS to bb", Process::Type::background, kCyan-3,
			      allfiles, num_gs_bb30>=2.));

  recos["gs_jet_frac"] = vector<shared_ptr<Process> >();
  recos["gs_jet_frac"].push_back(Process::MakeShared<Baby_full>
				 ("Has no b-jets from GS", Process::Type::background, kOrange,
				  allfiles,  total_gs_numb==0.));
  recos["gs_jet_frac"].push_back(Process::MakeShared<Baby_full>
				 ("Has one b-jet from GS", Process::Type::background, kGreen-3,
				  allfiles, total_gs_numb==1.));
  recos["gs_jet_frac"].push_back(Process::MakeShared<Baby_full>
				 ("Has two or more b-jets from GS", Process::Type::background, kCyan-3,
				  allfiles, total_gs_numb>=2.));

  recos["non_gs_jet"] = vector<shared_ptr<Process> >();
  recos["non_gs_jet"].push_back(Process::MakeShared<Baby_full>
				("Has 0 non-GS true b-jets", Process::Type::background, kOrange,
				 allfiles,  numNonGSTruebjets==0.));
  recos["non_gs_jet"].push_back(Process::MakeShared<Baby_full>
				("Has 1 non-GS true b-jets", Process::Type::background, kGreen-3,
				 allfiles, numNonGSTruebjets==1.));
  recos["non_gs_jet"].push_back(Process::MakeShared<Baby_full>
				("Has 2 non-GS true b-jets", Process::Type::background, kCyan-3,
				 allfiles, numNonGSTruebjets==2.));

  recos["non_gs_jet"].push_back(Process::MakeShared<Baby_full>
				("Has 3 non-GS true b-jets", Process::Type::background, kRed-4,
				 allfiles, numNonGSTruebjets==3.));

  recos["non_gs_jet"].push_back(Process::MakeShared<Baby_full>
				("Has 4 non-GS true b-jets", Process::Type::background, kMagenta+2,
				 allfiles,numNonGSTruebjets>=4.));


  recos["num_HS_and_gs"] = vector<shared_ptr<Process> >();
  recos["num_HS_and_gs"].push_back(Process::MakeShared<Baby_full>
				   ("Has no GS or additional HS b's", Process::Type::background, kOrange,
				    allfiles,  total_gs_numb_bb30==0.&&num_HS_b<=2.));
  recos["num_HS_and_gs"].push_back(Process::MakeShared<Baby_full>
				   ("Has 1 additional HS b and no GS", Process::Type::background, kGreen-3,
				    allfiles,  total_gs_numb_bb30==0.&&num_HS_b==3.));
  recos["num_HS_and_gs"].push_back(Process::MakeShared<Baby_full>
				   ("Has 2 additional HS b's and no GS", Process::Type::background, kCyan-3,
				    allfiles, total_gs_numb_bb30==0.&&num_HS_b>=4.));

  recos["num_HS_and_gs"].push_back(Process::MakeShared<Baby_full>
				   ("Has GS (30 GeV, bb) but no additional HS b's", Process::Type::background, kRed-4,
				    allfiles, total_gs_numb_bb30>0.&&num_HS_b<=2.));

  recos["num_HS_and_gs"].push_back(Process::MakeShared<Baby_full>
				   ("Has GS (30 GeV, bb) and additional HS b's", Process::Type::background, kMagenta+2,
				    allfiles, total_gs_numb_bb30>0.&&num_HS_b>=3.));





  for(auto &ipr: recos) pm.Push<Table>("chart_"+ipr.first,  reco_cats, ipr.second, true, true, true, false);
    
   
  pm.min_print_ = true;
  pm.MakePlots(lumi);

}




//Generate vector of GS final state categories
NamedFunc::VectorType gs_cats(const Baby &b){
  vector<double> cats;
  for(int imc = 0; imc < static_cast<int>(b.mc_pt()->size()); ++imc){
    if(abs(b.mc_id()->at(imc)) == 21 && b.mc_gs()->at(imc)){
      int jetmatch = b.mc_gs_dau_jetmatch()->at(imc);
      if(jetmatch<=1) cats.push_back(1.); //1: 1 or more quark lost
      else if(jetmatch==2) cats.push_back(2.); //2: both quarks in same jet
      else{//each quark matched to different jet; now check if they are btagged
	int jetidx1=-1;
	int jetidx2=-1;
	for(size_t idau = 0; idau < b.mc_pt()->size(); ++idau){
	  if(b.mc_momidx()->at(idau) == imc && (abs(b.mc_id()->at(idau))==4 || abs(b.mc_id()->at(idau))==5) ){
	    if(jetidx1<0) jetidx1=b.mc_jetidx()->at(idau);
	    else if (jetidx2<0) jetidx2=b.mc_jetidx()->at(idau);
	  }
	}
	if(jetidx1>=0 && jetidx2>=0){
	  if(b.jets_pt()->at(jetidx1)>30. && (!b.jets_islep()->at(jetidx1)) && b.jets_csv()->at(jetidx1) > CSVM && 
	     b.jets_pt()->at(jetidx2)>30. && (!b.jets_islep()->at(jetidx2)) && b.jets_csv()->at(jetidx2) > CSVM) cats.push_back(4.); //4: both jets b-tagged
	  else if(b.jets_islep()->at(jetidx1) || b.jets_islep()->at(jetidx2)) cats.push_back(5.);//5: a jet is eaten by lepton
	  else cats.push_back(3.); //3: One or both jets are not btagged (but neither are leptons)

	}
	else cats.push_back(1.); //In theory this is reached literally never, but just in case..
      }
    }
  }
  return cats;
}

//Generate vector of GS flavors
NamedFunc::VectorType gs_flavs(const Baby &b){
  vector<double> flavs;
  for(int imc = 0; imc < static_cast<int>(b.mc_pt()->size()); ++imc){
    if(abs(b.mc_id()->at(imc)) == 21 && b.mc_gs()->at(imc)){
      double flav=0;
      for(size_t idau = 0; idau < b.mc_pt()->size(); ++idau){
	//Take flavor of first daughter
	if(b.mc_momidx()->at(idau) == imc){ flav = abs(b.mc_id()->at(idau)); break; }
      }
      flavs.push_back(flav);
    }
  }
  return flavs;
}

//Generate vector of GS pT
NamedFunc::VectorType gs_pts(const Baby &b){
  vector<double> pts;
  for(int imc = 0; imc < static_cast<int>(b.mc_pt()->size()); ++imc){
    if(abs(b.mc_id()->at(imc)) == 21 && b.mc_gs()->at(imc)){
      pts.push_back(b.mc_pt()->at(imc));
    }
  }
  return pts;
}

//Generate vector of num b-jets per gluon
NamedFunc::VectorType gs_num_bjets(const Baby &b){
  vector<double> num_bjets;
  for(int imc = 0; imc < static_cast<int>(b.mc_pt()->size()); ++imc){
    if(abs(b.mc_id()->at(imc)) == 21 && b.mc_gs()->at(imc)){
      int jetmatch = b.mc_gs_dau_jetmatch()->at(imc);
      int numb=0;
      if(jetmatch>0){// 1 or both jets are matched to jets
	
	int jetidx1=-1;
	int jetidx2=-1;
	for(size_t idau = 0; idau < b.mc_pt()->size(); ++idau){
	  if(b.mc_momidx()->at(idau) == imc && (abs(b.mc_id()->at(idau))==4 || abs(b.mc_id()->at(idau))==5) ){
	    if(jetidx1<0) jetidx1=b.mc_jetidx()->at(idau);
	    else if (jetidx2<0) jetidx2=b.mc_jetidx()->at(idau);
	  }
	}

	if(jetidx1>=0){
	  if(b.jets_pt()->at(jetidx1)>30. && (!b.jets_islep()->at(jetidx1)) && b.jets_csv()->at(jetidx1) > CSVM) numb++;
	}
	if(jetidx2>=0 && jetidx2!=jetidx1){
	  if(b.jets_pt()->at(jetidx2)>30. && (!b.jets_islep()->at(jetidx2)) && b.jets_csv()->at(jetidx2) > CSVM) numb++;
	}

      }
      num_bjets.push_back(numb);
    }
  }
  return num_bjets;
}







NamedFunc::ScalarType num_gluon_split(const Baby &b, int flavor, int ptrange, int status){
  int number=0;
  for(int imc = 0; imc < static_cast<int>(b.mc_pt()->size()); ++imc){
    if(abs(b.mc_id()->at(imc)) == 21 && b.mc_gs()->at(imc)){
      if(status==1 && (b.mc_status()->at(imc)>29 || b.mc_status()->at(imc) < 21)) continue; //HS gluons only
      if(status==2 && (b.mc_status()->at(imc)<=29 && b.mc_status()->at(imc) >= 21)) continue; //PS gluons only
      if(ptrange>0){
	float pt = b.mc_pt()->at(imc);
	if(ptrange==1 && pt>30) continue;
	if(ptrange==2 && (pt<=30 || pt>100)) continue;
	if(ptrange==3 && (pt<=100 || pt>250)) continue;
	if(ptrange==4 && (pt<=250 || pt>500)) continue;
	if(ptrange==5 && (pt<=500)) continue;
	if(ptrange==6 && (pt<=30)) continue;
	if(ptrange==7 && (pt<=60)) continue;
      }


      //Need to check flavor if flag is set
      if(flavor > 3){
	bool flav=false;
	for(size_t idau = 0; idau < b.mc_pt()->size(); ++idau){
	  if(b.mc_momidx()->at(idau) == imc && abs(b.mc_id()->at(idau))==flavor) flav=true;
	}
	if(flav==false) continue;
      }

      number++;
      

    }
  }
  return number;
}


NamedFunc::ScalarType num_nonGS_true_bjets(const Baby &b){
  int number=0;
  for(unsigned int ijet = 0; ijet < b.jets_pt()->size(); ++ijet){
    if(b.jets_pt()->at(ijet)>30. && (!b.jets_islep()->at(ijet)) && b.jets_csv()->at(ijet) > CSVM && abs(b.jets_hflavor()->at(ijet))==5 && b.jets_gs_index()->at(ijet)<0) number++;
  }
  return number;
}


NamedFunc::ScalarType numHSb(const Baby &b){
  int number=0;
  for(int imc = 0; imc < static_cast<int>(b.mc_pt()->size()); ++imc){
    if(b.mc_status()->at(imc)<=29 && b.mc_status()->at(imc) >= 21 && abs(b.mc_id()->at(imc))==5) number++; //Only consider hard scatter
  }
  return number;
}


NamedFunc::VectorType num_bjets_from_gluon(const Baby &b, int flavor, int ptrange, int status, int scalar){
  vector<double> numbjets;

  for(int imc = 0; imc < static_cast<int>(b.mc_pt()->size()); ++imc){
    if(abs(b.mc_id()->at(imc)) == 21 && b.mc_gs()->at(imc)){
      if(status==1 && (b.mc_status()->at(imc)>29 || b.mc_status()->at(imc) < 21)) continue; //HS gluons only
      if(status==2 && (b.mc_status()->at(imc)<=29 && b.mc_status()->at(imc) >= 21)) continue; //PS gluons only
      if(ptrange>0){
	float pt = b.mc_pt()->at(imc);
	if(ptrange==1 && pt>30) continue;
	if(ptrange==2 && (pt<=30 || pt>100)) continue;
	if(ptrange==3 && (pt<=100 || pt>250)) continue;
	if(ptrange==4 && (pt<=250 || pt>500)) continue;
	if(ptrange==5 && (pt<=500)) continue;
	if(ptrange==6 && (pt<=30)) continue;
      }


      //Need to check flavor if flag is set
      if(flavor > 3){
	bool flav=false;
	for(size_t idau = 0; idau < b.mc_pt()->size(); ++idau){
	  if(b.mc_momidx()->at(idau) == imc && abs(b.mc_id()->at(idau))==flavor) flav=true;
	}
	if(flav==false) continue;
      }


      int jetmatch = b.mc_gs_dau_jetmatch()->at(imc);
      int numb=0;
      if(jetmatch>0){// 1 or both jets are matched to jets
	
	int jetidx1=-1;
	int jetidx2=-1;
	for(size_t idau = 0; idau < b.mc_pt()->size(); ++idau){
	  if(b.mc_momidx()->at(idau) == imc && (abs(b.mc_id()->at(idau))==4 || abs(b.mc_id()->at(idau))==5) ){
	    if(jetidx1<0) jetidx1=b.mc_jetidx()->at(idau);
	    else if (jetidx2<0) jetidx2=b.mc_jetidx()->at(idau);
	  }
	}

	if(jetidx1>=0){
	  if(b.jets_pt()->at(jetidx1)>30. && (!b.jets_islep()->at(jetidx1)) && b.jets_csv()->at(jetidx1) > CSVM) numb++;
	}
	if(jetidx2>=0 && jetidx2!=jetidx1){
	  if(b.jets_pt()->at(jetidx2)>30. && (!b.jets_islep()->at(jetidx2)) && b.jets_csv()->at(jetidx2) > CSVM) numb++;
	}


      }
      numbjets.push_back(numb);
    }
  }
  if(scalar == 0 )return numbjets;
  else{
    double num_b = 0;
    for(unsigned int i=0; i<numbjets.size();i++) {
      num_b+=numbjets.at(i);
    }
    vector<double> scal_numb;
    scal_numb.push_back(num_b);
    return scal_numb;
  }
}




