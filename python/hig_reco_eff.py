#! /usr/bin/env python

import argparse
from ROOT import *
import array

gStyle.SetOptStat(0)
gStyle.SetOptTitle(0)
gStyle.SetPaintTextFormat("6.2f")

pt_cuts = range(50,1001,50)

cuts = [#'all-nj4',
        # 'all-nj',
        'all-good_reco',
        # 'all_noisr-good_reco',
        # 'nj-dr',
        # 'nj-dm',
        # 'nj-am',
        # 'nj_dr-am',
        # 'nj_dm-am',
        # 'nj_dr_dm-am',
        # 'nj_am-dr',
        # 'nj_am-dm',
        'good_reco-am']
        # 'good_reco-dr',
        # 'good_reco-dm',
        # 'good_reco_noisr-am']

def getEfficiency(mass, outname, quick):

    TH1.SetDefaultSumw2()
    TH2.SetDefaultSumw2()

    folder = "/cms2r0/babymaker/babies/2017_03_17/TChiHH/merged_higmc_unskimmed/"

    c = TChain("tree", "tree")
    if (mass=="all"):
        c.Add(folder+"*SMS-TChiHH_mGluino-*mLSP-1_*.root")
    else:
        c.Add(folder+"*SMS-TChiHH_mGluino-"+mass+"*mLSP-1_*.root")

    num,den = {},{}
    for cut in cuts:                        
        num[cut] = TH2D("num_"+cut, "mass:"+mass+"num_"+cut, len(pt_cuts)-1, array.array('d', pt_cuts),
                              len(pt_cuts)-1, array.array('d', pt_cuts))

        den[cut] = num[cut].Clone(num[cut].GetName().replace( "num","den"))
    

    entry = 0
    num_entries = c.GetEntries()
    if quick: num_entries = 10000
    print "Running over ", num_entries
    maxentry = num_entries
    for event in c:
        if entry % 10000 == 0:
            print "Completed", '{:.2f}'.format(100.*entry/maxentry)+"%"
        entry = entry + 1
        if entry==maxentry: break

        # are all the b's within acceptance
        nb_in_acc = 0
        h1pt, h2pt = -999, -999
        for imc in xrange(len(c.mc_pt)):
            if abs(c.mc_id[imc])==25:
                if c.mc_pt[imc]>h1pt: 
                    h2pt = h1pt
                    h1pt = c.mc_pt[imc]
                else:
                    h2pt = c.mc_pt[imc]                   
            if abs(c.mc_id[imc])==5 and abs(c.mc_mom[imc])==25:
                if c.mc_pt[imc]>30 and abs(c.mc_eta[imc])<2.4:
                    nb_in_acc += 1

        if nb_in_acc<4: 
            continue
        elif nb_in_acc>4: 
            raise Exception("Bad Nb in acc.: {}".format(nb_in_acc))

        wgt = 1 #c.weight/c.w_btag*c.w_bhig_deep 

        # cut definitions
        pass_nj = c.njets>=4 and c.njets<=5
        pass_dr = c.higd_drmax<2.2
        pass_dm = c.higd_dm<40
        pass_am = c.higd_am>100 and c.higd_am<140
        pass_noisr = c.isr_tru_pt<20
        pass_goodreco = c.njets>=4
        for ijet in xrange(len(c.jets_pt)):
            if c.jets_h1d[ijet] or c.jets_h2d[ijet]:
                if not c.jets_hflavor[ijet]==5:
                    pass_goodreco = False
        
        # fill denom and num for each cut combination
        # den['all-nj4'].Fill(h1pt, h2pt, wgt)
        # if c.njets>=4: num['all-nj4'].Fill(h1pt, h2pt, wgt)

        # den['all-nj'].Fill(h1pt, h2pt, wgt)
        # if pass_nj: num['all-nj'].Fill(h1pt, h2pt, wgt)

        den['all-good_reco'].Fill(h1pt, h2pt, wgt)
        if pass_goodreco: num['all-good_reco'].Fill(h1pt, h2pt, wgt)

        # if pass_noisr:
        #     den['all_noisr-good_reco'].Fill(h1pt, h2pt, wgt)
        #     if pass_goodreco: num['all_noisr-good_reco'].Fill(h1pt, h2pt, wgt)

        # if pass_nj:
        #     den['nj-dr'].Fill(h1pt, h2pt, wgt)
        #     if pass_dr: num['nj-dr'].Fill(h1pt, h2pt, wgt)
        
        #     den['nj-dm'].Fill(h1pt, h2pt, wgt)
        #     if pass_dm: num['nj-dm'].Fill(h1pt, h2pt, wgt)
        
        #     den['nj-am'].Fill(h1pt, h2pt, wgt)
        #     if pass_am: num['nj-am'].Fill(h1pt, h2pt, wgt)
        
        # if pass_nj and pass_dr: 
        #     den['nj_dr-am'].Fill(h1pt, h2pt, wgt)
        #     if pass_am: num['nj_dr-am'].Fill(h1pt, h2pt, wgt)

        # if pass_nj and pass_dm: 
        #     den['nj_dm-am'].Fill(h1pt, h2pt, wgt)
        #     if pass_am: num['nj_dm-am'].Fill(h1pt, h2pt, wgt)

        # if pass_nj and pass_dr and pass_dm: 
        #     den['nj_dr_dm-am'].Fill(h1pt, h2pt, wgt)
        #     if pass_am: num['nj_dr_dm-am'].Fill(h1pt, h2pt, wgt)

        # if pass_am: 
        #     den['nj_am-dm'].Fill(h1pt, h2pt, wgt)
        #     if pass_dm: num['nj_am-dm'].Fill(h1pt, h2pt, wgt)

        #     den['nj_am-dr'].Fill(h1pt, h2pt, wgt)
        #     if pass_dr: num['nj_am-dr'].Fill(h1pt, h2pt, wgt)

        if pass_goodreco: 
            den['good_reco-am'].Fill(h1pt, h2pt, wgt)
            if pass_am: num['good_reco-am'].Fill(h1pt, h2pt, wgt)

        #     den['good_reco-dr'].Fill(h1pt, h2pt, wgt)
        #     if pass_dr: num['good_reco-dr'].Fill(h1pt, h2pt, wgt)

        #     den['good_reco-dm'].Fill(h1pt, h2pt, wgt)
        #     if pass_dm: num['good_reco-dm'].Fill(h1pt, h2pt, wgt)

        # if pass_goodreco and pass_noisr: 
        #     den['good_reco_noisr-am'].Fill(h1pt, h2pt, wgt)
        #     if pass_am: num['good_reco_noisr-am'].Fill(h1pt, h2pt, wgt)

    eff ={}
    for cut in cuts:
        eff[cut] = num[cut].Clone(num[cut].GetName().replace( "num","eff"))
        eff[cut].SetTitle(num[cut].GetTitle().replace( "num_",", denom: ").replace( "-",", num cut: "))
        eff[cut].Divide(den[cut])

    for cut in cuts:
        outfile = TFile(cut+'_'+outname, "recreate")
        # num[cut].Write()
        # den[cut].Write()
        eff[cut].Write()
        outfile.Close()

def makePlot(fname):
    gStyle.SetPalette(kBird)

    for cut in cuts:
        file = TFile(cut+'_'+fname, "read")

        for hname in ['eff']:#'num','den','eff']:
            hname+='_'+cut
            
            tMargin = 0.1
            bMargin = 0.12
            lMargin = 0.11
            rMargin = 0.15
            can = TCanvas()
            can.SetBottomMargin(bMargin)
            can.SetLeftMargin(lMargin)
            can.SetRightMargin(rMargin)
            hist = file.Get(hname)

            hist.GetXaxis().SetTitleSize(0.045)
            hist.GetXaxis().SetTitleOffset(1.15)
            hist.GetXaxis().SetLabelSize(0.04)
            hist.GetXaxis().SetLabelOffset(0.015)

            hist.GetYaxis().SetTitleSize(0.045)
            hist.GetYaxis().SetLabelSize(0.04)
            hist.GetYaxis().SetTitleOffset(1.26)

            hist.GetZaxis().SetRangeUser(0.,1.)
            hist.GetZaxis().SetTitleSize(0.045)
            hist.GetZaxis().CenterTitle(True)

            hist.GetXaxis().SetTitle("Leading Higgs boson p_{T} [GeV]")
            hist.GetYaxis().SetTitle("Subleading Higgs boson p_{T} [GeV]")
            hist.GetZaxis().SetTitle("Efficiency")
            hist.Draw("COLZTEXT")

            for i in range(len(pt_cuts)-2):
                a = TLine()
                a.SetLineWidth(1)
                a.SetLineStyle(3)
                a.SetLineColor(kBlack)
                a.DrawLine(pt_cuts[0],  pt_cuts[i+1], pt_cuts[-1],  pt_cuts[i+1])
                a.DrawLine(pt_cuts[i+1],pt_cuts[0],   pt_cuts[i+1], pt_cuts[-1])


            # title_label = TLatex()
            # title_label.DrawLatexNDC(0.1, 0.915,"#font[42]{"+hist.GetTitle()+"}")

            cmslabel = TLatex()
            cmslabel.DrawLatexNDC(lMargin, 0.915,"#font[62]{CMS}#scale[0.76]{#font[52]{ Simulation Preliminary}}")

            lumilabel = TLatex()
            lumilabel.SetTextAlign(31)
            lumilabel.DrawLatexNDC(1-rMargin, 0.915,"#font[42]{35.9 fb^{-1} (13 TeV)}")

            Xmin, Xmax, Ymin, Ymax = 85, 730, 780, 975
            box = TBox()
            box.SetFillColor(0) 
            box.SetFillStyle(1001)
            box.SetLineColor(1) 
            box.SetLineWidth(2) 
            box.SetLineStyle(1)
            box.DrawBox(Xmin, Ymin, Xmax, Ymax)
            box.SetFillColor(0) 
            box.SetFillStyle(0)
            box.SetLineColor(1) 
            box.SetLineWidth(2) 
            box.SetLineStyle(1)
            box.DrawBox(Xmin, Ymin, Xmax, Ymax)

            num = cut.split("-")[0]
            denom = cut.split("-")[1]
            caption = TLatex()
            caption.SetTextAlign(13)
            caption.SetTextSize(0.03)
            if cut=="all-good_reco":
                caption.DrawLatex(Xmin+13, Ymax-25, "#font[72]{Den:}#font[42]{ all b-quarks from hh#rightarrow4b have p_{T}>30 GeV, |#eta|<2.4}")
                caption.DrawLatex(Xmin+13, Ymax-80, "#font[72]{Num:}#font[42]{ Den. + all jets in di-Higgs reconstruction}")
                caption.DrawLatex(Xmin+13, Ymax-135,"#font[72]{        }#font[42]{ matched to B-hadrons}")
            if cut=="good_reco-am":
                caption.DrawLatex(Xmin+13, Ymax-25, "#font[72]{Den:}#font[42]{ all b-quarks from hh#rightarrow4b have p_{T}>30 GeV, |#eta|<2.4}")
                caption.DrawLatex(Xmin+13, Ymax-80, "#font[42]{ + all jets in di-Higgs reconstruction matched to B-hadrons}")
                caption.DrawLatex(Xmin+13, Ymax-135, "#font[72]{Num:}#font[42]{ Den. + 100 < #LTm#GT < 140 GeV}")


            can.Print("plots/"+hname+fname.replace(".root",".pdf"))

if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="Computes Higgs reco efficiency",
                                     formatter_class=argparse.ArgumentDefaultsHelpFormatter)
    parser.add_argument("-m", "--mass", default="7",
                        help="Higgsino mass")
    parser.add_argument("-q", "--quick", action="store_true", 
                        help="Do only one cut")
    parser.add_argument("-o", "--overwrite", action="store_true", 
                        help="Recalculate efficiency and recreate root file.")
    args = parser.parse_args()

    outname = "hig_reco_eff_"+args.mass+".root"
    if (args.overwrite): getEfficiency(args.mass, outname, args.quick)
    makePlot(outname)
