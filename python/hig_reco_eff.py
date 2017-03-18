#! /usr/bin/env python

import argparse
from ROOT import *
import array

gStyle.SetOptStat(0)
gStyle.SetOptTitle(0)
gStyle.SetPaintTextFormat("6.2f")

pt_cuts = range(50,651,50)

cuts = ['all-nj4',
        'all-nj',
        'nj-dr',
        'nj-dm',
        'nj-am',
        'nj_dr-am',
        'nj_dm-am',
        'nj_dr_dm-am']

def getEfficiency(mass, outname, quick):

    TH1.SetDefaultSumw2()
    TH2.SetDefaultSumw2()

    folder = "/cms2r0/babymaker/babies/2017_02_26/TChiHH/merged_higmc_unskimmed/"

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
            if abs(c.mc_mom[imc])==25:
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

        pass_nj = c.njets>=4 and c.njets<=5
        pass_dr = c.higd_drmax<2.2
        pass_dm = c.higd_dm<40
        pass_hig = c.higd_am>100 and c.higd_am<140
        
        # fill denom and num for each cut
        den['all-nj4'].Fill(h1pt, h2pt, wgt)
        if c.njets>=4: num['all-nj4'].Fill(h1pt, h2pt, wgt)

        den['all-nj'].Fill(h1pt, h2pt, wgt)
        if pass_nj: num['all-nj'].Fill(h1pt, h2pt, wgt)

        if pass_nj:
            den['nj-dr'].Fill(h1pt, h2pt, wgt)
            if pass_dr: num['nj-dr'].Fill(h1pt, h2pt, wgt)
        
            den['nj-dm'].Fill(h1pt, h2pt, wgt)
            if pass_dm: num['nj-dm'].Fill(h1pt, h2pt, wgt)
        
            den['nj-am'].Fill(h1pt, h2pt, wgt)
            if pass_hig: num['nj-am'].Fill(h1pt, h2pt, wgt)
        
        if pass_nj and pass_dr: 
            den['nj_dr-am'].Fill(h1pt, h2pt, wgt)
            if pass_hig: num['nj_dr-am'].Fill(h1pt, h2pt, wgt)

        if pass_nj and pass_dm: 
            den['nj_dm-am'].Fill(h1pt, h2pt, wgt)
            if pass_hig: num['nj_dm-am'].Fill(h1pt, h2pt, wgt)

        if pass_nj and pass_dr and pass_dm: 
            den['nj_dr_dm-am'].Fill(h1pt, h2pt, wgt)
            if pass_hig: num['nj_dr_dm-am'].Fill(h1pt, h2pt, wgt)

    eff ={}
    for cut in cuts:
        eff[cut] = num[cut].Clone(num[cut].GetName().replace( "num","eff"))
        eff[cut].SetTitle(num[cut].GetTitle().replace( "num_",", denom: ").replace( "-",", num cut: "))
        eff[cut].Divide(den[cut])

    outfile = TFile(outname, "recreate")
    for cut in cuts:
        num[cut].Write()
        den[cut].Write()
        eff[cut].Write()
    outfile.Close()

def makePlot(fname):
    gStyle.SetPalette(kBird)

    file = TFile(fname, "read")

    for cut in cuts:
        for hname in ['num','den','eff']:
            hname+='_'+cut
            can = TCanvas()
            can.SetBottomMargin(0.11)
            can.SetRightMargin(0.15)
            hist = file.Get(hname)

            hist.GetXaxis().SetTitle("Leading Higgs boson p_{T} [GeV]")
            hist.GetXaxis().SetTitleSize(0.045)
            hist.GetXaxis().SetTitleOffset(1.06)
            hist.GetXaxis().SetLabelSize(0.04)
            hist.GetXaxis().SetLabelOffset(0.015)
            hist.GetYaxis().SetTitle("Subleading Higgs boson p_{T} [GeV]")
            hist.GetYaxis().SetTitleSize(0.045)
            hist.GetYaxis().SetLabelSize(0.04)
            hist.GetYaxis().SetTitleOffset(1.05)
            hist.GetZaxis().SetRangeUser(0.,1.)
            hist.GetZaxis().SetTitle("Higgs boson reconstruction efficiency")
            hist.GetZaxis().SetTitleSize(0.045)
            hist.GetZaxis().SetLabelSize(0.04)
            hist.Draw("COLZTEXTe")

            for i in range(len(pt_cuts)-2):
                a = TLine()
                a.SetLineWidth(1)
                a.SetLineStyle(3)
                a.SetLineColor(kBlack)
                a.DrawLine(pt_cuts[0],  pt_cuts[i+1], pt_cuts[-1],  pt_cuts[i+1])
                a.DrawLine(pt_cuts[i+1],pt_cuts[0],   pt_cuts[i+1], pt_cuts[-1])


            title_label = TLatex()
            title_label.DrawLatexNDC(0.1, 0.915,"#font[42]{"+hist.GetTitle()+"}")

            # cmslabel = TLatex()
            # cmslabel.DrawLatexNDC(0.1, 0.915,"#font[62]{CMS}#scale[0.76]{#font[52]{ Preliminary}}")

            # lumilabel = TLatex()
            # lumilabel.SetTextAlign(31)
            # lumilabel.DrawLatexNDC(0.9, 0.915,"#font[42]{35.9 fb^{-1} (13 TeV)}")

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