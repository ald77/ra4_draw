#!/usr/bin/env python

from ROOT import *
from array import array

gStyle.SetOptStat(0)
gStyle.SetOptTitle(0)

def trig_eff(ht, met):
    if(ht>   0 and ht<= 200 and met> 150 and met<= 155): return [0.532, 0.013, 0.013]
    elif (ht> 200 and ht<= 600 and met> 150 and met<= 155): return [0.612, 0.005, 0.005]
    elif (ht> 600 and ht<= 800 and met> 150 and met<= 155): return [0.589, 0.023, 0.024]
    elif (ht> 800 and ht<=1000 and met> 150 and met<= 155): return [0.515, 0.042, 0.042]
    elif (ht>1000 and ht<=9999 and met> 150 and met<= 155): return [0.588, 0.052, 0.054]
    elif (ht>   0 and ht<= 200 and met> 155 and met<= 160): return [0.591, 0.014, 0.014]
    elif (ht> 200 and ht<= 600 and met> 155 and met<= 160): return [0.684, 0.005, 0.005]
    elif (ht> 600 and ht<= 800 and met> 155 and met<= 160): return [0.678, 0.022, 0.022]
    elif (ht> 800 and ht<=1000 and met> 155 and met<= 160): return [0.537, 0.042, 0.042]
    elif (ht>1000 and ht<=9999 and met> 155 and met<= 160): return [0.511, 0.057, 0.057]
    elif (ht>   0 and ht<= 200 and met> 160 and met<= 165): return [0.619, 0.016, 0.016]
    elif (ht> 200 and ht<= 600 and met> 160 and met<= 165): return [0.727, 0.006, 0.006]
    elif (ht> 600 and ht<= 800 and met> 160 and met<= 165): return [0.699, 0.022, 0.023]
    elif (ht> 800 and ht<=1000 and met> 160 and met<= 165): return [0.690, 0.039, 0.041]
    elif (ht>1000 and ht<=9999 and met> 160 and met<= 165): return [0.568, 0.048, 0.049]
    elif (ht>   0 and ht<= 200 and met> 165 and met<= 170): return [0.678, 0.017, 0.018]
    elif (ht> 200 and ht<= 600 and met> 165 and met<= 170): return [0.769, 0.006, 0.006]
    elif (ht> 600 and ht<= 800 and met> 165 and met<= 170): return [0.732, 0.022, 0.023]
    elif (ht> 800 and ht<=1000 and met> 165 and met<= 170): return [0.609, 0.042, 0.044]
    elif (ht>1000 and ht<=9999 and met> 165 and met<= 170): return [0.685, 0.058, 0.064]
    elif (ht>   0 and ht<= 200 and met> 170 and met<= 175): return [0.670, 0.019, 0.020]
    elif (ht> 200 and ht<= 600 and met> 170 and met<= 175): return [0.811, 0.005, 0.006]
    elif (ht> 600 and ht<= 800 and met> 170 and met<= 175): return [0.779, 0.022, 0.024]
    elif (ht> 800 and ht<=1000 and met> 170 and met<= 175): return [0.736, 0.041, 0.045]
    elif (ht>1000 and ht<=9999 and met> 170 and met<= 175): return [0.663, 0.056, 0.061]
    elif (ht>   0 and ht<= 200 and met> 175 and met<= 180): return [0.730, 0.020, 0.021]
    elif (ht> 200 and ht<= 600 and met> 175 and met<= 180): return [0.838, 0.005, 0.006]
    elif (ht> 600 and ht<= 800 and met> 175 and met<= 180): return [0.820, 0.022, 0.024]
    elif (ht> 800 and ht<=1000 and met> 175 and met<= 180): return [0.819, 0.037, 0.043]
    elif (ht>1000 and ht<=9999 and met> 175 and met<= 180): return [0.736, 0.055, 0.062]
    elif (ht>   0 and ht<= 200 and met> 180 and met<= 185): return [0.745, 0.023, 0.024]
    elif (ht> 200 and ht<= 600 and met> 180 and met<= 185): return [0.874, 0.005, 0.005]
    elif (ht> 600 and ht<= 800 and met> 180 and met<= 185): return [0.848, 0.021, 0.023]
    elif (ht> 800 and ht<=1000 and met> 180 and met<= 185): return [0.869, 0.038, 0.048]
    elif (ht>1000 and ht<=9999 and met> 180 and met<= 185): return [0.759, 0.048, 0.055]
    elif (ht>   0 and ht<= 200 and met> 185 and met<= 190): return [0.777, 0.024, 0.026]
    elif (ht> 200 and ht<= 600 and met> 185 and met<= 190): return [0.903, 0.005, 0.005]
    elif (ht> 600 and ht<= 800 and met> 185 and met<= 190): return [0.850, 0.021, 0.023]
    elif (ht> 800 and ht<=1000 and met> 185 and met<= 190): return [0.839, 0.041, 0.049]
    elif (ht>1000 and ht<=9999 and met> 185 and met<= 190): return [0.847, 0.044, 0.055]
    elif (ht>   0 and ht<= 200 and met> 190 and met<= 195): return [0.792, 0.026, 0.028]
    elif (ht> 200 and ht<= 600 and met> 190 and met<= 195): return [0.907, 0.005, 0.005]
    elif (ht> 600 and ht<= 800 and met> 190 and met<= 195): return [0.884, 0.020, 0.023]
    elif (ht> 800 and ht<=1000 and met> 190 and met<= 195): return [0.870, 0.036, 0.045]
    elif (ht>1000 and ht<=9999 and met> 190 and met<= 195): return [0.781, 0.051, 0.059]
    elif (ht>   0 and ht<= 200 and met> 195 and met<= 200): return [0.757, 0.033, 0.036]
    elif (ht> 200 and ht<= 600 and met> 195 and met<= 200): return [0.924, 0.005, 0.005]
    elif (ht> 600 and ht<= 800 and met> 195 and met<= 200): return [0.921, 0.016, 0.020]
    elif (ht> 800 and ht<=1000 and met> 195 and met<= 200): return [0.936, 0.027, 0.041]
    elif (ht>1000 and ht<=9999 and met> 195 and met<= 200): return [0.803, 0.049, 0.059]
    elif (ht>   0 and ht<= 200 and met> 200 and met<= 210): return [0.841, 0.022, 0.025]
    elif (ht> 200 and ht<= 600 and met> 200 and met<= 210): return [0.949, 0.003, 0.003]
    elif (ht> 600 and ht<= 800 and met> 200 and met<= 210): return [0.927, 0.013, 0.015]
    elif (ht> 800 and ht<=1000 and met> 200 and met<= 210): return [0.894, 0.023, 0.027]
    elif (ht>1000 and ht<=9999 and met> 200 and met<= 210): return [0.839, 0.036, 0.042]
    elif (ht>   0 and ht<= 200 and met> 210 and met<= 220): return [0.850, 0.028, 0.032]
    elif (ht> 200 and ht<= 600 and met> 210 and met<= 220): return [0.966, 0.003, 0.003]
    elif (ht> 600 and ht<= 800 and met> 210 and met<= 220): return [0.952, 0.011, 0.013]
    elif (ht> 800 and ht<=1000 and met> 210 and met<= 220): return [0.919, 0.024, 0.031]
    elif (ht>1000 and ht<=9999 and met> 210 and met<= 220): return [0.959, 0.018, 0.027]
    elif (ht>   0 and ht<= 200 and met> 220 and met<= 230): return [0.896, 0.029, 0.037]
    elif (ht> 200 and ht<= 600 and met> 220 and met<= 230): return [0.973, 0.003, 0.003]
    elif (ht> 600 and ht<= 800 and met> 220 and met<= 230): return [0.979, 0.008, 0.011]
    elif (ht> 800 and ht<=1000 and met> 220 and met<= 230): return [0.956, 0.017, 0.025]
    elif (ht>1000 and ht<=9999 and met> 220 and met<= 230): return [0.971, 0.016, 0.027]
    elif (ht>   0 and ht<= 200 and met> 230 and met<= 240): return [0.844, 0.042, 0.053]
    elif (ht> 200 and ht<= 600 and met> 230 and met<= 240): return [0.983, 0.003, 0.003]
    elif (ht> 600 and ht<= 800 and met> 230 and met<= 240): return [0.976, 0.009, 0.013]
    elif (ht> 800 and ht<=1000 and met> 230 and met<= 240): return [0.983, 0.011, 0.022]
    elif (ht>1000 and ht<=9999 and met> 230 and met<= 240): return [0.942, 0.028, 0.043]
    elif (ht>   0 and ht<= 200 and met> 240 and met<= 250): return [0.880, 0.047, 0.065]
    elif (ht> 200 and ht<= 600 and met> 240 and met<= 250): return [0.985, 0.003, 0.003]
    elif (ht> 600 and ht<= 800 and met> 240 and met<= 250): return [0.992, 0.005, 0.010]
    elif (ht> 800 and ht<=1000 and met> 240 and met<= 250): return [0.989, 0.010, 0.026]
    elif (ht>1000 and ht<=9999 and met> 240 and met<= 250): return [0.931, 0.030, 0.044]
    elif (ht>   0 and ht<= 200 and met> 250 and met<= 275): return [0.915, 0.033, 0.047]
    elif (ht> 200 and ht<= 600 and met> 250 and met<= 275): return [0.989, 0.002, 0.002]
    elif (ht> 600 and ht<= 800 and met> 250 and met<= 275): return [0.992, 0.004, 0.006]
    elif (ht> 800 and ht<=1000 and met> 250 and met<= 275): return [0.984, 0.009, 0.016]
    elif (ht>1000 and ht<=9999 and met> 250 and met<= 275): return [0.965, 0.015, 0.023]
    elif (ht>   0 and ht<= 200 and met> 275 and met<= 300): return [0.862, 0.065, 0.096]
    elif (ht> 200 and ht<= 600 and met> 275 and met<= 300): return [0.992, 0.002, 0.003]
    elif (ht> 600 and ht<= 800 and met> 275 and met<= 300): return [0.989, 0.005, 0.008]
    elif (ht> 800 and ht<=1000 and met> 275 and met<= 300): return [0.963, 0.016, 0.024]
    elif (ht>1000 and ht<=9999 and met> 275 and met<= 300): return [0.991, 0.007, 0.020]
    elif (ht>   0 and ht<= 200 and met> 300 and met<=9999): return [0.744, 0.074, 0.089]
    elif (ht> 200 and ht<= 600 and met> 300 and met<=9999): return [0.994, 0.001, 0.002]
    elif (ht> 600 and ht<= 800 and met> 300 and met<=9999): return [0.996, 0.002, 0.003]
    elif (ht> 800 and ht<=1000 and met> 300 and met<=9999): return [1.000, 0.000, 0.003]
    elif (ht>1000 and ht<=9999 and met> 300 and met<=9999): return [0.987, 0.005, 0.008]
#-------------------------------------------------------
# num = 2
# bands = 255
# colors = []
# stops = [0., 1.]
# red =   [207/255., 73/255.]
# green = [224/255., 97/255.]
# blue =  [254/255., 175/255.]
# fi = TColor.CreateGradientColorTable(num, array('d',stops), array('d',red), array('d',green), array('d',blue), bands)
# for ib in range(bands):
#     colors.append(fi+ib)  
# gStyle.SetNumberContours(bands)

# NRGBs = 5
# NCont = 255
# stops = array("d",[0.00, 0.34, 0.61, 0.84, 1.00])
# red= array("d",[0.50, 0.50, 1.00, 1.00, 1.00])
# green = array("d",[ 0.50, 1.00, 1.00, 0.60, 0.50])
# blue = array("d",[1.00, 1.00, 0.50, 0.40, 0.50])
# TColor.CreateGradientColorTable(NRGBs, stops, red, green, blue, NCont)
# gStyle.SetNumberContours(NCont)

gStyle.SetPalette(kBird)

htbins = [0,200,600,800,1000,1400]
metbins = range(150,200,5) + range(200,250,10) + range(250,301,25) + [400]

can = TCanvas()
can.SetBottomMargin(0.11)
hist = TH2D('hist','hist',len(metbins)-1, array('d', metbins),len(htbins)-1, array('d', htbins))

for i, imet in enumerate(metbins):
    for j, iht in enumerate(htbins):
        hist.SetBinContent(i+1, j+1, trig_eff(iht+1, imet+1)[0])

hist.GetYaxis().SetTitle("H_{T} [GeV]")
hist.GetYaxis().SetTitleSize(0.045)
hist.GetYaxis().SetLabelSize(0.04)
hist.GetYaxis().SetTitleOffset(1.1)
hist.GetXaxis().SetTitle("p_{T}^{miss} [GeV]")
hist.GetXaxis().SetTitleSize(0.045)
hist.GetXaxis().SetTitleOffset(1.04)
hist.GetXaxis().SetLabelSize(0.04)
hist.GetXaxis().SetLabelOffset(0.015)
hist.GetZaxis().SetRangeUser(0.,1.)
hist.Draw("COLZ")

for i in range(len(htbins)-2):
    a = TLine()
    a.SetLineWidth(1)
    a.SetLineStyle(3)
    a.SetLineColor(kBlack)
    a.DrawLine(metbins[0], htbins[i+1],metbins[-1],htbins[i+1])

for i in range(len(metbins)-2):
    a = TLine()
    a.SetLineWidth(1)
    a.SetLineStyle(3)
    a.SetLineColor(kBlack)
    a.DrawLine(metbins[i+1],htbins[0],metbins[i+1],htbins[-1])

cmslabel = TLatex()
cmslabel.DrawLatexNDC(0.1, 0.915,"#font[62]{CMS}#scale[0.76]{#font[52]{ Preliminary}}")

lumilabel = TLatex()
lumilabel.SetTextAlign(31)
lumilabel.DrawLatexNDC(0.9, 0.915,"#font[42]{35.9 fb^{-1} (13 TeV)}")

can.Print("hig_trigeff.pdf")

f = TFile("SUS-16-044_trig_eff.root","RECREATE")
hist.Write()
f.Close()


