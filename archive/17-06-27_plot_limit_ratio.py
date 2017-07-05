#! /usr/bin/env python

import os, numpy
import collections
from array import array
from ROOT import *

index = 4
index = 7

gr = {}
can = TCanvas()
leg = TLegend(0.5,0.2,0.6,0.4)
leg.SetLineWidth(0)
leg.SetTextSize(0.05)

for index in [4,7]:
    wspaces = {}
    with open("txt/limits/limits_TChiHH_lumi35p9_datacards.txt") as f:
        for line in f.readlines():
            mass = int(line.split()[0])
            if (mass<175): continue
            wspaces[mass] = float(line.split()[index])

    denom = collections.OrderedDict(sorted(wspaces.items()))
    

    files = {
    "lnU": "txt/limits/limits_TChiHH_lumi35p9_corrected.txt",
    }

    colors = [kGreen+1, kAzure,kRed+1,  kGreen+1]

    i=0
    for type in files.keys():
        cards = {}
        with open(files[type]) as f:
            for line in f.readlines():
                mass = int(line.split()[0])
                if (mass<175): continue
                cards[mass] = float(line.split()[index])

        num = collections.OrderedDict(sorted(cards.items()))

        ratios = array('d', [num[x]/denom[x] for x in denom.keys()])
        masses = array('d', denom.keys())

        key = type+str(index)
        gr[key] = TGraph(len(masses), masses, ratios)
        gr[key].SetLineColor(colors[i])
        gr[key].SetLineWidth(2)
        gr[key].SetMarkerStyle(20+i)
        if (index==7): 
            gr[key].SetLineStyle(2)
            gr[key].SetMarkerStyle(24+i)
        gr[key].SetMarkerColor(colors[i])
        gr[key].SetTitle("")
        gr[key].GetXaxis().SetTitle('Higgsino mass [GeV]')
        gr[key].GetYaxis().SetTitle('corrected limit / datacard limit')
        gr[key].GetYaxis().SetTitleOffset(1.2)
        gr[key].GetYaxis().SetRangeUser(0.7, 1.2)
        if (index==7): 
            gr[key].Draw("CP same")
            leg.AddEntry(gr[key], "Expected", "LP")
        else: 
            gr[key].Draw('ACP')
            leg.AddEntry(gr[key], "Observed", "LP")
        i+=1

leg.Draw()
can.Print("limits_ratio.pdf")
