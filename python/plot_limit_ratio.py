#! /usr/bin/env python

import os, numpy
import collections
from array import array
from ROOT import *

cards = {}
with open("txt/limits/limits_TChiHH_lumi35p9_datacards.txt") as f:
    for line in f.readlines():
        cards[int(line.split()[0])] = float(line.split()[4])

num = collections.OrderedDict(sorted(cards.items()))

wspaces = {}
with open("txt/limits/limits_TChiHH_lumi35p9_paper.txt") as f:
    for line in f.readlines():
        wspaces[int(line.split()[0])] = float(line.split()[4])

denom = collections.OrderedDict(sorted(wspaces.items()))

ratios = array('d', [num[x]/denom[x]-1. for x in denom.keys()])
masses = array('d', denom.keys())

# for i in range(len(masses)):
#     print masses[i], ratios[i]

can = TCanvas()
gr = TGraph(len(masses), masses, ratios)
gr.Draw()
gr.SetLineColor( kRed+1 )
gr.SetLineWidth( 2 )
gr.SetMarkerColor( kRed+1 )
gr.SetMarkerStyle( 20 )
gr.SetTitle( 'Datacard-based limit / workspace-based limit' )
gr.GetXaxis().SetTitle( 'Higgsino mass [GeV]' )
gr.GetYaxis().SetTitle( 'Limit ratio' )
gr.Draw( 'ACP' )
can.Print("limits_ratio.pdf")