#!/usr/bin/env python

###### Script to find all the closure tables
import os, sys, subprocess
import glob
import string

do_unblind = True
full_lumi = True

## Aggregate bins
methods = ["agg_himet", "agg_mixed", "agg_himult", "agg_1b"]

## 815 ipb tests
methods = ['m2lveto', 'm2lveto_el', 'm2lveto_mu', 'm2lveto_lonj', 'm2lveto_hinj', 'm2l', 'mveto', 'm2lvetomet150', 'm2lmet150', 
           'met200', 'met500', 'm5j', 'm1lmet150']

## Dilepton tests
methods = ['m2lveto', 'm2lveto_2015', 'm2lveto_lonj', 'm2lveto_hinj', 'm2l', 'mveto', 
           'm2lmet150', 'm2lvetomet150', 'mvetomet150', 'm2lmet150',
           'm2lveto_el', 'm2lveto_mu', 'm2lveto_el_2015', 'm2lveto_mu_2015', 
           'm2l_2015', 'mveto_2015', 'm2lmet150_2015', 'm2lvetomet150_2015', 
           'mvetomet150_2015', 'm2lmet150_2015']

methods = ['m2lvetomet150_lonj', 'm2lvetomet150_hinj']


os.system("./compile.sh")

options = ""
if full_lumi:  options = options + " -f"
if do_unblind: options = options + " -u"

goodtex = []
missing_files = []
for method in methods:
    cmd = "./run/table_all_preds.exe "+options+" -m "+method
    print "\n\n"+cmd+"\n"
    os.system(cmd)

    tag = "lumi"
    if full_lumi:
        if "2015" in method: tag = tag+"2p3"
        else: tag = tag+"2p6"
    else: tag = tag+"0p815"
    if do_unblind: tag = "unblind_"+tag

    texfile = "txt/table_predictions_"+tag+"_"+method+".tex"
    if os.path.isfile(texfile):  
        os.system("pdflatex "+texfile+" > /dev/null")
        goodtex.append(texfile)
    else: missing_files.append([texfile, cmd])

for mfile in missing_files:
    print "Could not find "+mfile[0]+". Command is  "+mfile[1]

sys.exit(0)
