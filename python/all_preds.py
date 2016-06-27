#!/usr/bin/env python

###### Script to find all the closure tables
import os, sys, subprocess
import glob
import string

os.system("./compile.sh")

##  -m   Method to run on (if you just want one)
##  -d   Debug: prints yields and cuts used
##  -n   Does not print signal columns

##  -b   Prints Other, tt1l, tt2l contributions
##  -l   Does tables for e/mu/emu as well
##  -u   Unblinds R4/D4
##  -f   Uses all data (does not apply nonblind)
##  -2   Makes tables only for dilepton tests
##  -s   Which skim to use: standard, met150, 2015 data


## 815 ipb
os.system("./run/table_preds.exe -u -l -b")
os.system("./run/table_preds.exe -u -l -b -s met150")
## Full lumi
os.system("./run/table_preds.exe -u -l -b -f -2")
os.system("./run/table_preds.exe -u -l -b -f -2 -s met150")
os.system("./run/table_preds.exe -u -l -b -f -2 -s 2015 ")

## Making pdfs and organizing them
os.system("./python/tex_all.py")
folders = ["lumi0p815", "lumi2p6", "2015"]
for folder in folders:
    os.sytem("mkdir out/"+folder)
    os.system("mv out/*"+folder+"*pdf out/"+folder)


sys.exit(0)
