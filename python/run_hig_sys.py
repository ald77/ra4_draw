#! /usr/bin/env python
import os

masses = ['127']
masses.extend([str(i) for i in range(150,1001,25)])

infolder = "/cms2r0/babymaker/babies/2017_02_26/TChiHH/merged_higsys_higsys/"

os.system("./compile.py")
for mass in masses:
    print 100*"="
    infile = "mergedbaby__SMS-TChiHH_mGluino-"+mass+"_mLSP-1_Tune_skim_higsys_higsys_nfiles_1.root"
    cmd = "./run/hig/syscalc_hig.exe -i "+infolder+" -f "+infile+" -o $PWD -l 35.9"
    os.system(cmd)