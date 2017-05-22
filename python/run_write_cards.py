#! /usr/bin/env python
import os

masses = ['127']
masses.extend([str(i) for i in range(150,1001,25)])

infolder = "/cms2r0/babymaker/babies/2017_03_17/TChiHH/merged_higsys_higsys/"

os.system("./compile.py")
for mass in masses:
    print 100*"="
    infile = "*SMS-TChiHH_mGluino-"+mass+"_mLSP-1_*.root"
    cmd = "./run/hig/write_datacards.exe -i "+infolder+" -f "+infile+" -o out/ -l 35.9"
    os.system(cmd)


