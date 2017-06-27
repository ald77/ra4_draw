#! /usr/bin/env python
import os

os.system("./compile.py")

masses = ['127']
masses.extend([str(i) for i in range(150,1001,25)])

infolder = "/cms2r0/babymaker/babies/2017_03_17/TChiHH/merged_higsys_higsys/"
os.system("./compile.py")
for mass in masses:
    print 100*"="
    infile = "*SMS-TChiHH*_mGluino-"+mass+"_mLSP-1_*.root"
    os.system("./run/hig/write_datacards.exe -i "+infolder+" -f "+infile+" -o out/ -l 35.9 --bf 1.")
    for bf in [str(x*.1) for x in range(0,10,1)]:
        os.system("./run/hig/write_datacards.exe -i "+infolder+" -f "+infile+" -o out/ -l 35.9 --bf "+bf+" --incl_nonhh")


#---------- For contamination studies
# infolder = "/cms2r0/babymaker/babies/2017_06_01/TChiHZ/merged_higsys_higsys/"
# for mass in masses:
#     print 100*"="
#     infile = "*SMS-TChiHZ*_mGluino-"+mass+"_mLSP-1_*.root"
#     os.system("./run/hig/write_datacards.exe -i "+infolder+" -f "+infile+" -o out/ -l 35.9 --bf 1.")
#     os.system("./run/hig/write_datacards.exe -i "+infolder+" -f "+infile+" -o out/ -l 35.9 --bf 1. --incl_nonbb")
#     for bf in [str(x*.1) for x in range(0,10,1)]:
#         os.system("./run/hig/write_datacards.exe -i "+infolder+" -f "+infile+" -o out/ -l 35.9 --bf "+bf)
#         os.system("./run/hig/write_datacards.exe -i "+infolder+" -f "+infile+" -o out/ -l 35.9 --bf "+bf+" --incl_nonhh --incl_nonbb")
