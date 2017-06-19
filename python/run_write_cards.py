#! /usr/bin/env python
import os

# masses = ['127']
# masses.extend([str(i) for i in range(150,1001,25)])

# infolder = "/cms2r0/babymaker/babies/2017_03_17/TChiHH/merged_higsys_higsys/"

# os.system("./compile.py")
# for mass in masses:
#     print 100*"="
#     infile = "*SMS-TChiHH_mGluino-"+mass+"_mLSP-1_*.root"
#     cmd = "./run/hig/write_datacards.exe -i "+infolder+" -f "+infile+" -o out/ -l 35.9"
#     os.system(cmd)


# masses = ['400']
# masses = ['225','700']

# infolder = "/cms2r0/babymaker/babies/2017_06_01/TChiHZ/merged_higsys_higsys/"

# os.system("./compile.py")
# for mass in masses:
#     print 100*"="
#     infile = "*SMS-TChiHZ*_mGluino-"+mass+"_mLSP-1_*.root"
#     os.system("./run/hig/write_datacards.exe -i "+infolder+" -f "+infile+" -o out/ -l 35.9 --bf 1.")
#     os.system("./run/hig/write_datacards.exe -i "+infolder+" -f "+infile+" -o out/ -l 35.9 --bf 1. --incl_nonbb")
#     os.system("./run/hig/write_datacards.exe -i "+infolder+" -f "+infile+" -o out/ -l 35.9 --bf 0.9")
#     os.system("./run/hig/write_datacards.exe -i "+infolder+" -f "+infile+" -o out/ -l 35.9 --bf 0.9 --incl_nonhh --incl_nonbb")
#     os.system("./run/hig/write_datacards.exe -i "+infolder+" -f "+infile+" -o out/ -l 35.9 --bf 0.8")
#     os.system("./run/hig/write_datacards.exe -i "+infolder+" -f "+infile+" -o out/ -l 35.9 --bf 0.8 --incl_nonhh --incl_nonbb")
#     os.system("./run/hig/write_datacards.exe -i "+infolder+" -f "+infile+" -o out/ -l 35.9 --bf 0.7")
#     os.system("./run/hig/write_datacards.exe -i "+infolder+" -f "+infile+" -o out/ -l 35.9 --bf 0.7 --incl_nonhh --incl_nonbb")
#     os.system("./run/hig/write_datacards.exe -i "+infolder+" -f "+infile+" -o out/ -l 35.9 --bf 0.6")
#     os.system("./run/hig/write_datacards.exe -i "+infolder+" -f "+infile+" -o out/ -l 35.9 --bf 0.6 --incl_nonhh --incl_nonbb")


masses = ['127']
masses.extend([str(i) for i in range(150,1001,25)])

infolder = "/cms2r0/babymaker/babies/2017_03_17/TChiHH/merged_higsys_higsys/"

os.system("./compile.py")
for mass in masses:
    print 100*"="
    infile = "*SMS-TChiHH*_mGluino-"+mass+"_mLSP-1_*.root"
    # os.system("./run/hig/write_datacards.exe -i "+infolder+" -f "+infile+" -o out_old/ -l 35.9 --bf 1. --old")
    os.system("./run/hig/write_datacards.exe -i "+infolder+" -f "+infile+" -o out_old/ -l 35.9 --bf 0.9 --incl_nonhh --old")
    os.system("./run/hig/write_datacards.exe -i "+infolder+" -f "+infile+" -o out_old/ -l 35.9 --bf 0.8 --incl_nonhh --old")
    os.system("./run/hig/write_datacards.exe -i "+infolder+" -f "+infile+" -o out_old/ -l 35.9 --bf 0.7 --incl_nonhh --old")
    os.system("./run/hig/write_datacards.exe -i "+infolder+" -f "+infile+" -o out_old/ -l 35.9 --bf 0.6 --incl_nonhh --old")
    os.system("./run/hig/write_datacards.exe -i "+infolder+" -f "+infile+" -o out_old/ -l 35.9 --bf 0.5 --incl_nonhh --old")
    os.system("./run/hig/write_datacards.exe -i "+infolder+" -f "+infile+" -o out_old/ -l 35.9 --bf 0.4 --incl_nonhh --old")
    os.system("./run/hig/write_datacards.exe -i "+infolder+" -f "+infile+" -o out_old/ -l 35.9 --bf 0.3 --incl_nonhh --old")
    os.system("./run/hig/write_datacards.exe -i "+infolder+" -f "+infile+" -o out_old/ -l 35.9 --bf 0.2 --incl_nonhh --old")
    os.system("./run/hig/write_datacards.exe -i "+infolder+" -f "+infile+" -o out_old/ -l 35.9 --bf 0.1 --incl_nonhh --old")
    os.system("./run/hig/write_datacards.exe -i "+infolder+" -f "+infile+" -o out_old/ -l 35.9 --bf 0. --incl_nonhh --old")

