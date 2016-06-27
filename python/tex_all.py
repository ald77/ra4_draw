#!/usr/bin/env python

###### Script to texify all .tex files in a folder with a tag in the name
import os, sys, subprocess
import glob
import string
import argparse

parser = argparse.ArgumentParser(description='Texify folder.')
parser.add_argument("-i", "--infolder",  help="Folder to find tex files in", default="tables/")
parser.add_argument("-o", "--outfolder", help="Folder to write pdfs to", default="out/")
parser.add_argument("-t", "--tag",       help="Tag the .tex files need to have in the name", default="full")
args = parser.parse_args()

os.system("mkdir -p "+args.outfolder)

print
files = glob.glob(args.infolder+"/*"+args.tag+"*.tex")
for file in files:
    os.system("pdflatex -output-directory="+args.outfolder+" "+file+"  > /dev/null")
    os.system("rm "+args.outfolder+"*.log")
    os.system("rm "+args.outfolder+"*.aux")
    outname = file.replace(args.infolder, args.outfolder).replace(".tex", ".pdf")
    print "Written "+outname

print
