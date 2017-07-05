#! /usr/bin/env python

name1 = "eventlists/hh4b.txt"
# name2 = "eventlists/razor.txt"
name2 = "eventlists/multilep.txt"

runs1 = []
lumis1 = []
evts1 = []

with open(name1) as f1:
    for line in f1:
        info = [int(i.strip()) for i in line.split()]
        runs1.append(info[0])
        lumis1.append(info[1])
        evts1.append(info[2])


with open(name2) as f2:
    for line in f2:
        info = [i.strip() for i in line.split()]
        run = int(info[0])
        lumi = int(info[1])
        evt = int(info[2])
        for ent in range(len(lumis1)):
            if lumis1[ent]==lumi and runs1[ent]==run and evts1[ent]==evt:
                print '{:>10}'.format(run), '{:>10}'.format(lumi), '{:>10}'.format(evt)



