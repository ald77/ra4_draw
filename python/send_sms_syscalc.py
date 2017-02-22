#! /usr/bin/env python

from __future__ import print_function

import argparse
import os
import subprocess
import numpy

def fullPath(path):
    return os.path.realpath(os.path.abspath(os.path.expanduser(path)))

def ensureDir(path):
    try:
        os.makedirs(path)
    except OSError:
        if not os.path.isdir(path):
            raise

def sendSysCalc(in_dir, out_dir, num_jobs, fake_PU):
    in_dir = fullPath(in_dir)
    out_dir = fullPath(out_dir)
    run_dir = os.path.join(out_dir, "run")

    ensureDir(out_dir)
    ensureDir(run_dir)

    ra4_draw_dir = os.path.dirname(os.path.dirname(__file__))
    exe_path = fullPath(os.path.join(ra4_draw_dir,"run","ra4","syscalc_scan.exe"))

    cmssw_dir = os.path.join(os.environ["CMSSW_BASE"],"src")

    if num_jobs < 1:
        num_jobs = 1

    in_files = [ os.path.join(in_dir,f) for f in os.listdir(in_dir)
                 if os.path.isfile(os.path.join(in_dir,f))
                 and os.path.splitext(f)[1] == ".root" ]
    in_files = numpy.array_split(numpy.array(in_files), num_jobs)

    num_submitted = 0

    for sublist in in_files:
        if len(sublist) == 0:
            continue
        job_files = sublist.tolist()
        run_path = os.path.join(run_dir,"syscalc_scan_{}.sh".format(num_submitted)) 
        with open(run_path, "w") as run_file:
            os.fchmod(run_file.fileno(), 0755)
            print("#! /bin/bash", file=run_file)
            print("", file=run_file)
            print("DIRECTORY=`pwd`", file=run_file)
            print("cd {}".format(cmssw_dir), file=run_file)
            print(". /net/cms2/cms2r0/babymaker/cmsset_default.sh", file=run_file)
            print("eval `scramv1 runtime -sh`", file=run_file)
            print("cd $DIRECTORY", file=run_file)
            for i in range(len(job_files)):
                f = job_files[i]
                command = "{} -i {} -f {} -o {}".format(exe_path,os.path.dirname(f),os.path.basename(f),out_dir)
                if fake_PU:
                    command += " --fake_PU"
                print("", file=run_file)
                print("echo Starting to process file {} of {}".format(i+1, len(job_files)), file=run_file)
                print(command, file=run_file)
        subprocess.call(["JobSubmit.csh",run_path])
        num_submitted += 1
        
    print("\nSubmitted {} jobs.".format(num_submitted))
    print("Text systematics files sent to {}.".format(out_dir))
    print("Shell scripts sent to {}.".format(run_dir))

if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="Submits batch jobs to compute systematics in SMS mass plane",
                                     formatter_class=argparse.ArgumentDefaultsHelpFormatter)
    parser.add_argument("-i","--in_dir", default="/net/cms29/cms29r0/babymaker/babies/2017_02_13_grooming/T1tttt/renormed/",
                        help="Directory containing signal ntuples")
    parser.add_argument("-o","--out_dir", default="/net/cms2/cms2r0/babymaker/sys/2017_02_21/T1tttt/",
                        help="Directory in which to store text systematics files")
    parser.add_argument("-n","--njobs", type=int, default=50, help="Number of jobs to submit")
    parser.add_argument("--fake_PU", action="store_true", help="Use dummy 10/15 percent PU systematic")
    args = parser.parse_args()

    sendSysCalc(args.in_dir, args.out_dir, args.njobs, args.fake_PU)
