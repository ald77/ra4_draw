#! /usr/bin/env python

import argparse
import os
import tempfile
import shutil
import subprocess

if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="Compiles all .tex documents in a directory")
    parser.add_argument("-i", "--input", default="tables", help="Directory with .tex files to compile")
    parser.add_argument("-o", "--output", default="", help="Directory in which to put resulting .pdf files")
    args = parser.parse_args()

    orig_dir = os.getcwd()

    in_dir = os.path.abspath(os.path.expanduser(args.input))
    out_dir = ""
    if args.output == "":
        out_dir = in_dir
    else:
        out_dir = os.path.abspath(os.path.expanduser(args.output))

    tmp_dir = tempfile.mkdtemp()

    files = [f for f in os.listdir(in_dir) if f.endswith(".tex")]
    for f in files:
        os.symlink(os.path.join(in_dir, f), os.path.join(tmp_dir, f))
        os.chdir(tmp_dir)
        subprocess.call(["pdflatex","--shell-escape",f])
        subprocess.call(["pdflatex","--shell-escape",f])
        subprocess.call(["pdflatex","--shell-escape",f])
        os.chdir(orig_dir)
        the_pdf = f.replace(".tex",".pdf")
        os.rename(os.path.join(tmp_dir, the_pdf), os.path.join(out_dir, the_pdf))

    shutil.rmtree(tmp_dir)
