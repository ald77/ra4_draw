#! /usr/bin/env python

import argparse
import os
import tempfile
import shutil
import subprocess

def Touch(file_path):
    with open(file_path, 'a'):
        os.utime(file_path, None)

def FullPath(path):
    return os.path.abspath(os.path.expanduser(path))

def Texify(input_dirs, output_dir, tag, keep_logs = False):
    orig_dir = os.getcwd()

    in_dirs = frozenset([FullPath(d) for d in input_dirs if os.path.exists(d)])

    for raw_in_dir in in_dirs:
        in_dir = FullPath(raw_in_dir)
        out_dir = ""
        if output_dir == None:
            out_dir = in_dir
        else:
            out_dir = FullPath(output_dir)
            if not os.path.exists(out_dir):
                os.makedirs(out_dir)
        tmp_dir = FullPath(tempfile.mkdtemp(prefix="tmp_texify_", dir=out_dir))

        files = [f for f in os.listdir(in_dir) if f.endswith(".tex") and (tag==None or tag in f)]
        for f in files:
            the_pdf = f.replace(".tex",".pdf")
            the_log = f.replace(".tex",".log")
            tmp_pdf = os.path.join(tmp_dir, the_pdf)
            out_pdf = os.path.join(out_dir, the_pdf)
            tmp_log = os.path.join(tmp_dir, the_log)
            out_log = os.path.join(out_dir, the_log)
            in_tex = os.path.join(in_dir, f)

            if os.path.exists(out_pdf) and os.path.getmtime(in_tex) < os.path.getmtime(out_pdf):
                print("Kept pre-existing "+out_pdf)
                continue
            elif os.path.exists(out_log) and os.path.getmtime(in_tex) < os.path.getmtime(out_log):
                print("Ignoring uncompilable "+in_tex)
                continue

            os.symlink(os.path.join(in_dir, f), os.path.join(tmp_dir, f))
            os.chdir(tmp_dir)
            command = ["pdflatex","--shell-escape","--interaction","batchmode",f]
            null_file = open(os.devnull, 'w')
            subprocess.call(command, stdout=null_file)
            if os.path.exists(tmp_pdf):
                subprocess.call(command, stdout=null_file)
                subprocess.call(command, stdout=null_file)
                os.rename(tmp_pdf, out_pdf)
                if keep_logs:
                    os.rename(tmp_log, out_log)
                print("Produced "+out_pdf)
            else:
                os.rename(tmp_log, out_log)
                print("Failed to compile "+in_tex)

            os.chdir(orig_dir)
        shutil.rmtree(tmp_dir)

if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="Compiles all .tex documents in a set of directories",
                                     formatter_class=argparse.ArgumentDefaultsHelpFormatter)
    parser.add_argument("-o", "--output", default=None, metavar="OUTPUT_DIR",
                        help="Directory in which to put resulting .pdf files. Uses input directories if omitted.")
    parser.add_argument("-t", "--tag", default=None, metavar="TAG",
                        help = "Only process files containing %(metavar)s in their names")
    parser.add_argument("-k", "--keep_logs", action='store_true',
                        help = "Keep log file even in case of successful compilation")
    parser.add_argument("input", nargs="*", default=["tables"], metavar="INPUT_DIR",
                        help="List of directories with .tex files to compile")
    args = parser.parse_args()

    Texify(args.input, args.output, args.tag, args.keep_logs)
