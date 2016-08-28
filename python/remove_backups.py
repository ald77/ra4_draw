#! /usr/bin/env python

from __future__ import print_function

import argparse
import os

def FullPath(path):
    return os.path.abspath(os.path.expanduser(path))

def RemoveBackups(file_dir):
    full_path = FullPath(file_dir)
    if os.path.isdir(full_path):
        contents = []
        try:
            contents = os.listdir(full_path)
        except OSError as err:
            print("Cannot access directory "+full_path+" (OSError "+str(err.errno)+": "+err.strerror+")")
        for fd in contents:
            RemoveBackups(os.path.join(full_path, fd))
    elif os.path.isfile(full_path):
        base_name = os.path.basename(full_path)
        if base_name.endswith('~') or (base_name.startswith('#') and base_name.endswith('#')):
            try:
                os.remove(full_path)
            except OSError as err:
                print("Cannot remove file "+full_path+" (OSError "+str(err.errno)+": "+err.strerror+")")

if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="Removes all emacs backup files from given directories",
                                     formatter_class=argparse.ArgumentDefaultsHelpFormatter)
    parser.add_argument("input_dir", nargs="*", default=["."], metavar="INPUT_DIR",
                        help="Directories from which to remove backup files")
    args = parser.parse_args()

    for d in args.input_dir:
        RemoveBackups(d)
