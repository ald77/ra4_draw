#! /usr/bin/env python

from __future__ import print_function

import argparse
import os

def fullPath(path):
    return os.path.realpath(os.path.abspath(os.path.expanduser(path)))

def ensureDir(path):
    try:
        os.makedirs(path)
    except OSError:
        if not os.path.isdir(path):
            raise

def ensureSubdirs(base, subdirs):
    for subdir in subdirs:
        ensureDir(os.path.join(base, subdir))

def getSubdirs(base, subdirs):
    base = fullPath(base)
    for root, dirs, files in os.walk(base):
        if root == base: continue
        subdirs.add(os.path.relpath(root, base))

def writeMakefile(subdirs, srcdir, objdir, makedir, exedir):
    ra4_draw = os.path.dirname(fullPath(__file__))
    with open(os.path.join(ra4_draw, ".subdirs.mk"), "w") as f:
        for subdir in subdirs:
            make = os.path.join("$(MAKEDIR)", os.path.join(subdir, "%.d"))
            cpp = os.path.join("$(SRCDIR)", os.path.join(subdir, "%.cpp"))
            cxx = os.path.join("$(SRCDIR)", os.path.join(subdir, "%.cxx"))
            obj = os.path.join("$(OBJDIR)", os.path.join(subdir, "%.o"))
            exe = os.path.join("$(EXEDIR)", os.path.join(subdir, "%.exe"))
            
            f.write("".join((make,": ",cpp,"\n")))
            f.write("\t$(GET_DEPS)\n\n")

            f.write("".join((make,": ",cxx,"\n")))
            f.write("\t$(GET_DEPS)\n\n")

            f.write("".join((obj,": ",cpp,"\n")))
            f.write("\t$(COMPILE)\n\n")

            f.write("".join((obj,": ",cxx,"\n")))
            f.write("\t$(COMPILE)\n\n")

            f.write("".join((exe,": ",obj," $(LIBFILE)\n")))
            f.write("\t$(LINK)\n\n")
            
def genSubdirMake(srcdir, incdir, objdir, makedir, exedir):
    print("Checking source code subdirectories")
    srcdir = fullPath(srcdir)
    incdir = fullPath(incdir)
    objdir = fullPath(objdir)
    makedir = fullPath(makedir)
    exedir = fullPath(exedir)

    subdirs = set()
    getSubdirs(srcdir, subdirs)
    getSubdirs(incdir, subdirs)

    ensureSubdirs(srcdir, subdirs)
    ensureSubdirs(incdir, subdirs)
    ensureSubdirs(objdir, subdirs)
    ensureSubdirs(makedir, subdirs)
    ensureSubdirs(exedir, subdirs)

    writeMakefile(subdirs, srcdir, objdir, makedir, exedir)

if __name__ == "__main__":
    parser = argparse.ArgumentParser(description = "Generates makefile with necessary rules for building subdirectories",
                                     formatter_class = argparse.ArgumentDefaultsHelpFormatter)
    parser.add_argument("srcdir", help = "Directory containing .cpp and .cxx files")
    parser.add_argument("incdir", help = "Directory containing .hpp files")
    parser.add_argument("objdir", help = "Directory in which to place .o files")
    parser.add_argument("makedir", help = "Directory in which to place .d files")
    parser.add_argument("exedir", help = "Directory in which to place .exe files")
    args = parser.parse_args()

    genSubdirMake(args.srcdir, args.incdir, args.objdir, args.makedir, args.exedir)
