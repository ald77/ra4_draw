#! /usr/bin/env python

from __future__ import print_function

import argparse
import shutil
import errno
import os
import subprocess
import multiprocessing
import sys
import fnmatch

class DirStructure:
    def __init__(self, src, inc, make, obj, exe):
        self.src = src
        self.inc = inc
        self.make = make
        self.obj = obj
        self.exe = exe

class Term(object):
    BLACK = '\033[30m'
    RED = '\033[31m'
    GREEN = '\033[32m'
    YELLOW = '\033[33m'
    BLUE = '\033[34m'
    MAGENTA = '\033[35m'
    CYAN = '\033[36m'
    WHITE = '\033[37m'

    BOLD = '\033[1m'
    UNDERLINE = '\033[4m'

    END = '\033[0m'

def fullPath(path):
    return os.path.realpath(os.path.abspath(os.path.expanduser(path)))

def ensureDir(path):
    try:
        os.makedirs(path)
    except OSError:
        if not os.path.isdir(path):
            raise

def tryRemove(directory, pattern):
    if pattern == None:
        try: shutil.rmtree(directory)
        except OSError as e:
            if e.errno != errno.ENOENT:
                raise
        finally:
            return

    for root, dirs, files in os.walk(fullPath(directory), topdown=False):
        for f in files:
            if fnmatch.fnmatch(f, pattern):
                os.remove(os.path.join(root, f))
        try: os.rmdir(root)
        except OSError as e:
            if e.errno != errno.ENOTEMPTY:
                raise
        
def clean(dirs):
    tryRemove(".", "*~")
    tryRemove(".", "*#")
    tryRemove(dirs.exe, "*.exe")
    tryRemove(dirs.make, "*.d")
    tryRemove(dirs.obj, "*.o")
    tryRemove(dirs.obj, "*.a")
    tryRemove(dirs.inc, "baby*.hpp")
    tryRemove(dirs.src, "baby*.cpp")
    tryRemove(".", ".subdirs.mk")
    pass

def ensureSubdirs(base, subdirs):
    for subdir in subdirs:
        ensureDir(os.path.join(base, subdir))

def getSubdirs(base, subdirs):
    base = fullPath(base)
    for root, dirs, files in os.walk(base):
        if root == base: continue
        subdirs.add(os.path.relpath(root, base))
            
def genSubdirs(dirs):
    ra4_draw = os.path.dirname(fullPath(__file__))
    srcdir = fullPath(os.path.join(ra4_draw, dirs.src))
    incdir = fullPath(os.path.join(ra4_draw, dirs.inc))
    objdir = fullPath(os.path.join(ra4_draw, dirs.obj))
    makedir = fullPath(os.path.join(ra4_draw, dirs.make))
    exedir = fullPath(os.path.join(ra4_draw, dirs.exe))

    subdirs = set()
    getSubdirs(srcdir, subdirs)
    getSubdirs(incdir, subdirs)

    ensureSubdirs(srcdir, subdirs)
    ensureSubdirs(incdir, subdirs)
    ensureSubdirs(objdir, subdirs)
    ensureSubdirs(makedir, subdirs)
    ensureSubdirs(exedir, subdirs)

    return subdirs

def writeMakefile(subdirs):
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

def genSubdirMake(dirs):
    subdirs = genSubdirs(dirs)
    writeMakefile(subdirs)

def build(dirs, verbosity):
    genSubdirs(dirs)
    
    command = ["make","-j",str(multiprocessing.cpu_count()),"-k","-r","-R",
               "SRCDIR="+dirs.src, "INCDIR="+dirs.inc,
               "MAKEDIR="+dirs.make, "OBJDIR="+dirs.obj, "EXEDIR="+dirs.exe]
    if verbosity < 1:
        command.append("--silent")
    elif verbosity > 1:
        command.append("--debug")
    p = subprocess.Popen(command, stderr=subprocess.PIPE)
    err_msg = p.communicate()[1]
    if p.returncode == 0:
        print("\n\n"+Term.GREEN+Term.BOLD
              +"Compilation succeeded!"
              +Term.END+"\n")
    else:
        print("\n\n"+Term.RED+Term.BOLD
              +"################ ERRORS AND WARNINGS ################"
              +Term.END+Term.END+"\n", file=sys.stderr)
        print(err_msg.decode("utf-8"), file=sys.stderr)
        print("\n\n"+Term.RED+Term.BOLD
              +"Compilation failed."
              +Term.END+"\n", file=sys.stderr)
        sys.exit(p.returncode)
    
def compile(mode, verbosity, dirs):
    if mode == "build":
        build(dirs, verbosity)
    elif mode == "clean":
        clean(dirs)
    elif mode == "set_dirs":
        genSubdirMake(dirs)
    elif mode== "print_vars":
        subprocess.check_call(["make","test","-r","-R","--silent",
                               "SRCDIR="+dirs.src, "INCDIR="+dirs.inc,
                               "MAKEDIR="+dirs.make, "OBJDIR="+dirs.obj, "EXEDIR="+dirs.exe, "print_vars"])
    else:
        raise Exception("Unrecognized option: "+mode)

if __name__ == "__main__":
    parser = argparse.ArgumentParser(description = "Compiles ra4_draw code",
                                     formatter_class = argparse.ArgumentDefaultsHelpFormatter)
    parser.add_argument("mode", nargs="?", default="build", choices=["build","clean","set_dirs","print_vars"],
                        help = "Selects which action to perform")
    parser.add_argument("-v","--verbosity", type=int, default=1, choices=[0,1,2],
                        help = "Set verbosity. Lower = less printing.")
    parser.add_argument("--src_dir", default = "src",
                        help = "Directory containing .cpp and .cxx files")
    parser.add_argument("--inc_dir", default = "inc",
                        help = "Directory containing .hpp files")
    parser.add_argument("--make_dir", default = "bin",
                        help = "Directory in which to store .d files")
    parser.add_argument("--obj_dir", default = "bin",
                        help = "Directory in which to place .o files")
    parser.add_argument("--exe_dir", default = "run",
                        help = "Directory in which to store .exe files")
    args = parser.parse_args()

    compile(args.mode, args.verbosity,
            DirStructure(args.src_dir, args.inc_dir, args.make_dir, args.obj_dir, args.exe_dir))
