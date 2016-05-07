#! /usr/bin/env python

from tempfile import mkdtemp
from subprocess import call
from shutil import rmtree, copytree

temp_dir = mkdtemp()

call(["git","clone","git@github.com:manuelfs/babymaker",temp_dir])
rmtree("txt/variables")
copytree(temp_dir+"/variables","txt/variables")
rmtree(temp_dir)
