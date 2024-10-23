import os
import stat
#
# Coded by Marcin Modrzejewski in June 2017
#
# Generate a python script for compiling files in parallel.
# FileList if a list of lists which represent batches of files.
# Files belonging to a single batch are compiled in parallel.
# If a file in the n-th batch is changed, all files in
# the subsequent batches n+1, n+2,...
# are assumed to depend on that change and recompiled.
#
def BuildScript(FileList, FilePath):
    Batches = ""
    for batch in FileList:
        if len(batch) > 0:
            s = "BatchList.append(["
            for File in batch[:-1]:
                s += '"{}", \n'.format(File)
            s += '"{}"])\n'.format(batch[-1])
            Batches += s
    Code = """#!/usr/bin/env python
from subprocess import Popen
import datetime
import argparse
from os import path
import os
import shutil
import sys

parser = argparse.ArgumentParser()
parser.add_argument("CompilerOptions", nargs="?", default="default", help="Name of compiler options file in ./CompilerFlags/")
parser.add_argument("-np", help="Number of concurrent compiler processes", type=int, default=4)
parser.add_argument("--clean", action="store_true", help="Remove all .obj and .mod files")
args = parser.parse_args()
CompilerOptions = args.CompilerOptions
nproc = args.np
datetime_start = datetime.datetime.utcnow()
RootDir = path.dirname(path.realpath(__file__))
ObjDir = path.join(RootDir, "../obj")
if not path.exists(ObjDir):
   os.makedirs(ObjDir)
OptionsDir = path.join(RootDir, "CompilerFlags")
OptionsFile = path.join(OptionsDir, CompilerOptions)
DefaultOptionsFile = path.join(OptionsDir, "default")
CompilerCmd = open(OptionsFile).readline().rstrip()

def move_files(BaseName, DstDir):
    ObjFileFrom = BaseName + ".o"
    ObjFileTo = path.join(DstDir, BaseName + ".o")
    if path.exists(ObjFileFrom):
        shutil.move(ObjFileFrom, ObjFileTo)
    ModFileFrom = BaseName + ".mod"
    ModFileTo = path.join(DstDir, BaseName + ".mod")
    if path.exists(ModFileFrom):
        shutil.move(ModFileFrom, ModFileTo)

def test_recomp(SrcFile):
    SrcCreated = path.getmtime(SrcFile)
    BaseName, Ext = path.splitext(path.split(SrcFile)[1])
    ObjFile = path.join(ObjDir, BaseName + ".o")
    if path.exists(ObjFile):
        ObjCreated = path.getmtime(ObjFile)
        return SrcCreated >= ObjCreated
    else:
        return True

def remove_files(RelSrcPath):
    BaseName, Ext = path.splitext(path.split(RelSrcPath)[1])
    remove_list = []
    remove_list.append(path.join(ObjDir, BaseName + ".o"))
    remove_list.append(path.join(ObjDir, BaseName + ".mod"))
    remove_list.append(path.join(ObjDir, BaseName.lower() + ".mod"))
    for p in remove_list:
        if path.exists(p):
            os.remove(p)
#
# The compilation process is carried out in parallel based on the following assumptions:
# 1. Files belonging to the same batch do not depend on each other. The program will attempt to use
# the maximum amount of concurrent processes to compile the batch of files in parallel.
# 2. Files belonging to the given batch depend on all previous batches. If a file
# in n-th batch is changed, all files in batches n+1, n+2, ... are recompiled.
#
BatchList = []
{}

if args.clean:
    for batch in BatchList:
        for File in batch:
            remove_files(File)
    sys.exit()
#
# The current compiler options will be the default options
# for the next compilation. (The default options are used
# if the corresponding positional argument is skipped.)
#
if args.CompilerOptions != "default":
    shutil.copyfile(OptionsFile, DefaultOptionsFile)

RecompBatch = False
for batch in BatchList:
    processes = []
    outfiles = []
    nstarted = 0
    ncompiled = 0
    nunchanged = 0
    while ncompiled+nunchanged < len(batch):
        if nstarted - ncompiled < nproc and nunchanged+nstarted < len(batch):
            RelSrcPath = batch[nunchanged+nstarted]
            AbsSrcPath = path.join(RootDir, RelSrcPath)
            if test_recomp(AbsSrcPath) or RecompBatch:
                BaseName, Ext = path.splitext(path.split(AbsSrcPath)[1])
                Display = CompilerCmd + " " + RelSrcPath
                Command = CompilerCmd + " " + AbsSrcPath
                print(Display)
                processes.append(Popen(Command, shell=True))
                outfiles.append(BaseName)
                nstarted += 1
            else:
                nunchanged += 1
        else:
            error_code = processes[ncompiled].wait()
            if error_code != 0:
                sys.exit(error_code)
            move_files(outfiles[ncompiled], ObjDir)
            #
            # When the camel case convention is used for module names,
            # the compiler may produce lower case names for the mod files
            #
            move_files(outfiles[ncompiled].lower(), ObjDir)
            ncompiled += 1
    if ncompiled > 0:
        RecompBatch = True

datetime_finish = datetime.datetime.utcnow()
td = datetime_finish - datetime_start
dateformat = "%a %b %d %H:%M:%S UTC %Y"
datestr_finish = datetime_finish.strftime(dateformat)
print("% Compilation completed at {{}}".format(datestr_finish))
print("% Employed {{}} concurrent processes".format(nproc))
print("% Total wall clock time [s]: {{walltime:.1f}}".format(walltime=td.total_seconds()))
""".format(Batches)

    f = open(FilePath, "w")
    f.write(Code)
    f.close()
    #
    # Make the build script executable
    #
    st = os.stat(FilePath)
    os.chmod(FilePath, st.st_mode | stat.S_IEXEC)


    
