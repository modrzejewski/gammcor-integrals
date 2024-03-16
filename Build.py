#!/usr/bin/env python3
from subprocess import Popen
import datetime
import argparse
from os import path
import os
import shutil
import sys

parser = argparse.ArgumentParser()
parser.add_argument("CompilerOptions", default="", nargs="?", help="Name of a subdirectory in ./CompilerFlags/")
parser.add_argument("-np", help="Number of concurrent compiler processes", type=int, default=4)
parser.add_argument("--clean", action="store_true", help="Remove all .obj and .mod files")
args = parser.parse_args()
CompilerOptions = args.CompilerOptions
nproc = args.np
datetime_start = datetime.datetime.utcnow()
RootDir = path.dirname(path.realpath(__file__))
ObjDir = path.join(RootDir, "obj")
IncludeDir = path.join(RootDir, "include")
if not path.exists(ObjDir):
   os.makedirs(ObjDir)
if not path.exists(IncludeDir):
    os.makedirs(IncludeDir)
LibDir = path.join(RootDir, "lib")
if not path.exists(LibDir):
    os.makedirs(LibDir)

if not args.clean:
    OptionsDir = path.join(RootDir, "CompilerFlags")
    CompilerFile = path.join(OptionsDir, CompilerOptions, "compiler")
    LinkerFile = path.join(OptionsDir, CompilerOptions, "linker")
    ArchiveFile = path.join(OptionsDir, CompilerOptions, "static-library")
    CompilerCmd = open(CompilerFile).readline().format(Include=IncludeDir).rstrip()
    LinkerCmd = open(LinkerFile).readline().format(Include=IncludeDir).rstrip()
    ArchiveCmd = open(ArchiveFile).readline().rstrip()

StaticLibFile = path.join(LibDir, "cholesky.a")
ExeFile = "test"

def move_files(BaseName, LibDir, IncludeDir):
    ObjFileFrom = BaseName + ".o"
    ObjFileTo = path.join(LibDir, BaseName + ".o")
    if path.exists(ObjFileFrom):
        shutil.move(ObjFileFrom, ObjFileTo)
    ModFileFrom = BaseName + ".mod"
    ModFileTo = path.join(IncludeDir, BaseName + ".mod")
    if path.exists(ModFileFrom):
        shutil.move(ModFileFrom, ModFileTo)

def test_recomp(SrcFile):
    SrcCreated = path.getmtime(SrcFile)
    BaseName, Ext = path.splitext(path.split(SrcFile)[1])
    ObjFile = path.join(ObjDir, BaseName + ".o")
    ModFile = path.join(IncludeDir, BaseName + ".mod")
    ModFile2 = path.join(IncludeDir, BaseName.lower() + ".mod")
    if path.exists(ObjFile) and (path.exists(ModFile) or path.exists(ModFile2)):
        ObjCreated = path.getmtime(ObjFile)
        return SrcCreated >= ObjCreated
    else:
        return True

def remove_files(RelSrcPath):
    BaseName, Ext = path.splitext(path.split(RelSrcPath)[1])
    remove_list = []
    remove_list.append(path.join(ObjDir, BaseName + ".o"))
    remove_list.append(path.join(IncludeDir, BaseName + ".mod"))
    remove_list.append(path.join(IncludeDir, BaseName.lower() + ".mod"))
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
BatchList.append(["src/common/arithmetic.f90"])
BatchList.append(["src/common/math_constants.f90"])
BatchList.append(["src/common/string.f90"])
BatchList.append(["src/common/clock.f90"])
BatchList.append(["src/common/display.f90"])
BatchList.append(["src/common/io.f90",
                  "src/integrals/spherh.f90"])
BatchList.append(["src/integrals/chebinterp.f90"])
BatchList.append(["src/integrals/boys.f90"])
BatchList.append(["src/integrals/grid/grid_definitions.f90",
                  "src/integrals/grid/LebedevGrid.f90"])
BatchList.append(["src/integrals/THC/thc_definitions.f90"])
BatchList.append(["src/integrals/Auto2e/src/auto2e_BraTransform.f90", 
"src/integrals/Auto2e/src/auto2e_Hermite.f90", 
"src/integrals/Auto2e/src/auto2e_KetTransfer.f90", 
"src/integrals/Auto2e/src/auto2e_KetTransform_10_1.f90", 
"src/integrals/Auto2e/src/auto2e_KetTransform_10_10.f90", 
"src/integrals/Auto2e/src/auto2e_KetTransform_10_2.f90", 
"src/integrals/Auto2e/src/auto2e_KetTransform_10_3.f90", 
"src/integrals/Auto2e/src/auto2e_KetTransform_10_4.f90", 
"src/integrals/Auto2e/src/auto2e_KetTransform_10_5.f90", 
"src/integrals/Auto2e/src/auto2e_KetTransform_10_6.f90", 
"src/integrals/Auto2e/src/auto2e_KetTransform_10_7.f90", 
"src/integrals/Auto2e/src/auto2e_KetTransform_10_8.f90", 
"src/integrals/Auto2e/src/auto2e_KetTransform_10_9.f90", 
"src/integrals/Auto2e/src/auto2e_KetTransform_2_1.f90", 
"src/integrals/Auto2e/src/auto2e_KetTransform_2_2.f90", 
"src/integrals/Auto2e/src/auto2e_KetTransform_3_1.f90", 
"src/integrals/Auto2e/src/auto2e_KetTransform_3_2.f90", 
"src/integrals/Auto2e/src/auto2e_KetTransform_3_3.f90", 
"src/integrals/Auto2e/src/auto2e_KetTransform_4_1.f90", 
"src/integrals/Auto2e/src/auto2e_KetTransform_4_2.f90", 
"src/integrals/Auto2e/src/auto2e_KetTransform_4_3.f90", 
"src/integrals/Auto2e/src/auto2e_KetTransform_4_4.f90", 
"src/integrals/Auto2e/src/auto2e_KetTransform_5_1.f90", 
"src/integrals/Auto2e/src/auto2e_KetTransform_5_2.f90", 
"src/integrals/Auto2e/src/auto2e_KetTransform_5_3.f90", 
"src/integrals/Auto2e/src/auto2e_KetTransform_5_4.f90", 
"src/integrals/Auto2e/src/auto2e_KetTransform_5_5.f90", 
"src/integrals/Auto2e/src/auto2e_KetTransform_6_1.f90", 
"src/integrals/Auto2e/src/auto2e_KetTransform_6_2.f90", 
"src/integrals/Auto2e/src/auto2e_KetTransform_6_3.f90", 
"src/integrals/Auto2e/src/auto2e_KetTransform_6_4.f90", 
"src/integrals/Auto2e/src/auto2e_KetTransform_6_5_part1.f90", 
"src/integrals/Auto2e/src/auto2e_KetTransform_6_5_part2.f90", 
"src/integrals/Auto2e/src/auto2e_KetTransform_6_6_part1.f90", 
"src/integrals/Auto2e/src/auto2e_KetTransform_6_6_part2.f90", 
"src/integrals/Auto2e/src/auto2e_KetTransform_6_6_part3.f90", 
"src/integrals/Auto2e/src/auto2e_KetTransform_7_1.f90", 
"src/integrals/Auto2e/src/auto2e_KetTransform_7_2.f90", 
"src/integrals/Auto2e/src/auto2e_KetTransform_7_3.f90", 
"src/integrals/Auto2e/src/auto2e_KetTransform_7_4.f90", 
"src/integrals/Auto2e/src/auto2e_KetTransform_7_5.f90", 
"src/integrals/Auto2e/src/auto2e_KetTransform_7_6.f90", 
"src/integrals/Auto2e/src/auto2e_KetTransform_7_7.f90", 
"src/integrals/Auto2e/src/auto2e_KetTransform_8_1.f90", 
"src/integrals/Auto2e/src/auto2e_KetTransform_8_2.f90", 
"src/integrals/Auto2e/src/auto2e_KetTransform_8_3.f90", 
"src/integrals/Auto2e/src/auto2e_KetTransform_8_4.f90", 
"src/integrals/Auto2e/src/auto2e_KetTransform_8_5.f90", 
"src/integrals/Auto2e/src/auto2e_KetTransform_8_6.f90", 
"src/integrals/Auto2e/src/auto2e_KetTransform_8_7.f90", 
"src/integrals/Auto2e/src/auto2e_KetTransform_8_8.f90", 
"src/integrals/Auto2e/src/auto2e_KetTransform_9_1.f90", 
"src/integrals/Auto2e/src/auto2e_KetTransform_9_2.f90", 
"src/integrals/Auto2e/src/auto2e_KetTransform_9_3.f90", 
"src/integrals/Auto2e/src/auto2e_KetTransform_9_4.f90", 
"src/integrals/Auto2e/src/auto2e_KetTransform_9_5.f90", 
"src/integrals/Auto2e/src/auto2e_KetTransform_9_6.f90", 
"src/integrals/Auto2e/src/auto2e_KetTransform_9_7.f90", 
"src/integrals/Auto2e/src/auto2e_KetTransform_9_8.f90", 
"src/integrals/Auto2e/src/auto2e_KetTransform_9_9.f90", 
"src/integrals/Auto2e/src/auto2e_SpherTransf.f90"])
BatchList.append(["src/integrals/Auto2e/src/auto2e_WMatrix_part1.f90", 
"src/integrals/Auto2e/src/auto2e_WMatrix_part2.f90", 
"src/integrals/Auto2e/src/auto2e_WMatrix_part3.f90"])
BatchList.append(["src/integrals/Auto2e/src/auto2e_eri_dddd.f90", 
"src/integrals/Auto2e/src/auto2e_eri_dddp.f90", 
"src/integrals/Auto2e/src/auto2e_eri_ddds.f90", 
"src/integrals/Auto2e/src/auto2e_eri_ddpp.f90", 
"src/integrals/Auto2e/src/auto2e_eri_ddps.f90", 
"src/integrals/Auto2e/src/auto2e_eri_ddss.f90", 
"src/integrals/Auto2e/src/auto2e_eri_dpdp.f90", 
"src/integrals/Auto2e/src/auto2e_eri_dpds.f90", 
"src/integrals/Auto2e/src/auto2e_eri_dppp.f90", 
"src/integrals/Auto2e/src/auto2e_eri_dpps.f90", 
"src/integrals/Auto2e/src/auto2e_eri_dpss.f90", 
"src/integrals/Auto2e/src/auto2e_eri_dsds.f90", 
"src/integrals/Auto2e/src/auto2e_eri_dspp.f90", 
"src/integrals/Auto2e/src/auto2e_eri_dsps.f90", 
"src/integrals/Auto2e/src/auto2e_eri_dsss.f90", 
"src/integrals/Auto2e/src/auto2e_eri_fddd.f90", 
"src/integrals/Auto2e/src/auto2e_eri_fddp.f90", 
"src/integrals/Auto2e/src/auto2e_eri_fdds.f90", 
"src/integrals/Auto2e/src/auto2e_eri_fdfd.f90", 
"src/integrals/Auto2e/src/auto2e_eri_fdfp.f90", 
"src/integrals/Auto2e/src/auto2e_eri_fdfs.f90", 
"src/integrals/Auto2e/src/auto2e_eri_fdpp.f90", 
"src/integrals/Auto2e/src/auto2e_eri_fdps.f90", 
"src/integrals/Auto2e/src/auto2e_eri_fdss.f90", 
"src/integrals/Auto2e/src/auto2e_eri_ffdd.f90", 
"src/integrals/Auto2e/src/auto2e_eri_ffdp.f90", 
"src/integrals/Auto2e/src/auto2e_eri_ffds.f90", 
"src/integrals/Auto2e/src/auto2e_eri_fffd.f90", 
"src/integrals/Auto2e/src/auto2e_eri_ffff.f90", 
"src/integrals/Auto2e/src/auto2e_eri_fffp.f90", 
"src/integrals/Auto2e/src/auto2e_eri_fffs.f90", 
"src/integrals/Auto2e/src/auto2e_eri_ffpp.f90", 
"src/integrals/Auto2e/src/auto2e_eri_ffps.f90", 
"src/integrals/Auto2e/src/auto2e_eri_ffss.f90", 
"src/integrals/Auto2e/src/auto2e_eri_fpdd.f90", 
"src/integrals/Auto2e/src/auto2e_eri_fpdp.f90", 
"src/integrals/Auto2e/src/auto2e_eri_fpds.f90", 
"src/integrals/Auto2e/src/auto2e_eri_fpfp.f90", 
"src/integrals/Auto2e/src/auto2e_eri_fpfs.f90", 
"src/integrals/Auto2e/src/auto2e_eri_fppp.f90", 
"src/integrals/Auto2e/src/auto2e_eri_fpps.f90", 
"src/integrals/Auto2e/src/auto2e_eri_fpss.f90", 
"src/integrals/Auto2e/src/auto2e_eri_fsdd.f90", 
"src/integrals/Auto2e/src/auto2e_eri_fsdp.f90", 
"src/integrals/Auto2e/src/auto2e_eri_fsds.f90", 
"src/integrals/Auto2e/src/auto2e_eri_fsfs.f90", 
"src/integrals/Auto2e/src/auto2e_eri_fspp.f90", 
"src/integrals/Auto2e/src/auto2e_eri_fsps.f90", 
"src/integrals/Auto2e/src/auto2e_eri_fsss.f90", 
"src/integrals/Auto2e/src/auto2e_eri_gddd.f90", 
"src/integrals/Auto2e/src/auto2e_eri_gddp.f90", 
"src/integrals/Auto2e/src/auto2e_eri_gdds.f90", 
"src/integrals/Auto2e/src/auto2e_eri_gdfd.f90", 
"src/integrals/Auto2e/src/auto2e_eri_gdff.f90", 
"src/integrals/Auto2e/src/auto2e_eri_gdfp.f90", 
"src/integrals/Auto2e/src/auto2e_eri_gdfs.f90", 
"src/integrals/Auto2e/src/auto2e_eri_gdgd.f90", 
"src/integrals/Auto2e/src/auto2e_eri_gdgp.f90", 
"src/integrals/Auto2e/src/auto2e_eri_gdgs.f90", 
"src/integrals/Auto2e/src/auto2e_eri_gdpp.f90", 
"src/integrals/Auto2e/src/auto2e_eri_gdps.f90", 
"src/integrals/Auto2e/src/auto2e_eri_gdss.f90", 
"src/integrals/Auto2e/src/auto2e_eri_gfdd.f90", 
"src/integrals/Auto2e/src/auto2e_eri_gfdp.f90", 
"src/integrals/Auto2e/src/auto2e_eri_gfds.f90", 
"src/integrals/Auto2e/src/auto2e_eri_gffd.f90", 
"src/integrals/Auto2e/src/auto2e_eri_gfff.f90", 
"src/integrals/Auto2e/src/auto2e_eri_gffp.f90", 
"src/integrals/Auto2e/src/auto2e_eri_gffs.f90", 
"src/integrals/Auto2e/src/auto2e_eri_gfgd.f90", 
"src/integrals/Auto2e/src/auto2e_eri_gfgf.f90", 
"src/integrals/Auto2e/src/auto2e_eri_gfgp.f90", 
"src/integrals/Auto2e/src/auto2e_eri_gfgs.f90", 
"src/integrals/Auto2e/src/auto2e_eri_gfpp.f90", 
"src/integrals/Auto2e/src/auto2e_eri_gfps.f90", 
"src/integrals/Auto2e/src/auto2e_eri_gfss.f90", 
"src/integrals/Auto2e/src/auto2e_eri_ggdd.f90", 
"src/integrals/Auto2e/src/auto2e_eri_ggdp.f90", 
"src/integrals/Auto2e/src/auto2e_eri_ggds.f90", 
"src/integrals/Auto2e/src/auto2e_eri_ggfd.f90", 
"src/integrals/Auto2e/src/auto2e_eri_ggff.f90", 
"src/integrals/Auto2e/src/auto2e_eri_ggfp.f90", 
"src/integrals/Auto2e/src/auto2e_eri_ggfs.f90", 
"src/integrals/Auto2e/src/auto2e_eri_gggd.f90", 
"src/integrals/Auto2e/src/auto2e_eri_gggf.f90", 
"src/integrals/Auto2e/src/auto2e_eri_gggg.f90", 
"src/integrals/Auto2e/src/auto2e_eri_gggp.f90", 
"src/integrals/Auto2e/src/auto2e_eri_gggs.f90", 
"src/integrals/Auto2e/src/auto2e_eri_ggpp.f90", 
"src/integrals/Auto2e/src/auto2e_eri_ggps.f90", 
"src/integrals/Auto2e/src/auto2e_eri_ggss.f90", 
"src/integrals/Auto2e/src/auto2e_eri_gpdd.f90", 
"src/integrals/Auto2e/src/auto2e_eri_gpdp.f90", 
"src/integrals/Auto2e/src/auto2e_eri_gpds.f90", 
"src/integrals/Auto2e/src/auto2e_eri_gpfd.f90", 
"src/integrals/Auto2e/src/auto2e_eri_gpff.f90", 
"src/integrals/Auto2e/src/auto2e_eri_gpfp.f90", 
"src/integrals/Auto2e/src/auto2e_eri_gpfs.f90", 
"src/integrals/Auto2e/src/auto2e_eri_gpgp.f90", 
"src/integrals/Auto2e/src/auto2e_eri_gpgs.f90", 
"src/integrals/Auto2e/src/auto2e_eri_gppp.f90", 
"src/integrals/Auto2e/src/auto2e_eri_gpps.f90", 
"src/integrals/Auto2e/src/auto2e_eri_gpss.f90", 
"src/integrals/Auto2e/src/auto2e_eri_gsdd.f90", 
"src/integrals/Auto2e/src/auto2e_eri_gsdp.f90", 
"src/integrals/Auto2e/src/auto2e_eri_gsds.f90", 
"src/integrals/Auto2e/src/auto2e_eri_gsfd.f90", 
"src/integrals/Auto2e/src/auto2e_eri_gsff.f90", 
"src/integrals/Auto2e/src/auto2e_eri_gsfp.f90", 
"src/integrals/Auto2e/src/auto2e_eri_gsfs.f90", 
"src/integrals/Auto2e/src/auto2e_eri_gsgs.f90", 
"src/integrals/Auto2e/src/auto2e_eri_gspp.f90", 
"src/integrals/Auto2e/src/auto2e_eri_gsps.f90", 
"src/integrals/Auto2e/src/auto2e_eri_gsss.f90", 
"src/integrals/Auto2e/src/auto2e_eri_hddd.f90", 
"src/integrals/Auto2e/src/auto2e_eri_hddp.f90", 
"src/integrals/Auto2e/src/auto2e_eri_hdds.f90", 
"src/integrals/Auto2e/src/auto2e_eri_hdfd.f90", 
"src/integrals/Auto2e/src/auto2e_eri_hdff.f90", 
"src/integrals/Auto2e/src/auto2e_eri_hdfp.f90", 
"src/integrals/Auto2e/src/auto2e_eri_hdfs.f90", 
"src/integrals/Auto2e/src/auto2e_eri_hdgd.f90", 
"src/integrals/Auto2e/src/auto2e_eri_hdgf.f90", 
"src/integrals/Auto2e/src/auto2e_eri_hdgg.f90", 
"src/integrals/Auto2e/src/auto2e_eri_hdgp.f90", 
"src/integrals/Auto2e/src/auto2e_eri_hdgs.f90", 
"src/integrals/Auto2e/src/auto2e_eri_hdhd.f90", 
"src/integrals/Auto2e/src/auto2e_eri_hdhp.f90", 
"src/integrals/Auto2e/src/auto2e_eri_hdhs.f90", 
"src/integrals/Auto2e/src/auto2e_eri_hdpp.f90", 
"src/integrals/Auto2e/src/auto2e_eri_hdps.f90", 
"src/integrals/Auto2e/src/auto2e_eri_hdss.f90", 
"src/integrals/Auto2e/src/auto2e_eri_hfdd.f90", 
"src/integrals/Auto2e/src/auto2e_eri_hfdp.f90", 
"src/integrals/Auto2e/src/auto2e_eri_hfds.f90", 
"src/integrals/Auto2e/src/auto2e_eri_hffd.f90", 
"src/integrals/Auto2e/src/auto2e_eri_hfff.f90", 
"src/integrals/Auto2e/src/auto2e_eri_hffp.f90", 
"src/integrals/Auto2e/src/auto2e_eri_hffs.f90", 
"src/integrals/Auto2e/src/auto2e_eri_hfgd.f90", 
"src/integrals/Auto2e/src/auto2e_eri_hfgf.f90", 
"src/integrals/Auto2e/src/auto2e_eri_hfgg.f90", 
"src/integrals/Auto2e/src/auto2e_eri_hfgp.f90", 
"src/integrals/Auto2e/src/auto2e_eri_hfgs.f90", 
"src/integrals/Auto2e/src/auto2e_eri_hfhd.f90", 
"src/integrals/Auto2e/src/auto2e_eri_hfhf.f90", 
"src/integrals/Auto2e/src/auto2e_eri_hfhp.f90", 
"src/integrals/Auto2e/src/auto2e_eri_hfhs.f90", 
"src/integrals/Auto2e/src/auto2e_eri_hfpp.f90", 
"src/integrals/Auto2e/src/auto2e_eri_hfps.f90", 
"src/integrals/Auto2e/src/auto2e_eri_hfss.f90", 
"src/integrals/Auto2e/src/auto2e_eri_hgdd.f90", 
"src/integrals/Auto2e/src/auto2e_eri_hgdp.f90", 
"src/integrals/Auto2e/src/auto2e_eri_hgds.f90", 
"src/integrals/Auto2e/src/auto2e_eri_hgfd.f90", 
"src/integrals/Auto2e/src/auto2e_eri_hgff.f90", 
"src/integrals/Auto2e/src/auto2e_eri_hgfp.f90", 
"src/integrals/Auto2e/src/auto2e_eri_hgfs.f90", 
"src/integrals/Auto2e/src/auto2e_eri_hggd.f90", 
"src/integrals/Auto2e/src/auto2e_eri_hggf.f90", 
"src/integrals/Auto2e/src/auto2e_eri_hggg.f90", 
"src/integrals/Auto2e/src/auto2e_eri_hggp.f90", 
"src/integrals/Auto2e/src/auto2e_eri_hggs.f90", 
"src/integrals/Auto2e/src/auto2e_eri_hghd.f90", 
"src/integrals/Auto2e/src/auto2e_eri_hghf.f90", 
"src/integrals/Auto2e/src/auto2e_eri_hghg.f90", 
"src/integrals/Auto2e/src/auto2e_eri_hghp.f90", 
"src/integrals/Auto2e/src/auto2e_eri_hghs.f90", 
"src/integrals/Auto2e/src/auto2e_eri_hgpp.f90", 
"src/integrals/Auto2e/src/auto2e_eri_hgps.f90", 
"src/integrals/Auto2e/src/auto2e_eri_hgss.f90", 
"src/integrals/Auto2e/src/auto2e_eri_hhdd.f90", 
"src/integrals/Auto2e/src/auto2e_eri_hhdp.f90", 
"src/integrals/Auto2e/src/auto2e_eri_hhds.f90", 
"src/integrals/Auto2e/src/auto2e_eri_hhfd.f90", 
"src/integrals/Auto2e/src/auto2e_eri_hhff.f90", 
"src/integrals/Auto2e/src/auto2e_eri_hhfp.f90", 
"src/integrals/Auto2e/src/auto2e_eri_hhfs.f90", 
"src/integrals/Auto2e/src/auto2e_eri_hhgd.f90", 
"src/integrals/Auto2e/src/auto2e_eri_hhgf.f90", 
"src/integrals/Auto2e/src/auto2e_eri_hhgg.f90", 
"src/integrals/Auto2e/src/auto2e_eri_hhgp.f90", 
"src/integrals/Auto2e/src/auto2e_eri_hhgs.f90", 
"src/integrals/Auto2e/src/auto2e_eri_hhhd.f90", 
"src/integrals/Auto2e/src/auto2e_eri_hhhf.f90", 
"src/integrals/Auto2e/src/auto2e_eri_hhhg.f90", 
"src/integrals/Auto2e/src/auto2e_eri_hhhh.f90", 
"src/integrals/Auto2e/src/auto2e_eri_hhhp.f90", 
"src/integrals/Auto2e/src/auto2e_eri_hhhs.f90", 
"src/integrals/Auto2e/src/auto2e_eri_hhpp.f90", 
"src/integrals/Auto2e/src/auto2e_eri_hhps.f90", 
"src/integrals/Auto2e/src/auto2e_eri_hhss.f90", 
"src/integrals/Auto2e/src/auto2e_eri_hpdd.f90", 
"src/integrals/Auto2e/src/auto2e_eri_hpdp.f90", 
"src/integrals/Auto2e/src/auto2e_eri_hpds.f90", 
"src/integrals/Auto2e/src/auto2e_eri_hpfd.f90", 
"src/integrals/Auto2e/src/auto2e_eri_hpff.f90", 
"src/integrals/Auto2e/src/auto2e_eri_hpfp.f90", 
"src/integrals/Auto2e/src/auto2e_eri_hpfs.f90", 
"src/integrals/Auto2e/src/auto2e_eri_hpgd.f90", 
"src/integrals/Auto2e/src/auto2e_eri_hpgf.f90", 
"src/integrals/Auto2e/src/auto2e_eri_hpgg.f90", 
"src/integrals/Auto2e/src/auto2e_eri_hpgp.f90", 
"src/integrals/Auto2e/src/auto2e_eri_hpgs.f90", 
"src/integrals/Auto2e/src/auto2e_eri_hphp.f90", 
"src/integrals/Auto2e/src/auto2e_eri_hphs.f90", 
"src/integrals/Auto2e/src/auto2e_eri_hppp.f90", 
"src/integrals/Auto2e/src/auto2e_eri_hpps.f90", 
"src/integrals/Auto2e/src/auto2e_eri_hpss.f90", 
"src/integrals/Auto2e/src/auto2e_eri_hsdd.f90", 
"src/integrals/Auto2e/src/auto2e_eri_hsdp.f90", 
"src/integrals/Auto2e/src/auto2e_eri_hsds.f90", 
"src/integrals/Auto2e/src/auto2e_eri_hsfd.f90", 
"src/integrals/Auto2e/src/auto2e_eri_hsff.f90", 
"src/integrals/Auto2e/src/auto2e_eri_hsfp.f90", 
"src/integrals/Auto2e/src/auto2e_eri_hsfs.f90", 
"src/integrals/Auto2e/src/auto2e_eri_hsgd.f90", 
"src/integrals/Auto2e/src/auto2e_eri_hsgf.f90", 
"src/integrals/Auto2e/src/auto2e_eri_hsgg.f90", 
"src/integrals/Auto2e/src/auto2e_eri_hsgp.f90", 
"src/integrals/Auto2e/src/auto2e_eri_hsgs.f90", 
"src/integrals/Auto2e/src/auto2e_eri_hshs.f90", 
"src/integrals/Auto2e/src/auto2e_eri_hspp.f90", 
"src/integrals/Auto2e/src/auto2e_eri_hsps.f90", 
"src/integrals/Auto2e/src/auto2e_eri_hsss.f90", 
"src/integrals/Auto2e/src/auto2e_eri_pppp.f90", 
"src/integrals/Auto2e/src/auto2e_eri_ppps.f90", 
"src/integrals/Auto2e/src/auto2e_eri_ppss.f90", 
"src/integrals/Auto2e/src/auto2e_eri_psps.f90", 
"src/integrals/Auto2e/src/auto2e_eri_psss.f90", 
"src/integrals/Auto2e/src/auto2e_eri_ssss.f90"])
BatchList.append(["src/integrals/Auto2e/src/auto2e.f90"])
BatchList.append(["src/integrals/sphergto.f90",
                  "src/common/periodic.f90",
                  "src/common/sort.f90"])
BatchList.append(["src/common/sys_definitions.f90",
                  "src/integrals/Cholesky/chol_definitions.f90",
                  "src/integrals/Cholesky/TwoStepCholesky_definitions.f90"])
BatchList.append(["src/integrals/basis_sets.f90"])
BatchList.append(["src/integrals/Cholesky/sorter_Cholesky.f90"])
BatchList.append(["src/Main/linalg.f90",
                  "src/linear-algebra/real_linalg.f90"])
BatchList.append(["src/integrals/Cholesky/TwoStepCholesky_Step1.f90",
                  "src/integrals/Cholesky/TwoStepCholesky_Step2.f90"])
BatchList.append(["src/integrals/Cholesky/TwoStepCholesky.f90"])
BatchList.append(["src/integrals/Cholesky/Cholesky.f90",
                  "src/integrals/Cholesky/CholeskyOTF.f90"])
BatchList.append(["src/integrals/OneElectronInts.f90"])
BatchList.append(["src/Main/CholeskyExchange.f90",
                  "src/Main/CholeskyCoulomb.f90"])
BatchList.append(["src/Main/CholeskyFock.f90"])
BatchList.append(["src/integrals/grid/GridFunctions.f90"])
BatchList.append(["src/integrals/grid/BeckeGrid.f90"])
BatchList.append(["src/integrals/THC/TensorHypercontraction.f90"])
BatchList.append(["src/gammcor-interface/Auto2eInterface.f90"])
BatchList.append(["src/gammcor-interface/THC_Gammcor.f90"])
BatchList.append(["src/Main/Cholesky_driver.f90"])
BatchList.append(["src/Main/test.f90"])
MainProgram = "src/Main/test.f90"


if args.clean:
    for batch in BatchList:
        for File in batch:
            remove_files(File)
    sys.exit()


ObjList = []
LibraryObjList = []
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
            BaseName, Ext = path.splitext(path.split(AbsSrcPath)[1])
            ObjFile = path.join(ObjDir, BaseName + ".o")
            ObjList.append(ObjFile)
            if RelSrcPath != MainProgram:
                LibraryObjList.append(ObjFile)
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
            move_files(outfiles[ncompiled], ObjDir, IncludeDir)
            #
            # When the camel case convention is used for module names,
            # the compiler may produce lower case names for the mod files
            #
            move_files(outfiles[ncompiled].lower(), ObjDir, IncludeDir)
            ncompiled += 1
    if ncompiled > 0:
        RecompBatch = True


# ############################################################################
#                                   Linker
# ############################################################################
LinkerName, LinkerOptions = LinkerCmd.split(maxsplit=1)
Command = LinkerName + " " + " ".join(ObjList) + " " + LinkerOptions + " " + ExeFile
Display = LinkerName + " " + "[object files in {}]".format(ObjDir) + " " + LinkerOptions + " " + ExeFile
print(Display)
p = Popen(Command, shell=True)
error_code = p.wait()
if error_code != 0:
    sys.exit(error_code)

# ############################################################################
#                               Static library
# ############################################################################
print("Creating a static library")
Display = ArchiveCmd + " " + StaticLibFile + " " + "[object files in {}]".format(ObjDir)
Command = ArchiveCmd + " " + StaticLibFile + " " + " ".join(LibraryObjList)
print(Display)
p = Popen(Command, shell=True)
error_code = p.wait()
if error_code != 0:
    sys.exit(error_code)

datetime_finish = datetime.datetime.utcnow()
td = datetime_finish - datetime_start
dateformat = "%a %b %d %H:%M:%S UTC %Y"
datestr_finish = datetime_finish.strftime(dateformat)
print("% Compilation completed at {}".format(datestr_finish))
print("% Employed {} concurrent processes".format(nproc))
print("% Total wall clock time [s]: {walltime:.1f}".format(walltime=td.total_seconds()))
