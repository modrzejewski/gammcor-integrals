#!/usr/bin/env python
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
BatchList.append(["src/auto2e_BraTransform.f90", 
"src/auto2e_Hermite.f90", 
"src/auto2e_KetTransfer.f90", 
"src/auto2e_KetTransform_10_1.f90", 
"src/auto2e_KetTransform_10_10.f90", 
"src/auto2e_KetTransform_10_2.f90", 
"src/auto2e_KetTransform_10_3.f90", 
"src/auto2e_KetTransform_10_4.f90", 
"src/auto2e_KetTransform_10_5.f90", 
"src/auto2e_KetTransform_10_6.f90", 
"src/auto2e_KetTransform_10_7.f90", 
"src/auto2e_KetTransform_10_8.f90", 
"src/auto2e_KetTransform_10_9.f90", 
"src/auto2e_KetTransform_11_1.f90", 
"src/auto2e_KetTransform_11_10.f90", 
"src/auto2e_KetTransform_11_11.f90", 
"src/auto2e_KetTransform_11_2.f90", 
"src/auto2e_KetTransform_11_3.f90", 
"src/auto2e_KetTransform_11_4.f90", 
"src/auto2e_KetTransform_11_5.f90", 
"src/auto2e_KetTransform_11_6.f90", 
"src/auto2e_KetTransform_11_7.f90", 
"src/auto2e_KetTransform_11_8.f90", 
"src/auto2e_KetTransform_11_9.f90", 
"src/auto2e_KetTransform_12_1.f90", 
"src/auto2e_KetTransform_12_10.f90", 
"src/auto2e_KetTransform_12_11.f90", 
"src/auto2e_KetTransform_12_12.f90", 
"src/auto2e_KetTransform_12_2.f90", 
"src/auto2e_KetTransform_12_3.f90", 
"src/auto2e_KetTransform_12_4.f90", 
"src/auto2e_KetTransform_12_5.f90", 
"src/auto2e_KetTransform_12_6.f90", 
"src/auto2e_KetTransform_12_7.f90", 
"src/auto2e_KetTransform_12_8.f90", 
"src/auto2e_KetTransform_12_9.f90", 
"src/auto2e_KetTransform_2_1.f90", 
"src/auto2e_KetTransform_2_2.f90", 
"src/auto2e_KetTransform_3_1.f90", 
"src/auto2e_KetTransform_3_2.f90", 
"src/auto2e_KetTransform_3_3.f90", 
"src/auto2e_KetTransform_4_1.f90", 
"src/auto2e_KetTransform_4_2.f90", 
"src/auto2e_KetTransform_4_3.f90", 
"src/auto2e_KetTransform_4_4.f90", 
"src/auto2e_KetTransform_5_1.f90", 
"src/auto2e_KetTransform_5_2.f90", 
"src/auto2e_KetTransform_5_3.f90", 
"src/auto2e_KetTransform_5_4.f90", 
"src/auto2e_KetTransform_5_5.f90", 
"src/auto2e_KetTransform_6_1.f90", 
"src/auto2e_KetTransform_6_2.f90", 
"src/auto2e_KetTransform_6_3.f90", 
"src/auto2e_KetTransform_6_4.f90", 
"src/auto2e_KetTransform_6_5_part1.f90", 
"src/auto2e_KetTransform_6_5_part2.f90", 
"src/auto2e_KetTransform_6_6_part1.f90", 
"src/auto2e_KetTransform_6_6_part2.f90", 
"src/auto2e_KetTransform_6_6_part3.f90", 
"src/auto2e_KetTransform_7_1.f90", 
"src/auto2e_KetTransform_7_2.f90", 
"src/auto2e_KetTransform_7_3.f90", 
"src/auto2e_KetTransform_7_4.f90", 
"src/auto2e_KetTransform_7_5.f90", 
"src/auto2e_KetTransform_7_6.f90", 
"src/auto2e_KetTransform_7_7.f90", 
"src/auto2e_KetTransform_8_1.f90", 
"src/auto2e_KetTransform_8_2.f90", 
"src/auto2e_KetTransform_8_3.f90", 
"src/auto2e_KetTransform_8_4.f90", 
"src/auto2e_KetTransform_8_5.f90", 
"src/auto2e_KetTransform_8_6.f90", 
"src/auto2e_KetTransform_8_7.f90", 
"src/auto2e_KetTransform_8_8.f90", 
"src/auto2e_KetTransform_9_1.f90", 
"src/auto2e_KetTransform_9_2.f90", 
"src/auto2e_KetTransform_9_3.f90", 
"src/auto2e_KetTransform_9_4.f90", 
"src/auto2e_KetTransform_9_5.f90", 
"src/auto2e_KetTransform_9_6.f90", 
"src/auto2e_KetTransform_9_7.f90", 
"src/auto2e_KetTransform_9_8.f90", 
"src/auto2e_KetTransform_9_9.f90", 
"src/auto2e_SpherTransf.f90"])
BatchList.append(["src/auto2e_WMatrix_part1.f90", 
"src/auto2e_WMatrix_part2.f90", 
"src/auto2e_WMatrix_part3.f90", 
"src/auto2e_WMatrix_part4.f90"])
BatchList.append(["src/auto2e_eri_dddd.f90", 
"src/auto2e_eri_dddp.f90", 
"src/auto2e_eri_ddds.f90", 
"src/auto2e_eri_ddpp.f90", 
"src/auto2e_eri_ddps.f90", 
"src/auto2e_eri_ddss.f90", 
"src/auto2e_eri_dpdp.f90", 
"src/auto2e_eri_dpds.f90", 
"src/auto2e_eri_dppp.f90", 
"src/auto2e_eri_dpps.f90", 
"src/auto2e_eri_dpss.f90", 
"src/auto2e_eri_dsds.f90", 
"src/auto2e_eri_dspp.f90", 
"src/auto2e_eri_dsps.f90", 
"src/auto2e_eri_dsss.f90", 
"src/auto2e_eri_fddd.f90", 
"src/auto2e_eri_fddp.f90", 
"src/auto2e_eri_fdds.f90", 
"src/auto2e_eri_fdfd.f90", 
"src/auto2e_eri_fdfp.f90", 
"src/auto2e_eri_fdfs.f90", 
"src/auto2e_eri_fdpp.f90", 
"src/auto2e_eri_fdps.f90", 
"src/auto2e_eri_fdss.f90", 
"src/auto2e_eri_ffdd.f90", 
"src/auto2e_eri_ffdp.f90", 
"src/auto2e_eri_ffds.f90", 
"src/auto2e_eri_fffd.f90", 
"src/auto2e_eri_ffff.f90", 
"src/auto2e_eri_fffp.f90", 
"src/auto2e_eri_fffs.f90", 
"src/auto2e_eri_ffpp.f90", 
"src/auto2e_eri_ffps.f90", 
"src/auto2e_eri_ffss.f90", 
"src/auto2e_eri_fpdd.f90", 
"src/auto2e_eri_fpdp.f90", 
"src/auto2e_eri_fpds.f90", 
"src/auto2e_eri_fpfp.f90", 
"src/auto2e_eri_fpfs.f90", 
"src/auto2e_eri_fppp.f90", 
"src/auto2e_eri_fpps.f90", 
"src/auto2e_eri_fpss.f90", 
"src/auto2e_eri_fsdd.f90", 
"src/auto2e_eri_fsdp.f90", 
"src/auto2e_eri_fsds.f90", 
"src/auto2e_eri_fsfs.f90", 
"src/auto2e_eri_fspp.f90", 
"src/auto2e_eri_fsps.f90", 
"src/auto2e_eri_fsss.f90", 
"src/auto2e_eri_gddd.f90", 
"src/auto2e_eri_gddp.f90", 
"src/auto2e_eri_gdds.f90", 
"src/auto2e_eri_gdfd.f90", 
"src/auto2e_eri_gdff.f90", 
"src/auto2e_eri_gdfp.f90", 
"src/auto2e_eri_gdfs.f90", 
"src/auto2e_eri_gdgd.f90", 
"src/auto2e_eri_gdgp.f90", 
"src/auto2e_eri_gdgs.f90", 
"src/auto2e_eri_gdpp.f90", 
"src/auto2e_eri_gdps.f90", 
"src/auto2e_eri_gdss.f90", 
"src/auto2e_eri_gfdd.f90", 
"src/auto2e_eri_gfdp.f90", 
"src/auto2e_eri_gfds.f90", 
"src/auto2e_eri_gffd.f90", 
"src/auto2e_eri_gfff.f90", 
"src/auto2e_eri_gffp.f90", 
"src/auto2e_eri_gffs.f90", 
"src/auto2e_eri_gfgd.f90", 
"src/auto2e_eri_gfgf.f90", 
"src/auto2e_eri_gfgp.f90", 
"src/auto2e_eri_gfgs.f90", 
"src/auto2e_eri_gfpp.f90", 
"src/auto2e_eri_gfps.f90", 
"src/auto2e_eri_gfss.f90", 
"src/auto2e_eri_ggdd.f90", 
"src/auto2e_eri_ggdp.f90", 
"src/auto2e_eri_ggds.f90", 
"src/auto2e_eri_ggfd.f90", 
"src/auto2e_eri_ggff.f90", 
"src/auto2e_eri_ggfp.f90", 
"src/auto2e_eri_ggfs.f90", 
"src/auto2e_eri_gggd.f90", 
"src/auto2e_eri_gggf.f90", 
"src/auto2e_eri_gggg.f90", 
"src/auto2e_eri_gggp.f90", 
"src/auto2e_eri_gggs.f90", 
"src/auto2e_eri_ggpp.f90", 
"src/auto2e_eri_ggps.f90", 
"src/auto2e_eri_ggss.f90", 
"src/auto2e_eri_gpdd.f90", 
"src/auto2e_eri_gpdp.f90", 
"src/auto2e_eri_gpds.f90", 
"src/auto2e_eri_gpfd.f90", 
"src/auto2e_eri_gpff.f90", 
"src/auto2e_eri_gpfp.f90", 
"src/auto2e_eri_gpfs.f90", 
"src/auto2e_eri_gpgp.f90", 
"src/auto2e_eri_gpgs.f90", 
"src/auto2e_eri_gppp.f90", 
"src/auto2e_eri_gpps.f90", 
"src/auto2e_eri_gpss.f90", 
"src/auto2e_eri_gsdd.f90", 
"src/auto2e_eri_gsdp.f90", 
"src/auto2e_eri_gsds.f90", 
"src/auto2e_eri_gsfd.f90", 
"src/auto2e_eri_gsff.f90", 
"src/auto2e_eri_gsfp.f90", 
"src/auto2e_eri_gsfs.f90", 
"src/auto2e_eri_gsgs.f90", 
"src/auto2e_eri_gspp.f90", 
"src/auto2e_eri_gsps.f90", 
"src/auto2e_eri_gsss.f90", 
"src/auto2e_eri_hddd.f90", 
"src/auto2e_eri_hddp.f90", 
"src/auto2e_eri_hdds.f90", 
"src/auto2e_eri_hdfd.f90", 
"src/auto2e_eri_hdff.f90", 
"src/auto2e_eri_hdfp.f90", 
"src/auto2e_eri_hdfs.f90", 
"src/auto2e_eri_hdgd.f90", 
"src/auto2e_eri_hdgf.f90", 
"src/auto2e_eri_hdgg.f90", 
"src/auto2e_eri_hdgp.f90", 
"src/auto2e_eri_hdgs.f90", 
"src/auto2e_eri_hdhd.f90", 
"src/auto2e_eri_hdhp.f90", 
"src/auto2e_eri_hdhs.f90", 
"src/auto2e_eri_hdpp.f90", 
"src/auto2e_eri_hdps.f90", 
"src/auto2e_eri_hdss.f90", 
"src/auto2e_eri_hfdd.f90", 
"src/auto2e_eri_hfdp.f90", 
"src/auto2e_eri_hfds.f90", 
"src/auto2e_eri_hffd.f90", 
"src/auto2e_eri_hfff.f90", 
"src/auto2e_eri_hffp.f90", 
"src/auto2e_eri_hffs.f90", 
"src/auto2e_eri_hfgd.f90", 
"src/auto2e_eri_hfgf.f90", 
"src/auto2e_eri_hfgg.f90", 
"src/auto2e_eri_hfgp.f90", 
"src/auto2e_eri_hfgs.f90", 
"src/auto2e_eri_hfhd.f90", 
"src/auto2e_eri_hfhf.f90", 
"src/auto2e_eri_hfhp.f90", 
"src/auto2e_eri_hfhs.f90", 
"src/auto2e_eri_hfpp.f90", 
"src/auto2e_eri_hfps.f90", 
"src/auto2e_eri_hfss.f90", 
"src/auto2e_eri_hgdd.f90", 
"src/auto2e_eri_hgdp.f90", 
"src/auto2e_eri_hgds.f90", 
"src/auto2e_eri_hgfd.f90", 
"src/auto2e_eri_hgff.f90", 
"src/auto2e_eri_hgfp.f90", 
"src/auto2e_eri_hgfs.f90", 
"src/auto2e_eri_hggd.f90", 
"src/auto2e_eri_hggf.f90", 
"src/auto2e_eri_hggg.f90", 
"src/auto2e_eri_hggp.f90", 
"src/auto2e_eri_hggs.f90", 
"src/auto2e_eri_hghd.f90", 
"src/auto2e_eri_hghf.f90", 
"src/auto2e_eri_hghg.f90", 
"src/auto2e_eri_hghp.f90", 
"src/auto2e_eri_hghs.f90", 
"src/auto2e_eri_hgpp.f90", 
"src/auto2e_eri_hgps.f90", 
"src/auto2e_eri_hgss.f90", 
"src/auto2e_eri_hhdd.f90", 
"src/auto2e_eri_hhdp.f90", 
"src/auto2e_eri_hhds.f90", 
"src/auto2e_eri_hhfd.f90", 
"src/auto2e_eri_hhff.f90", 
"src/auto2e_eri_hhfp.f90", 
"src/auto2e_eri_hhfs.f90", 
"src/auto2e_eri_hhgd.f90", 
"src/auto2e_eri_hhgf.f90", 
"src/auto2e_eri_hhgg.f90", 
"src/auto2e_eri_hhgp.f90", 
"src/auto2e_eri_hhgs.f90", 
"src/auto2e_eri_hhhd.f90", 
"src/auto2e_eri_hhhf.f90", 
"src/auto2e_eri_hhhg.f90", 
"src/auto2e_eri_hhhh.f90", 
"src/auto2e_eri_hhhp.f90", 
"src/auto2e_eri_hhhs.f90", 
"src/auto2e_eri_hhpp.f90", 
"src/auto2e_eri_hhps.f90", 
"src/auto2e_eri_hhss.f90", 
"src/auto2e_eri_hpdd.f90", 
"src/auto2e_eri_hpdp.f90", 
"src/auto2e_eri_hpds.f90", 
"src/auto2e_eri_hpfd.f90", 
"src/auto2e_eri_hpff.f90", 
"src/auto2e_eri_hpfp.f90", 
"src/auto2e_eri_hpfs.f90", 
"src/auto2e_eri_hpgd.f90", 
"src/auto2e_eri_hpgf.f90", 
"src/auto2e_eri_hpgg.f90", 
"src/auto2e_eri_hpgp.f90", 
"src/auto2e_eri_hpgs.f90", 
"src/auto2e_eri_hphp.f90", 
"src/auto2e_eri_hphs.f90", 
"src/auto2e_eri_hppp.f90", 
"src/auto2e_eri_hpps.f90", 
"src/auto2e_eri_hpss.f90", 
"src/auto2e_eri_hsdd.f90", 
"src/auto2e_eri_hsdp.f90", 
"src/auto2e_eri_hsds.f90", 
"src/auto2e_eri_hsfd.f90", 
"src/auto2e_eri_hsff.f90", 
"src/auto2e_eri_hsfp.f90", 
"src/auto2e_eri_hsfs.f90", 
"src/auto2e_eri_hsgd.f90", 
"src/auto2e_eri_hsgf.f90", 
"src/auto2e_eri_hsgg.f90", 
"src/auto2e_eri_hsgp.f90", 
"src/auto2e_eri_hsgs.f90", 
"src/auto2e_eri_hshs.f90", 
"src/auto2e_eri_hspp.f90", 
"src/auto2e_eri_hsps.f90", 
"src/auto2e_eri_hsss.f90", 
"src/auto2e_eri_iddd.f90", 
"src/auto2e_eri_iddp.f90", 
"src/auto2e_eri_idds.f90", 
"src/auto2e_eri_idfd.f90", 
"src/auto2e_eri_idff.f90", 
"src/auto2e_eri_idfp.f90", 
"src/auto2e_eri_idfs.f90", 
"src/auto2e_eri_idgd.f90", 
"src/auto2e_eri_idgf.f90", 
"src/auto2e_eri_idgg.f90", 
"src/auto2e_eri_idgp.f90", 
"src/auto2e_eri_idgs.f90", 
"src/auto2e_eri_idhd.f90", 
"src/auto2e_eri_idhf.f90", 
"src/auto2e_eri_idhg.f90", 
"src/auto2e_eri_idhh.f90", 
"src/auto2e_eri_idhp.f90", 
"src/auto2e_eri_idhs.f90", 
"src/auto2e_eri_idid.f90", 
"src/auto2e_eri_idip.f90", 
"src/auto2e_eri_idis.f90", 
"src/auto2e_eri_idpp.f90", 
"src/auto2e_eri_idps.f90", 
"src/auto2e_eri_idss.f90", 
"src/auto2e_eri_ifdd.f90", 
"src/auto2e_eri_ifdp.f90", 
"src/auto2e_eri_ifds.f90", 
"src/auto2e_eri_iffd.f90", 
"src/auto2e_eri_ifff.f90", 
"src/auto2e_eri_iffp.f90", 
"src/auto2e_eri_iffs.f90", 
"src/auto2e_eri_ifgd.f90", 
"src/auto2e_eri_ifgf.f90", 
"src/auto2e_eri_ifgg.f90", 
"src/auto2e_eri_ifgp.f90", 
"src/auto2e_eri_ifgs.f90", 
"src/auto2e_eri_ifhd.f90", 
"src/auto2e_eri_ifhf.f90", 
"src/auto2e_eri_ifhg.f90", 
"src/auto2e_eri_ifhh.f90", 
"src/auto2e_eri_ifhp.f90", 
"src/auto2e_eri_ifhs.f90", 
"src/auto2e_eri_ifid.f90", 
"src/auto2e_eri_ifif.f90", 
"src/auto2e_eri_ifip.f90", 
"src/auto2e_eri_ifis.f90", 
"src/auto2e_eri_ifpp.f90", 
"src/auto2e_eri_ifps.f90", 
"src/auto2e_eri_ifss.f90", 
"src/auto2e_eri_igdd.f90", 
"src/auto2e_eri_igdp.f90", 
"src/auto2e_eri_igds.f90", 
"src/auto2e_eri_igfd.f90", 
"src/auto2e_eri_igff.f90", 
"src/auto2e_eri_igfp.f90", 
"src/auto2e_eri_igfs.f90", 
"src/auto2e_eri_iggd.f90", 
"src/auto2e_eri_iggf.f90", 
"src/auto2e_eri_iggg.f90", 
"src/auto2e_eri_iggp.f90", 
"src/auto2e_eri_iggs.f90", 
"src/auto2e_eri_ighd.f90", 
"src/auto2e_eri_ighf.f90", 
"src/auto2e_eri_ighg.f90", 
"src/auto2e_eri_ighh.f90", 
"src/auto2e_eri_ighp.f90", 
"src/auto2e_eri_ighs.f90", 
"src/auto2e_eri_igid.f90", 
"src/auto2e_eri_igif.f90", 
"src/auto2e_eri_igig.f90", 
"src/auto2e_eri_igip.f90", 
"src/auto2e_eri_igis.f90", 
"src/auto2e_eri_igpp.f90", 
"src/auto2e_eri_igps.f90", 
"src/auto2e_eri_igss.f90", 
"src/auto2e_eri_ihdd.f90", 
"src/auto2e_eri_ihdp.f90", 
"src/auto2e_eri_ihds.f90", 
"src/auto2e_eri_ihfd.f90", 
"src/auto2e_eri_ihff.f90", 
"src/auto2e_eri_ihfp.f90", 
"src/auto2e_eri_ihfs.f90", 
"src/auto2e_eri_ihgd.f90", 
"src/auto2e_eri_ihgf.f90", 
"src/auto2e_eri_ihgg.f90", 
"src/auto2e_eri_ihgp.f90", 
"src/auto2e_eri_ihgs.f90", 
"src/auto2e_eri_ihhd.f90", 
"src/auto2e_eri_ihhf.f90", 
"src/auto2e_eri_ihhg.f90", 
"src/auto2e_eri_ihhh.f90", 
"src/auto2e_eri_ihhp.f90", 
"src/auto2e_eri_ihhs.f90", 
"src/auto2e_eri_ihid.f90", 
"src/auto2e_eri_ihif.f90", 
"src/auto2e_eri_ihig.f90", 
"src/auto2e_eri_ihih.f90", 
"src/auto2e_eri_ihip.f90", 
"src/auto2e_eri_ihis.f90", 
"src/auto2e_eri_ihpp.f90", 
"src/auto2e_eri_ihps.f90", 
"src/auto2e_eri_ihss.f90", 
"src/auto2e_eri_iidd.f90", 
"src/auto2e_eri_iidp.f90", 
"src/auto2e_eri_iids.f90", 
"src/auto2e_eri_iifd.f90", 
"src/auto2e_eri_iiff.f90", 
"src/auto2e_eri_iifp.f90", 
"src/auto2e_eri_iifs.f90", 
"src/auto2e_eri_iigd.f90", 
"src/auto2e_eri_iigf.f90", 
"src/auto2e_eri_iigg.f90", 
"src/auto2e_eri_iigp.f90", 
"src/auto2e_eri_iigs.f90", 
"src/auto2e_eri_iihd.f90", 
"src/auto2e_eri_iihf.f90", 
"src/auto2e_eri_iihg.f90", 
"src/auto2e_eri_iihh.f90", 
"src/auto2e_eri_iihp.f90", 
"src/auto2e_eri_iihs.f90", 
"src/auto2e_eri_iiid.f90", 
"src/auto2e_eri_iiif.f90", 
"src/auto2e_eri_iiig.f90", 
"src/auto2e_eri_iiih.f90", 
"src/auto2e_eri_iiii.f90", 
"src/auto2e_eri_iiip.f90", 
"src/auto2e_eri_iiis.f90", 
"src/auto2e_eri_iipp.f90", 
"src/auto2e_eri_iips.f90", 
"src/auto2e_eri_iiss.f90", 
"src/auto2e_eri_ipdd.f90", 
"src/auto2e_eri_ipdp.f90", 
"src/auto2e_eri_ipds.f90", 
"src/auto2e_eri_ipfd.f90", 
"src/auto2e_eri_ipff.f90", 
"src/auto2e_eri_ipfp.f90", 
"src/auto2e_eri_ipfs.f90", 
"src/auto2e_eri_ipgd.f90", 
"src/auto2e_eri_ipgf.f90", 
"src/auto2e_eri_ipgg.f90", 
"src/auto2e_eri_ipgp.f90", 
"src/auto2e_eri_ipgs.f90", 
"src/auto2e_eri_iphd.f90", 
"src/auto2e_eri_iphf.f90", 
"src/auto2e_eri_iphg.f90", 
"src/auto2e_eri_iphh.f90", 
"src/auto2e_eri_iphp.f90", 
"src/auto2e_eri_iphs.f90", 
"src/auto2e_eri_ipip.f90", 
"src/auto2e_eri_ipis.f90", 
"src/auto2e_eri_ippp.f90", 
"src/auto2e_eri_ipps.f90", 
"src/auto2e_eri_ipss.f90", 
"src/auto2e_eri_isdd.f90", 
"src/auto2e_eri_isdp.f90", 
"src/auto2e_eri_isds.f90", 
"src/auto2e_eri_isfd.f90", 
"src/auto2e_eri_isff.f90", 
"src/auto2e_eri_isfp.f90", 
"src/auto2e_eri_isfs.f90", 
"src/auto2e_eri_isgd.f90", 
"src/auto2e_eri_isgf.f90", 
"src/auto2e_eri_isgg.f90", 
"src/auto2e_eri_isgp.f90", 
"src/auto2e_eri_isgs.f90", 
"src/auto2e_eri_ishd.f90", 
"src/auto2e_eri_ishf.f90", 
"src/auto2e_eri_ishg.f90", 
"src/auto2e_eri_ishh.f90", 
"src/auto2e_eri_ishp.f90", 
"src/auto2e_eri_ishs.f90", 
"src/auto2e_eri_isis.f90", 
"src/auto2e_eri_ispp.f90", 
"src/auto2e_eri_isps.f90", 
"src/auto2e_eri_isss.f90", 
"src/auto2e_eri_pppp.f90", 
"src/auto2e_eri_ppps.f90", 
"src/auto2e_eri_ppss.f90", 
"src/auto2e_eri_psps.f90", 
"src/auto2e_eri_psss.f90", 
"src/auto2e_eri_ssss.f90"])
BatchList.append(["src/auto2e.f90"])


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
print("% Compilation completed at {}".format(datestr_finish))
print("% Employed {} concurrent processes".format(nproc))
print("% Total wall clock time [s]: {walltime:.1f}".format(walltime=td.total_seconds()))
