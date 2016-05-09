#!/usr/bin/python
# vim: tabstop=8 expandtab shiftwidth=4 softtabstop=4
import datetime
import gzip
import os
import shutil
import subprocess
import sys
import time
import pybedtools

PD_DIST = 10
allelicPrimExe = "/qgen/home/xuc/software/vcflib/bin/vcfallelicprimitives"
vcfComparatorCmd = "java -jar -Xmx22g /mnt/fdkbio05/rvijaya/software/USeq_8.7.6/Apps/VCFComparator"
#nistGIBVcfFile = "/qgen/home/rvijaya/misc/NA12878_NIST/NISTIntegratedCalls_14datasets_131103_allcall_UGHapMerge_HetHomVarPASS_VQSRv2.18_all.primitives.vcf"
#nistGIBBedFile = "/qgen/home/rvijaya/misc/NA12878_NIST/union13callableMQonlymerged_addcert_nouncert_excludesimplerep_excludesegdups_excludedecoy_excludeRepSeqSTRs_noCNVs_v2.18_2mindatasets_5minYesNoRatio.bed"
nistGIBVcfFile = "/qgen/home/xuc/GIAB/NA12878/NISTIntegratedCalls_14datasets_131103_allcall_UGHapMerge_HetHomVarPASS_VQSRv2.19_2mindatasets_5minYesNoRatio_all_nouncert_excludesimplerep_excludesegdups_excludedecoy_excludeRepSeqSTRs_noCNVs.vcf"
nistGIBBedFile = "/qgen/home/xuc/GIAB/NA12878/union13callableMQonlymerged_addcert_nouncert_excludesimplerep_excludesegdups_excludedecoy_excludeRepSeqSTRs_noCNVs_v2.19_2mindatasets_5minYesNoRatio.bed"

#-------------------------------------------------------------------------------------
# function to call subprocess.check_call
#-------------------------------------------------------------------------------------
def runAndTimeCmd(cmdToRun):
   timeStart = datetime.datetime.now()
   #print("Running cmd:" + cmdToRun)
   subprocess.check_call(cmdToRun, shell=True)
   timeEnd = datetime.datetime.now()
   #print("Time taken: " + str(timeEnd-timeStart))
   
#------------------------------------------------------------------------------------
# filter function to filter features present in hash. customized for mutect output. 
#-----------------------------------------------------------------------------------
def presentFilter_1(feature, vcfDict):
   if (feature[0], feature[1], feature[3]) in vcfDict:
      return False
   else:
      return True

#------------------------------------------------------------------------------------
# filter function to filter features present in hash. customized for variant calling output. 
#-----------------------------------------------------------------------------------
def presentFilter_2(feature, vcfDict):
   if (feature[0], feature[1], feature[2]) in vcfDict:
      return False
   else:
      return True

#------------------------------------------------------------------------------------
# subtract vcfB calls from vcfA calls: subtractBed is not able to get this right, customized for mutect output.
#------------------------------------------------------------------------------------
def subtractCalls_1(vcfA, vcfB, fileName):
   bCalls = {}
   # load calls from vcfB
   for feature in vcfB:
      bCalls[(feature[0], feature[1], feature[3])] = feature[4]
   outVcf = vcfA.filter(presentFilter_1, vcfDict=bCalls)
   outVcf.saveas(fileName)

#------------------------------------------------------------------------------------
# subtract vcfB calls from vcfA calls: subtractBed is not able to get this right, customized for variant calling output.
#------------------------------------------------------------------------------------
def subtractCalls_2(vcfA, vcfB, fileName):
   bCalls = {}
   # load calls from vcfB
   for feature in vcfB:
      bCalls[(feature[0], feature[1], feature[3])] = feature[4]
   outVcf = vcfA.filter(presentFilter_2, vcfDict=bCalls)
   outVcf.saveas(fileName)

#-------------------------------------------------------------------------------------
# limit candidate variants to panel and nist bed region
#-------------------------------------------------------------------------------------
def cleanVcf_2(runPath, bedName, vcfName, RunOutVcf, sharedRegion, dilutionVcfFile = None):
   os.chdir(runPath)
   sharedBed = pybedtools.BedTool(sharedRegion)

   # convert vcf to allelic primitives
   rootName = os.path.splitext(vcfName)[0]
   extName = os.path.splitext(vcfName)[1]
   allelePVcfName = rootName + '.primitive' + extName
   runAndTimeCmd(allelicPrimExe + " " + vcfName + " > " + allelePVcfName)

   # create vcf with no chr
   allelePVcfWOChr = 'nochr.vcf'
   runAndTimeCmd("cat " + allelePVcfName + " | sed 's/chr//' > " + allelePVcfWOChr)

   # load vcfs and limit them to the shared regions
   runVcf = pybedtools.BedTool(allelePVcfWOChr).intersect(sharedBed, f=1.0, u=True)

   # if dilution vcf is provided, subtract those mutations out. Note: dilution vcf must have been intersected with the shared region BED. 
   if dilutionVcfFile != None:
      dilutionVcf =  pybedtools.BedTool(dilutionVcfFile)
      dilutionCalls = len(dilutionVcf)
      subtractCalls_1(runVcf, dilutionVcf, RunOutVcf)
   else:
      runVcf.saveas(RunOutVcf)

   os.remove('nochr.vcf')

#-------------------------------------------------------------------------------------
# main program for running from shell
#-------------------------------------------------------------------------------------
if __name__ == "__main__":
   runPath = sys.argv[1]
   bedName = sys.argv[2]
   vcfName = sys.argv[3]
   RunOutVcf = sys.argv[4]
   sharedRegion = sys.argv[5]

   if len(sys.argv) > 6 and sys.argv[6] != "None":
      dilutionVcfFile = sys.argv[6]
   else:
      dilutionVcfFile = None
   cleanVcf_2(runPath, bedName, vcfName, RunOutVcf, sharedRegion, dilutionVcfFile)
