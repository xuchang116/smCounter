#!/usr/bin/python
# vim: tabstop=9 expandtab shiftwidth=3 softtabstop=3
# Downsample MTs
# Chang Xu. 07DEC2015

import os
import sys
import datetime
import subprocess
import shutil
import time
import math
import pysam
from collections import defaultdict
import argparse
import random

bamHeader = "/qgen/home/xuc/VariantCallingPaper/N0030/bamHeader.txt"

#--------------------------------------------------------------------------------------
# main function
#--------------------------------------------------------------------------------------
def main(args):
   random.seed(args.seed)
   # change working directory to runDir and make output directories 
   os.chdir(args.runPath)

   timeStart = datetime.datetime.now()
   #print("started at " + str(timeStart))

   bcDict = defaultdict(list)
   samfile = pysam.AlignmentFile(args.inBam, "rb")
   dsBam = pysam.AlignmentFile(args.outBam, "wb", template=samfile) # this does not work sometimes due to BAM header error
   #dsBam = pysam.AlignmentFile(args.outBam, "wb", text=bamHeader) # error


   selectedMTs = set()

   print('reading in bam file')
   for read in samfile.fetch():
      # read ID
      qname = read.query_name
      # Barcode
      BC = qname.strip().split(':')[-2]
      if qname not in bcDict[BC]:
         bcDict[BC].append(qname)

   # select or drop MTs 
   #print('Random sampling')
   for bc in bcDict.keys():
      r = random.random()
      if r <= args.pct:
         selectedMTs.add(bc)

   # write to file
   print('Writing to BAM file')
   for read in samfile.fetch():
      # read ID
      qname = read.query_name
      # Barcode
      BC = qname.strip().split(':')[-2]
      if BC in selectedMTs:  
         dsBam.write(read)

   samfile.close()
   dsBam.close()

   # log run completion
   timeEnd = datetime.datetime.now()
   #print("completed running at " + str(timeEnd) + "\n")
   #print("total time: "+ str(timeEnd-timeStart) + "\n")   


#--------------------------------------------------------------------------------------
# Run 
#--------------------------------------------------------------------------------------
def run():
   parser = argparse.ArgumentParser(description='Downsample MTs')
   parser.add_argument('--runPath', default=None, help='path to working directory')
   parser.add_argument('--inBam', default=None, help='Input BAM file')
   parser.add_argument('--outBam', default=None, help='Output BAM file')
   parser.add_argument('--pct', type=float, default=0.5, help='Percent of MTs kept')
   parser.add_argument('--seed', type=int, default=1234567, help='Seed for random number generation')
   args = parser.parse_args()
   main(args)

#----------------------------------------------------------------------------------------------
#pythonism to run from the command line
#----------------------------------------------------------------------------------------------
if __name__ == "__main__":
    run()





