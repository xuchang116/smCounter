#!/usr/bin/python
# vim: tabstop=9 expandtab shiftwidth=3 softtabstop=3
# Downsample reads, while maintaining MT count
# Chang Xu. 23NOV2015

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
   dsBam = pysam.AlignmentFile(args.outBam, "wb", template=samfile)

   selectedReads = set()
   droppedReads = set()

   #print('reading in bam file')
   for read in samfile.fetch():
      # read ID
      qname = read.query_name
      # Barcode
      BC = qname.strip().split(':')[-2]
      if qname not in bcDict[BC]:
         bcDict[BC].append(qname)

   # Count the number of one-read MT and multi-read MT
   #print('Counting MTs')
   oneReadMtCnt = 0
   multiReadMtCnt = 0
   multiReadMtReadCnt = 0
   for bc in bcDict.values():
      if len(bc) == 1:
         oneReadMtCnt += 1
      else:
         multiReadMtCnt += 1
         multiReadMtReadCnt += len(bc)

   # calculate probability to keep reads in multi-read MT, after keeping 1 read
   probKeep = 1.0 * (args.rpb - 1.0) * (oneReadMtCnt + multiReadMtCnt) / (multiReadMtReadCnt - multiReadMtCnt)

   # select or drop the reads
   #print('Random sampling')
   for bc in bcDict.values():
      # always keep one read MT
      if len(bc) == 1:
         selectedReads.add(bc[0])
      # keep the first read in multi-read MT, then random sample
      else:
         selectedReads.add(bc[0])
         for rid in bc[1:]:
            r = random.random()
            if r <= probKeep:
               selectedReads.add(rid)

   # write to file
   #print('Writing to BAM file')
   for read in samfile.fetch():
      # read ID
      qname = read.query_name
      if qname in selectedReads:  
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
   parser.add_argument('--rpb', type=float, default=1.0, help='target reads per MT')
   parser.add_argument('--seed', type=int, default=1234567, help='Seed for random number generation')
   args = parser.parse_args()
   main(args)

#----------------------------------------------------------------------------------------------
#pythonism to run from the command line
#----------------------------------------------------------------------------------------------
if __name__ == "__main__":
    run()





