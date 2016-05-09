#!/usr/bin/python
# vim: tabstop=9 expandtab shiftwidth=3 softtabstop=3
# Updates from v6.6.1: 
# correct the bias against INDEL in model
# added simple repeat, low complexity, satellite region filter based on repeat masker definition
# relaxed mismatch filter
# modified strand bias filter -- only for SNPs with AF lower than 50%
# Chang Xu. 12Mar2016

import os
import sys
import datetime
import subprocess
import shutil
import time
import gzip
import runLog
import re
import math
import operator
import pysam
import multiprocessing as mp
from collections import defaultdict
import argparse
import scipy.stats
import random

refg = '/qgen/home/rvijaya/downloads/alt_hap_masked_ref/ucsc.hg19.fasta'
repBed = '/qgen/home/xuc/UCSC/simpleRepeat.bed'
srBed = '/qgen/home/xuc/UCSC/SR_LC_SL.nochr.bed'
pcr_error = 1e-6
pcr_no_error = 1.0 - 3e-5
atgc = ['A', 'T', 'G', 'C']
#cutoff = [40, 45, 50, 55, 60, 65, 70, 75, 80, 85, 90, 95, 100, 110, 120, 130, 140, 150, 160, 170] 
cutoff = [10, 15, 20, 25, 30, 35, 40, 45, 50, 55, 60, 65, 70, 75, 80, 85, 90, 95, 100, 125] 

#-------------------------------------------------------------------------------------
# function to calculate posterior probability for each barcode. 
#-------------------------------------------------------------------------------------
def calProb(oneBC, MTdrop):
    outDict = defaultdict(float)
    if len(oneBC) <= MTdrop:
        outDict['A'] = 0.0
        outDict['T'] = 0.0
        outDict['G'] = 0.0
        outDict['C'] = 0.0
    else:
        prodP = defaultdict(float)
        cnt = defaultdict(int)
        tmpOut = defaultdict(float)
        rightP = 1.0
        sumP = 0.0
        pcrP = defaultdict(float)
        
        # set ATGC count = 0
        for char in atgc:
           cnt[char] = 0

        # get unique bases. Make sure uniqBaseList contains 4 members, unless the barcode already contains more than or equal to 4 bases/indels
        # NOTE: existBase contains only the alleles, including indels, with at least 1 read in the MT. uniqBase may contain more. 
        existBase = set([info[0][0] for info in oneBC.values()])
        uniqBase = set([info[0][0] for info in oneBC.values()])
        if len(uniqBase) < 4:
           for b in atgc:
              if b not in uniqBase:
                 uniqBase.add(b)
                 if len(uniqBase) == 4:
                    break

        uniqBaseList = list(uniqBase)

        # set initial value in prodP to be 1.0
        for b in uniqBaseList:
           prodP[b] = 1.0

        for info in oneBC.values():
           base = info[0][0]
           prob = info[0][1]
           pairOrder = info[0][2]
           if pairOrder != 'Paired':
               prob = 0.1
           prodP[base] *= 1.0 - prob
           cnt[base] += 1

           for char in list(uniqBase - set([base])):
               prodP[char] *= prob

           rightP *= 1.0 - prob 
            
        for char in uniqBaseList:
           ratio = (cnt[char] + 0.5) / (len(oneBC) + 0.5 * len(uniqBaseList))
           pcrP[char] = 10.0 ** (-6.0 * ratio)

        for key in prodP.keys():
            if key in existBase:
               tmpOut[key] = pcr_no_error * prodP[key] + rightP * min([pcrP[char] for char in pcrP.keys() if char != key])
            else:
               tmpOut[key] = rightP
               for char in existBase:
                  if char != key:
                     tmpOut[key] *= pcrP[char] 

            sumP += tmpOut[key]
        
        for key in prodP.keys():
            outDict[key] = tmpOut[key] / sumP
        
    return outDict

#-------------------------------------------------------------------------------------
# check if a locus is within or flanked by homopolymer region and/or low complexity region
#-------------------------------------------------------------------------------------
def isHPorLowComp(chrom, pos, length, refb, altb):
   # get reference base
   refs = pysam.FastaFile(refg)
   # ref sequence of [pos-length, pos+length] interval
   Lseq = refs.fetch(reference=chrom, start=int(pos)-1-length, end=int(pos)-1).upper()
   Rseq_ref = refs.fetch(reference=chrom, start=int(pos)-1+len(refb), end=int(pos)-1+len(refb)+length).upper()
   refSeq = Lseq + refb + Rseq_ref
   # alt sequence
   Rseq_alt = refs.fetch(reference=chrom, start=int(pos)-1+len(altb), end=int(pos)-1+len(altb)+length).upper()
   altSeq = Lseq + altb + Rseq_alt
   # check homopolymer
   homoA = True if refSeq.find('A'*length) >= 0 or altSeq.find('A'*length) >= 0 else False
   homoT = True if refSeq.find('T'*length) >= 0 or altSeq.find('T'*length) >= 0 else False
   homoG = True if refSeq.find('G'*length) >= 0 or altSeq.find('G'*length) >= 0 else False
   homoC = True if refSeq.find('C'*length) >= 0 or altSeq.find('C'*length) >= 0 else False
   homop = True if homoA or homoT or homoG or homoC else False

   # check low complexity -- window length is 2 * homopolymer region. If any 2 nucleotide >= 99% 
   len2 = int(2 * length)
   LseqLC = refs.fetch(reference=chrom, start=int(pos)-1-len2, end=int(pos)-1).upper()
   # ref seq
   Rseq_refLC = refs.fetch(reference=chrom, start=int(pos)-1+len(refb), end=int(pos)-1+len(refb)+len2).upper()
   refSeqLC = LseqLC + refb + Rseq_refLC
   # alt seq
   Rseq_altLC = refs.fetch(reference=chrom, start=int(pos)-1+len(altb), end=int(pos)-1+len(altb)+len2).upper()
   altSeqLC = LseqLC + altb + Rseq_altLC

   lowcomp = False

   # Ref seq
   totalLen = len(refSeqLC)
   for i in range(totalLen-len2):
      subseq = refSeqLC[i:(i+len2)]
      countA = subseq.count('A')
      countT = subseq.count('T')
      countG = subseq.count('G')
      countC = subseq.count('C')
      sortedCounts = sorted([countA, countT, countG, countC], reverse=True)
      top2Freq = 1.0 * (sortedCounts[0] + sortedCounts[1]) / len2
      if top2Freq >= 0.99:
         lowcomp = True
         break
      
   # If ref seq is not LC, check alt seq
   if not lowcomp:
      totalLen = len(altSeqLC)
      for i in range(totalLen-len2):
         subseq = altSeqLC[i:(i+len2)]
         countA = subseq.count('A')
         countT = subseq.count('T')
         countG = subseq.count('G')
         countC = subseq.count('C')
         sortedCounts = sorted([countA, countT, countG, countC], reverse=True)
         top2Freq = 1.0 * (sortedCounts[0] + sortedCounts[1]) / len2
         if top2Freq >= 0.99:
            lowcomp = True
            break

   return [homop, lowcomp]

#-------------------------------------------------------------------------------------
# function to call variants
#-------------------------------------------------------------------------------------
def vc(bamName, chrom, pos, minBQ, minMQ, smt, length, mismatchThr, MTdrop, ds):
   samfile = pysam.AlignmentFile(bamName, 'rb')
   idx = 0
   cvg = 0
   bcDict = defaultdict(lambda: defaultdict(list)) 
   allBcDict = defaultdict(list) 
   alleleCnt = defaultdict(int)
   MTCnt = defaultdict(int)
   r1EndPos = defaultdict(list)
   r2EndPos = defaultdict(list)
   MT3Cnt = 0
   MT5Cnt = 0
   MT7Cnt = 0
   MT10Cnt = 0
   strongMTCnt = defaultdict(int)
   weakMTCnt = defaultdict(int)
   predIndex = defaultdict(lambda: defaultdict(float))
   finalDict = defaultdict(float)
   r1Cnt = defaultdict(int)
   r2Cnt = defaultdict(int)
   forwardCnt = defaultdict(int)
   reverseCnt = defaultdict(int)
   concordPairCnt = defaultdict(int)
   discordPairCnt = defaultdict(int)
   mismatchCnt = defaultdict(float)
   bqSum = defaultdict(int)

   # get reference base
   refseq = pysam.FastaFile(refg)
   origRef = refseq.fetch(reference=chrom, start=int(pos)-1, end=int(pos))
   origRef = origRef.upper()

   # pile up reads
   for read in samfile.pileup(region = chrom + ':' + pos + ':' + pos, truncate=True, max_depth=1000000, stepper='nofilter'):
      for pileupRead in read.pileups:
         # read ID
         qname = pileupRead.alignment.query_name
         readid = ':'.join(qname.strip().split(':')[:-2])
         # barcode sequence
         BC = qname.strip().split(':')[-2]
         # mapping quality
         mq = pileupRead.alignment.mapping_quality
         # get NM tag 
         NM = 0 
         allTags = pileupRead.alignment.tags
         for (tag, value) in allTags:
            if tag == 'NM':
               NM = value
               break
         # count number of INDELs in the read sequence
         nIndel = 0
         cigar = pileupRead.alignment.cigar
         cigarOrder = 1
         leftSP = 0  # soft clipped bases on the left
         rightSP = 0  # soft clipped bases on the right
         for (op, value) in cigar:
            # 1 for insertion
            if op == 1 or op == 2:
               nIndel += value 
            if cigarOrder == 1 and op == 4:
               leftSP = value
            if cigarOrder > 1 and op == 4:
               rightSP += value
            cigarOrder += 1

         # Number of mismatches except INDEL, including softcilpped sequences 
         mismatch = max(0, NM - nIndel)
         # read length, including softclip
         readLen = pileupRead.alignment.query_length
         # calculate mismatch per 100 bases
         mismatchPer100b = 100.0 * mismatch / readLen if readLen > 0 else 0.0

         # paired read
         if pileupRead.alignment.is_read1:
            pairOrder = 'R1'
         if pileupRead.alignment.is_read2:
            pairOrder = 'R2'

         # +/- strand
         strand = 'Reverse' if pileupRead.alignment.is_reverse else 'Forward'

         # coverage -- read, not fragment
         cvg += 1

         # check if the site is the beginning of insertion
         if pileupRead.indel > 0:
            site = pileupRead.alignment.query_sequence[pileupRead.query_position]
            inserted = pileupRead.alignment.query_sequence[(pileupRead.query_position + 1) : (pileupRead.query_position + 1 +  pileupRead.indel)]
            base = 'INS|' + site + '|' + site + inserted
            bq = pileupRead.alignment.query_qualities[pileupRead.query_position]
            bqSum[base] += bq
            # inclusion condition
            incCond = bq >= minBQ and mq >= minMQ and mismatchPer100b <= mismatchThr
            alleleCnt[base] += 1
            mismatchCnt[base] += mismatchPer100b
            if pairOrder == 'R1':
               r1Cnt[base] += 1
            if pairOrder == 'R2':
               r2Cnt[base] += 1

            if strand == 'Reverse':
               reverseCnt[base] += 1
            else:
               forwardCnt[base] += 1
            
         # check if the site is the beginning of deletion
         elif pileupRead.indel < 0:
            site = pileupRead.alignment.query_sequence[pileupRead.query_position]
            deleted = refseq.fetch(reference=chrom, start=int(pos), end=int(pos)+abs(pileupRead.indel))
            deleted = deleted.upper()
            base = 'DEL|' + site + deleted + '|' + site
            bq = pileupRead.alignment.query_qualities[pileupRead.query_position]
            bqSum[base] += bq
            # inclusion condition
            incCond = bq >= minBQ and mq >= minMQ and mismatchPer100b <= mismatchThr
            alleleCnt[base] += 1
            mismatchCnt[base] += mismatchPer100b
            if pairOrder == 'R1':
               r1Cnt[base] += 1
            if pairOrder == 'R2':
               r2Cnt[base] += 1

            if strand == 'Reverse':
               reverseCnt[base] += 1
            else:
               forwardCnt[base] += 1

            
         # site is not beginning of any INDEL
         else:
            # If the site ifself is a deletion, set quality = minBQ 
            if pileupRead.is_del:
               base = 'DEL'
               bq = minBQ
               bqSum[base] += bq
               # inclusion condition
               incCond = bq >= minBQ and mq >= minMQ and mismatchPer100b <= mismatchThr
            # if the site is a SNP, 
            else: 
               base = pileupRead.alignment.query_sequence[pileupRead.query_position] # note: query_sequence includes soft clipped bases
               bq = pileupRead.alignment.query_qualities[pileupRead.query_position]
               bqSum[base] += bq
               # inclusion condition
               incCond = bq >= minBQ and mq >= minMQ and mismatchPer100b <= mismatchThr
               if pairOrder == 'R1':
                  # distance to the starting end in R1; 
                  if pileupRead.alignment.is_reverse:
                     distToEnd = pileupRead.alignment.query_alignment_length - (pileupRead.query_position - leftSP)
                  else:
                     distToEnd = pileupRead.query_position - leftSP
                  if incCond:
                     r1EndPos[base].append(distToEnd)
                  r1Cnt[base] += 1
               if pairOrder == 'R2':
                  # distance to the random end in R2. Different cases for forward and reverse strand
                  if pileupRead.alignment.is_reverse:
                     distToEnd = pileupRead.query_position - leftSP
                  else:
                     distToEnd = pileupRead.alignment.query_alignment_length - (pileupRead.query_position - leftSP)
                  if incCond:
                     r2EndPos[base].append(distToEnd)
                  r2Cnt[base] += 1

               if strand == 'Reverse':
                  reverseCnt[base] += 1
               else:
                  forwardCnt[base] += 1

            alleleCnt[base] += 1
            mismatchCnt[base] += mismatchPer100b

         # count total number of fragments and MTs
         if readid not in allBcDict[BC]:
            allBcDict[BC].append(readid)

         # decide which read goes into analysis
         if incCond:
            if readid not in bcDict[BC]:
                prob = pow(10.0, -bq / 10.0)
                readinfo = [base, prob, pairOrder]
                bcDict[BC][readid].append(readinfo)
            elif base == bcDict[BC][readid][0][0] or base in ['N', '*']:
                bcDict[BC][readid][0][1] = max((pow(10.0, -bq / 10.0) , bcDict[BC][readid][0][1]))
                bcDict[BC][readid][0][2] = 'Paired'
                if base == bcDict[BC][readid][0][0]:
                   concordPairCnt[base] += 1
            else:
                del bcDict[BC][readid]
                discordPairCnt[base] += 1
         
   # total number of MT, fragments, reads, including those dropped from analysis
   allMT = len(allBcDict)
   allFrag = sum([len(allBcDict[bc]) for bc in allBcDict])

   # downsampling MTs (not dropped) to args.ds
   usedMT = min(ds, len(bcDict))
   if len(bcDict) > ds:
      bcKeys = random.sample(bcDict.keys(), ds)
   else:
      bcKeys = bcDict.keys()
   usedFrag = sum([len(bcDict[bc]) for bc in bcKeys])

   totalR1 = sum(r1Cnt.values())
   totalR2 = sum(r2Cnt.values())
   if usedMT == 0:
      out_long = '\t'.join([chrom, pos, origRef] + ['']*(46+len(cutoff)) + ['Zero_Coverage']) + '\n'
   else:
      for bc in bcKeys:
         bcProb = calProb(bcDict[bc], MTdrop)
         for char in bcProb.keys():
             log10P = -math.log10(1.0-bcProb[char])
             predIndex[bc][char] = log10P
             finalDict[char] += log10P

         max_base = [x for x in predIndex[bc].keys() if  predIndex[bc][x] == max(predIndex[bc].values())]
         if len(max_base) == 1:
             cons = max_base[0]
             MTCnt[cons] += 1
             if predIndex[bc][cons] > smt:
                 strongMTCnt[cons] += 1           
             if bcProb[cons] > 0.3 and bcProb[cons] < 0.85:
                 weakMTCnt[cons] += 1           
         # Tie in max predIndex is most likely due to single read MT. 
         elif len(bcDict[bc]) == 1:
             cons = bcDict[bc].values()[0][0][0]
             MTCnt[cons] += 1

         if len(bcDict[bc]) >= 3:
             MT3Cnt += 1
         if len(bcDict[bc]) >= 5:
             MT5Cnt += 1
         if len(bcDict[bc]) >= 7:
             MT7Cnt += 1
         if len(bcDict[bc]) >= 10:
             MT10Cnt += 1

      sortedList = sorted(finalDict.items(), key=operator.itemgetter(1), reverse=True)
      maxBase = sortedList[0][0]
      maxPI = sortedList[0][1]
      secondMaxBase = sortedList[1][0]
      secondMaxPI = sortedList[1][1]
      
      # call variants at different cutoffs and write to output
      inits = []
      finals = []            
      fltr = '|'
      origAlt = secondMaxBase if maxBase == origRef else maxBase
      altPI = secondMaxPI if maxBase == origRef else maxPI

      # reset variant type, reference base, variant base 
      vtype = '.'
      ref = origRef
      alt = origAlt

      if len(origAlt) == 1:
         vtype = 'SNP'
      elif origAlt == 'DEL':
         vtype = 'SDEL'
      elif origAlt.split('|')[0] in ['DEL', 'INS']:
         vtype = 'INDEL'
         ref = origAlt.split('|')[1]
         alt = origAlt.split('|')[2]

      # initial call comparing pred index with cutoff
      for cut in cutoff:
         if altPI > cut:
             finals.append('Y')
         else:
             finals.append('N')

      # if PI is higher than the minimum cutoff, and locus not in a deletion, apply filters
      if altPI > min(cutoff) and vtype in ['SNP', 'INDEL']:
         # low coverage filter
         if usedMT < 5:
            fltr += 'LM|' 

         # low number of strong MTs filter
         if strongMTCnt[origAlt] < 2 :
            fltr += 'LSM|'

         # homopolymer filter 
         if isHPorLowComp(chrom, pos, length, ref, alt)[0] and 1.0 * MTCnt[origAlt] / usedMT < 0.99:
            fltr += 'HP|'

         # low complexity filter
         if isHPorLowComp(chrom, pos, length, ref, alt)[1] and 1.0 * MTCnt[origAlt] / usedMT < 0.99:
            fltr += 'LC|'

         # strand bias and discordant pairs filter
         af_alt = 100.0 * alleleCnt[origAlt] / cvg
         pairs = discordPairCnt[origAlt] + concordPairCnt[origAlt] # total number of paired reads covering the pos
         if pairs >= 1000 and 1.0 * discordPairCnt[origAlt] / pairs >= 0.5:
            fltr += 'DP|'
         elif af_alt <= 60.0:
            refR = reverseCnt[origRef]
            refF = forwardCnt[origRef]
            altR = reverseCnt[origAlt]
            altF = forwardCnt[origAlt]
            fisher = scipy.stats.fisher_exact([[refR, refF], [altR, altF]])
            oddsRatio = fisher[0]
            pvalue = fisher[1]
            if pvalue < 0.00001 and (oddsRatio >= 50 or oddsRatio <= 1.0/50):
               fltr += 'SB|'

         # allele fraction ranking filter. Reject if AF of ALT is not top 2. Do not use. 
         #sortedAllele = sorted(alleleCnt.items(), key=operator.itemgetter(1), reverse=True)
         #maxAllele = sortedAllele[0][0]
         #secondAllele = sortedAllele[1][0]
         #if origAlt not in [maxAllele, secondAllele]:
         #   fltr += 'AFRank|'

         # base quality filter. Reject of Mean base quality < 22. Do not use. 
         # if origAlt in alleleCnt.keys():
         #    bqAlt =  1.0 * bqSum[origAlt] / alleleCnt[origAlt] 
         # else:
         #    bqAlt = 0.0
         # if bqAlt < 22.0:
         #    fltr += 'LowQ|'


         # mismatch filter. If average mismatches in the variant reads per 100 bases >= threshold, reject. Do not use. 
         #if alleleCnt[origAlt] > 0:
         #   if mismatchCnt[origAlt] / alleleCnt[origAlt] >= mismatchThr:
         #      fltr += 'MM|'
         
         # clustered position filter
         if vtype=='SNP':
            # R1
            refLe10 = sum(d <= 10 for d in r1EndPos[origRef])  # number of REF R2 reads with distance <= 10
            refGt10 = len(r1EndPos[origRef]) - refLe10         # number of REF R2 reads with distance > 10
            altLe10 = sum(d <= 10 for d in r1EndPos[origAlt])  # number of ALT R2 reads with distance <= 10
            altGt10 = len(r1EndPos[origAlt]) - altLe10         # number of ALT R2 reads with distance > 10

            fisher = scipy.stats.fisher_exact([[refLe10, refGt10], [altLe10, altGt10]])
            oddsRatio = fisher[0]
            pvalue = fisher[1]
            if pvalue < 0.001 and oddsRatio < 1.0/20:
               fltr += 'R1CP|'

            # R2
            refLe10 = sum(d <= 10 for d in r2EndPos[origRef])  # number of REF R2 reads with distance <= 10
            refGt10 = len(r2EndPos[origRef]) - refLe10         # number of REF R2 reads with distance > 10
            altLe10 = sum(d <= 10 for d in r2EndPos[origAlt])  # number of ALT R2 reads with distance <= 10
            altGt10 = len(r2EndPos[origAlt]) - altLe10         # number of ALT R2 reads with distance > 10

            fisher = scipy.stats.fisher_exact([[refLe10, refGt10], [altLe10, altGt10]])
            oddsRatio = fisher[0]
            pvalue = fisher[1]
            if pvalue < 0.001 and oddsRatio < 1.0/20:
               fltr += 'R2CP|'

      # write detailed output
      frac_alt = str(round((100.0 * alleleCnt[origAlt] / cvg),1)) 
      frac_A = str(round((100.0 * alleleCnt['A'] / cvg),1)) 
      frac_T = str(round((100.0 * alleleCnt['T'] / cvg),1)) 
      frac_G = str(round((100.0 * alleleCnt['G'] / cvg),1)) 
      frac_C = str(round((100.0 * alleleCnt['C'] / cvg),1))
      fracs = [str(alleleCnt['A']), str(alleleCnt['T']), str(alleleCnt['G']), str(alleleCnt['C']), frac_A, frac_T, frac_G, frac_C]
      
      MT_f_alt = str(round((100.0 * MTCnt[origAlt] / usedMT),1))
      MT_f_A = str(round((100.0 * MTCnt['A'] / usedMT),1))
      MT_f_T = str(round((100.0 * MTCnt['T'] / usedMT),1))
      MT_f_G = str(round((100.0 * MTCnt['G'] / usedMT),1))
      MT_f_C = str(round((100.0 * MTCnt['C'] / usedMT),1))
      MTs = [str(MT3Cnt), str(MT5Cnt), str(MT7Cnt), str(MT10Cnt), str(MTCnt['A']), str(MTCnt['T']), str(MTCnt['G']), str(MTCnt['C']), MT_f_A, MT_f_T, MT_f_G, MT_f_C]
      
      strongMT = [str(strongMTCnt['A']), str(strongMTCnt['T']), str(strongMTCnt['G']), str(strongMTCnt['C'])] 
      weakMT = [str(weakMTCnt['A']), str(weakMTCnt['T']), str(weakMTCnt['G']), str(weakMTCnt['C'])] 
      predIdx = [str(round(finalDict['A'], 2)), str(round(finalDict['T'], 2)), str(round(finalDict['G'], 2)), str(round(finalDict['C'], 2))]

      out_long = '\t'.join([chrom, pos, ref, alt, vtype, str(cvg), str(allFrag), str(allMT), str(usedFrag), str(usedMT), str(round(finalDict[origAlt], 2)), str(alleleCnt[origAlt]), frac_alt, str(MTCnt[origAlt]), MT_f_alt, str(strongMTCnt[origAlt]), str(weakMTCnt[origAlt])] + fracs + MTs + strongMT + weakMT + predIdx + finals + [fltr]) + '\n'
      
   return out_long

#--------------------------------------------------------------------------------------
# main function
#--------------------------------------------------------------------------------------
def main(args):
   # change working directory to runDir and make output directories 
   os.chdir(args.runPath)
   # make /VCF directory to keep all VCF files
   if not os.path.exists('vcf'):
      os.makedirs('vcf')
   # initialize logger (param is whether or not to print to stdout)
   runLog.init(False, args.logFile)
   # log run start
   timeStart = datetime.datetime.now()
   print("started at " + str(timeStart))

   outfile_long = open(args.outLong, 'w')
   outfile_short = open(args.outShort, 'w')

   init_names = ['InitialCall_' + str(abs(x)) for x in cutoff]
   final_names = ['FinalCall_' + str(abs(x)) for x in cutoff]
   cat_names = ['Cat_' + str(abs(x)) for x in cutoff]

   header_1 = ['CHROM', 'POS', 'REF', 'ALT', 'TYPE', 'DP', 'FR', 'MT', 'UFR', 'UMT', 'PI', 'VDP', 'VAF', 'VMT', 'VMF', 'VSM', 'VWM', 'DP_A', 'DP_T', 'DP_G', 'DP_C', 'AF_A', 'AF_T', 'AF_G', 'AF_C', 'MT_3RPM', 'MT_5RPM', 'MT_7RPM', 'MT_10RPM',  'UMT_A', 'UMT_T', 'UMT_G', 'UMT_C', 'UMF_A', 'UMF_T', 'UMF_G', 'UMF_C', 'VSM_A', 'VSM_T', 'VSM_G', 'VSM_C', 'VWM_A', 'VWM_T', 'VWM_G', 'VWM_C', 'PI_A', 'PI_T', 'PI_G', 'PI_C']
   lenHeaderBeforeCutoff = len(header_1)
   header_long = '\t'.join(header_1 + final_names + ['FILTER']) + '\n'

   outfile_long.write(header_long)
   outfile_short.write(header_long)

   # read in tandem repeat list
   repRegion = defaultdict(list)
   for line in open(repBed, 'r'):
      lineList = line.strip().split()
      chrom = 'chr' + lineList[0]
      regionStart = lineList[1]
      regionEnd = lineList[2]
      repRegion[chrom].append([int(regionStart), int(regionEnd)])

   # read in simple repeat, low complexity, satelite list
   srRegion = defaultdict(list)
   for line in open(srBed, 'r'):
      lineList = line.strip().split()
      chrom = 'chr' + lineList[0]
      regionStart = lineList[1]
      regionEnd = lineList[2]
      if lineList[3] == 'Simple_repeat':
         repType = 'SR'
      elif lineList[3] == 'Low_complexity':
         repType = 'LowC'
      elif lineList[3] == 'Satellite':
         repType = 'SL'
      else:
         repType = 'Other_Repeat'
      srRegion[chrom].append([regionStart, regionEnd, repType])

   locList = []
   for line in open(args.bedName, 'r'):
       lineList = line.strip().split('\t')
       chrom = lineList[0]
       regionStart = lineList[1]
       regionEnd = lineList[2]

       pos = int(regionStart)
       lineEnd = False

       while not lineEnd:
          locList.append((chrom, str(pos)))
          if str(pos) == regionEnd:
             lineEnd = True
          else:
             pos += 1

   pool = mp.Pool(processes=args.nCPU)

   results = [pool.apply_async(vc, args=(args.bamName, x[0], x[1], args.minBQ, args.minMQ, args.StrongMtThr, args.hpLen, args.mismatchThr, args.MTdrop, args.ds)) for x in locList]
   output = [p.get() for p in results]
   
   # ALL repeats filter. If MT fraction < 40% and the variant is inside the tandem repeat region, reject. 
   for i in range(len(output)):
      outline = output[i]
      lineList = outline.strip().split('\t')
      chromTr = lineList[header_1.index('CHROM')]
      altTr = lineList[header_1.index('ALT')]
      lowestCut = lineList[lenHeaderBeforeCutoff]
      try:
         posTr = int(lineList[header_1.index('POS')])
      except ValueError:
         continue
      try:
         altMtFracTr = float(lineList[header_1.index('VMF')])
      except ValueError:
         continue
      if lowestCut == 'Y' and altTr != 'DEL':
         # check tandem repeat if MT fraction < 40%
         if altMtFracTr < 40:
            for values in repRegion[chromTr]:
               if posTr >= values[0] and posTr <= values[1]:
                  lineList[-1] += 'TR|'
                  output[i] = '\t'.join(lineList) + '\n'
                  break

         # check simple repeat, lc, sl
         for values in srRegion[chromTr]:
            if posTr >= int(values[0]) and posTr <= int(values[1]):
               lineList[-1] += values[2] + '|'
               output[i] = '\t'.join(lineList) + '\n'
               break


   # create a vcf file at each cutoff
   base = os.path.basename(args.bamName)
   sampleName = os.path.splitext(base)[0]
   header_vcf = '##fileformat=VCFv4.2' + '\n' + \
         '##FILTER=<ID=LM,Description="Less than 5 MTs">' + '\n' + \
         '##FILTER=<ID=TR,Description="Tandem repeat region, as defined in TFR">' + '\n' + \
         '##FILTER=<ID=SR,Description="Simple repeat region, as defined in Repeat Masker, repClass=Simple Repeat">' + '\n' + \
         '##FILTER=<ID=LowC,Description="Tandem repeat region, as defined in Repeat Masker, repClass=Low_Complexity">' + '\n' + \
         '##FILTER=<ID=SL,Description="Satelite region, as defined in Repeat Masker, repClass=Satelite">' + '\n' + \
         '##FILTER=<ID=LSM,Description="Less than 2 strong MTs">' + '\n' + \
         '##FILTER=<ID=SB,Description="Strand Bias">' + '\n' + \
         '##FILTER=<ID=DP,Description="Two many discordant pairs">' + '\n' + \
         '##FILTER=<ID=HP,Description="Inside or flanked by homopolymer regions">' + '\n' + \
         '##FILTER=<ID=LC,Description="Inside or flanked by low complexity regions (any two nucleotides combine for 99% frequency near the variant)">' + '\n' + \
         '##FILTER=<ID=MM,Description="Too many mismatches in a read. Default threshold is 6.5 per 100 bases">' + '\n' + \
         '##FILTER=<ID=R1CP,Description="Variant are clustered at the end of R1 reads">' + '\n' + \
         '##FILTER=<ID=R2CP,Description="Variant are clustered at the end of R2 reads">' + '\n' + \
         '##FILTER=<ID=AFRank,Description="Variant is not the highest AF in non reference">' + '\n' + \
         '##FILTER=<ID=LowQ,Description="Low base quality (mean base Q < 22)">' + '\n' + \
         '##INFO=<ID=SAMPLE,Number=1,Type=String,Description="Sample name">' + '\n' + \
         '##INFO=<ID=TYPE,Number=1,Type=String,Description="Variant Type: SNV Insertion Deletion">' + '\n' + \
         '##INFO=<ID=DP,Number=1,Type=Integer,Description="Total Depth">' + '\n' + \
         '##INFO=<ID=VDP,Number=1,Type=Integer,Description="Variant Depth">' + '\n' + \
         '##INFO=<ID=VAF,Number=1,Type=Integer,Description="Allele Frequency">' + '\n' + \
         '##INFO=<ID=MT,Number=1,Type=Integer,Description="Total MT Depth">' + '\n' + \
         '##INFO=<ID=VMT,Number=1,Type=Integer,Description="Variant MT Depth">' + '\n' + \
         '##INFO=<ID=VMF,Number=1,Type=Integer,Description="MT Allele Frequency">' + '\n' + \
         '##INFO=<ID=VSM,Number=1,Type=Integer,Description="Strong MT">' + '\n' + \
         '##FORMAT=<ID=GT,Number=1,Type=String,Description="Genotype">' + '\n' + \
         '##FORMAT=<ID=DP,Number=1,Type=Integer,Description="Total Depth">' + '\n' + \
         '##FORMAT=<ID=VDP,Number=1,Type=Integer,Description="Variant Depth">' + '\n' + \
         '##FORMAT=<ID=VAF,Number=1,Type=Integer,Description="Allele Frequency">' + '\n' + \
         '\t'.join(['#CHROM', 'POS', 'ID', 'REF', 'ALT', 'QUAL', 'FILTER', 'INFO', 'FORMAT'] + [sampleName]) + '\n'

   for i in cutoff:
      vcfName = 'vcf/' + args.vcfNamePrefix + '.' + str(i) + '.vcf'
      outvcf = open(vcfName, 'w')
      outvcf.write(header_vcf)
      outvcf.close()

   for line in output:
      outfile_long.write(line)
      # This is the result for the least conservative cutoff. Write to file if only 'Y'
      fields = line.strip().split('\t')
      lowestCutoff = fields[lenHeaderBeforeCutoff]
      if lowestCutoff == 'Y':
         outfile_short.write(line)
      # output to vcf format
      CHROM = fields[header_1.index('CHROM')]
      POS = fields[header_1.index('POS')]
      ID = '.'
      REF = fields[header_1.index('REF')]
      ALT = fields[header_1.index('ALT')]
      TYPE = fields[header_1.index('TYPE')]
      DP = fields[header_1.index('DP')]
      MT = fields[header_1.index('MT')]

      try:
         QUAL = str(int(float(fields[header_1.index('PI')]))) 
      except ValueError:
         QUAL = '0'

      try:
         VAF = str(float(fields[header_1.index('VAF')])/100)
      except ValueError:
         VAF = '-1'

      try:
         VMF = str(float(fields[header_1.index('VMF')])/100)
      except ValueError:
         VMF = '-1'

      VMT = fields[header_1.index('VMT')]
      VDP = fields[header_1.index('VDP')]
      VSM = fields[header_1.index('VSM')]

      FILTER = 'PASS' if fields[-1]=='|' else fields[-1].strip('|')
      INFO = ';'.join(['SAMPLE=' + sampleName, 'TYPE=' + TYPE, 'DP=' + DP, 'VDP=' + VDP, 'VAF=' + VAF, 'MT=' + MT, 'VMT=' + VMT, 'VMF=' + VMF, 'VSM=' + VSM]) 
      FORMAT = 'GT:DP:VDP:VAF'
      SAMPLE = ':'.join(['0/1', DP, VDP, VAF])

      for i in range(len(cutoff)):
         vcfName = 'vcf/' + args.vcfNamePrefix + '.' + str(cutoff[i]) + '.vcf' 
         outvcf = open(vcfName, 'a')
         vcfLine = '\t'.join([CHROM, POS, ID, REF, ALT, QUAL, FILTER, INFO, FORMAT, SAMPLE]) + '\n'
         if fields[lenHeaderBeforeCutoff+i] == 'Y' and ALT != 'DEL': # if PI > cutoff, write to vcf (regardless of filters)
            outvcf.write(vcfLine)
         outvcf.close()


   outfile_long.close()
   outfile_short.close()

   # log run completion
   timeEnd = datetime.datetime.now()
   print("completed running at " + str(timeEnd) + "\n")
   print("total time: "+ str(timeEnd-timeStart) + "\n")   


#--------------------------------------------------------------------------------------
# Run 
#--------------------------------------------------------------------------------------
def run():
   parser = argparse.ArgumentParser(description='Perform SNV and INDEL calling with barcoded reads')
   parser.add_argument('--runPath', default=None, help='path to working directory')
   parser.add_argument('--bedName', default=None, help='BED file with target region')
   parser.add_argument('--bamName', default=None, help='BAM file')
   parser.add_argument('--logFile', default=None, help='log file')
   parser.add_argument('--vcfNamePrefix', default=None, help='Prefix of output VCF file name')
   parser.add_argument('--outLong', default=None, help='long version of output including all positions in the target region')
   parser.add_argument('--outShort', default=None, help='short version of output including only the positions called at the most lenient cutoff')
   parser.add_argument('--nCPU', type=int, default=1, help='number of CPU used in parallel')
   parser.add_argument('--minBQ', type=int, default=20, help='minimum base quality allowed for analysis')
   parser.add_argument('--minMQ', type=int, default=30, help='minimum mapping quality allowed for analysis')
   parser.add_argument('--StrongMtThr', type=float, default=4.0, help='Threshold for strong MT')
   parser.add_argument('--hpLen', type=int, default=8, help='Minimum length for homopolymers. 2*hpLen will be low complexity length')
   parser.add_argument('--mismatchThr', type=float, default=6.0, help='average number of mismatches per 100 bases allowed')
   parser.add_argument('--MTdrop', type=int, default=0, help='Drop MTs with lower than or equal to X reads')
   parser.add_argument('--ds', type=int, default=10000, help='Randomly downsample to X MTs (max number of MTs at any position)')
   args = parser.parse_args()
   main(args)

#----------------------------------------------------------------------------------------------
#pythonism to run from the command line
#----------------------------------------------------------------------------------------------
if __name__ == "__main__":
    run()





