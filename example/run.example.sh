#!/bin/bash
# vim: tabstop=9 expandtab shiftwidth=3 softtabstop=3

python ../smCounter.py \
   --outPrefix example \
   --bamFile example.bam \
   --bedTarget example.bed \
   --mtDepth 3612 \
   --rpb 8.6 \
   --nCPU 10 \
   --minBQ 20 \
   --minMQ 30 \
   --hpLen 8 \
   --mismatchThr 6.0 \
   --mtDrop 1 \
   --maxMT 0 \
   --primerDist 2 \
   --threshold 0 \
   --refGenome /qgen/home/rvijaya/downloads/alt_hap_masked_ref/ucsc.hg19.fasta \
   --bedTandemRepeats ../simpleRepeat.bed \
   --bedRepeatMaskerSubset ../SR_LC_SL.nochr.bed \
   --bedtoolsPath /qgen/bin/bedtools-2.25.0/bin/ \
   --runPath ./ \
   --logFile example





