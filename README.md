This repository contains scripts and data files supporting the manuscript "Detecting very low allele fraction variants using targeted DNA sequencing with molecular barcodes". 

# File description 
  * smCounter.py -- Python script for smCounter, a barcode aware somatic variant caller that integrates molecular barcode information into the variant calling algorithm. The script was developed and tested under Python v2.7.3. Python modules required: pysam, math, scipy, random, multiprocessing. Samtools v0.1.19 and Bedtools are also required. 

  * smCounter.command.sh -- An example command for running smCounter. The parameters in the example are the default settings. 

  * ds.mt.py -- Python script for downsampling barcode over the entire target region. 

  * ds.reads.withinMT.py -- Python script for downsampling reads within barcodes. 

  * ds.allele.fraction.py -- Python script for reducing the variant allele fraction at given variant loci. 

  * primers.NA12878-194-genes-63-indels.10867.coding.bed -- Bed file for the target region of N0030 panel


