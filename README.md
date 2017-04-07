<<<<<<< HEAD
This repository contains scripts and data files supporting smCounter, a versatile UMI-aware variant caller that detects both somatic and germline SNVs and indels with high sensitivity and specificity. An example of running smCounter is included. The algorithm and validation results were published in "Detecting very low allele fraction variants using targeted DNA sequencing and a novel molecular barcode-aware variant caller", BMC Genomics, 2017 18:5. https://bmcgenomics.biomedcentral.com/articles/10.1186/s12864-016-3425-4
=======
This repository contains scripts and data files supporting the article "Detecting very low allele fraction variants using targeted DNA sequencing and a novel molecular barcode-aware variant caller", BMC Genomics, 2017 18:5. https://bmcgenomics.biomedcentral.com/articles/10.1186/s12864-016-3425-4
>>>>>>> 8eb8cf01029b9b8fbecb5d3326628c747fc1ccad

# File description 
  * smCounter.py -- Python script for smCounter, a barcode aware somatic variant caller that integrates molecular barcode information into the variant calling algorithm. The script was developed and tested under Python v2.7.3. Python modules required: pysam, math, scipy, random, multiprocessing. Samtools v0.1.19 and Bedtools are also required. 
  * run_log.py -- custom python script to direct stdout to log file
  * ds.mt.py -- Python script for downsampling barcode over the entire target region. 
  * ds.reads.withinMT.py -- Python script for downsampling reads within barcodes. 
  * ds.allele.fraction.py -- Python script for reducing the variant allele fraction at given variant loci. 
  * primers.NA12878-194-genes-63-indels.10867.coding.bed -- Bed file for the target region of N0030 panel
  * SR_LC_SL.nochr.bed -- simple repeat, low complexity, satellite region
  * simpleRepeat.bed -- tandem repeat region
<<<<<<< HEAD
  * example -- The folder contains an example of running smCounter on a BAM file with UMIs. The BAM file is a subset of N0030 data (see BMC Genomics paper) that covers part of BRCA1 gene
=======
>>>>>>> 8eb8cf01029b9b8fbecb5d3326628c747fc1ccad


