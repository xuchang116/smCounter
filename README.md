# DeepCall
DeepCall is a variant caller that integrates molecular barcode information into the inherent statistical model. DeepCall has demonstrated very good accuracy in detecting SNVs and short INDELs at very low allele fractions. 

# Requirements
* Python (developed and tested with v2.7.3, with the following modules installed: math, pysam, scipy, random, multiprocessing) 
* samtools v0.1.19
* bedtools 

# Input/output files
* Input: BAM file with index, BED file (sorted and merged, with 'chr')
* Output: Variant call set in VCF format, details of each target locus in text file, log file.  
