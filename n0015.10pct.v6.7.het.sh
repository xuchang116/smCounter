#!/bin/bash
# vim: tabstop=9 expandtab shiftwidth=3 softtabstop=3
# 1% variants only
# Chang Xu, 09APR2016 

RefSDF=/qgen/home/nr/ref/ucsc.hg19.sdf
#bedFile=/qgen/home/xuc/VariantCallingPaper/N0015/primers.NA12878.3587.target.coding.bed   # non HC region included. has 'chr'
runPath=/qgen/home/xuc/VariantCallingPaper/N0015/10pct/

cd /qgen/home/xuc/VariantCallingPaper/N0015

# get heterozygous variants (0/1) only
grep 0/1 GIB.NA12878.hc.dil.additional.noheader.sorted.vcf > GIB.NA12878.hc.dil.het.additional.noheader.sorted.vcf

cat VCFheader.txt GIB.NA12878.hc.dil.het.additional.noheader.sorted.vcf > GIB.NA12878.hc.dil.het.additional.header.sorted.vcf
# move 1/1 and multi-allelic variants to background
grep -v 0/1 GIB.NA12878.hc.dil.additional.noheader.sorted.vcf > misc/not.het.vcf
cat HG002-multiall-fullcombine.additional.N0015.vcf misc/not.het.vcf | bedtools sort -header -i > HG002-multiall-fullcombine.additional.N0015.nohet.vcf

#########################################################
### whole region 
#########################################################
caller=v6.7.het
cd /qgen/home/xuc/VariantCallingPaper/N0015
# Process GIB vcf and create ground truth set
awk '{if (/^#/) {print} else {if (/^chr/) {print} else {print "chr"$0}}}' GIB.NA12878.hc.dil.het.additional.header.sorted.vcf | bedtools sort -header -i | bgzip -c > GIB.NA12878.hc.dil.het.additional.header.sorted.vcf.gz
tabix -f -p vcf GIB.NA12878.hc.dil.het.additional.header.sorted.vcf.gz > /dev/null
GroundTruth=/qgen/home/xuc/VariantCallingPaper/N0015/GIB.NA12878.hc.dil.het.additional.header.sorted.vcf.gz

cd $runPath
Rscript /qgen/home/xuc/VariantCallingPaper/code/getSummary.withRTG.R     $runPath   /qgen/home/xuc/VariantCallingPaper/N0015/GIB.NA12878.hc.dil.het.additional.noheader.sorted.vcf  result/v6.7/Table.txt result/v6.7/Table-baseline.txt        395303  $caller   roc.$caller.n0015.10pct.png    summary.$caller.n0015.10pct.csv  1>/dev/null 2>$runPath/rtg.log	

#########################################################
### coding region 
#########################################################
cd /qgen/home/xuc/VariantCallingPaper/N0015/
caller=v6.7.coding.het
bedFile=primers.NA12878.3587.target.coding.bed
# Process GIB vcf
bedtools intersect -a GIB.NA12878.hc.dil.het.additional.header.sorted.vcf.gz -b $bedFile -header > GIB.NA12878.hc.dil.het.coding.additional.header.sorted.vcf
awk '{if (/^#/) {print} else {if (/^chr/) {print} else {print "chr"$0}}}' GIB.NA12878.hc.dil.het.coding.additional.header.sorted.vcf | bedtools sort -header -i | bgzip -c > GIB.NA12878.hc.dil.het.coding.additional.header.sorted.vcf.gz
tabix -f -p vcf GIB.NA12878.hc.dil.het.coding.additional.header.sorted.vcf.gz > /dev/null
GroundTruth=/qgen/home/xuc/VariantCallingPaper/N0015/GIB.NA12878.hc.dil.het.coding.additional.header.sorted.vcf.gz   # NOTE: it's already in HC region

cd $runPath
sed '/^#/ d' /qgen/home/xuc/VariantCallingPaper/N0015/GIB.NA12878.hc.dil.het.coding.additional.header.sorted.vcf > /qgen/home/xuc/VariantCallingPaper/N0015/GIB.NA12878.hc.dil.het.coding.additional.noheader.sorted.vcf 
sed -i 's/chr//g' /qgen/home/xuc/VariantCallingPaper/N0015/GIB.NA12878.hc.dil.het.coding.additional.noheader.sorted.vcf
Rscript /qgen/home/xuc/VariantCallingPaper/code/getSummary.withRTG.R     $runPath   /qgen/home/xuc/VariantCallingPaper/N0015/GIB.NA12878.hc.dil.het.coding.additional.noheader.sorted.vcf  result/v6.7.coding/Table.txt result/v6.7.coding/Table-baseline.txt        42339  $caller   roc.$caller.n0015.10pct.png    summary.$caller.n0015.10pct.csv  1>/dev/null 2>$runPath/rtg.log	


#########################################################
### noncoding region 
#########################################################
cd /qgen/home/xuc/VariantCallingPaper/N0015/
caller=v6.7.noncoding.het
bedFile=primers.NA12878.3587.target.noncoding.bed
# Process GIB vcf
bedtools intersect -a GIB.NA12878.hc.dil.het.additional.header.sorted.vcf.gz -b $bedFile -header > GIB.NA12878.hc.dil.het.noncoding.additional.header.sorted.vcf
awk '{if (/^#/) {print} else {if (/^chr/) {print} else {print "chr"$0}}}' GIB.NA12878.hc.dil.het.noncoding.additional.header.sorted.vcf | bedtools sort -header -i | bgzip -c > GIB.NA12878.hc.dil.het.noncoding.additional.header.sorted.vcf.gz
tabix -f -p vcf GIB.NA12878.hc.dil.het.noncoding.additional.header.sorted.vcf.gz > /dev/null
GroundTruth=/qgen/home/xuc/VariantCallingPaper/N0015/GIB.NA12878.hc.dil.het.noncoding.additional.header.sorted.vcf.gz   # NOTE: it's already in HC region

cd $runPath
# get summary results
sed '/^#/ d' /qgen/home/xuc/VariantCallingPaper/N0015/GIB.NA12878.hc.dil.het.noncoding.additional.header.sorted.vcf > /qgen/home/xuc/VariantCallingPaper/N0015/GIB.NA12878.hc.dil.het.noncoding.additional.noheader.sorted.vcf
sed -i 's/chr//g' /qgen/home/xuc/VariantCallingPaper/N0015/GIB.NA12878.hc.dil.het.noncoding.additional.noheader.sorted.vcf

Rscript /qgen/home/xuc/VariantCallingPaper/code/getSummary.withRTG.R     $runPath   /qgen/home/xuc/VariantCallingPaper/N0015/GIB.NA12878.hc.dil.het.noncoding.additional.noheader.sorted.vcf  result/v6.7.noncoding/Table.txt result/v6.7.noncoding/Table-baseline.txt        352964  $caller   roc.$caller.n0015.10pct.png    summary.$caller.n0015.10pct.csv  1>/dev/null 2>$runPath/rtg.log	



