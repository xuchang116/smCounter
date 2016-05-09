#!/bin/bash
# vim: tabstop=9 expandtab shiftwidth=3 softtabstop=3
# variant calling on 2% mix with V6.7, coding and noncoding region
# Chang Xu, 10APR2016 

RefSDF=/qgen/home/nr/ref/ucsc.hg19.sdf
#bedFile=/qgen/home/xuc/VariantCallingPaper/N0015/primers.NA12878.3587.target.coding.bed   # non HC region included. has 'chr'
runPath=/qgen/home/xuc/VariantCallingPaper/N0015/2pct/

cd /qgen/home/xuc/VariantCallingPaper/N0015

# create coding and non-coding bed files
bedtools intersect -a primers.NA12878.3587.target.bed -b step99.output.NCBI-GTF-coding-GRCh37.bed | sortBed -i |  bedtools merge -i > primers.NA12878.3587.target.coding.bed
bedtools subtract -a primers.NA12878.3587.target.bed -b step99.output.NCBI-GTF-coding-GRCh37.bed | sortBed -i |   bedtools merge -i > primers.NA12878.3587.target.noncoding.bed

#########################################################
### coding region 
#########################################################

caller=v6.7.coding
bedFile=primers.NA12878.3587.target.coding.bed
# Process GIB vcf
bedtools intersect -a GIB.NA12878.hc.dil.additional.header.sorted.vcf.gz -b $bedFile -header > GIB.NA12878.hc.dil.coding.additional.header.sorted.vcf
awk '{if (/^#/) {print} else {if (/^chr/) {print} else {print "chr"$0}}}' GIB.NA12878.hc.dil.coding.additional.header.sorted.vcf | bedtools sort -header -i | bgzip -c > GIB.NA12878.hc.dil.coding.additional.header.sorted.vcf.gz
tabix -f -p vcf GIB.NA12878.hc.dil.coding.additional.header.sorted.vcf.gz > /dev/null
GroundTruth=/qgen/home/xuc/VariantCallingPaper/N0015/GIB.NA12878.hc.dil.coding.additional.header.sorted.vcf.gz   # NOTE: it's already in HC region

# create HC bed file, 42339 bp
sed 's/chr//g' $bedFile > misc/primers.NA12878.3587.target.coding.nochr.bed
bedtools intersect -a misc/primers.NA12878.3587.target.coding.nochr.bed -b /qgen/home/xuc/GIAB/NA12878/union13callableMQonlymerged_addcert_nouncert_excludesimplerep_excludesegdups_excludedecoy_excludeRepSeqSTRs_noCNVs_v2.19_2mindatasets_5minYesNoRatio.bed | sortBed -i | bedtools merge -i > misc/shared.region.coding.merged.bed

# create NA24385 variant list in the shared region (32 variants)
bedtools intersect -a HG002-multiall-fullcombine.additional.N0015.vcf -b misc/shared.region.coding.merged.bed > misc/HG002-multiall-fullcombine.shared.region.coding.vcf

cd $runPath

# create directory to keep processed vcfs
if [ ! -s processed ]; then
   mkdir processed
fi
# create results directory
if [ ! -s result ]; then
   mkdir result 
fi
# create caller directory
if [ ! -s result/$caller ]; then
   mkdir result/$caller
fi
# create tmp directory
if [ ! -s tmp ]; then
   mkdir tmp 
fi

echo "RTG for coding region"
# limit variants in GIB HC region and dilute NA24385 variants
python /qgen/home/xuc/VariantCallingPaper/code/HCandDilute.fast.alt.py  $runPath	$bedFile v6.7.VariantList.short.txt	$caller.VariantList.short.hc.dil.txt  /qgen/home/xuc/VariantCallingPaper/N0015/misc/shared.region.coding.merged.bed   /qgen/home/xuc/VariantCallingPaper/N0015/misc/HG002-multiall-fullcombine.shared.region.coding.vcf

# Process output vcfs
for cutoff in 10 15 20 25 30 35 40 45 50 55 60 65 70 75 80 85 90 95 100 125; do
   # create results cutoff directory 
   if [ ! -s result/$caller/$cutoff ]; then
      mkdir result/$caller/$cutoff 
   fi

   # limit existing processed vcf from whole region to coding region
   file1=vcf/v6.7.$cutoff.hc.dil.header.vcf
   file2=vcf/$caller.$cutoff.hc.dil.header.vcf
   bedtools intersect -a $file1 -b /qgen/home/xuc/VariantCallingPaper/N0015/misc/shared.region.coding.merged.bed -header > $file2
   awk '{if (/^#/) {print} else {if (/^chr/) {print} else {print "chr"$0}}}' $file2 > processed/${file2:4}
   bedtools sort -header -i processed/${file2:4} | bgzip -c > processed/${file2:4}.gz
   tabix -f -p vcf processed/${file2:4}.gz > /dev/null

   # Run RTG evaluation
   file=processed/$caller.$cutoff.hc.dil.header.vcf.gz
   # check if RTG results not already exist
   if [ -s result/$caller/$cutoff/filter.applied ]; then
      echo "RTG tools evaluation results already exist for ${file:10}"
   else
      /qgen/home/xuc/software/RTG/rtg-tools-3.5/rtg vcfeval --all-records -T 1 -b $GroundTruth -c $file --squash-ploidy --no-gzip -t $RefSDF -f "QUAL" -o result/$caller/$cutoff/filter.ignored &>> $runPath/rtg.log 
      /qgen/home/xuc/software/RTG/rtg-tools-3.5/rtg vcfeval --baseline-tp -T 1 -b $GroundTruth -c $file --squash-ploidy --no-gzip -t $RefSDF -f "QUAL" -o result/$caller/$cutoff/filter.applied &>> $runPath/rtg.log
   fi
   
   # create tables by merging info from TP, FP, and FN vcf files
   awk '!/^#/ {OFS="\t"; if ($7 == "PASS") {Pred="TP"} else {Pred="FN"}; print $1":"$2";"$4"/"$5,Pred}' result/$caller/$cutoff/filter.ignored/tp.vcf > tmp/tp.txt 2>>$runPath/rtg.log
   awk '!/^#/ {OFS="\t"; if ($7 == "PASS") {Pred="FP"} else {Pred="TN"}; print $1":"$2";"$4"/"$5,Pred}' result/$caller/$cutoff/filter.ignored/fp.vcf > tmp/fp.txt 2>>$runPath/rtg.log
   awk '!/^#/ {OFS="\t"; print $1":"$2";"$4"/"$5,"FN"}' result/$caller/$cutoff/filter.ignored/fn.vcf > tmp/fn.txt 2>>$runPath/rtg.log
   awk '!/^#/ {OFS="\t"; print $1":"$2";"$4"/"$5,"TP"}' result/$caller/$cutoff/filter.applied/tp-baseline.vcf > result/$caller/$cutoff/tp-baseline.txt 2>>$runPath/rtg.log
   
   # merge tables to one for one cutoff and remove the duplicates (caused by RTG tools limitations/bugs)
   cat tmp/tp.txt tmp/fp.txt tmp/fn.txt | sort | uniq | awk -v duplicate=0 '{OFS="\t"; indb=ind; typeb=type; ind=$1; type=$2; if (ind==indb) {if (typeb ~ /T/) {print indb,typeb} else if (type ~ /T/) {print ind,type}; duplicate=1} else if (duplicate==0) {print indb,typeb} else {duplicate=0}}' | awk 'NF>0' > result/$caller/$cutoff/table.txt 2>>$runPath/rtg.log
done    # end of cutoff loop
   
# integrate tables to one Table for each group -- this is at rpb level, after looping all cutoffs
if [ ! -s result/$caller ]; then
   echo "warning: no RTG tools results found for group: "$caller 
else
   # initialize
   count=1
   OutField="0"
   OutHeader="#Variant"
   
   # folders for each cutoff
   cutOffs=$(ls result/$caller)
   for cutOff in $cutOffs; do
      if [ $count -eq 1 ]; then
         cp result/$caller/${cutOff}/table.txt result/$caller/Table.txt
         cp result/$caller/${cutOff}/tp-baseline.txt result/$caller/Table-baseline.txt			
      else			
         OutField=${OutField}",1."$count
         #
         # join tables from all cutoff values
         join -a1 -a2 -1 1 -2 1 -e "TN" -o ${OutField}",2.2" -t $'\t' <(sort -t $'\t' -k1,1 result/$caller/Table.txt | uniq) <(sort -t $'\t' -k1,1 result/$caller/${cutOff}/table.txt | uniq) > result/$caller/Table.tmp
         mv result/$caller/Table.tmp result/$caller/Table.txt
         #
         # same for baseline
         join -a1 -a2 -1 1 -2 1 -e "FN" -o ${OutField}",2.2" -t $'\t' <(sort -t $'\t' -k1,1 result/$caller/Table-baseline.txt | uniq) <(sort -t $'\t' -k1,1 result/$caller/${cutOff}/tp-baseline.txt | uniq) > result/$caller/Table.tmp
         mv result/$caller/Table.tmp result/$caller/Table-baseline.txt				
      fi
      #
      count=$((count+1))
      OutHeader=${OutHeader}$'\t'"Thr="$cutOff
   done
   echo "$OutHeader" > result/$caller/Header.txt
fi
   
# add header info (cutoff values)
sort result/$caller/Table.txt | uniq > result/$caller/Table.tmp
cat result/$caller/Header.txt result/$caller/Table.tmp > result/$caller/Table.txt

sort result/$caller/Table-baseline.txt | uniq > result/$caller/Table-baseline.tmp
cat result/$caller/Header.txt result/$caller/Table-baseline.tmp > result/$caller/Table-baseline.txt

rm result/$caller/Header.txt result/$caller/Table.tmp result/$caller/Table-baseline.tmp
   
# get summary results
sed '/^#/ d' /qgen/home/xuc/VariantCallingPaper/N0015/GIB.NA12878.hc.dil.coding.additional.header.sorted.vcf > /qgen/home/xuc/VariantCallingPaper/N0015/GIB.NA12878.hc.dil.coding.additional.noheader.sorted.vcf 
sed -i 's/chr//g' /qgen/home/xuc/VariantCallingPaper/N0015/GIB.NA12878.hc.dil.coding.additional.noheader.sorted.vcf
Rscript /qgen/home/xuc/VariantCallingPaper/code/validateWithGIB.withRTG.v6.7.R     $runPath	$caller.VariantList.short.hc.dil.txt    /qgen/home/xuc/VariantCallingPaper/N0015/GIB.NA12878.hc.dil.coding.additional.noheader.sorted.vcf  result/$caller/Table.txt result/$caller/Table-baseline.txt        42339  roc.$caller.n0015.2pct.png    summary.$caller.n0015.2pct.csv        details.$caller.n0015.2pct.csv   $caller 1>/dev/null 2>$runPath/rtg.log	



#########################################################
### noncoding region 
#########################################################

cd /qgen/home/xuc/VariantCallingPaper/N0015
caller=v6.7.noncoding
bedFile=primers.NA12878.3587.target.noncoding.bed
# Process GIB vcf
bedtools intersect -a GIB.NA12878.hc.dil.additional.header.sorted.vcf.gz -b $bedFile -header > GIB.NA12878.hc.dil.noncoding.additional.header.sorted.vcf
awk '{if (/^#/) {print} else {if (/^chr/) {print} else {print "chr"$0}}}' GIB.NA12878.hc.dil.noncoding.additional.header.sorted.vcf | bedtools sort -header -i | bgzip -c > GIB.NA12878.hc.dil.noncoding.additional.header.sorted.vcf.gz
tabix -f -p vcf GIB.NA12878.hc.dil.noncoding.additional.header.sorted.vcf.gz > /dev/null
GroundTruth=/qgen/home/xuc/VariantCallingPaper/N0015/GIB.NA12878.hc.dil.noncoding.additional.header.sorted.vcf.gz   # NOTE: it's already in HC region

# create HC bed file, 352964 bp
sed 's/chr//g' $bedFile > misc/primers.NA12878.3587.target.noncoding.nochr.bed
bedtools intersect -a misc/primers.NA12878.3587.target.noncoding.nochr.bed -b /qgen/home/xuc/GIAB/NA12878/union13callableMQonlymerged_addcert_nouncert_excludesimplerep_excludesegdups_excludedecoy_excludeRepSeqSTRs_noCNVs_v2.19_2mindatasets_5minYesNoRatio.bed | sortBed -i | bedtools merge -i > misc/shared.region.noncoding.merged.bed

# create NA24385 variant list in the shared region (32 variants)
bedtools intersect -a HG002-multiall-fullcombine.additional.N0015.vcf -b misc/shared.region.noncoding.merged.bed > misc/HG002-multiall-fullcombine.shared.region.noncoding.vcf

cd $runPath

# create directory to keep processed vcfs
if [ ! -s processed ]; then
   mkdir processed
fi
# create results directory
if [ ! -s result ]; then
   mkdir result 
fi
# create caller directory
if [ ! -s result/$caller ]; then
   mkdir result/$caller
fi
# create tmp directory
if [ ! -s tmp ]; then
   mkdir tmp 
fi

echo "RTG for noncoding region"
# limit variants in GIB HC region and dilute NA24385 variants
python /qgen/home/xuc/VariantCallingPaper/code/HCandDilute.fast.alt.py  $runPath	$bedFile v6.7.VariantList.short.txt	$caller.VariantList.short.hc.dil.txt  /qgen/home/xuc/VariantCallingPaper/N0015/misc/shared.region.noncoding.merged.bed   /qgen/home/xuc/VariantCallingPaper/N0015/misc/HG002-multiall-fullcombine.shared.region.noncoding.vcf

# Process output vcfs
for cutoff in 10 15 20 25 30 35 40 45 50 55 60 65 70 75 80 85 90 95 100 125; do
   # create results cutoff directory 
   if [ ! -s result/$caller/$cutoff ]; then
      mkdir result/$caller/$cutoff 
   fi

   # limit existing processed vcf from whole region to noncoding region
   file1=vcf/v6.7.$cutoff.hc.dil.header.vcf
   file2=vcf/$caller.$cutoff.hc.dil.header.vcf
   bedtools intersect -a $file1 -b /qgen/home/xuc/VariantCallingPaper/N0015/misc/shared.region.noncoding.merged.bed -header > $file2
   awk '{if (/^#/) {print} else {if (/^chr/) {print} else {print "chr"$0}}}' $file2 > processed/${file2:4}
   bedtools sort -header -i processed/${file2:4} | bgzip -c > processed/${file2:4}.gz
   tabix -f -p vcf processed/${file2:4}.gz > /dev/null

   # Run RTG evaluation
   file=processed/$caller.$cutoff.hc.dil.header.vcf.gz
   # check if RTG results not already exist
   if [ -s result/$caller/$cutoff/filter.applied ]; then
      echo "RTG tools evaluation results already exist for ${file:10}"
   else
      /qgen/home/xuc/software/RTG/rtg-tools-3.5/rtg vcfeval --all-records -T 1 -b $GroundTruth -c $file --squash-ploidy --no-gzip -t $RefSDF -f "QUAL" -o result/$caller/$cutoff/filter.ignored &>> $runPath/rtg.log 
      /qgen/home/xuc/software/RTG/rtg-tools-3.5/rtg vcfeval --baseline-tp -T 1 -b $GroundTruth -c $file --squash-ploidy --no-gzip -t $RefSDF -f "QUAL" -o result/$caller/$cutoff/filter.applied &>> $runPath/rtg.log
   fi
   
   # create tables by merging info from TP, FP, and FN vcf files
   awk '!/^#/ {OFS="\t"; if ($7 == "PASS") {Pred="TP"} else {Pred="FN"}; print $1":"$2";"$4"/"$5,Pred}' result/$caller/$cutoff/filter.ignored/tp.vcf > tmp/tp.txt 2>>$runPath/rtg.log
   awk '!/^#/ {OFS="\t"; if ($7 == "PASS") {Pred="FP"} else {Pred="TN"}; print $1":"$2";"$4"/"$5,Pred}' result/$caller/$cutoff/filter.ignored/fp.vcf > tmp/fp.txt 2>>$runPath/rtg.log
   awk '!/^#/ {OFS="\t"; print $1":"$2";"$4"/"$5,"FN"}' result/$caller/$cutoff/filter.ignored/fn.vcf > tmp/fn.txt 2>>$runPath/rtg.log
   awk '!/^#/ {OFS="\t"; print $1":"$2";"$4"/"$5,"TP"}' result/$caller/$cutoff/filter.applied/tp-baseline.vcf > result/$caller/$cutoff/tp-baseline.txt 2>>$runPath/rtg.log
   
   # merge tables to one for one cutoff and remove the duplicates (caused by RTG tools limitations/bugs)
   cat tmp/tp.txt tmp/fp.txt tmp/fn.txt | sort | uniq | awk -v duplicate=0 '{OFS="\t"; indb=ind; typeb=type; ind=$1; type=$2; if (ind==indb) {if (typeb ~ /T/) {print indb,typeb} else if (type ~ /T/) {print ind,type}; duplicate=1} else if (duplicate==0) {print indb,typeb} else {duplicate=0}}' | awk 'NF>0' > result/$caller/$cutoff/table.txt 2>>$runPath/rtg.log
done    # end of cutoff loop
   
# integrate tables to one Table for each group -- this is at rpb level, after looping all cutoffs
if [ ! -s result/$caller ]; then
   echo "warning: no RTG tools results found for group: "$caller 
else
   # initialize
   count=1
   OutField="0"
   OutHeader="#Variant"
   
   # folders for each cutoff
   cutOffs=$(ls result/$caller)
   for cutOff in $cutOffs; do
      if [ $count -eq 1 ]; then
         cp result/$caller/${cutOff}/table.txt result/$caller/Table.txt
         cp result/$caller/${cutOff}/tp-baseline.txt result/$caller/Table-baseline.txt			
      else			
         OutField=${OutField}",1."$count
         #
         # join tables from all cutoff values
         join -a1 -a2 -1 1 -2 1 -e "TN" -o ${OutField}",2.2" -t $'\t' <(sort -t $'\t' -k1,1 result/$caller/Table.txt | uniq) <(sort -t $'\t' -k1,1 result/$caller/${cutOff}/table.txt | uniq) > result/$caller/Table.tmp
         mv result/$caller/Table.tmp result/$caller/Table.txt
         #
         # same for baseline
         join -a1 -a2 -1 1 -2 1 -e "FN" -o ${OutField}",2.2" -t $'\t' <(sort -t $'\t' -k1,1 result/$caller/Table-baseline.txt | uniq) <(sort -t $'\t' -k1,1 result/$caller/${cutOff}/tp-baseline.txt | uniq) > result/$caller/Table.tmp
         mv result/$caller/Table.tmp result/$caller/Table-baseline.txt				
      fi
      #
      count=$((count+1))
      OutHeader=${OutHeader}$'\t'"Thr="$cutOff
   done
   echo "$OutHeader" > result/$caller/Header.txt
fi
   
# add header info (cutoff values)
sort result/$caller/Table.txt | uniq > result/$caller/Table.tmp
cat result/$caller/Header.txt result/$caller/Table.tmp > result/$caller/Table.txt

sort result/$caller/Table-baseline.txt | uniq > result/$caller/Table-baseline.tmp
cat result/$caller/Header.txt result/$caller/Table-baseline.tmp > result/$caller/Table-baseline.txt

rm result/$caller/Header.txt result/$caller/Table.tmp result/$caller/Table-baseline.tmp
   
# get summary results
sed '/^#/ d' /qgen/home/xuc/VariantCallingPaper/N0015/GIB.NA12878.hc.dil.noncoding.additional.header.sorted.vcf > /qgen/home/xuc/VariantCallingPaper/N0015/GIB.NA12878.hc.dil.noncoding.additional.noheader.sorted.vcf 
sed -i 's/chr//g' /qgen/home/xuc/VariantCallingPaper/N0015/GIB.NA12878.hc.dil.noncoding.additional.noheader.sorted.vcf
Rscript /qgen/home/xuc/VariantCallingPaper/code/validateWithGIB.withRTG.v6.7.R     $runPath	$caller.VariantList.short.hc.dil.txt    /qgen/home/xuc/VariantCallingPaper/N0015/GIB.NA12878.hc.dil.noncoding.additional.noheader.sorted.vcf  result/$caller/Table.txt result/$caller/Table-baseline.txt        352964  roc.$caller.n0015.2pct.png    summary.$caller.n0015.2pct.csv        details.$caller.n0015.2pct.csv   $caller 1>/dev/null 2>$runPath/rtg.log	

