#!/bin/bash
# vim: tabstop=9 expandtab shiftwidth=3 softtabstop=3
# call variants with MuTect and VarDict
# Chang Xu, 11APR2016

RefGenome=/qgen/home/rvijaya/downloads/alt_hap_masked_ref/ucsc.hg19.fasta
RefSDF=/qgen/home/nr/ref/ucsc.hg19.sdf
VarDictJava_Path=/qgen/home/nr/tools/VarDictJava
BedFile=/qgen/home/xuc/VariantCallingPaper/N0015/primers.NA12878.3587.target.bed
bamFile0=/qgen/home/jdicarlo/spe/spe29/29/N0015.Next0015-customR1-primer_S1.oligoClip.bam
bamFile=/qgen/home/xuc/VariantCallingPaper/N0015/N0015.Next0015-customR1-primer_S1.oligoClip.RG.bam
runPath=/qgen/home/xuc/VariantCallingPaper/N0015/10pct

# create read groups (RG) 
#echo "Creating read groups..."
#java -Xmx2g -jar /qgen/home/xuc/software/Picard/picard-tools-1.119/AddOrReplaceReadGroups.jar RGLB=lib01 RGPL=illumina RGPU=run RGSM=N0015 I=$bamFile0 O=$bamFile SORT_ORDER=coordinate CREATE_INDEX=TRUE VALIDATION_STRINGENCY=LENIENT > RG.log 2>&1
#samtools index $bamFile

cd $runPath

# ground truth
GroundTruth=/qgen/home/xuc/VariantCallingPaper/N0015/GIB.NA12878.hc.dil.additional.header.sorted.vcf.gz

# create directory to keep processed vcfs
if [ ! -s processed ]; then
   mkdir processed
fi
# create results directory
if [ ! -s result ]; then
   mkdir result 
fi
# create tmp directory
if [ ! -s tmp ]; then
   mkdir tmp 
fi

###############################################
########### MuTect 
###############################################

caller=MuTect

if [ ! -s result/$caller ]; then
   mkdir result/$caller
fi

echo "Calling with MuTect..."
python /qgen/home/xuc/VariantCallingPaper/code/runMutect.n0015.py   $runPath/vcf $bamFile  $BedFile 
            
# Process output vcfs
echo "Evaluating peformance..."
for filter in strict partial relaxed; do
   # create results cutoff directory 
   if [ ! -s result/$caller/$filter ]; then
      mkdir result/$caller/$filter 
   fi
   for cutoff in 1.0 5.0 6.3 10.0 20.0 30.0 40.0 50.0 75.0 100.0 150.0 200.0 250.0 300.0 350.0 400.0 450.0 500.0; do
      # create results cutoff directory 
      if [ ! -s result/$caller/$filter/$cutoff ]; then
         mkdir result/$caller/$filter/$cutoff 
      fi

      python /qgen/home/xuc/VariantCallingPaper/code/HCandDilute.fast.py  $runPath $BedFile vcf/$caller.$filter.$cutoff.vcf vcf/$caller.$filter.$cutoff.hc.dil.vcf   /qgen/home/xuc/VariantCallingPaper/N0015/misc/shared.region.merged.bed /qgen/home/xuc/VariantCallingPaper/N0015/misc/HG002-multiall-fullcombine.shared.region.vcf
      cat /qgen/home/xuc/VariantCallingPaper/N0015/muTectHeader.txt vcf/$caller.$filter.$cutoff.hc.dil.vcf > vcf/$caller.$filter.$cutoff.hc.dil.header.vcf
      # add GT field to MuTect
      sed -i 's/NAD/GT/g' vcf/$caller.$filter.$cutoff.hc.dil.header.vcf
      sed -i 's/0,0/0\/1/g' vcf/$caller.$filter.$cutoff.hc.dil.header.vcf

      file=vcf/$caller.$filter.$cutoff.hc.dil.header.vcf
      awk '{if (/^#/) {print} else {if (/^chr/) {print} else {print "chr"$0}}}' $file > processed/${file:4}
      bedtools sort -header -i processed/${file:4} | bgzip -c > processed/${file:4}.gz
      tabix -f -p vcf processed/${file:4}.gz > /dev/null

      # Run RTG evaluation
      file=processed/$caller.$filter.$cutoff.hc.dil.header.vcf.gz
      # check if RTG results not already exist
      if [ -s result/$caller/$filter/$cutoff/filter.applied ]; then
         echo "RTG tools evaluation results already exist for ${file:10}"
      else
         /qgen/home/xuc/software/RTG/rtg-tools-3.5/rtg vcfeval --all-records -T 1 -b $GroundTruth -c $file --squash-ploidy --no-gzip -t $RefSDF -f "QUAL" -o result/$caller/$filter/$cutoff/filter.ignored &>> rtg.log 
         /qgen/home/xuc/software/RTG/rtg-tools-3.5/rtg vcfeval --baseline-tp -T 1 -b $GroundTruth -c $file --squash-ploidy --no-gzip -t $RefSDF -f "QUAL" -o result/$caller/$filter/$cutoff/filter.applied &>> rtg.log
      fi
      
      # create tables by merging info from TP, FP, and FN vcf files
      awk '!/^#/ {OFS="\t"; if ($7 == "PASS") {Pred="TP"} else {Pred="FN"}; print $1":"$2";"$4"/"$5,Pred}' result/$caller/$filter/$cutoff/filter.ignored/tp.vcf > tmp/tp.txt 2>>rtg.log
      awk '!/^#/ {OFS="\t"; if ($7 == "PASS") {Pred="FP"} else {Pred="TN"}; print $1":"$2";"$4"/"$5,Pred}' result/$caller/$filter/$cutoff/filter.ignored/fp.vcf > tmp/fp.txt 2>>rtg.log
      awk '!/^#/ {OFS="\t"; print $1":"$2";"$4"/"$5,"FN"}' result/$caller/$filter/$cutoff/filter.ignored/fn.vcf > tmp/fn.txt 2>>rtg.log
      awk '!/^#/ {OFS="\t"; print $1":"$2";"$4"/"$5,"TP"}' result/$caller/$filter/$cutoff/filter.applied/tp-baseline.vcf > result/$caller/$filter/$cutoff/tp-baseline.txt 2>>rtg.log
      
      # merge tables to one for one cutoff and remove the duplicates (caused by RTG tools limitations/bugs)
      cat tmp/tp.txt tmp/fp.txt tmp/fn.txt | sort | uniq | awk -v duplicate=0 '{OFS="\t"; indb=ind; typeb=type; ind=$1; type=$2; if (ind==indb) {if (typeb ~ /T/) {print indb,typeb} else if (type ~ /T/) {print ind,type}; duplicate=1} else if (duplicate==0) {print indb,typeb} else {duplicate=0}}' | awk 'NF>0' > result/$caller/$filter/$cutoff/table.txt 2>>rtg.log
   done # end of cutoff loop

   # integrate tables to one Table for each group -- this is at rpb level, after looping all cutoffs
   if [ ! -s result/$caller/$filter ]; then
      echo "warning: no RTG tools results found for group: "$caller $filter
   else
      # initialize
      count=1
      OutField="0"
      OutHeader="#Variant"
      
      # folders for each cutoff
      cutOffs=$(ls result/$caller/$filter)
      for cutOff in $cutOffs; do
         if [ $count -eq 1 ]; then
            cp result/$caller/$filter/$cutOff/table.txt result/$caller/$filter/Table.txt
            cp result/$caller/$filter/$cutOff/tp-baseline.txt result/$caller/$filter/Table-baseline.txt			
         else			
            OutField=${OutField}",1."$count
            #
            # join tables from all cutoff values
            join -a1 -a2 -1 1 -2 1 -e "TN" -o ${OutField}",2.2" -t $'\t' <(sort -t $'\t' -k1,1 result/$caller/$filter/Table.txt | uniq) <(sort -t $'\t' -k1,1 result/$caller/$filter/$cutOff/table.txt | uniq) > result/$caller/$filter/Table.tmp
            mv result/$caller/$filter/Table.tmp result/$caller/$filter/Table.txt
            #
            # same for baseline
            join -a1 -a2 -1 1 -2 1 -e "FN" -o ${OutField}",2.2" -t $'\t' <(sort -t $'\t' -k1,1 result/$caller/$filter/Table-baseline.txt | uniq) <(sort -t $'\t' -k1,1 result/$caller/$filter/$cutOff/tp-baseline.txt | uniq) > result/$caller/$filter/Table.tmp
            mv result/$caller/$filter/Table.tmp result/$caller/$filter/Table-baseline.txt				
         fi
         #
         count=$((count+1))
         OutHeader=${OutHeader}$'\t'"Thr="$cutOff
      done
      echo "$OutHeader" > result/$caller/$filter/Header.txt
   fi
      
   # add header info (cutoff values)
   sort result/$caller/$filter/Table.txt | uniq > result/$caller/$filter/Table.tmp
   cat result/$caller/$filter/Header.txt result/$caller/$filter/Table.tmp > result/$caller/$filter/Table.txt
   sort result/$caller/$filter/Table-baseline.txt | uniq > result/$caller/$filter/Table-baseline.tmp
   cat result/$caller/$filter/Header.txt result/$caller/$filter/Table-baseline.tmp > result/$caller/$filter/Table-baseline.txt
   rm result/$caller/$filter/Header.txt result/$caller/$filter/Table.tmp result/$caller/$filter/Table-baseline.tmp
   rm tmp/*.txt
      
   # get summary results
   Rscript /qgen/home/xuc/VariantCallingPaper/code/getSummary.withRTG.R     $runPath  /qgen/home/xuc/VariantCallingPaper/N0015/GIB.NA12878.hc.dil.additional.noheader.sorted.vcf  result/$caller/$filter/Table.txt result/$caller/$filter/Table-baseline.txt        395303  $caller:$filter   roc.$caller.$filter.n0015.10pct.png    summary.$caller.$filter.n0015.10pct.csv   1>/dev/null 2>>rtg.log	

done # end of filter loop


###############################################
########### VarDict 
###############################################
alleleThrLow=(0.002 0.004 0.006 0.008 0.01 0.02)
caller=VarDict 
alleleThr=${alleleThrLow[*]}
for cutoff in $alleleThr; do
#echo "Calling with VarDict at cutoff="$cutoff"..."
   if [ ! -s result/$caller ]; then
      mkdir result/$caller
   fi
   # create results cutoff directory 
   if [ ! -s result/$caller/$cutoff ]; then
      mkdir result/$caller/$cutoff 
   fi

   # run VarDict JAVA version.
   ${VarDictJava_Path}/build/install/VarDict/bin/VarDict -th 20 -F 0 -q 20 -G $RefGenome -f $cutoff -N N0015 -b $bamFile  -z -c 1 -S 2 -E 3 $runPath/$BedFile | ${VarDictJava_Path}/VarDict/teststrandbias.R | ${VarDictJava_Path}/VarDict/var2vcf_valid.pl -N N0015 -E -f $cutoff > vcf/$caller.$cutoff.vcf 2> error.log

# Process output vcfs
echo "evaluating VCF with RTG..."
   python /qgen/home/xuc/VariantCallingPaper/code/HCandDilute.fast.py  $runPath $BedFile vcf/$caller.$cutoff.vcf vcf/$caller.$cutoff.hc.dil.vcf   /qgen/home/xuc/VariantCallingPaper/N0015/misc/shared.region.merged.bed /qgen/home/xuc/VariantCallingPaper/N0015/misc/HG002-multiall-fullcombine.shared.region.vcf
   sort -u -k1,4 vcf/$caller.$cutoff.hc.dil.vcf -o vcf/$caller.$cutoff.hc.dil.vcf
   cat /qgen/home/xuc/VariantCallingPaper/N0015/VarDictHeader.txt vcf/$caller.$cutoff.hc.dil.vcf > vcf/$caller.$cutoff.hc.dil.header.vcf


   file=vcf/$caller.$cutoff.hc.dil.header.vcf
   awk '{if (/^#/) {print} else {if (/^chr/) {print} else {print "chr"$0}}}' $file > processed/${file:4}
   bedtools sort -header -i processed/${file:4} | bgzip -c > processed/${file:4}.gz
   tabix -f -p vcf processed/${file:4}.gz > /dev/null

   # Run RTG evaluation
   file=processed/$caller.$cutoff.hc.dil.header.vcf.gz
   # check if RTG results not already exist
   if [ -s result/$caller/$cutoff/filter.applied ]; then
      echo "RTG tools evaluation results already exist for ${file:10}"
   else
      /qgen/home/xuc/software/RTG/rtg-tools-3.5/rtg vcfeval --all-records -T 1 -b $GroundTruth -c $file --squash-ploidy --no-gzip -t $RefSDF -f "QUAL" -o result/$caller/$cutoff/filter.ignored &>> rtg.log 
      /qgen/home/xuc/software/RTG/rtg-tools-3.5/rtg vcfeval --baseline-tp -T 1 -b $GroundTruth -c $file --squash-ploidy --no-gzip -t $RefSDF -f "QUAL" -o result/$caller/$cutoff/filter.applied &>> rtg.log
   fi
   
   # create tables by merging info from TP, FP, and FN vcf files
   awk '!/^#/ {OFS="\t"; if ($7 == "PASS") {Pred="TP"} else {Pred="FN"}; print $1":"$2";"$4"/"$5,Pred}' result/$caller/$cutoff/filter.ignored/tp.vcf > tmp/tp.txt 2>>rtg.log
   awk '!/^#/ {OFS="\t"; if ($7 == "PASS") {Pred="FP"} else {Pred="TN"}; print $1":"$2";"$4"/"$5,Pred}' result/$caller/$cutoff/filter.ignored/fp.vcf > tmp/fp.txt 2>>rtg.log
   awk '!/^#/ {OFS="\t"; print $1":"$2";"$4"/"$5,"FN"}' result/$caller/$cutoff/filter.ignored/fn.vcf > tmp/fn.txt 2>>rtg.log
   awk '!/^#/ {OFS="\t"; print $1":"$2";"$4"/"$5,"TP"}' result/$caller/$cutoff/filter.applied/tp-baseline.vcf > result/$caller/$cutoff/tp-baseline.txt 2>>rtg.log
   
   # merge tables to one for one cutoff and remove the duplicates (caused by RTG tools limitations/bugs)
   cat tmp/tp.txt tmp/fp.txt tmp/fn.txt | sort | uniq | awk -v duplicate=0 '{OFS="\t"; indb=ind; typeb=type; ind=$1; type=$2; if (ind==indb) {if (typeb ~ /T/) {print indb,typeb} else if (type ~ /T/) {print ind,type}; duplicate=1} else if (duplicate==0) {print indb,typeb} else {duplicate=0}}' | awk 'NF>0' > result/$caller/$cutoff/table.txt 2>>/qgen/home/xuc/N0015/ds.mt.all/rtg.log
done # end of cutoff loop

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
         cp result/$caller/$cutOff/table.txt result/$caller/Table.txt
         cp result/$caller/$cutOff/tp-baseline.txt result/$caller/Table-baseline.txt			
      else			
         OutField=${OutField}",1."$count
         #
         # join tables from all cutoff values
         join -a1 -a2 -1 1 -2 1 -e "TN" -o ${OutField}",2.2" -t $'\t' <(sort -t $'\t' -k1,1 result/$caller/Table.txt | uniq) <(sort -t $'\t' -k1,1 result/$caller/$cutOff/table.txt | uniq) > result/$caller/Table.tmp
         mv result/$caller/Table.tmp result/$caller/Table.txt
         #
         # same for baseline
         join -a1 -a2 -1 1 -2 1 -e "FN" -o ${OutField}",2.2" -t $'\t' <(sort -t $'\t' -k1,1 result/$caller/Table-baseline.txt | uniq) <(sort -t $'\t' -k1,1 result/$caller/$cutOff/tp-baseline.txt | uniq) > result/$caller/Table.tmp
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
rm tmp/*.txt
   
# get summary results
Rscript /qgen/home/xuc/VariantCallingPaper/code/getSummary.withRTG.R     $runPath  /qgen/home/xuc/VariantCallingPaper/N0015/GIB.NA12878.hc.dil.additional.noheader.sorted.vcf  result/$caller/Table.txt result/$caller/Table-baseline.txt        395303  $caller   roc.$caller.n0015.10pct.png    summary.$caller.n0015.10pct.csv   1>/dev/null 2>>rtg.log	
