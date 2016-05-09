#!/bin/bash
# vim: tabstop=9 expandtab shiftwidth=3 softtabstop=3
# benchmark using RTG tools for 1% (heterozygous) variants only! 
# Chang Xu, 28APR2016


RefSDF=/qgen/home/nr/ref/ucsc.hg19.sdf
BedFile=primers.NA12878-194-genes-63-indels.10867.coding.bed
runPath=/qgen/home/xuc/VariantCallingPaper/N0030
caller=v6.7.1.ds.coding

# Process GIB vcf
GroundTruth=$runPath/GIB.NA12878.hc.dil.coding.header.sorted.vcf.gz  
GroundTruth_het_vcf=$runPath/GIB.NA12878.hc.dil.het.coding.noheader.sorted.vcf

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

for mt_pct in 80; do 
#for mt_pct in 10 20 40 60 100; do 
   # create MT percent directory
   if [ ! -s result/$caller/$mt_pct ]; then
      mkdir result/$caller/$mt_pct 
   fi

   for ds_rpb in 1.1 1.5 2.0 4.0 6.0 8.6; do

      # skip the full data 
      if [ $mt_pct -eq 100 ] && [ $ds_rpb -eq 8.6 ]; then
         continue
      fi

      # create results rpb directory 
      if [ ! -s result/$caller/$mt_pct/$ds_rpb ]; then
         mkdir result/$caller/$mt_pct/$ds_rpb 
      fi

      echo "----------------------"$mt_pct"% MT, "$ds_rpb "rpb, evaluating VCF with RTG-------------------------"

      # limit variants in GIB HC region and dilute NA24385 variants
      python /qgen/home/xuc/VariantCallingPaper/code/HCandDilute.fast.alt.py  $runPath	$BedFile	intermediate.output/$caller.$mt_pct.$ds_rpb.VariantList.short.txt   intermediate.output/$caller.$mt_pct.$ds_rpb.VariantList.short.hc.dil.txt  shared.region.coding.merged.bed HG002-multiall-fullcombine.additional.n0030.shared.region.coding.vcf

      # Process output vcfs
      for cutoff in 10 15 20 25 30 35 40 45 50 55 60 65 70 75 80 85 90 95 100 125; do
         # create results cutoff directory 
         if [ ! -s result/$caller/$mt_pct/$ds_rpb/$cutoff ]; then
            mkdir result/$caller/$mt_pct/$ds_rpb/$cutoff 
         fi

         python /qgen/home/xuc/VariantCallingPaper/code/HCandDilute.fast.py  $runPath $BedFile vcf/$caller.$mt_pct.$ds_rpb.$cutoff.vcf vcf/$caller.$mt_pct.$ds_rpb.$cutoff.hc.dil.vcf   shared.region.coding.merged.bed HG002-multiall-fullcombine.additional.n0030.shared.region.coding.vcf
         cat VCFheader.txt vcf/$caller.$mt_pct.$ds_rpb.$cutoff.hc.dil.vcf > vcf/$caller.$mt_pct.$ds_rpb.$cutoff.hc.dil.header.vcf 

         file=vcf/$caller.$mt_pct.$ds_rpb.$cutoff.hc.dil.header.vcf
         awk '{if (/^#/) {print} else {if (/^chr/) {print} else {print "chr"$0}}}' $file > processed/${file:4}
         bedtools sort -header -i processed/${file:4} | bgzip -c > processed/${file:4}.gz
         tabix -f -p vcf processed/${file:4}.gz > /dev/null

         # Run RTG evaluation
         file=processed/$caller.$mt_pct.$ds_rpb.$cutoff.hc.dil.header.vcf.gz
         # check if RTG results not already exist
         if [ -s result/$caller/$mt_pct/$ds_rpb/$cutoff/filter.applied ]; then
            echo "RTG tools evaluation results already exist for ${file:10}"
         else
            /qgen/home/xuc/software/RTG/rtg-tools-3.5/rtg vcfeval --all-records -T 10 -b $GroundTruth -c $file --squash-ploidy --no-gzip -t $RefSDF -f "QUAL" -o result/$caller/$mt_pct/$ds_rpb/$cutoff/filter.ignored &>> $runPath/rtg.log 
            /qgen/home/xuc/software/RTG/rtg-tools-3.5/rtg vcfeval --baseline-tp -T 10 -b $GroundTruth -c $file --squash-ploidy --no-gzip -t $RefSDF -f "QUAL" -o result/$caller/$mt_pct/$ds_rpb/$cutoff/filter.applied &>> $runPath/rtg.log
         fi

         # create tables by merging info from TP, FP, and FN vcf files
         awk '!/^#/ {OFS="\t"; if ($7 == "PASS") {Pred="TP"} else {Pred="FN"}; print $1":"$2";"$4"/"$5,Pred}' result/$caller/$mt_pct/$ds_rpb/$cutoff/filter.ignored/tp.vcf > tmp/tp.txt 2>>$runPath/rtg.log
         awk '!/^#/ {OFS="\t"; if ($7 == "PASS") {Pred="FP"} else {Pred="TN"}; print $1":"$2";"$4"/"$5,Pred}' result/$caller/$mt_pct/$ds_rpb/$cutoff/filter.ignored/fp.vcf > tmp/fp.txt 2>>$runPath/rtg.log
         awk '!/^#/ {OFS="\t"; print $1":"$2";"$4"/"$5,"FN"}' result/$caller/$mt_pct/$ds_rpb/$cutoff/filter.ignored/fn.vcf > tmp/fn.txt 2>>$runPath/rtg.log
         awk '!/^#/ {OFS="\t"; print $1":"$2";"$4"/"$5,"TP"}' result/$caller/$mt_pct/$ds_rpb/$cutoff/filter.applied/tp-baseline.vcf > result/$caller/$mt_pct/$ds_rpb/$cutoff/tp-baseline.txt 2>>$runPath/rtg.log
         
         # merge tables to one for one cutoff and remove the duplicates (caused by RTG tools limitations/bugs)
         cat tmp/tp.txt tmp/fp.txt tmp/fn.txt | sort | uniq | awk -v duplicate=0 '{OFS="\t"; indb=ind; typeb=type; ind=$1; type=$2; if (ind==indb) {if (typeb ~ /T/) {print indb,typeb} else if (type ~ /T/) {print ind,type}; duplicate=1} else if (duplicate==0) {print indb,typeb} else {duplicate=0}}' | awk 'NF>0' > result/$caller/$mt_pct/$ds_rpb/$cutoff/table.txt 2>>$runPath/rtg.log
      done    # end of cutoff loop
         
      # integrate tables to one Table for each group -- this is at rpb level, after looping all cutoffs
      if [ ! -s result/$caller/$mt_pct/$ds_rpb ]; then
         echo "warning: no RTG tools results found for group: "$caller $mt_pct $ds_rpb
      else
         # initialize
         count=1
         OutField="0"
         OutHeader="#Variant"
         
         # folders for each cutoff
         cutOffs=$(ls result/$caller/$mt_pct/$ds_rpb)
         for cutOff in $cutOffs; do
            if [ $count -eq 1 ]; then
               cp result/$caller/$mt_pct/$ds_rpb/${cutOff}/table.txt result/$caller/$mt_pct/$ds_rpb/Table.txt
               cp result/$caller/$mt_pct/$ds_rpb/${cutOff}/tp-baseline.txt result/$caller/$mt_pct/$ds_rpb/Table-baseline.txt			
            else			
               OutField=${OutField}",1."$count
               #
               # join tables from all cutoff values
               join -a1 -a2 -1 1 -2 1 -e "TN" -o ${OutField}",2.2" -t $'\t' <(sort -t $'\t' -k1,1 result/$caller/$mt_pct/$ds_rpb/Table.txt | uniq) <(sort -t $'\t' -k1,1 result/$caller/$mt_pct/$ds_rpb/${cutOff}/table.txt | uniq) > result/$caller/$mt_pct/$ds_rpb/Table.tmp
               mv result/$caller/$mt_pct/$ds_rpb/Table.tmp result/$caller/$mt_pct/$ds_rpb/Table.txt
               #
               # same for baseline
               join -a1 -a2 -1 1 -2 1 -e "FN" -o ${OutField}",2.2" -t $'\t' <(sort -t $'\t' -k1,1 result/$caller/$mt_pct/$ds_rpb/Table-baseline.txt | uniq) <(sort -t $'\t' -k1,1 result/$caller/$mt_pct/$ds_rpb/${cutOff}/tp-baseline.txt | uniq) > result/$caller/$mt_pct/$ds_rpb/Table.tmp
               mv result/$caller/$mt_pct/$ds_rpb/Table.tmp result/$caller/$mt_pct/$ds_rpb/Table-baseline.txt				
            fi
            #
            count=$((count+1))
            OutHeader=${OutHeader}$'\t'"Thr="$cutOff
         done
         echo "$OutHeader" > result/$caller/$mt_pct/$ds_rpb/Header.txt
      fi
         
      # add header info (cutoff values)
      sort result/$caller/$mt_pct/$ds_rpb/Table.txt | uniq > result/$caller/$mt_pct/$ds_rpb/Table.tmp
      cat result/$caller/$mt_pct/$ds_rpb/Header.txt result/$caller/$mt_pct/$ds_rpb/Table.tmp > result/$caller/$mt_pct/$ds_rpb/Table.txt
      
      sort result/$caller/$mt_pct/$ds_rpb/Table-baseline.txt | uniq > result/$caller/$mt_pct/$ds_rpb/Table-baseline.tmp
      cat result/$caller/$mt_pct/$ds_rpb/Header.txt result/$caller/$mt_pct/$ds_rpb/Table-baseline.tmp > result/$caller/$mt_pct/$ds_rpb/Table-baseline.txt
      
      rm result/$caller/$mt_pct/$ds_rpb/Header.txt result/$caller/$mt_pct/$ds_rpb/Table.tmp result/$caller/$mt_pct/$ds_rpb/Table-baseline.tmp

      # get summary results
      Rscript /qgen/home/xuc/VariantCallingPaper/code/getSummary.withRTG.het.R     $runPath   $GroundTruth_het_vcf  result/$caller/$mt_pct/$ds_rpb/Table.txt result/$caller/$mt_pct/$ds_rpb/Table-baseline.txt        681980  $caller.het.$mt_pct.$ds_rpb   roc.$caller.het.$mt_pct.$ds_rpb.n0030.png    summary.$caller.het.$mt_pct.$ds_rpb.n0030.csv  1>/dev/null 2>$runPath/rtg.log	
         
      #Rscript /qgen/home/xuc/VariantCallingPaper/code/validateWithGIB.withRTG.v6.7.R   $runPath $caller.$mt_pct.$ds_rpb.VariantList.short.hc.dil.txt  GIB.NA12878.hc.dil.het.coding.noheader.sorted.vcf result/$caller/$mt_pct/$ds_rpb/Table.txt result/$caller/$mt_pct/$ds_rpb/Table-baseline.txt   681980  roc.$caller.$mt_pct.$ds_rpb.n0030.png    summary.$caller.$mt_pct.$ds_rpb.n0030.csv        details.$caller.$mt_pct.$ds_rpb.n0030.csv   $caller:$mt_pct:$ds_rpb 1>/dev/null 2>$runPath/rtg.log	
   done
   rm tmp/*.*
done



