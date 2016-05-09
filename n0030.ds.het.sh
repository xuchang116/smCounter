#!/bin/bash
# vim: tabstop=9 expandtab shiftwidth=3 softtabstop=3
# 1% variants only
# Chang Xu, 09APR2016 

runPath=/qgen/home/xuc/VariantCallingPaper/N0030  
cd $runPath
bedFile=primers.NA12878-194-genes-63-indels.10867.coding.bed
GroundTruth=$runPath/GIB.NA12878.hc.dil.het.coding.noheader.sorted.vcf

#########################################################
### coding region 
#########################################################
#for mt_pct in 20 40 60; do 
for mt_pct in 60; do 
   #for ds_rpb in 1.1 1.5 2.0 4.0 6.0 8.6; do 
   for ds_rpb in 8.6; do 

      echo $mt_pct $ds_rpb
      caller=v6.7.1.coding.het.$mt_pct.$ds_rpb
      Rscript /qgen/home/xuc/VariantCallingPaper/code/getSummary.withRTG.het.R     $runPath   $GroundTruth  result/v6.7.1.ds.coding/$mt_pct/$ds_rpb/Table.txt result/v6.7.1.ds.coding/$mt_pct/$ds_rpb/Table-baseline.txt        681980  $caller   roc.$caller.n0030.png    summary.$caller.n0030.csv  1>/dev/null 2>$runPath/rtg.log	

   done
done

