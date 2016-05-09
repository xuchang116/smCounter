#!/usr/bin/python
# vim: tabstop=8 expandtab shiftwidth=4 softtabstop=4
import os
import sys
import datetime
import subprocess
import mutect2vcf
import shutil
import time
import gzip
import pybedtools

#---------------------------------------------------------------------------------------------
# main function
#----------------------------------------------------------------------------------------------
def runMutect(runPath, bamName, bedName):
   os.chdir(runPath)
   # this is the threshold for proximal gap filter, set at 30% of the variant reads with INDEL within 11-base window
   gapThr = str(int(0.01 * 51764 * 0.3))
   
   # ### run without dbsnp prior and cosmic, --artifact_detection_mode doesn't work well cmp to std and keep everything in out.stats file
   mutectCmd = "/qgen/bin/java/jre1.6.0_45/bin/java -Xmx50g -jar /qgen/home/rvijaya/software/muTect-1.1.4/muTect-1.1.4.jar \
   --analysis_type MuTect \
   --reference_sequence /qgen/home/rvijaya/downloads/alt_hap_masked_ref/ucsc.hg19.fasta \
   --intervals " + bedName + " \
   --input_file:tumor " + bamName + " \
   --out mutect_extended_stats.out \
   --coverage_file mutect_extended.wig.txt \
   --vcf mutect_extended.vcf --enable_extended_output \
   --fraction_contamination 0.003 \
   --minimum_mutation_cell_fraction 0.004 \
   --heavily_clipped_read_fraction 0.75 \
   --min_qscore 20 \
   --num_threads 20 \
   --gap_events_threshold " + gapThr + "\
   --downsample_to_coverage 200000 2>mutect.log 1>/dev/null"

   subprocess.check_call(mutectCmd,shell=True)

   lods = [10, 50, 63, 100, 200, 300, 400, 500, 750, 1000, 1500, 2000, 2500, 3000, 3500, 4000, 4500, 5000]
   for lod in lods:
       lod_num = str(round(lod/10.0, 1))

       mutect2vcf.run("mutect_extended_stats.out", 'MuTect.strict.' + lod_num + ".vcf", tlodThresh=float(lod_num), ignoredbSNP=True, ignoreClusteredPos=True, ignoreStrandArtifact=False, ignoreContamination=False, ignoreFstarLOD=False, ignorePoorMapping=False, nearByGapEvents=False)
       mutect2vcf.run("mutect_extended_stats.out", 'MuTect.partial.' + lod_num + ".vcf", tlodThresh=float(lod_num), ignoredbSNP=True, ignoreClusteredPos=True, ignoreStrandArtifact=False, ignoreContamination=True, ignoreFstarLOD=False, ignorePoorMapping=False, nearByGapEvents=False)
       mutect2vcf.run("mutect_extended_stats.out", 'MuTect.relaxed.' + lod_num + ".vcf", tlodThresh=float(lod_num), ignoredbSNP=True, ignoreClusteredPos=True, ignoreStrandArtifact=False, ignoreContamination=True, ignoreFstarLOD=False, ignorePoorMapping=False, nearByGapEvents=True)


#----------------------------------------------------------------------------------------------
#pythonism to run from the command line
#----------------------------------------------------------------------------------------------
if __name__ == "__main__":
    runPath = sys.argv[1] 
    bamName = sys.argv[2]
    bedName = sys.argv[3]
    runMutect(runPath, bamName, bedName)





