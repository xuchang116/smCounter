#!/usr/bin/python
import string
import sys

#-------------------------------------------------------------------------------
# main function
#-------------------------------------------------------------------------------
def run(inMutectFile, outVCFFile, tlodThresh=-1000, minDepth=1, minAlt=1, ignoredbSNP=False, ignoreClusteredPos =True, ignoreStrandArtifact=False, \
	ignoreContamination=False, ignoreFstarLOD=False, ignorePoorMapping=False, nearByGapEvents=False, ignoreTriallelicSite=False):

   # open output vcf file
   fileOut = open(outVCFFile, "w")
	#write standard header lines
   fileOut.write("##fileformat=VCFv4.1\n")
   fileOut.write("#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\tSAMPLE\n")
  
   filters_to_ignore = []
   filterHash = { "clustered_read_position":ignoreClusteredPos,
			"DBSNP Site":ignoredbSNP,
			"strand_artifact":ignoreStrandArtifact,
         "possible_contamination":ignoreContamination,
         "poor_mapping_region_alternate_allele_mapq":ignorePoorMapping,
         "nearby_gap_events":nearByGapEvents,
			"fstar_tumor_lod":ignoreFstarLOD, 
			"triallelic_site":ignoreTriallelicSite}

   for filterName in filterHash:
      if filterHash[filterName] == True:
         filters_to_ignore.append(filterName)

   header = ''

   filein = open(inMutectFile, "r")
   lineCnt = 0	
   fields = {}
   for line in filein:
      lineCnt += 1
      if lineCnt == 1:
         continue
      elif lineCnt == 2:
         colHeaderElements = line.replace('#','').rstrip('\n').split()
         continue
      cols = line.rstrip().split("\t")
      for i in range(len(colHeaderElements)):
         fields[colHeaderElements[i]]= cols[i]
      m_filter = "PASS"
      fail_reasons = fields["failure_reasons"].split(",")
      if fields["judgement"] != "KEEP":
         for reason in fail_reasons:
             if reason not in filters_to_ignore:
                m_filter = "FAIL"
                break
             
      tumorDP = str(int(fields["t_ref_count"]) + int(fields["t_alt_count"]) + int(fields["t_ins_count"]) + int(fields["t_alt_count"]))
      AD = fields["t_ref_count"]+","+fields["t_alt_count"]
      NAD = fields["n_ref_count"]+","+fields["n_alt_count"]

      # apply filter on t_lod_fstar
      if float(fields["t_lod_fstar"]) < tlodThresh:
         m_filter= "FAIL"

      if (int(fields["t_ref_count"]) + int(fields["t_alt_count"])) < minDepth:
          m_filter= "FAIL"

      if int(fields["t_alt_count"]) < minAlt:
          m_filter= "FAIL"

      # adjust t_lod_fstar
      if float(fields["t_lod_fstar"]) < 0:
         fields["t_lod_fstar"] = 0
      fields["t_lod_fstar"] = str(int(round(float(fields["t_lod_fstar"]), 0))) 

      if m_filter == "PASS":
#         fileOut.write("\t".join((fields["contig"], fields["position"], ".", fields["ref_allele"], fields["alt_allele"], "100", m_filter, "DP="+tumorDP, "AD:DP:NAD", AD+":"+tumorDP+":"+NAD))+"\n")
         fileOut.write("\t".join((fields["contig"], fields["position"], ".", fields["ref_allele"], fields["alt_allele"], fields["t_lod_fstar"], m_filter, "DP="+tumorDP, "AD:DP:NAD", AD+":"+tumorDP+":"+NAD))+"\n")
   filein.close()
   fileOut.close()

#-------------------------------------------------------------------------------
# weird Python-ism for running directly from command line
#-------------------------------------------------------------------------------
if __name__ == "__main__":
   if len(sys.argv) == 3:
      run(sys.argv[1], sys.argv[2])
   else:		
      run(sys.argv[1], sys.argv[2], sys.argv[3])
