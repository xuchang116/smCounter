#!/usr/bin/Rscript
# combine RTGtools results with variant caller output. Use with V6.7
# vim: tabstop=8 expandtab shiftwidth=2 softtabstop=2
# Chang Xu, 12MAR2016 

args <- commandArgs(TRUE)
library(ggplot2)
library(plyr)
setwd(args[1])

# Read in variant caller output data -- short version
var.list <- read.delim(args[2], header=F, colClass='character')
#cutoff <- c(50, 55, 60, 65, 70, 75, 80, 85, 90, 95, 100, 110, 120, 130, 140, 150, 160, 170, 175, 180)
cutoff <- c(10, 15, 20, 25, 30, 35, 40, 45, 50, 55, 60, 65, 70, 75, 80, 85, 90, 95, 100, 125)
cutoff.names <- paste0('FinalCall_', cutoff)
colnames(var.list) <- c('CHROM', 'POS', 'REF', 'ALT', 'TYPE', 'DP', 'FR', 'MT', 'UFR', 'UMT', 'PI', 'VDP', 'VAF', 'VMT', 'VMF', 'VSM', 'VWM', 'DP_A', 'DP_T', 'DP_G', 'DP_C', 'AF_A', 'AF_T', 'AF_G', 'AF_C', 'MT_3RPM', 'MT_5RPM', 'MT_7RPM', 'MT_10RPM',  'UMT_A', 'UMT_T', 'UMT_G', 'UMT_C', 'UMF_A', 'UMF_T', 'UMF_G', 'UMF_C', 'VSM_A', 'VSM_T', 'VSM_G', 'VSM_C', 'VWM_A', 'VWM_T', 'VWM_G', 'VWM_C', 'PI_A', 'PI_T', 'PI_G', 'PI_C', cutoff.names, 'FILTER') 
var.list$mergeID <- paste0('chr', var.list$CHROM, ':', var.list$POS, ';', var.list$REF, '/', var.list$ALT)
sub.var.list <- subset(var.list, select=c(mergeID, CHROM, POS, REF, ALT, TYPE, DP, FR, UFR, MT, UMT, PI, VDP, VAF, VMT, VMF, VSM, VWM, FILTER))

# read in GIB ground truth vcf
GIBvcf <- read.delim(args[3], header=F, colClass='character')
colnames(GIBvcf) <- c('CHROM', 'POS', 'ID', 'REF', 'ALT', 'Qual', 'FILTER_GIB', 'Info', 'Var1', 'Var2')
GIBvcf$VarType  <- ifelse(nchar(GIBvcf$REF)==1 & nchar(GIBvcf$ALT)==1, 'SNP', ifelse(nchar(GIBvcf$REF) > 1 & nchar(GIBvcf$ALT) > 1 | grepl(',', GIBvcf$ALT) | grepl(',', GIBvcf$REF), 'Other', 'INDEL'))
GIBvcf$mergeID <- paste0('chr', GIBvcf$CHROM, ':', GIBvcf$POS, ';', GIBvcf$REF, '/', GIBvcf$ALT)
n.snp <- sum(GIBvcf$VarType=='SNP')
n.indel <- sum(GIBvcf$VarType=='INDEL')
n.other <- sum(GIBvcf$VarType=='Other')

# read in results -- variant caller
res.vc <- unique(read.delim(args[4], header=T, colClass='character'))
colnames(res.vc)[1] <- 'mergeID'
# read in results -- baseline
res.base <- unique(read.delim(args[5], header=T, colClass='character'))
colnames(res.base)[1] <- 'mergeID'
# get shared region size (could be diluted)
regionSize <- as.numeric(args[6]) / 1e6

# merge results with variant caller output
dat1 <- merge(sub.var.list, res.vc, by='mergeID', all=T)
# merge baseline results with GIB vcf
dat2 <- merge(GIBvcf, res.base, by='mergeID', all=T)
for(cut in cutoff){
  dat2[,paste0('Thr.', cut)] <- ifelse(is.na(dat2[,paste0('Thr.', cut)]), 'FN', dat2[,paste0('Thr.', cut)])
}

# count true positives, regarding to baseline
TP <- apply(dat2[,(ncol(dat2)-length(cutoff)+1):ncol(dat2)], 2, function(x) {sum(x=='TP', na.rm=T)})
TP.snp <- apply(dat2[dat2$VarType=='SNP',(ncol(dat2)-length(cutoff)+1):ncol(dat2)], 2, function(x) {sum(x=='TP', na.rm=T)})
TP.indel <- apply(dat2[dat2$VarType=='INDEL',(ncol(dat2)-length(cutoff)+1):ncol(dat2)], 2, function(x) {sum(x=='TP', na.rm=T)})

# false negatives = total variants - TP
FN <- nrow(GIBvcf) - TP
FN.snp <- n.snp - TP.snp
FN.indel <- n.indel - TP.indel

# count false positives
FP <- apply(dat1[,(ncol(dat1)-length(cutoff)+1):ncol(dat1)], 2, function(x) {sum(x=='FP', na.rm=T)})
FP.snp <- apply(dat1[dat1$TYPE=='SNP',(ncol(dat1)-length(cutoff)+1):ncol(dat1)], 2, function(x) {sum(x=='FP', na.rm=T)})
FP.indel <- apply(dat1[dat1$TYPE=='INDEL',(ncol(dat1)-length(cutoff)+1):ncol(dat1)], 2, function(x) {sum(x=='FP', na.rm=T)})

thrs <- colnames(dat1)[(ncol(dat1)-length(cutoff)+1):ncol(dat1)]
thr <- as.numeric(substr(thrs, 5, nchar(thrs)))

summary <- data.frame(thr, TP, FP, FN)
colnames(summary) <- c('Cutoff', 'TP', 'FP', 'FN')
summary$Sensitivity <- summary$TP / nrow(GIBvcf)
summary$FP_per_Mbp <- summary$FP / regionSize

# Summary on SNPs
summary.snp <- data.frame(thr, TP.snp, FP.snp, FN.snp)
colnames(summary.snp) <- c('Cutoff', 'TP', 'FP', 'FN')
summary.snp$Sensitivity <- summary.snp$TP / n.snp
summary.snp$FP_per_Mbp <- summary.snp$FP / regionSize

# Summary on INDELs
summary.indel <- data.frame(thr, TP.indel, FP.indel, FN.indel)
colnames(summary.indel) <- c('Cutoff', 'TP', 'FP', 'FN')
summary.indel$Sensitivity <- summary.indel$TP / n.indel
summary.indel$FP_per_Mbp <- summary.indel$FP / regionSize

# make ROC curve
summary$VarType <- 'Overall'
summary.snp$VarType <- 'SNP'
summary.indel$VarType <- 'INDEL'
summary.long <- rbind(summary, summary.snp, summary.indel)
summary.long$Run <- args[10]

summary.long <- subset(summary.long, select=c(Run, Cutoff, VarType, TP, FP, FN, Sensitivity, FP_per_Mbp))
summary.long$VarType <- factor(summary.long$VarType, levels=c('SNP', 'INDEL', 'Overall'))
summary.long <- arrange(summary.long, Run, Cutoff, VarType)

ggplot(summary.long, aes(FP_per_Mbp, Sensitivity, group=VarType, colour=VarType)) + geom_point(size=3) + geom_path(size=1) + xlab('FP/Mbp') + ylab('Sensitivity') + xlim(0,200) + scale_y_continuous(minor_breaks = seq(0 , 1, .025), breaks = seq(0, 1, .05), limits=c(0,1), expand=c(0,0))
ggsave(args[7])

# save outputs
sub.dat1 <- dat1[is.na(dat1$TYPE) | dat1$TYPE == 'SNP' | dat1$TYPE == 'INDEL', ]
dat1.sorted <- arrange(sub.dat1, CHROM, POS)
write.csv(summary.long, args[8], row.names=F)
write.csv(dat1.sorted, args[9], row.names=F)

 


