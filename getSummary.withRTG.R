#!/usr/bin/Rscript
# combine RTGtools results with variant caller output
# vim: tabstop=8 expandtab shiftwidth=2 softtabstop=2
# get summary for 1% variants only
# Chang Xu, 10APR2016

args <- commandArgs(TRUE)
library(ggplot2)
library(plyr)
setwd(args[1])

#cutoff <- c(10, 15, 20, 25, 50, 75, 100, 125, 150, 175, 200, 225, 250, 275, 300, 400)
# read in GIB ground truth vcf
GIBvcf <- read.delim(args[2], header=F, colClass='character')
colnames(GIBvcf) <- c('Chrom', 'Position', 'ID', 'Ref', 'Alt', 'Qual', 'Filter_GIB', 'Info', 'Var1', 'Var2')
GIBvcf$VarType  <- ifelse(nchar(GIBvcf$Ref)==1 & nchar(GIBvcf$Alt)==1, 'SNP', ifelse(nchar(GIBvcf$Ref) > 1 & nchar(GIBvcf$Alt) > 1 | grepl(',', GIBvcf$Alt) | grepl(',', GIBvcf$Ref), 'Other', 'INDEL'))
GIBvcf$mergeID <- paste0('chr', GIBvcf$Chrom, ':', GIBvcf$Position, ';', GIBvcf$Ref, '/', GIBvcf$Alt)
n.snp <- sum(GIBvcf$VarType=='SNP')
n.indel <- sum(GIBvcf$VarType=='INDEL')
n.other <- sum(GIBvcf$VarType=='Other')

# read in background sample GIB vcf. 
#bkgvcf <- read.delim(args[3], header=F, colClass='character')
#colnames(bkgvcf) <- c('Chrom', 'Position', 'ID', 'Ref', 'Alt', 'Qual', 'Filter_GIB', 'Info', 'Var1', 'Var2')
#bkgvcf$VarType  <- ifelse(nchar(bkgvcf$Ref)==1 & nchar(bkgvcf$Alt)==1, 'SNP', ifelse(nchar(bkgvcf$Ref) > 1 & nchar(bkgvcf$Alt) > 1 | grepl(',', bkgvcf$Alt) | grepl(',', bkgvcf$Ref), 'Other', 'INDEL'))
#bkgvcf$mergeID <- paste0('chr', bkgvcf$Chrom, ':', bkgvcf$Position, ';', bkgvcf$Ref, '/', bkgvcf$Alt)

# read in results -- variant caller. NOTE these are both het and homo variants 
res.vc <- unique(read.delim(args[3], header=T, colClass='character'))
colnames(res.vc)[1] <- 'mergeID'
cutoff.names <- colnames(res.vc)[-1]
cutoff <- as.numeric(substr(cutoff.names, 5, nchar(cutoff.names)))

res.vc$GT <- sapply(res.vc$mergeID, function(x) {unlist(strsplit(x, ';'))[2]})
res.vc$Ref <- sapply(res.vc$GT, function(x) {unlist(strsplit(x, '/'))[1]})
res.vc$Alt <- sapply(res.vc$GT, function(x) {unlist(strsplit(x, '/'))[2]})
res.vc$Type <- ifelse(nchar(res.vc$Ref) == 1 & nchar(res.vc$Alt) == 1, 'SNP', ifelse(nchar(res.vc$Ref) > 1 & nchar(res.vc$Alt) > 1, 'COMPLEX', 'INDEL'))

# read in results -- baseline. NOTE these are both het and homo variants 
res.base <- unique(read.delim(args[4], header=T, colClass='character'))
colnames(res.base)[1] <- 'mergeID'
# get shared region size (could be diluted)
regionSize <- as.numeric(args[5]) / 1e6

# merge baseline results with GIB vcf, include only the het variants in GIB vcf 
dat1 <- merge(GIBvcf, res.base, by='mergeID', all.y=F, all.x=T)
for(cut in cutoff.names){
  dat1[,cut] <- ifelse(is.na(dat1[, cut]), 'FN', dat1[,cut])
}

# count true positives, regarding to baseline
TP <- apply(dat1[,(ncol(dat1)-length(cutoff)+1):ncol(dat1)], 2, function(x) {sum(x=='TP', na.rm=T)})
TP.snp <- apply(dat1[dat1$VarType=='SNP',(ncol(dat1)-length(cutoff)+1):ncol(dat1)], 2, function(x) {sum(x=='TP', na.rm=T)})
TP.indel <- apply(dat1[dat1$VarType=='INDEL',(ncol(dat1)-length(cutoff)+1):ncol(dat1)], 2, function(x) {sum(x=='TP', na.rm=T)})

# false negatives = total variants - TP
FN <- nrow(GIBvcf) - TP
FN.snp <- n.snp - TP.snp
FN.indel <- n.indel - TP.indel

# count false positives

FP <- apply(res.vc[,(ncol(res.vc)-length(cutoff)-3):(ncol(res.vc)-4)], 2, function(x) {sum(x=='FP', na.rm=T)})
FP.snp <- apply(res.vc[res.vc$Type=='SNP',(ncol(res.vc)-length(cutoff)-3):(ncol(res.vc)-4)], 2, function(x) {sum(x=='FP', na.rm=T)})
FP.indel <- apply(res.vc[res.vc$Type=='INDEL',(ncol(res.vc)-length(cutoff)-3):(ncol(res.vc)-4)], 2, function(x) {sum(x=='FP', na.rm=T)})

summary <- data.frame(cutoff, TP, FP, FN)
colnames(summary) <- c('Cutoff', 'TP', 'FP', 'FN')
summary$Sensitivity <- summary$TP / nrow(GIBvcf)
summary$FP_per_Mbp <- summary$FP / regionSize

# Summary on SNPs
summary.snp <- data.frame(cutoff, TP.snp, FP.snp, FN.snp)
colnames(summary.snp) <- c('Cutoff', 'TP', 'FP', 'FN')
summary.snp$Sensitivity <- summary.snp$TP / n.snp
summary.snp$FP_per_Mbp <- summary.snp$FP / regionSize

# Summary on INDELs
summary.indel <- data.frame(cutoff, TP.indel, FP.indel, FN.indel)
colnames(summary.indel) <- c('Cutoff', 'TP', 'FP', 'FN')
summary.indel$Sensitivity <- summary.indel$TP / n.indel
summary.indel$FP_per_Mbp <- summary.indel$FP / regionSize

# make ROC curve
summary$VarType <- 'Overall'
summary.snp$VarType <- 'SNP'
summary.indel$VarType <- 'INDEL'
summary.long <- rbind(summary, summary.snp, summary.indel)
summary.long$Run <- args[6] 

summary.long <- subset(summary.long, select=c(Run, Cutoff, VarType, TP, FP, FN, Sensitivity, FP_per_Mbp))
summary.long$VarType <- factor(summary.long$VarType, levels=c('SNP', 'INDEL', 'Overall'))
summary.long <- arrange(summary.long, Run, Cutoff, VarType)

ggplot(summary.long, aes(FP_per_Mbp, Sensitivity, group=VarType, colour=VarType)) + geom_point(size=3) + geom_path(size=1) + xlab('FP/Mbp') + ylab('Sensitivity') + xlim(0,200) + scale_y_continuous(minor_breaks = seq(0 , 1, .025), breaks = seq(0, 1, .05), limits=c(0,1), expand=c(0,0))
ggsave(args[7])

# save outputs
write.csv(summary.long, args[8], row.names=F)

 


