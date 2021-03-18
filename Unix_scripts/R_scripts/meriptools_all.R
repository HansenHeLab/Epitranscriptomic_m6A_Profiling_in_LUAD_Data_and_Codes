
args <- commandArgs(trailingOnly = TRUE)

sample_txt <- args[1];
thread <- as.numeric(args[2]);     ## peak caller needs numeric input

library(RADAR)
library(MeRIPtools)

## sample list to vector of sample
IN <- read.table(sample_txt, as.is = T);
sample <- IN$V1;
sample_group <-IN$V2

## load bin read counts from RADAR output
count <- readRDS("all_sample_bin_readsCount_MeRIP.rds")

####################################
## call peaks start from class MeRIP 
####################################

if(TRUE)
{

#################
## binominla test
#################
if(TRUE){
print("calling peaks by binomial test")
peak_bi <- callPeakBinomial(MeRIP = count, min_counts = 15, peak_cutoff_fdr = 0.05, peak_cutoff_oddRatio = 1, threads = thread)
peak_bi <- reportJointPeak(MeRIPdata = peak_bi, joint_threshold = 2, threads = thread)

## joint peak read counts and normalization 
peak_bi <- jointPeakCount(peak_bi)
peak_bi <- normalizeLibrary(peak_bi)
peak_bi <- adjustExprLevel(peak_bi, adjustBy = "geneSum")
variable(peak_bi)  <- data.frame(group = sample_group)

## differential analysis:  QNB test 
peak_bi <- QNBtest(peak_bi)
saveRDS(peak_bi@test.est, file = "all_sample_binomialTest_peaks_QNB_diff.rds") 

#########
## output
#########
## all results
saveRDS(peak_bi, file="all_sample_binomialTest_peaks.rds")

## joint peaks BED12
write.table(peak_bi@jointPeaks, file ="all_sample_binomialTest_joint_peaks.bed", col.names = F, row.names = F, quote = F, sep = "\t")

## Joint peak count
norm_input <- t( t(peak_bi@jointPeak_input) / peak_bi@sizeFactor$input )
norm_ip <-  peak_bi@norm.jointPeak_ip
adj_ip <- peak_bi@jointPeak_adjExpr
save(norm_input, norm_ip, adj_ip, file = "all_sample_binomialTest_joint_peaks_norm_INPUT_IP_adjIP.Rdata")

## normlized geneSum
saveRDS(peak_bi@geneSum, file = "all_sample_normalized_geneSum.rds")

## metagene by Guitar
pdf("Binomial_jointPeak_metagene.pdf", width = 7, height = 3.5)
MetaGene(peak_bi)
dev.off()

}


##############
## fisher test
############## 
if(TRUE){
print("calling peaks by fisher test")
peak_fi <- callPeakFisher(MeRIP = count, min_counts = 15, peak_cutoff_fdr = 0.05, peak_cutoff_oddRatio = 1, threads = thread)
peak_fi <- reportJointPeak(MeRIPdata = peak_fi, joint_threshold = 2, threads = thread)

## joint peak read counts and normalizations
peak_fi <- jointPeakCount(peak_fi)
peak_fi <- normalizeLibrary(peak_fi)
peak_fi <- adjustExprLevel(peak_fi, adjustBy = "geneSum")
variable(peak_fi)  <- data.frame(group = sample_group)

## differential analysis:  QNB test
peak_fi <- QNBtest(peak_fi)
saveRDS(peak_fi@test.est, file = "all_sample_fisherTest_peaks_QNB_diff.rds")

#########
## output
#########
## all results
saveRDS(peak_fi, file="all_sample_fisherTest_peaks.rds")

## joint peaks BED12
write.table(peak_fi@jointPeaks, file ="all_sample_fisherTest_joint_peaks.bed", col.names = F, row.names = F, quote = F, sep = "\t")

## Joint peak count 
norm_input <- t( t(peak_fi@jointPeak_input) / peak_fi@sizeFactor$input ) 
norm_ip <-  peak_fi@norm.jointPeak_ip
adj_ip <- peak_fi@jointPeak_adjExpr 
save(norm_input, norm_ip, adj_ip, file = "all_sample_fisherTest_joint_peaks_norm_INPUT_IP_adjIP.Rdata")

## metagene by Guitar
pdf("Fisher_jointPeak_metagene.pdf", width = 7, height = 3.5)
MetaGene(peak_fi)
dev.off()
}

}



