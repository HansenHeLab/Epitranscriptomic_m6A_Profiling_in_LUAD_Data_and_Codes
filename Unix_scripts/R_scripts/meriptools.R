
args <- commandArgs(trailingOnly = TRUE)

sample <- args[1];
GTF <- args[2];
bam_dir <- args[3];
out_dir <- args[4];
strand <- args[5];
thread <- args[6];

library(MeRIPtools)

if(FALSE)
{
## Bin read counts
MeRIP <- countReads(samplenames = sample,
                    gtf = GTF,
                    bamFolder = bam_dir,
                    outputDir = out_dir,
                    modification = "m6A",
                    strandToKeep = strand,
                    fragmentLength = 200,
                    binSize = 50,
                    threads = thread,
                    saveOutput = F,
                    )

saveRDS(MeRIP, file=paste(sample, "_bin_readsCount.rds", sep = ""));
}

MeRIP <- readRDS(paste(sample, "_bin_readsCount.rds", sep = ""))

if(TRUE)
{
###############################
## call peak 
## binominla test
print("calling peaks by binominal test")
peak_bi <- callPeakBinomial(MeRIP = MeRIP, min_counts = 15, peak_cutoff_fdr = 0.05, peak_cutoff_oddRatio = 1, threads = thread)
#peak_bi <- reportJointPeak(MeRIPdata = peak_bi, joint_threshold = 2, threads = 10)
#peak_bi <- jointPeakCount(peak_bi)

saveRDS(peak_bi, file=paste(sample, "_binominal_test_peaks.rds", sep = ""))


## fisher test 
print("calling peaks by fisher test")
peak_fi <- callPeakFisher(MeRIP = MeRIP, min_counts = 15, peak_cutoff_fdr = 0.05, peak_cutoff_oddRatio = 1, threads = thread)
#peak_fi <- reportJointPeak(MeRIPdata = peak_fi, joint_threshold = 2, threads = 10)
#peak_fi <- jointPeakCount(peak_fi)

saveRDS(peak_fi, file=paste(sample, "_fisher_test_peaks.rds", sep = ""))

}



