
args <- commandArgs(trailingOnly = TRUE)

sample_txt <- args[1];
GTF <- args[2];
bam_dir <- args[3];
out_dir <- args[4];
strand <- args[5];
frag <- as.numeric(args[6])         ## fragment length 
bin <- as.numeric(args[7])	   ## bin size
thread <- as.numeric(args[8]);     ## peak caller needs numeric input

library(RADAR)

## sample list to vector of sample
IN <- read.table(sample_txt, as.is = T);
sample <- IN$V1;
sample_group <-IN$V2

###########################################
## Bin-wise reads count: calss MeRIP: count
###########################################

if(TRUE){
count <- countReads(samplenames = sample,
                    gtf = GTF,
                    bamFolder = bam_dir,
                    outputDir = out_dir,
                    modification = "m6A",
                    strandToKeep = strand,
                    fragmentLength = frag,
                    binSize = bin,
                    threads = thread,
                    saveOutput = F,
                    )

saveRDS(count, file="all_sample_bin_readsCount_MeRIP.rds");
getSlots("MeRIP")

###################################
###### MeRIP.RADAR class: count_adj
###################################

## INPUT/IP normalization and IP adjusment without filtering 
count_adj <- count

## nomalization
pdf("Normalization_before_after_IP_INPUT_boxplot.pdf")
count_adj <- normalizeLibrary(count_adj, boxPlot = T) 
dev.off()

## IP adjustment and filtering
count_adj <- adjustExprLevel(count_adj, adjustBy = "geneSum")		  ## adjusment of IP to measure m6A level 		
variable(count_adj)  <- data.frame(group = sample_group)
count_adj <- filterBins(count_adj ,minCountsCutOff = 15)

## for diff analysis
count_adj <- diffIP_parallel(count_adj, thread = thread)
count_adj <- reportResult(count_adj, cutoff = 0.1, Beta_cutoff = 0.5 , threads = thread)

pdf("DMRs_heatmap.pdf")
plotHeatMap(count_adj)
dev.off()

######################
## geneplot 
count_adj <- PrepCoveragePlot(count_adj)

## eg: BLVRA whole gene
pdf("BLVRA_signal.pdf")
plotGeneCov(count_adj, geneName = "ENSG00000106605.10", center = median, libraryType = "same")
dev.off()

pdf("BLVRA_3UTR_signal.pdf")
plotGeneCov(count_adj, geneName = "ENSG00000106605.10", center = mean, libraryType = "same",ZoomIn = c(43806969, 43807365 ), adjustExprLevel = F )
dev.off()


pdf("BLVRA_3UTR_adj_signal.pdf")
plotGeneCov(count_adj, geneName = "ENSG00000106605.10", center = mean, libraryType = "same",ZoomIn = c(43806969, 43807365 ), adjustExprLevel = T )
dev.off()


## PCA plots based on top 10000 normalized IP and INPUT
ip_top_bins <- extractIP(count_adj, filtered = T)[order(rowMeans( extractIP(count_adj,filtered = T) ),decreasing = T)[1:10000],]
input_top_bins <- extractInput(count_adj)[order(rowMeans( extractInput(count_adj) ),decreasing = T)[1:10000],]

pdf("Normalization_IP_top_10000_PCA.pdf")
plotPCAfromMatrix(ip_top_bins,group = unlist(variable(count_adj))) 
dev.off()

pdf("Normalization_Input_top_10000_PCA.pdf")
plotPCAfromMatrix(input_top_bins,group = unlist(variable(count_adj))) 
dev.off()

save(ip_top_bins, input_top_bins, file = "Top_10000_IP_Input_normalizaed_bin.Rdata")

## saveoutput 
saveRDS(count_adj@test.est, file = "all_sample_RADAR_diff_all_bins.rds")      ## all bin-based diff results
saveRDS(count_adj@mergedBins, file = "all_sample_RADAR_sig_diff_FDR-0.1_log2FC-0.5.rds")      ## significant diff bins (merged to peaks) 
saveRDS(count_adj, file="all_sample_bin_readsCount_MeRIP.RADAR.rds");
getSlots("MeRIP.RADAR")

}


