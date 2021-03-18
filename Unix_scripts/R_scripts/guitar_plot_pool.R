## run guitar plot for pooled samples' m6A summits
args <- commandArgs(trailingOnly = TRUE)
sample_list <- args[1];                        ## sample list for all samples' summit bed files
name_out <-args[2];
ref <- args[3]

library(Guitar)
library(regioneR)
load(ref)

## Bed6 summit to Granges list
sample_l <- read.table(sample_list, as.is = T)
L <- length(sample_l$V1)
bed_l <- list()
count_l <- list()

for (i in 1:L)
{
	tmp_bed <- toGRanges(sample_l$V1[i])
	#tmp_bed <- toGRanges(paste(sample_l$V1[i], "_peaks_summit.bed", sep = ""))
	tmp <- split(tmp_bed, f=tmp_bed$score)
	bed_l[[i]]  <- tmp
}

names(bed_l) <- sample_l$V1
count <- GuitarPlot(gfeatures = bed_l, returnCount = TRUE,  GuitarCoordsFromTxDb = guitar_txdb_hg38, rescaleComponent=T, saveToPDFprefix = name_out)
save(count, file = paste(name_out, ".Rdata", sep = ""))

## Bed12 to Granges list
#tmp1 <- BED12toGRangesList("CPC0344_ip_peaks_sorted_named.bed", header = FALSE)
#tmp2 <- BED12toGRangesList("CPC0267_ip_peaks_sorted_named.bed", header = FALSE)

