
## Rscript --vanilla /cluster/projects/hansengroup/Yong/PCa_m6a/script/MeRIP-SAP/run_guitar_plot.R $wd $sample_list $name_out
##  wd : for working directory 

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
	tmp <- split(tmp_bed, f=tmp_bed$circ_id)
	bed_l[[i]]  <- tmp 
        
	 
}

names(bed_l) <- c("ALL_3p", "ALL_5p", "Sahred_3p", "Shared_5p", "INPUT_3p", "INPUT_5p", "IP_3p", "IP_5p")
# names(bed_l) <- sample_l$V1

count <- GuitarPlot(gfeatures = bed_l, returnCount = TRUE,  GuitarCoordsFromTxDb = guitar_txdb_hg38, rescaleComponent=T, saveToPDFprefix = name_out) 

save(count, file = paste(name_out, ".Rdata", sep = ""))


## Bed12 to Granges list
#tmp1 <- BED12toGRangesList("CPC0344_ip_peaks_sorted_named.bed", header = FALSE)
#tmp2 <- BED12toGRangesList("CPC0267_ip_peaks_sorted_named.bed", header = FALSE)

