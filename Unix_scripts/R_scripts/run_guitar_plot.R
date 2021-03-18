rm(list = ls())

## Rscript --vanilla /cluster/projects/hansengroup/Yong/PCa_m6a/script/MeRIP-SAP/run_guitar_plot.R $wd $sample_list $name_out

##  wd : for working directory 
args <- commandArgs(trailingOnly = TRUE)
wd <- args[1];
sample_list <- args[2];
name_out <-args[3];

setwd(wd)

library(Guitar)
library(regioneR)

load("/cluster/projects/hansengroup/Yong/DB/hg38_ek12/Guitar_TxDb_hg38.Rdata")


## Bed6 summit to Granges list
sample_l <- read.table(sample_list, as.is = T)
L <- length(sample_l)
bed_l <- list()

for (i in 1:L)
{
	tmp_bed <- toGRanges(paste(sample_l[i], "_peaks_summit.bed", sep = ""))
	tmp <- split(tmp_bed, f=tmp_bed$socre)
	bed_l[[i]]  <- tmp	 
}

names(bed_l) <- sample_l

GuitarPlot(gfeatures = bed_l, GuitarCoordsFromTxDb = guitar_txdb_hg38, rescaleComponent=T, saveToPDFprefix = name_out) 



## Bed12 to Granges list
#tmp1 <- BED12toGRangesList("CPC0344_ip_peaks_sorted_named.bed", header = FALSE)
#tmp2 <- BED12toGRangesList("CPC0267_ip_peaks_sorted_named.bed", header = FALSE)

