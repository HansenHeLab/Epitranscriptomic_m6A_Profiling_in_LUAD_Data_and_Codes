## run guitar plot for each individual sample m6A summits
args <- commandArgs(trailingOnly = TRUE)
sample_summit <- args[1];                     # for one sample's called peakss
sample_name <-args[2];
ref <- args[3]  			## ref txdb

load(ref)
name_out <- paste(sample_name, "_summit_Guitar_plot", sep = "")

library(Guitar)
library(regioneR)
library(ggplot2)

print(sample_summit)
tmp_bed <- toGRanges(sample_summit)
sum_GR <- split(tmp_bed, f=tmp_bed$score)
names(sum_GR) <- sample_name
head(tmp_bed)
## plot for onesample

GuitarPlot(gfeatures = sum_GR, returnCount = TRUE,  GuitarCoordsFromTxDb = guitar_txdb_hg38, rescaleComponent=T, saveToPDFprefix = name_out)
