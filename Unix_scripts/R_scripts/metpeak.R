## Run MeTPeak
args <- commandArgs(trailingOnly = TRUE)

gtf <- args[1];
ip <- args[2];
input <-args[3];
name <-args[4];

library(MeTPeak)
metpeak(GENE_ANNO_GTF = gtf, IP_BAM = ip, INPUT_BAM = input, EXPERIMENT_NAME = name )



