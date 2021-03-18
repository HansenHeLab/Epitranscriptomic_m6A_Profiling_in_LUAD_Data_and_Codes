## convert gene_ID to gene symbol
args <- commandArgs(trailingOnly = TRUE)

input_file <- args[1];
gene_file <- args[2];
out_file <- args[3];
sample <- args[4]

input <-read.table(input_file);
gene_list <- read.table(gene_file)
idx <- match(input$V4, gene_list$V1);
input$V4 <- gene_list$V2[idx];
input$V7 <- sample

write.table(input, file=out_file, col.names=FALSE, row.name=FALSE, quote=FALSE, sep="\t")
