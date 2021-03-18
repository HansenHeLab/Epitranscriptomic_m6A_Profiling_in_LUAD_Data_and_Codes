
args <- commandArgs(trailingOnly = TRUE)

gene_file <- args[1];
peak_file <- args[2];
out_name  <- args[3];

gene <- read.table(gene_file)
peak <- read.table(peak_file)

idx <- !is.na(match(peak$V4, gene$V1))


# count number
all_gene_cnt <- length(levels(peak$V4));

query_gene_list <- as.factor(as.character(peak$V4[idx]));
query_gene_cnt <- length(unique(query_gene_list));
query_gene_tt <- as.data.frame(table(query_gene_list));
peaks_in_query <- sum(query_gene_tt[ , 2]);
peaks_per_gene_mean <- mean(query_gene_tt[ , 2]);
peaks_per_gene_mid <- median(query_gene_tt[ , 2]);
peaks_per_gene_max <- max(query_gene_tt[ , 2]);

exon <- peak$V10[idx]
peaks_all <- length(exon);
peaks_1exon_cnt <- sum(exon==1);

# length 
len <- peak$V11[idx]
len_l <- strsplit(as.character(len), ',', fixed=TRUE)
len_v <- vector();
for (i in 1:length(len_l))
{
	len_v[i] = sum(as.numeric(len_l[[i]]));
}
len_mean <- mean(len_v);
len_mid <- median(len_v)
len_min <- min(len_v)
len_max <- max(len_v)

out_v <- c(all_gene_cnt, query_gene_cnt, peaks_all, peaks_1exon_cnt, peaks_per_gene_mean, peaks_per_gene_mid, peaks_per_gene_max, len_mean, len_mid, len_min, len_max);
out_v <- t(out_v);

colnames(out_v) <- c("all_gene_cnt", "query_gene_cnt", "peaks_all", "peaks_1exon_cnt", "peaks_per_gene_mean", "peaks_per_gene_mid", "peaks_per_gene_max", "len_mean", "len_mid", "len_min", "len_max");
rownames(out_v) <- out_name;

write.table(out_v, file=out_name, col.names=TRUE, row.names=TRUE, sep="\t");



