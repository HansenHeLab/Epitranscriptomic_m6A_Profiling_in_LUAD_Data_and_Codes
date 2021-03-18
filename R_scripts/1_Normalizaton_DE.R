##########################################################################
## This script implements
## 1. IP and INPUT raw read counts normalizaiton
## 2. PCA plots before and after normalization for both IP and INPUT
## 3. MA plots between normal and tumor for both IP and INPUT
## 4. Density plots of RPKM for both IP and INPUT
## 5. Differential expression analysis between normal and tumor
## 6. TCGA_LUAD Differential expression analysis between normal and tumor
##########################################################################

rm(list=ls())
setwd("/working/dir/")

##################
## load raw counts
##################
{
# for all input
IN <- read.table("all_input.count")
N <- ncol(IN)
IN <- IN[, -seq(3, N , 2)]
rownames(IN) <- IN[, 1]
IN <- as.matrix(IN[, -1])
sample <- read.table("input_sample.list")
colnames(IN) <- sample$V1
input_count <- IN       # tumor46 & tumor 29 with replicates

IN <- read.table("all_input_ek12.count")
N <- ncol(IN)
IN <- IN[, -seq(3, N , 2)]
rownames(IN) <- IN[, 1]
IN <- as.matrix(IN[, -1])
sample <- read.table("input_sample_ek12.list")
colnames(IN) <- sample$V1
input_ek12_count <- IN       # tumor46 & tumor 29 with replicates

# for all ip
IN <- read.table("all_ip.count")
N <- ncol(IN)
IN <- IN[, -seq(3, N , 2)]
rownames(IN) <- IN[, 1]
IN <- as.matrix(IN[, -1])
sample <- read.table("ip_sample.list")
colnames(IN) <- sample$V1
ip_count <- IN       # normal11 with ip only

IN <- read.table("all_ip_ek12.count")
N <- ncol(IN)
IN <- IN[, -seq(3, N , 2)]
rownames(IN) <- IN[, 1]
IN <- as.matrix(IN[, -1])
sample <- read.table("ip_sample_ek12.list")
colnames(IN) <- sample$V1
ip_ek12_count <- IN       # normal11 with ip only

input_count_log2 <- log2(input_count + 1)
input_ek12_count_log2 <- log2(input_ek12_count + 1)
ip_count_log2 <- log2(ip_count + 1)
ip_ek12_count_log2 <- log2(ip_ek12_count + 1)
}

###################################
# ERCC log regression normalization
###################################
{
# for input
input_count_log2_ercc <- input_count_log2
N1 <- ncol(input_count_log2)

for(i in 2:N1)
{
	idx2 <- rowSums(input_ek12_count_log2[, c(1, i)] > 4) == 2	
	tmp2 <- lm(input_ek12_count_log2[idx2, i]~input_ek12_count_log2[idx2, 1])	
	idx_2 <- input_count_log2[, i] == 0
	input_count_log2_ercc[!idx_2, i] <- (input_count_log2[!idx_2, i] - tmp2$coefficients[1]) / tmp2$coefficients[2]	
}
input_count_log2_ercc[input_count_log2_ercc < 0] <- 0
input_count_ercc <- 2^input_count_log2_ercc -1

# for ip 
ip_count_log2_ercc <- ip_count_log2
N2 <- ncol(ip_count_log2)

for(i in 2:N2)
{
	idx1 <- rowSums(ip_ek12_count_log2[, c(1, i)] > 4) == 2
	tmp1 <- lm(ip_ek12_count_log2[idx1, i]~ip_ek12_count_log2[idx1, 1])
	idx_1 <- ip_count_log2[, i] == 0
	ip_count_log2_ercc[!idx_1, i] <- (ip_count_log2[!idx_1, i] - tmp1$coefficients[1]) / tmp1$coefficients[2]
}
ip_count_log2_ercc[ip_count_log2_ercc < 0] <- 0
ip_count_ercc <- 2^ip_count_log2_ercc -1
}

####################################
# DESeq  with ercc for normalization
####################################
{
# init set up    
library(DESeq2)
library("GenomicFeatures")

txdb <- makeTxDbFromGFF(file ="/Users/Yong/Yong/m6a_profiling/1_mapping_quant/5_RNA-seq/hg38_v25_ek12.gtf", format = "gtf")
#save(txdb, file = "hg38_v25_ek12_txdb.Rdata")
exonBygene <- exonsBy(txdb, by="gene")
rm(txdb)

all_gene <- read.table("/Users/Yong/Yong/m6a_profiling/1_mapping_quant/5_RNA-seq/all_gene_list.txt")
idx <- match(rownames(input_count_ercc), all_gene$V2)
tmp <- all_gene$V1[idx]
idx1 <- match(tmp, names(exonBygene))
exonBygene_f <- exonBygene[idx1]
names(exonBygene_f) <- all_gene$V2[idx]

#######################
## deseq2 normalization

# for input
expr <- round(input_count_ercc)
type <- factor(c(rep("Normal", 10), rep("Tumor", 56)))
dds <- DESeqDataSetFromMatrix(expr, DataFrame(type), ~type)
rowRanges(dds) <- exonBygene_f
dds <- DESeq(dds)
dds <- estimateSizeFactors(dds)
input_count_ercc_deseq <- counts(dds, normalized = T) 
input_count_ercc_deseq_rpkm <- fpkm(dds, robust = TRUE)

# for ip
expr <- round(ip_count_ercc)
type <- factor(c(rep("Normal", 11), rep("Tumor", 53)))
dds <- DESeqDataSetFromMatrix(expr, DataFrame(type), ~type)
rowRanges(dds) <- exonBygene_f
dds <- DESeq(dds)
dds <- estimateSizeFactors(dds)
ip_count_ercc_deseq <- counts(dds, normalized = T) 
ip_count_ercc_deseq_rpkm <- fpkm(dds, robust = TRUE)

save(input_count, ip_count, input_count_ercc, input_count_ercc_deseq, ip_count_ercc, ip_count_ercc_deseq, file = "all_count.Rdata")
save(input_count_ercc_deseq_rpkm, ip_count_ercc_deseq_rpkm, file = "all_rpkm.Rdata")
}

##########
# PCA plot
##########
{
    
library(devtools)
library(ggbiplot)
    
# extract 63 samples with both ip and input
idx_r1 <- c(30, 50, 51)  # remove tumor46 tumor29 rep
idx_r2 <- 11   # remove normal 11 ip

## raw read counts
input_count <- input_count[, -idx_r1]
ip_count <- ip_count[, -idx_r2]

## normalized counts
input_count_norm <- input_count_ercc_deseq[, -idx_r1]
ip_count_norm <- ip_count_ercc_deseq[, -idx_r2]

##########################################
## PCA plot before and after normalization
### with split screen
pca.plot <- function(count_l, name, file_n, shape, col_ip, col_input)
{
    pdf(file_n, height = 5.5 , width = 5.5)
    split.screen( figs = c(2, 2))
    for ( i in 1:length(count_l))
    {
        expr <- count_l[[i]]
        N = ncol(expr)
        idx <- rowSums(expr > 1) == ncol(expr)
        expr <- expr[idx, ]
        all_pca <- prcomp(t(expr), center = T, scale.=T)
        
        screen(i)
        if (i==2| i==4)
        {
            par(mar = c(3, 3, 2, 0.5), mgp  = c(2, 0.75, 0))
            plot(all_pca$x[, 1], all_pca$x[, 2], xlab = "PC1", ylab = "PC2", pch = shape, col = col_input, main = name[i], ylim = c(-200, 120), xlim = c(-180, 180))
            legend(x = 70, y = -130,c("Normal","Tumor"),pch=c(0,1), cex = 0.65)
            
        } else {
            par(mar = c(3, 3, 2, 0.5), mgp  = c(2, 0.75, 0))
            plot(all_pca$x[, 1], all_pca$x[, 2], xlab = "PC1", ylab = "PC2", pch = shape, col = col_ip, main = name[i], ylim = c(-100, 180), xlim = c(-100, 180))
            legend(x = 90, y = 180,c("Normal","Tumor"),pch=c(0,1), cex = 0.65)       
        }
        #text(x = all_pca$x[, 1], y = all_pca$x[, 2], labels = rownames(all_pca$x), pos = 3)	
        
    }
    close.screen(all = TRUE)
    dev.off()
}

## color and shape: 
name <- c( "IP_before_Normalization", "Input_before_Normalization",  "IP_after_Normalization", "Input_after_Normalization")
col_n <- c(rep("firebrick4", 5), rep("deepskyblue4", 5))       # batch1: normal01-05;   batch2: normal06-10
batch_c <- read.table("/Users/Yong/Yong/m6a_profiling/1_mapping_quant/5_RNA-seq/seq_batches_53.txt", header = T)
idx_s <- match(colnames(input_count)[11:63], rownames(batch_c))
col_tip <- as.character(batch_c$col_tip[idx_s])
col_tinput <- as.character(batch_c$col_tinput[idx_s])

col_ip <- c(col_n, col_tip)
col_input <- c(col_n, col_tinput)
shape <- c(rep(0, 10), rep(1, 53))

## tumor and normal together
count_all <- list( ip_count, input_count,  ip_count_norm, input_count_norm)
pca.plot(count_all, name, "normal_tumor_before_and_after_norm_PCA.pdf", shape, col_ip, col_input)
}

##########################
# MA plot : tumorvs normal
##########################
{
count_l <- list(ip_count, input_count,  ip_count_norm, input_count_norm)
name <- c( "ip_tumor_vs_normal", "input_tumor_vs_normal",  "ip_norm_tumor_vs_normal",  "input_norm_tumor_vs_normal")
L <- length(count_l)
N <- ncol(input_count)

pdf("normalized_before_after_MAplot.pdf", height = 14 , width = 14)
split.screen( figs = c(2, 2))
for ( i in 1:L)
{
    normal <- count_l[[i]][, 1:10] + 1
    tumor <- count_l[[i]][, 11:61] + 1
    M <- log2(rowMedians(tumor)/rowMedians(normal))
    A <- (1/2) * log2(rowMedians(tumor)*rowMedians(normal))
    
    screen(i)
    plot(A, M, xlab = "A", ylab = "M", main = name[i])
    
    idx1 <- is.na(A) | is.infinite(A)
    idx2 <-  is.na(M) | is.infinite(M)
    idx <- idx1&idx2
    A <- A[!idx]
    M <- M[!idx]
    
    abline(h = 0, lty = 2)
    reg <- lm(M~A)
    abline(a = reg$coefficients[1], b = reg$coefficients[2], col = "red")
}
close.screen(all = TRUE)
dev.off()
}

#############################
# log2(RPKM + 1) density plot
#############################
{
library("ggplot2")
input_rpkm <- input_count_ercc_deseq_rpkm[, -idx_r1]
ip_rpkm <- ip_count_ercc_deseq_rpkm[, -idx_r2]

rpkm_list <- list(ip_rpkm, input_rpkm)
L <- length(rpkm_list)
prefix <- c("ip", "input")

for(i in 1: L)
{
    expr_in <- rpkm_list[[i]]
    N <- ncol(expr_in)
    idx <- rowSums(expr_in > 1) == N
    expr_in <- log2(expr_in[idx, ])
    
    M <- nrow(expr_in)
    expr <- vector()
    group <- factor()
    for (i in 1: N )
    {
        expr <- c(expr, expr_in[, i])
        tag <- colnames(expr_in)[i]
        group  <- c(group, rep(tag, M))
    } 
    
    dat <- data.frame(expr, group)
    g <- ggplot(dat, aes(expr, col = group)) + geom_density(alpha = 0.2) +scale_x_continuous(name = "log2(RPKM + 1)") 
    g <- g + theme(legend.direction = 'horizontal', legend.position = 'top')+  theme_bw() 
    g <- g + theme(axis.line = element_line(colour = "black"),panel.grid.major = element_blank(),panel.grid.minor = element_blank(),panel.border =  element_blank(),   
                                                                                                                                                                                                             panel.background = element_blank())
    ggsave(filename=paste0(prefix[i], "_rpkm_density.pdf"), width = 6, heigh = 4.5)
}
}

#######################################
##  DE analysis for both IP and INPUT
#######################################
{
expr <- as.matrix(input_count)
type_cc <- factor(c(rep("Normal", 10), rep("Tumor", 53)))
dds_1 <- DESeqDataSetFromMatrix(expr, DataFrame(type_cc), ~type_cc)
dds_1 <- DESeq(dds_1)
res_1 <- results(dds_1)

# for ip_count
expr <- as.matrix(ip_count)
type_cc <- factor(c(rep("Normal", 10), rep("Tumor", 53)))
dds_2 <- DESeqDataSetFromMatrix(expr, DataFrame(type_cc), ~type_cc)
dds_2 <- DESeq(dds_2)
res_2 <- results(dds_2)

input_de <- res_1
ip_de <- res_2

save(input_de, ip_de, file = "input_ip_deseq_de.Rdata")
}


################
# GDC TCGA LUAD
# publisehd data
###############
{

##################################################   
## split published tcga luad into normal and tumor 
{
tcga <- read.table("/Users/Yong/Yong/m6a_profiling/1_mapping_quant/5_RNA-seq/TCGA/GDC_LUAD/TCGA-LUAD.htseq_counts.tsv", header = T)
#tcga <- read.table("/Users/Yong/Yong/m6a_profiling/8_RNA-seq/TCGA/GDC_LUAD/TCGA-LUAD.htseq_fpkm.tsv", header = T)
## have done log2( * + 1) transform

en_id <- as.character(tcga[, 1])
tcga <- tcga[, -1]
en_id_t <-strsplit(en_id, "\\.")
M <- length(en_id) 
en_id_o <- vector()
for (i in 1:M)
{
    en_id_o[i] <- en_id_t[[i]][1]
}

g_name <- read.table("/Users/Yong/Yong/m6a_profiling/1_mapping_quant/5_RNA-seq/TCGA/GDC_LUAD/GDC_gene_ID_with_name.txt")
idx_m <- match(g_name$V1, en_id_o)
tcga <- tcga[idx_m, ]

g_a <- levels(g_name$V2)
idx_u <- table(g_name$V2) == 1            # unique gene names
gene_u <- g_a[idx_u]

idx_k <- match(gene_u, g_name$V2)
tcga <- tcga[idx_k, ]
rownames(tcga) <- gene_u

type <- colnames(tcga)
type <- (strsplit(type, "\\."))
N <- length(type)
type_c <- vector()
for (i in 1:N)
{
    type_c[i] <- type[[i]][4]
}
idx_c <- type_c == "01A" | type_c == "02A" | type_c =="01B" | type_c == "01C"

tcga_type <- vector()
tcga_type[idx_c] <- "tumor"    # tumor
tcga_type[!idx_c] <- "normal"	# normal
col_c <- vector()
col_c[idx_c] <- "gray"    # tumor
col_c[!idx_c] <- "red"	# normal

save(tcga, tcga_type,  file="/Users/Yong/Yong/m6a_profiling/1_mapping_quant/5_RNA-seq/TCGA/GDC_LUAD/tcga_htseq_count_log2.Rdata")  
}

####################################  
# DEGseq2 normaliztion and DEGs
{
    library(DESeq2)
    tcga <- as.matrix(round(tcga))
    type_cc <- factor(col_c)
    dds_c <- DESeqDataSetFromMatrix(tcga, DataFrame(type_cc), ~type_cc)
    dds_c <- DESeq(dds_c)
    dds_c <- estimateSizeFactors(dds_c)
    tcga_deseq <- counts(dds_c, normalized = T)  
    res <- results(dds_c)
    res$log2FoldChange <- -res$log2FoldChange    ## switch to tumor vs normal
    write.csv(res, file = "tcga_sample_DESeq_DE.csv")
    
    ## PCA plot before and after normalization 
    {
    pdf("TCGA_normalized_before_after_pca.pdf", height = 14 , width = 7)
    split.screen( figs = c(2, 1))

    screen(1)
    expr <- tcga 
    N = ncol(expr)
    idx <- rowSums(expr > 1) == N
    expr <- expr[idx, ]
    all_pca <- prcomp(t(expr), center = T, scale.=T)
    plot(all_pca$x[, 1], all_pca$x[, 2], xlab = "PC1", ylab = "PC2", col = col_c, ylim = c(-150, 200), xlim = c(-150, 200))
    
    screen(2)
    expr <- tcga_deseq
    N = ncol(expr)
    idx <- rowSums(expr > 1) == N
    expr <- expr[idx, ]
    all_pca <- prcomp(t(expr), center = T, scale.=T)
    plot(all_pca$x[, 1], all_pca$x[, 2], xlab = "PC1", ylab = "PC2", col = col_c,  ylim = c(-150, 200), xlim = c(-150, 200))
    
    close.screen(all = TRUE)
    dev.off()
    }
    
}

}