########################################################################################
## This script implements:
## 1. immume infiltration estimation comparison
## 2. immume infiltration estimation comparison in TCGA LUAD
## 3. enriched GO terms andimmune related genes
## 4. immune related gene expression compariosn among sample methylated and unmethylated
########################################################################################

rm(list = ls())
setwd("/work/dir/")

library(heatmap.plus)
library(ggplot2)
library(gplots)

############################################
## immume infiltration estimation comparison
############################################
{
###################################################
## read-in immune celss estimated fraction by TIMER
{
IN <- read.csv(check.names = F, "/Users/Yong/Yong/m6a_profiling/2_peak_calling/1_peaks/merge/bed30_61/2_micro_envi/TIMER/6types_score_matrix.csv")
rownames(IN) <- IN$sampleID
frac<- as.matrix(IN[, 2: 7])

## heatmap  : need scale by cell type
hmcols<-colorRampPalette(c("blue","white","red"))(256)
png(file = "6_types_immune_cells_fraction_TIMER_heatmap.png", width = 8, height = 5, units = "in", res = 300)
heatmap.2(t(frac), dendrogram = "row", scale = "row", Rowv = T, Colv = F, col = hmcols, colsep = 10, lhei = c(1.5, 5), lwid = c(1, 4),
          key.title = NULL, key.xlab = "Scaled Cell fraction", density.info = "none", trace = "none", margins = c(5, 7), cexRow = 0.9)
dev.off()
}

###############################################
## test difference between the tumor and normal
{
p_val <- vector()
L <-  ncol(frac)
for (i in 1:L)
{
    tmp <- wilcox.test(frac[1:10, i], frac[11:61, i])
    p_val[i] <- tmp$p.value
}
idx_s <- p_val  < 1  #< 0.05              ### pvalue cutoff 0.01

## draw  boxplot
frac_v <- as.vector(frac)
g_1 <- rep(c(rep("normal",  10), rep("tumor", 51)), 6)
g_2 <- rep(colnames(frac), each = 61)

frac_go <- paste(g_2, g_1, sep = "_")  
frac_gou <- unique(frac_go)
# remian sort order 
frac_gf <- as.factor(frac_go)
frac_g <- factor(frac_gf, level = frac_gou)

## significant ones with color 
idx_ss = which(idx_s == TRUE)
col_s <- rep("gray", 12)
col_s[idx_ss*2] <- "coral2"        ## Tumor
col_s[idx_ss*2-1] <- "cadetblue3"  ## normal 

dat_g <- data.frame(frac_v, frac_g)

## boxplot
g <- ggplot(dat_g, aes(x = frac_g, y = frac_v)) + geom_boxplot(aes(col = frac_g), outlier.shape = NA) + theme(legend.position = "none") 
g <- g + geom_point(aes(col = frac_g), position = position_jitter(width = 0.3), size = 0.5)  + theme(legend.position = "none")  
g <- g + theme_bw() + theme(axis.line = element_line(colour = "black"),panel.grid.major = element_blank(),panel.grid.minor = element_blank(),panel.border =  element_blank(),  panel.background = element_blank()) 
g <- g + theme(axis.text.x = element_text(angle = 90)) + scale_color_manual(values= col_s) + theme(legend.position = "none") 
g <- g + xlab("Immune Cell types") + ylab("Franction")
ggsave("6_types_immune_cells_fraction_TIMER_boxplot.png",  width = 5, height = 5, units = "in", dpi = 600)

## for violin plots
g <- ggplot(dat_g, aes(x = frac_g, y = frac_v, col = frac_g, alpha = 1))  
g <- g + geom_violin(aes( fill=frac_g), trim=FALSE, scale = "width", lwd = 0.2) + geom_boxplot(width=0.5, lwd = 0.3) + theme_classic()
g <- g + scale_color_manual( values = col_s)  + scale_fill_manual( values = col_s) + theme(legend.position = "none")
ggsave("6_types_immune_cells_fraction_TIMER_violin_plot.pdf",  width = 4, height = 3)
}

}


#########################################################
## immume infiltration estimation comparison in TCGA LUAD
#########################################################
{
###################################################
## read-in immune celss estimated fraction by TIMER
{
luad <- read.table("/Users/Yong/Yong/m6a_profiling/2_peak_calling/1_peaks/merge/bed30_61/2_micro_envi/TIMER/immuneEstimation_TCGA.txt", header = T, row.names = 1)   ## all TCGA samples

## load LUAD normal  and tumor smaple list
load("/Users/Yong/Yong/m6a_profiling/1_mapping_quant/5_RNA-seq/TCGA/GDC_LUAD/TCGA_LUAD_normal_tumor_sample_ID.Rdata")
idx_n <- match(tcga_normal_sample, rownames(luad), nomatch = 0)
sum(idx_n != 0)
frac_n <- luad[idx_n, ]

idx_t <- match(tcga_tumor_sample, rownames(luad), nomatch = 0)
sum(idx_t != 0)
frac_t <- luad[idx_t, ]
frac <- as.matrix(rbind(frac_n, frac_t))

## heatmap  : need scale by cell type
hmcols<-colorRampPalette(c("blue","white","red"))(256)
png(file = "6_types_immune_cells_fraction_TIMER_heatmap_TCGA_LUAD.png", width = 8, height = 5, units = "in", res = 300)
heatmap.2(t(frac), labRow = NULL, dendrogram = "row", scale = "row", Rowv = T, Colv = F, col = hmcols, colsep = nrow(frac_n), lhei = c(1.5, 5), lwid = c(1, 4),
          key.title = NULL, key.xlab = "Scaled Cell fraction", density.info = "none", trace = "none", margins = c(5, 7), cexRow = 0.9)
dev.off()
}

#################################
## draw  boxplot and violin plots
{
p_val <- vector()
L <-  ncol(frac)
for (i in 1:L)
{
        tmp <- wilcox.test(frac[1:nrow(frac_n), i], frac[59:nrow(frac), i])
        p_val[i] <- tmp$p.value
    }
idx_s <- p_val < 1             ### pvalue cutoff 0.01    
    
frac_v <- as.vector(frac)
g_1 <- rep(c(rep("normal",  nrow(frac_n)), rep("tumor", nrow(frac_t))), 6)
g_2 <- rep(colnames(frac), each = nrow(frac))

frac_go <- paste(g_2, g_1, sep = "_")  
frac_gou <- unique(frac_go)
# remian sort order 
frac_gf <- as.factor(frac_go)
frac_g <- factor(frac_gf, level = frac_gou)

## significant ones with color 
idx_ss = which(idx_s == TRUE)
col_s <- rep("gray", 12)
col_s[idx_ss*2] <- "coral2"        ## Tumor
col_s[idx_ss*2-1] <- "cadetblue3"  ## normal 

dat_g <- data.frame(frac_v, frac_g)
g <- ggplot(dat_g, aes(x = frac_g, y = frac_v)) + geom_boxplot(aes(col = frac_g), outlier.shape = NA) + theme(legend.position = "none") 
g <- g + geom_point(aes(col = frac_g), position = position_jitter(width = 0.3), size = 0.5)  + theme(legend.position = "none")  
g <- g + theme_bw() + theme(axis.line = element_line(colour = "black"),panel.grid.major = element_blank(),panel.grid.minor = element_blank(),panel.border =  element_blank(),  panel.background = element_blank()) 
g <- g + theme(axis.text.x = element_text(angle = 90)) + scale_color_manual(values= col_s) + theme(legend.position = "none") 
g <- g + xlab("Immune Cell types") + ylab("Franction")

ggsave("6_types_immune_cells_fraction_TIMER_boxplot_TCGA_LUAD.png", width = 5, height = 5, units = "in", dpi = 600)

## for violin plots
g <- ggplot(dat_g, aes(x = frac_g, y = frac_v, col = frac_g, alpha = 1))  
g <- g + geom_violin(aes( fill=frac_g), trim=FALSE, scale = "width", lwd = 0.2) + geom_boxplot(width=0.5, lwd = 0.3) + theme_classic()
g <- g + scale_color_manual( values = col_s)  + scale_fill_manual( values = col_s) + theme(legend.position = "none")
ggsave("6_types_immune_cells_fraction_TIMER_boxplot_TCGA_LUAD.pdf",  width = 4, height = 3)
}

}


############################################
## enriched GO terms andimmune related genes
############################################
{
#################
## TOP2 MF genes
{
go_g <- read.table("/Users/Yong/Yong/m6a_profiling/2_peak_calling/1_peaks/merge/bed30_61/2_micro_envi/Top2_MF_terms_unique_genes.txt")

## read DE analysis resutls
de_g <- read.csv("/Users/Yong/Yong/m6a_profiling/1_mapping_quant/5_RNA-seq/htseq/count/61_sample_Input_DESeq_DE.csv")

idx <- match(go_g$V1, de_g$X)
go_gr <- de_g[idx, ] 

idx_de <-  abs(go_gr$log2FoldChange) > 1 & go_gr$padj < 0.05
table(idx_de)
col_de <- rep("gray", nrow(go_gr))
col_de[idx_de] <- "red"
fc <- go_gr$log2FoldChange
fdr_t <- -log10(go_gr$padj)
idx_c <- fdr_t > 10
fdr_t[idx_c] <- 10 

pdf("Top2_MF_tumor_only_genes_DE_plot.pdf", width = 3, height = 3)
par(mar = c(3, 3, 0.5, 0.5), mgp = c(2, 0.75, 0))
plot(fc, fdr_t, col = col_de, xlab = "log2FC", ylab = "-log10FDR")
abline(h = -log10(0.05), v = 0, lty = 2)
dev.off()
}


#############################################
## Tumor only genes related to the immune GO
{
go_g <- read.table("/Users/Yong/Yong/m6a_profiling/2_peak_calling/1_peaks/merge/bed30_61/2_micro_envi/tumor_only_methylated_gene_GO_immune_related_terms_unique_genes.txt")

idx <- match(go_g$V1, de_g$X)
go_gr <- de_g[idx, ] 

idx_de <-  abs(go_gr$log2FoldChange) > 1 & go_gr$padj < 0.05
table(idx_de)
col_de <- rep("gray", nrow(go_gr))
col_de[idx_de] <- "red"
fc <- go_gr$log2FoldChange
fdr_t <- -log10(go_gr$padj)
idx_c <- fdr_t > 10
fdr_t[idx_c] <- 10 

pdf("Tumor_only_immune_related_genes_DE_plot.pdf", width = 3, height = 3)
par(mar = c(3, 3, 0.5, 0.5), mgp = c(2, 0.75, 0))
plot(fc, fdr_t, col = col_de, xlab = "log2FC", ylab = "-log10FDR")
abline(h = -log10(0.05), v = 0,  lty = 2)
dev.off()
}
}


#####################################################################################
## immune related gene expression compariosn among sample methylated and unmethylated
#####################################################################################
{ 
###########################################################################################
## for tumor only : split tumor sample to groups with called peaks and without called peaks
{
    load("/Users/Yong/Yong/m6a_profiling/1_mapping_quant/5_RNA-seq/htseq/count/paired_61_rpkm.Rdata")
    gene_to <- read.table("/Users/Yong/Yong/m6a_profiling/2_peak_calling/1_peaks/merge/bed30_61/0_merged_peaks/tumor_only.bed", as.is = T)

    eg_g <- go_g$V1
    
    idx_g <- match(eg_g, gene_to$V4)
    idx_gs <- idx_g[!is.na(idx_g)]
    eg_gs <- eg_g[!is.na(idx_g)]
    gene_tos <- gene_to[idx_gs, ]
    gene_tos_sl <- strsplit(gene_tos$V5, ",")       ## samples burdern tumor only called peaks 
    
    ano_diff_p <- vector()              ##  ANOVA test
    kw_diff_p <- vector()               ##  Kruskal–Wallis test : nonparametic 
    
    L <- length(eg_gs)
    expr_to_mean <- matrix(0, L, 3)

    for(i in 1:L)
    {
        idx_gg <- match(eg_gs[i], rownames(input_rpkm))
        
        expr_n <- input_rpkm[idx_gg, 1:10]
        
        idx_tp <- match(gene_tos_sl[[i]], colnames(input_rpkm))          ## tumor samples with tumor only methylated genes
        expr_tp <- input_rpkm[idx_gg, idx_tp]            
        
        expr_tnp <- input_rpkm[idx_gg, -c(1:10, idx_tp)]
        
        expr_to_mean[i, ] <- c(log2(mean(expr_n + 1)), log2(mean(expr_tnp + 1)), log2(mean(expr_tp + 1)))
        
        expr <- c(expr_n, expr_tnp, expr_tp)
        group <- c(rep(1, length(expr_n)), rep(2, length(expr_tnp)), rep(3, length(expr_tp)))
        
        dat <- data.frame(expr, group)
        
        ## annova test 
        #res.aov <- aov(expr ~ group, data = dat)
        #ano_diff_p[i] <-  summary(res.aov)[[1]]["Pr(>F)"][[1]][1]
        
        ## Kruskal-Wallis test
        res.kw <- kruskal.test(expr ~group, data = dat)
        
        kw_diff_p[i] <-  res.kw$p.value
        
    }
    
    ## for annova
    # ano_diff_padj <- p.adjust(ano_diff_p, method = "BH")
    # sum(ano_diff_padj < 0.05, na.rm = T)
    # sum(ano_diff_padj >= 0.05, na.rm = T)
    
    ## for Kruskal-wallis test
    kw_diff_padj <- p.adjust(kw_diff_p, method = "BH")
    sum(kw_diff_padj < 0.05, na.rm = T)
    sum(kw_diff_padj >= 0.05, na.rm = T)
    
    ### based on Kruskal-wallis test
    if(TRUE){
        names(kw_diff_padj) <- eg_gs                              
        kw_diff_padj_s <- sort(kw_diff_padj, decreasing = F)      
        L1 <- length(kw_diff_padj_s)
        col_d <- rep("gray", L1)
        idx_r <- kw_diff_padj_s < 0.05
        col_d[idx_r] <- "red"
        
        ## B cell specifically expressed genes:  CD79B, BLK, PAX5;  CCR6 no expression   ##  T cell specific : CD3G
        # g_imm <- c("BLK", "CD79B", "PAX5", "CD3G")
        # g_imm <- c("B4GALT1", "CD79B", "PAX5", "MSN")        ## slected 2 non-diff
        # idx_imm <- c(1, match(g_imm, names(kw_diff_padj_s)), L1)      ## add the most sig and non-sig DE gene
        
        g_imm <- c("CD79B", "PAX5", "TNFSF9", "ISG15", "IGKV4-1", "HLA-DRA")     ## backup "HLA-DPB1"
        idx_imm <- match(g_imm, names(kw_diff_padj_s))
        pch_t <- rep(2, L1)
        pch_t[idx_imm] <- 17
        
        pdf(file = "Tumor_only_immune_related_genes_DE_plot_KW_test_new.pdf", width = 3, height = 3)
        par(mar = c(3, 3, 1, 1), mgp =c(2, 0.75, 0))
        plot(kw_diff_padj_s, col = col_d, cex = 0.5, lwd = 0.5, pch = pch_t, ylab = "FDR", xlab = "Ranked genes")
        abline(h = 0.05, lty = 2)
        text(x = idx_imm, y = kw_diff_padj_s[idx_imm], labels = names(kw_diff_padj_s)[idx_imm] , pos = 2, cex = 0.5)
        dev.off()
    }
    
    sort(kw_diff_padj_s[idx_imm])
}
 
#####################
## example gene plots
{
gene_eg <- c(names(kw_diff_padj_s)[sort(idx_imm)])   

{
    eg_g <- gene_eg
    idx_g <- match(eg_g, gene_to$V4)
    idx_gs <- idx_g[!is.na(idx_g)]
    eg_gs <- eg_g[!is.na(idx_g)]
    gene_tos <- gene_to[idx_gs, ]
    gene_tos_sl <- strsplit(gene_tos$V5, ",")       ## samples burdern tumor only called peaks 
    
    ano_diff_p <- vector()
    
    L <- length(eg_gs)
    
    expr_to_mean <- matrix(0, L, 3)
    expr_sg <- data.frame()
    
    for(i in 1:L)
    {
        idx_gg <- match(eg_gs[i], rownames(input_rpkm))
        expr_n <- input_rpkm[idx_gg, 1:10]
        idx_tp <- match(gene_tos_sl[[i]], colnames(input_rpkm))          ## tumor samples with tumor only methylated genes
        expr_tp <- input_rpkm[idx_gg, idx_tp]            
        
        expr_tnp <- input_rpkm[idx_gg, -c(1:10, idx_tp)]
        expr_to_mean[i, ] <- c(log2(mean(expr_n + 1)), log2(mean(expr_tnp + 1)), log2(mean(expr_tp + 1)))
        
        expr <- log2(c(expr_n, expr_tnp, expr_tp) + 1)
        group <- as.character(c(rep(1, length(expr_n)), rep(2, length(expr_tnp)), rep(3, length(expr_tp))))
        gene <- rep(eg_gs[i], length(expr))
        
        dat <- data.frame(expr, group, gene)   ## plot_boxplot with
        g <- ggplot(dat, aes(x = group, y = expr, col = group)) + geom_boxplot(outlier.shape = NA) 
        g <- g + geom_point(aes(col = group), position = position_jitter(width = 0.3)) 
        g <- g + theme_bw() + theme(axis.line = element_line(colour = "black"),panel.grid.major = element_blank(),panel.grid.minor = element_blank(),panel.border =  element_blank(),  panel.background = element_blank()) 
        g <- g + scale_color_manual( values = c("cadetblue3", "coral2", "coral3")) + theme(legend.position = "none")  + ggtitle(eg_gs[i])
        name <- paste("Tumor_only_immune_related_genes_DE_eg_plot_anova_",eg_gs[i],".pdf", sep = "")
        ggsave(name, width = 2, height = 3, units = "in")
        
        ## merged slected gene for facet plot
        expr_sg <- rbind(expr_sg, dat)
        
    }

}

## facet plot for genes
g <- ggplot(expr_sg, aes(x = group, y = expr, col = group)) + geom_boxplot(outlier.shape = NA) 
g <- g + geom_point(aes(col = group), position = position_jitter(width = 0.3), cex = 0.75) 
g <- g + scale_color_manual( values = c("cadetblue3", "coral2", "coral3")) + theme_bw()  + theme(legend.position = "top") 
g <- g + facet_grid(cols = vars(gene))
ggsave("Tumor_only_immune_related_genes_DE_6eg_facet_plot.pdf", width = 9, height = 3, units = "in")
}
}