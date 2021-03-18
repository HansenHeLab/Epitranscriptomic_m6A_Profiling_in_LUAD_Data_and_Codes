################################################################
## This script implements:
## 1. gene_m6A_level_and_filtering
## 2. gene m6A level vs H1ESC percentage m6A level
## 3. gene m6A level distribution and global diff 
## 4. association between tumor gene m6A level and clinical info
## 5. differentially m6A mehtylated genes
## 6. gene m6A level diff vs RNA level diff 
#################################################################

rm(list = ls())
setwd("/work/dir/")
               
library("matrixStats")

###############################
## gene_m6A_level_and_filtering
###############################
{
#### load data ######################################
# normalized by ERCC and DESeq with all duplications
{
    load("/Users/Yong/Yong/m6a_profiling/1_mapping_quant/5_RNA-seq/htseq/count/paired_61_count_ercc_deseq.Rdata") 
    load("/Users/Yong/Yong/m6a_profiling/1_mapping_quant/5_RNA-seq/htseq/count/paired_61_rpkm.Rdata") 
    
    ip_rpkm_norm <- ip_rpkm
    input_rpkm_norm <- input_rpkm
    gene_a <- read.table("~/Yong/m6a_profiling/2_peak_calling/1_peaks/merge/bed30_61/gene_with_m6a_peak_in_tumor_or_normal.txt")    # all gene with m6a Peak  
    pc_gene <- read.table("/Users/Yong/Yong/m6a_profiling/1_mapping_quant/5_RNA-seq/pc_gene_list.txt")
    lincRNA <- read.table("/Users/Yong/Yong/m6a_profiling/1_mapping_quant/5_RNA-seq/lincRNA_gene_list.txt")
    
    s_info <- read.table("/Users/Yong/Yong/m6a_profiling/2_peak_calling/1_peaks/merge/sample_info.txt", header = T)       # sample information 
}

##### RPKM based filtering ###############################################################
# gene with >= 1 RPKM in both ip or input for at least half of sample  && called m6a peaks
{
    #### for different RPKM threshold
    N <- ncol(ip_rpkm_norm)
    N_n <- 10                        #  normal sample size
    N_t <- 51                        #  tumor sample size
    gene_ll_1 <- vector()            #  length for genes pass the cutoff, and with m6A peaks as well  
    gene_ll_2 <- vector()            #  length for genes pass the cutoff, without detect m6a peaks
    
    for(cutoff in seq(0, 10, 1))
    {
        idx_ip <- ip_rpkm_norm > cutoff
        idx_input <- input_rpkm_norm > cutoff
        idx_m <- idx_ip & idx_input
        idx1 <- rowSums(idx_m[, 1:10]) >= round(N_n/2)
        idx2 <- rowSums(idx_m[, 11:61]) >= round(N_t/2)
        idx_g <- idx1 | idx2                                    # more than half of sample with ip and input > cutoff in tumor or normal separatly 
        gene_tmp <- rownames(ip_count_norm)[idx_g]
        
        idx_1 <- !is.na(match(gene_tmp, gene_a$V1))             # gene pass cut off and with m6A peaks
        gene_1 <- gene_tmp[idx_1]
        gene_ll_1 <- c(gene_ll_1, length(gene_1))
        
        idx_2 <- is.na(match(gene_tmp, gene_a$V1))             # gene pass cut off and without m6A peaks
        gene_2 <- gene_tmp[idx_2]
        gene_ll_2 <- c(gene_ll_2, length(gene_2))
        
    }
    
    #png(file = "genes_filtering_based_RPKM_and_peaks.png", width = 4, height = 4, units = "in", res = 300)
    pdf(file = "genes_filtering_based_RPKM_and_peaks.pdf", width = 3.5, height = 3)
    par(mar = c(4, 4, 0.5, 0.5), mgp =c(3, 0.7, 0))
    plot(seq(0, 10, 1), gene_ll_1, type = "b", lty = 1, col = "red", xaxt = "n", yaxt = "n", ylim = c(0, 16000), xlab = "RPKM threshold for both IP and Input", ylab = "Numbers of gene")
    axis(1, at = seq(0, 10, 1))
    axis(2, at = seq(0, 16000, 4000), las = 1, labels = c("0", "4,000", "8,000", "12,000", "16,000"))
    lines(seq(0, 10, 1), gene_ll_2, type = "b", lty = 1)
    lines(x = c(1, 1), y = c(1, gene_ll_1[2]), lty = 2)
    legend(6, 16000, legend = c("With m6A", "Without m6A"), lty = c(1, 2), col = c("red", "black"),  cex = 0.7)
    dev.off()
    
    #### selected RPKM threshold: RPKM > 1 in both IP and Input
    cutoff <- 1 
    idx_ip <- ip_rpkm_norm > cutoff
    idx_input <- input_rpkm_norm > cutoff
    idx_m <- idx_ip & idx_input
    idx1 <- rowSums(idx_m[, 1:10]) >= round(N_n/2)
    idx2 <- rowSums(idx_m[, 11:61]) >= round(N_t/2)
    idx_g <- idx1 | idx2                                    # more than half of sample with ip and input > cutoff in tumor or normal separatly 
    gene_tmp <- rownames(ip_count_norm)[idx_g]
    
    idx_1 <- !is.na(match(gene_tmp, gene_a$V1))             # gene pass cut off and with m6A peaks
    gene_m6a <- gene_tmp[idx_1]
    
    idx_2 <- is.na(match(gene_tmp, gene_a$V1))             # gene pass cut off and without m6A peaks
    gene_nm6a <- gene_tmp[idx_2]
    
}

#### genes types #####################################
# gene categories: ptrotein coding, lincRNA and other 
{
    gene.ty <- function(gene_se, pc_gene, lincRNA, o_name)
    {
        pc_cnt <- sum(!is.na(match(gene_se,pc_gene$V2)))
        li_cnt <- sum(!is.na(match(gene_se,lincRNA$V2)))
        other <-length(gene_se) - pc_cnt - li_cnt
        cnt <- c(pc_cnt, li_cnt, other)
        print(cnt) 
        pnt <- round((cnt/sum(cnt))* 100, 2)
        lab <- c("protein coding", "lincRNA", "Others")
        lab1 <- paste(lab, pnt, sep = " ")
        lab2 <- paste(lab1, "%", sep = "")
        
        file_name <- paste(o_name, ".pdf", sep = "")
        pdf(file = file_name, width = 3, height = 3)
        par(mar = c(0.5, 6, 0.5, 4), mgp =c(2, 0.7, 0))
        pie(cnt, labels = lab2,  col = c("orange", "violet", "spring green"))
        legend("top", legend = lab, fill = c("orange", "violet", "spring green"), bty ="n", cex = 0.7)
        dev.off()
    }
    gene.ty(gene_m6a, pc_gene, lincRNA, "Gene_with_m6A_RPKM>1")
    gene.ty(gene_nm6a, pc_gene, lincRNA, "Gene_without_m6A_RPKM>1")
}

####  m6a level ##########
## m6a level : IP / input 
{
idx_1 <- match(gene_m6a, rownames(ip_count_norm))
idx_2 <- match(gene_nm6a, rownames(ip_count_norm))

# based on normalized count 
m6a_r <- ip_count_norm / input_count_norm              ## based on counts 
# idx_i <- is.infinite(m6a_r) | (m6a_r == 0)           ## ratio to be zero or infinite
idx_i <- is.infinite(m6a_r)                            ## updated 20190329
m6a_r[idx_i] <- NA

###########################################
##  m6a level for all gene withoutfiltering 
m6a_rn_all <- m6a_r[, 1:10]
m6a_rt_all <- m6a_r[, 11:61]

## normal tumor separate for all genes
input_count_n_all <- input_count_norm[, 1:10]
input_count_t_all <- input_count_norm[, 11:61]
ip_count_n_all <- ip_count_norm[, 1:10]
ip_count_t_all <- ip_count_norm[, 11:61]
input_rpkm_n_all <- input_rpkm_norm[, 1:10]
input_rpkm_t_all <- input_rpkm_norm[, 11:61]
ip_rpkm_n_all <- ip_rpkm_norm[, 1:10]
ip_rpkm_t_all <- ip_rpkm_norm[, 11:61]  

save(m6a_rn_all, m6a_rt_all, input_count_n_all, input_count_t_all, ip_count_n_all, ip_count_t_all,
     input_rpkm_n_all, input_rpkm_t_all, ip_rpkm_n_all, ip_rpkm_t_all, file = "m6A_level_4gene_all.Rdata")

#######################################
##  m6a level for gene passed filtering 
m6a_rs <- m6a_r[idx_1, ]  
input_count_norm_s <- input_count_norm[idx_1, ]
ip_count_norm_s <- ip_count_norm[idx_1, ]
input_rpkm_norm_s <- input_rpkm_norm[idx_1, ]
ip_rpkm_norm_s <- ip_rpkm_norm[idx_1, ]

## normal tumor separate 
m6a_rn <- m6a_rs[, 1:10]
m6a_rt <- m6a_rs[, 11:61]
input_count_n <- input_count_norm_s[, 1:10]
input_count_t <- input_count_norm_s[, 11:61]
ip_count_n <- ip_count_norm_s[, 1:10]
ip_count_t <- ip_count_norm_s[, 11:61]
input_rpkm_n <- input_rpkm_norm_s[, 1:10]
input_rpkm_t <- input_rpkm_norm_s[, 11:61]
ip_rpkm_n <- ip_rpkm_norm_s[, 1:10]
ip_rpkm_t <- ip_rpkm_norm_s[, 11:61]   

save(m6a_rn, m6a_rt, input_count_n, input_count_t, ip_count_n, ip_count_t,
     input_rpkm_n, input_rpkm_t, ip_rpkm_n, ip_rpkm_t, file = "m6A_level_4gene_passed_filtering.Rdata")
}

############################################
## comparsion with the RPKM based IP / Input
{
m6a_rr <- ip_rpkm_norm / input_rpkm_norm
#idx_i <- is.infinite(m6a_rr) | (m6a_rr == 0)           ## ratio to be zero or infinite
idx_i <- is.infinite(m6a_rr) 
m6a_rr[idx_i] <- NA
m6a_rrs <- m6a_rr[idx_1, ]   ##  m6a level for gene passed filtering 

L <- nrow(m6a_rrs)
cor_cmp <- vector()
for(i in 1:L)
{
    cor_cmp[i] <- cor(m6a_rs[i, ], m6a_rrs[i, ])
}

# Conlusion: almost the same with the count based the m6a level
}

##################################################
## IP/Input for gene passed the filtering and not 
{
normal_r_m6a <- log2(rowMeans(m6a_r[idx_1, 1:10], na.rm = T))
tumor_r_m6a<- log2(rowMeans(m6a_r[idx_1, 11:61], na.rm = T))

normal_r_nm6a <- log2(rowMeans(m6a_r[idx_2, 1:10], na.rm = T))
tumor_r_nm6a<- log2(rowMeans(m6a_r[idx_2, 11:61], na.rm = T))
idx_nr <- !(is.infinite(normal_r_nm6a) | is.na(normal_r_nm6a))
idx_tr <- !(is.infinite(tumor_r_nm6a) | is.na(tumor_r_nm6a))




#png(file = "mean_m6a_level_log2_boxplot.png", width = 5, height = 4, units = "in", res = 300)

pdf(file = "mean_m6a_level_log2_boxplot.pdf", width = 4, height = 3)
par(mar = c(3, 3, 0.5, 0.5), mgp =c(2, 0.75, 0))
boxplot(list(normal_r_nm6a, normal_r_m6a, tumor_r_nm6a, tumor_r_m6a),  col = c("gray80", "cadetblue3", "gray80", "coral2"), outcol = c("gray80", "cadetblue3", "gray80", "coral2"), 
        axes = F, ylim = c(-6, 6), ylab = "log2(mean(IP/Input))", xlab = "Normal                          Tumor")
xl <- parse(text = paste(rep(c("Without_m", "With_m"), 2), "^6", "*A", sep=""))      # superscript    *A   
axis(1, at = seq(1, 4, 1), labels = xl )     # for x axis
axis(2, at = seq(-6, 6, 2))
dev.off()

library(BSDA)  ## for Z-test
z.test(normal_r_nm6a[idx_nr], normal_r_m6a, sigma.x = sd(normal_r_nm6a[idx_nr]), sigma.y = sd(normal_r_m6a), alternative = "less")
z.test(tumor_r_nm6a[idx_tr], tumor_r_m6a, sigma.x = sd(tumor_r_nm6a[idx_tr]), sigma.y = sd(tumor_r_m6a), alternative = "less")

t.test(normal_r_nm6a[idx_nr], normal_r_m6a, alternative = "less")
t.test(tumor_r_nm6a[idx_tr], tumor_r_m6a, alternative = "less")


}

}


################################################
##  gene m6A level vs H1ESC percentage m6A level 
################################################
{
load("/Users/Yong/Yong/m6a_profiling/3_m6a_level/m6a_laic/m6a-seq/m6a_level_cmp_all_log2.RData")

idx_s <- match(names(m6a_rff), names(tumor_r_m6a))
sum(!is.na(idx_s))
m6a_rffs <- m6a_rff[!is.na(idx_s)]
tumor_r_m6a_s <- tumor_r_m6a[idx_s[!is.na(idx_s)]]
normal_r_m6a_s <- normal_r_m6a[idx_s[!is.na(idx_s)]]

## scatter plot
cor.plot(normal_r_m6a_s, m6a_rffs, "normal_m6a_level_vs_H1sec_m6a_level", "Normal_log2(mean(IP/Input))",  "H1ESC_log2(IP/Input)")
cor.plot(tumor_r_m6a_s, m6a_rffs, "tumor_m6a_level_vs_H1sec_m6a_level", "Tumor_log2(mean(IP/Input))",  "H1ESC_log2(IP/Input)")

#########################
## ggplot hexagnoal heatmap
## for normal
dat <- data.frame(normal_r_m6a_s, m6a_rffs)
reg <- lm(m6a_rffs ~ normal_r_m6a_s)
g <- ggplot(dat, aes(normal_r_m6a_s, m6a_rffs))
g <- g + geom_hex(bins = 50) +  geom_abline(intercept = reg$coefficients[1], slope = reg$coefficients[2], color="red", linetype="dashed")
g <- g + theme_classic()
ggsave("normal_m6a_level_vs_H1sec_m6a_level_Hexagonal.pdf", width = 4, height = 3)

## for tumor
dat <- data.frame(tumor_r_m6a_s, m6a_rffs)
reg <- lm(m6a_rffs ~ tumor_r_m6a_s)
g <- ggplot(dat, aes(tumor_r_m6a_s, m6a_rffs))
g <- g + geom_hex(bins = 50) +  geom_abline(intercept = reg$coefficients[1], slope = reg$coefficients[2], color="red", linetype="dashed")
g <- g + theme_classic()
ggsave("tumor_m6a_level_vs_H1sec_m6a_level_Hexagonal.pdf", width = 4, height = 3)

##################
## density plot  
pdf(file = "mean_m6a_level_log2_density_with_h1esc_median.pdf", width = 4, height = 4)
par(mar = c(3, 3, 1, 1), mgp =c(2, 0.7, 0))
plot(density(tumor_r_m6a_s, na.rm = T), col = "coral2", xlab = "log2(mean(IP/Input))", main = "", xlim = c(-6, 6))
lines(density(normal_r_m6a_s, na.rm = T), col = "cadetblue3")
lines(density(m6a_rffs, na.rm = T), col = "gray", lty = 1)
legend(2, 0.32, legend = c("Tumor", "Normal", "H1ESC"), col = c("coral2", "cadetblue3", "gray"), lty = c(1, 1, 1), cex = 0.8)
abline(v = median(tumor_r_m6a_s, na.rm = T), col = "coral2", lty = 2)
abline(v = median(normal_r_m6a_s, na.rm = T), col = "cadetblue3", lty = 2)
abline(v = median(m6a_rffs, na.rm = T), col = "gray", lty = 2)
dev.off()
}

################################################
##  gene m6A level distribution and global diff 
################################################
{
    library("Hmisc")                 # for correlaiton with NA
    library("ggplot2")
    library("matrixStats")
    library("heatmap.plus")
    
    ## load m6a level, mrna level for selected genes 
    load("/Users/Yong/Yong/m6a_profiling/3_m6a_level/m6a_level_keep_dup/gene_with_peak/m6A_level_4gene_passed_filtering.Rdata")
    pc_gene <- read.table("/Users/Yong/Yong/m6a_profiling/1_mapping_quant/5_RNA-seq/pc_gene_list.txt") 
    s_info <- read.table("/Users/Yong/Yong/m6a_profiling/2_peak_calling/1_peaks/merge/sample_info.txt", header = T) 
    source("/Users/Yong/Yong/R/functions/cor.plot.R")
    
    ###############################
    ## mean and cv of m6A per gene 
    {
        ## box plot of m6A level per gene for each samples
        {
            m6a_rm  <- cbind(m6a_rn, m6a_rt)
            png(file = "m6a_level_log2_boxplot_per_sample.png", width = 10, height = 5, units = "in", res = 300)
            par(mar = c(5, 3, 2, 1), mgp =c(2, 0.7, 0))
            boxplot(log2(m6a_rm), outline = F, xaxt = "n", ylab = "log2(IP/Input)")
            axis(1, at = seq(1, 61, 1), labels = colnames(m6a_rm), las = 2, cex = 0.4)
            abline(v = 10.5, lty = 2 , col = "red")
            dev.off()
        }
        
        ## mean m6a level per slected genes 
        {
            normal_r_m6a <- log2(rowMeans(m6a_rn, na.rm = T))
            tumor_r_m6a<- log2(rowMeans(m6a_rt, na.rm = T))
            
            ## boxplot
            pdf(file = "mean_m6a_level_log2_boxplot_with_m6A.pdf", width = 3, height = 3)
            par(mar = c(3, 3, 0.5, 0.5), mgp =c(2, 0.75, 0))
            boxplot(list(normal_r_m6a,  tumor_r_m6a), col = c("cadetblue3", "coral2"), axes = F, ylim = c(-6, 6), ylab = "log2(mean(IP/Input))")
            xl <- c("Normal", "Tumor")       
            axis(1, at = seq(1, 2, 1), labels = xl )     # for x axis
            axis(2, at = seq(-6, 6, 3))
            dev.off()
            
            ## densidy plot for gene with m6A only
            png(file = "mean_m6a_level_log2_density_mean.png", width = 4, height = 4, units = "in", res = 300)
            par(mar = c(3, 3, 2, 1), mgp =c(2, 0.7, 0))
            plot(density(tumor_r_m6a, na.rm = T), col = "coral2", xlab = "log2(mean(IP/Input))", main = "Gene with m6A")
            lines(density(normal_r_m6a, na.rm = T), col = "cadetblue3")
            legend(2.5, 0.32, legend = c("Tumor", "Normal"), col = c("red", "blue"), lty = c(1, 1), cex = 0.8)
            abline(v = mean(tumor_r_m6a, na.rm = T), col = "coral2", lty = 2)
            abline(v = mean(normal_r_m6a, na.rm = T), col = "cadetblue", lty = 2)
            dev.off()
            
            z.test(normal_r_m6a, tumor_r_m6a, sigma.x = sd(normal_r_m6a), sigma.y = sd(tumor_r_m6a), alternative = "greater")
            #t.test(normal_r_m6a, tumor_r_m6a, alternative = "greater")
            #wilcox.test(normal_r_m6a, tumor_r_m6a, alternative = "greater")
            
            
            ####################################################################
            ##    interprate m6A level by comparing to H1ESC data (m6A level %)
            load("/Users/Yong/Yong/m6a_profiling/3_m6a_level/m6a_laic/m6a-seq/m6a_level_cmp_all_log2.RData")
            
            idx_s <- match(names(m6a_rff), names(tumor_r_m6a))
            sum(!is.na(idx_s))
            m6a_rffs <- m6a_rff[!is.na(idx_s)]
            tumor_r_m6a_s <- tumor_r_m6a[idx_s[!is.na(idx_s)]]
            normal_r_m6a_s <- normal_r_m6a[idx_s[!is.na(idx_s)]]
            
            ## scatter plot
            cor.plot(normal_r_m6a_s, m6a_rffs, "normal_m6a_level_vs_H1sec_m6a_level", "Normal_log2(mean(IP/Input))",  "H1ESC_log2(IP/Input)")
            cor.plot(tumor_r_m6a_s, m6a_rffs, "tumor_m6a_level_vs_H1sec_m6a_level", "Tumor_log2(mean(IP/Input))",  "H1ESC_log2(IP/Input)")
            
            #########################
            ## ggplot hexagnoal heatmap
            ## for normal
            dat <- data.frame(normal_r_m6a_s, m6a_rffs)
            reg <- lm(m6a_rffs ~ normal_r_m6a_s)
            g <- ggplot(dat, aes(normal_r_m6a_s, m6a_rffs))
            g <- g + geom_hex(bins = 50) +  geom_abline(intercept = reg$coefficients[1], slope = reg$coefficients[2], color="red", linetype="dashed")
            g <- g + theme_classic()
            ggsave("normal_m6a_level_vs_H1sec_m6a_level_Hexagonal.pdf", width = 4, height = 3)
            
            ## for tumor
            dat <- data.frame(tumor_r_m6a_s, m6a_rffs)
            reg <- lm(m6a_rffs ~ tumor_r_m6a_s)
            g <- ggplot(dat, aes(tumor_r_m6a_s, m6a_rffs))
            g <- g + geom_hex(bins = 50) +  geom_abline(intercept = reg$coefficients[1], slope = reg$coefficients[2], color="red", linetype="dashed")
            g <- g + theme_classic()
            ggsave("tumor_m6a_level_vs_H1sec_m6a_level_Hexagonal.pdf", width = 4, height = 3)
            
            ##################
            ## density plot  
            pdf(file = "mean_m6a_level_log2_density_with_h1esc_median.pdf", width = 4, height = 4)
            par(mar = c(3, 3, 1, 1), mgp =c(2, 0.7, 0))
            plot(density(tumor_r_m6a_s, na.rm = T), col = "coral2", xlab = "log2(mean(IP/Input))", main = "", xlim = c(-6, 6))
            lines(density(normal_r_m6a_s, na.rm = T), col = "cadetblue3")
            lines(density(m6a_rffs, na.rm = T), col = "gray", lty = 1)
            legend(2, 0.32, legend = c("Tumor", "Normal", "H1ESC"), col = c("coral2", "cadetblue3", "gray"), lty = c(1, 1, 1), cex = 0.8)
            abline(v = median(tumor_r_m6a_s, na.rm = T), col = "coral2", lty = 2)
            abline(v = median(normal_r_m6a_s, na.rm = T), col = "cadetblue3", lty = 2)
            abline(v = median(m6a_rffs, na.rm = T), col = "gray", lty = 2)
            dev.off()
        }
        
  
        ## m6a level variance per gene
        {
            normal_m6a_cv <- vector()
            tumor_m6a_cv <- vector()
            L <- nrow(m6a_rn)
            for (i in 1:L )
            {
                normal_m6a_cv[i] <- sd(m6a_rn[i, ]) / mean(m6a_rn[i, ])
                tumor_m6a_cv[i] <- sd(m6a_rt[i, ]) / mean(m6a_rt[i, ])
            }
            
            ks.test(normal_m6a_cv, tumor_m6a_cv, alternative = "greater")
            t.test(normal_m6a_cv, tumor_m6a_cv, alternative = "less")
            
            type <- c(rep("Normal", L), rep("Tumor", L))
            m6a_cv <- c(normal_m6a_cv, tumor_m6a_cv)
            dat <- data.frame(m6a_cv, type)
            
            g <- ggplot (dat, aes(x=m6a_cv, group = type, color = type)) + stat_ecdf(aes(linetype = type)) + theme_bw() + theme(axis.line = element_line(colour = "black"),panel.grid.major = element_blank(),panel.grid.minor = element_blank(),panel.border =  element_blank(),  panel.background = element_blank()) 
            g <- g + scale_color_manual( values = c("cadetblue3", "coral2")) + scale_linetype_manual(values=c("solid",  "solid"))
            g <- g + labs(x = "m6A level CV", y = "F(X)") + theme(legend.position = "none")  #theme(legend.position = c(0.85, 0.75))
            ggsave("m6a_level_cv_cdf.pdf",  width = 3, height = 3)
        }
        
    }
    
    #########################################
    ## samples m6A levels based correlation 
    {
        library("Hmisc")         
        library("corrplot")
        library("s2dverification")       # for color bar separately
        
        cor_m <- rcorr(m6a_rm)
        
        pdf("corrplot_based_on_tumor_normal_m6A.pdf", width = 6, height = 6)
        corrplot(cor_m$r,  order = "hclust", addrect = 20, method = "square")
        dev.off()
        
        ## color bar for sides after clustring 
        lims <- 1:62
        sample_c <- c(rep("coral2", 42), "cadetblue3", "coral2", "coral2", "cadetblue3", "cadetblue3", "coral2", "coral2", rep("cadetblue3", 7), rep("coral2", 5))
        
        pdf(file = "corrplot_based_on_tumor_normal_m6A_colorbar.pdf",  width = 3.5, heigh = 0.5)
        par(mar = c(0, 0, 0,0), mgp = c(0, 0, 0))
        ColorBar(lims, sample_c, vertical = F, draw_ticks = F, draw_separators = T, xaxt = "n", lwd = 0.5)
        dev.off()
    }
    
    #########################################################
    ## m6a level global comparison between tumor and normal
    ##  scatter plot for inter- and intra group comparion 
    {
        
        ## overlapping wiht % m6A level genes
        IN <- read.csv(file = "/Users/Yong/Yong/m6a_profiling/3_m6a_level/m6a_laic/m6a_laci_seq.csv")
        idx <- as.character(IN[, 8]) == "TRUE"
        m6a_p <- IN[idx, 4]                          # H1-ESC rep 1 
        names(m6a_p) <- IN[idx, 1]
        
        # normal vs tumor
        {
            cor.plot(normal_r_m6a, tumor_r_m6a , "Tumor_vs_Normal_m6A_level", "Mean m6A level (Normal)", "Mean m6A level (Tumor)")
            
            ## for highly methylated genes in both normal and tumor
            idx_hm <- normal_r_m6a > 2 & tumor_r_m6a > 2
            gene_hm <- names(normal_r_m6a)[idx_hm]
            idx_lm <- normal_r_m6a < -0 & tumor_r_m6a < -0
            gene_lm <- names(normal_r_m6a)[idx_lm]
            idx_mm <- normal_r_m6a < 2 & normal_r_m6a > 0 & tumor_r_m6a > 0 & tumor_r_m6a  < 2
            gene_mm <- names(normal_r_m6a)[idx_mm]
            
            idx_hm_p <- match(gene_hm , names(m6a_p))
            idx_lm_p <- match(gene_lm, names(m6a_p))
            idx_mm_p <- match(gene_mm, names(m6a_p))
            
            gene_hm_p <- m6a_p[idx_hm_p[!is.na(idx_hm_p)]]
            gene_lm_p <- m6a_p[idx_lm_p[!is.na(idx_lm_p)]]
            gene_mm_p <- m6a_p[idx_mm_p[!is.na(idx_mm_p)]]
            
            g_l <- c(length(gene_lm_p), length(gene_mm_p), length(gene_hm_p)) 
            g_l <- paste("N = ", g_l, sep = "")
            t_l <- c("m6A<0", "0<m6A<2", "m6A>2")
            
            #png(file = "m6a_percentage_level_for_3group_genes_in_both_tumor_and_normal.png", width = 4, height = 4, units = "in", res = 600)
            pdf(file = "m6a_percentage_level_for_3group_genes_in_both_tumor_and_normal.pdf", width = 3, height = 3)
            par(mar = c(3, 3, 0.5, 0.5), mgp = c(2, 0.75, 0))
            boxplot(list(gene_lm_p, gene_mm_p, gene_hm_p), ylab = "m6A Percentage(%)", xlab = "In both normal and tumor", xaxt = "n", ylim = c(0, 1.1),
                    col = c("steelblue4", "steelblue3", "steelblue1"), outcol = c("steelblue4", "steelblue3", "steelblue1"))
            axis(1, at = 1:3, labels = t_l) 
            text(x = 1:3, y = 1.05, labels = g_l)
            dev.off()
        }
        
        #################################
        # normal separated into two parts 
        M <- ncol(m6a_rn)
        M_m <- floor(M/2)
        n_mean1 <- log2(rowMeans(m6a_rn[, 1:M_m], na.rm = T))
        n_mean2 <- log2(rowMeans(m6a_rn[, (M_m + 1): M], na.rm = T))
        cor.plot(n_mean1, n_mean2, "Normal_intra_m6A_level", "Mean m6A level (Normal_p1)", "Mean m6A level (Normal_p2)")
        
        
        ## for highly methylated genes in normal only
        idx_hm <- normal_r_m6a > 2 
        gene_hm <- names(normal_r_m6a)[idx_hm]
        idx_lm <- normal_r_m6a < -0 
        gene_lm <- names(normal_r_m6a)[idx_lm]
        idx_mm <- normal_r_m6a < 2 & normal_r_m6a > 0 
        gene_mm <- names(normal_r_m6a)[idx_mm]
        
        idx_hm_p <- match(gene_hm , names(m6a_p))
        idx_lm_p <- match(gene_lm, names(m6a_p))
        idx_mm_p <- match(gene_mm, names(m6a_p))
        gene_hm_p <- m6a_p[idx_hm_p[!is.na(idx_hm_p)]]
        gene_lm_p <- m6a_p[idx_lm_p[!is.na(idx_lm_p)]]
        gene_mm_p <- m6a_p[idx_mm_p[!is.na(idx_mm_p)]]
        
        g_l <- c(length(gene_lm_p), length(gene_mm_p), length(gene_hm_p)) 
        g_l <- paste("N = ", g_l, sep = "")
        t_l <- c("m6A<0", "0<m6A<2", "m6A>2")
        #png(file = "m6a_percentage_level_for_3group_genes_in_normal.png", width = 4, height = 4, units = "in", res = 600)
        pdf(file = "m6a_percentage_level_for_3group_genes_in_normal.pdf", width = 3, height = 3)
        par(mar = c(3, 3, 0.5, 0.5), mgp = c(2, 0.75, 0))
        boxplot(list(gene_lm_p, gene_mm_p, gene_hm_p), ylab = "m6A Percentage(%)", xlab = "Normal", xaxt = "n", ylim = c(0, 1.1), 
                col = c("steelblue4", "steelblue3", "steelblue1"), outcol = c("steelblue4", "steelblue3", "steelblue1"))
        axis(1, at = 1:3, labels = t_l) 
        text(x = 1:3, y = 1.05, labels = g_l)
        dev.off()
        
        ##
        
        
        #################################
        # tumor separated into two parts
        M <- ncol(m6a_rt)
        M_m <- floor(M/2)
        t_mean1 <- log2(rowMeans(m6a_rt[, 1:M_m], na.rm = T))
        t_mean2 <- log2(rowMeans(m6a_rt[, (M_m + 1): M], na.rm = T))
        cor.plot(t_mean1, t_mean2, "Tumor_intra_m6A_level", "Mean m6A level (Tumor_p1)", "Mean m6A level (Tumor_p2)")
        
        idx_hm <- tumor_r_m6a > 2 
        gene_hm <- names(tumor_r_m6a)[idx_hm]
        idx_lm <- tumor_r_m6a < -0 
        gene_lm <- names(tumor_r_m6a)[idx_lm]
        idx_mm <- tumor_r_m6a < 2 & tumor_r_m6a > 0 
        gene_mm <- names(tumor_r_m6a)[idx_mm]
        
        idx_hm_p <- match(gene_hm , names(m6a_p))
        idx_lm_p <- match(gene_lm, names(m6a_p))
        idx_mm_p <- match(gene_mm, names(m6a_p))
        gene_hm_p <- m6a_p[idx_hm_p[!is.na(idx_hm_p)]]
        gene_lm_p <- m6a_p[idx_lm_p[!is.na(idx_lm_p)]]
        gene_mm_p <- m6a_p[idx_mm_p[!is.na(idx_mm_p)]]
        
        g_l <- c(length(gene_lm_p), length(gene_mm_p), length(gene_hm_p)) 
        g_l <- paste("N = ", g_l, sep = "")
        t_l <- c("m6A<0", "0<m6A<2", "m6A>2")
        #png(file = "m6a_percentage_level_for_3group_genes_in_tumor.png", width = 4, height = 4, units = "in", res = 600)
        pdf(file = "m6a_percentage_level_for_3group_genes_in_tumor.pdf", width = 3, height = 3)
        par(mar = c(3, 3, 0.5, 0.5), mgp = c(2, 0.75, 0))
        boxplot(list(gene_lm_p, gene_mm_p, gene_hm_p), ylab = "m6A Percentage(%)", xlab = "tumor", xaxt = "n", ylim = c(0, 1.1),
                col = c("steelblue4", "steelblue3", "steelblue1"), outcol = c("steelblue4", "steelblue3", "steelblue1"))
        axis(1, at = 1:3, labels = t_l) 
        text(x = 1:3, y = 1.05, labels = g_l)
        dev.off()    
        
    }
    
}


##############################################################
##  association between tumor gene m6A level and clinical info  
##############################################################
{
############################################
## tumor clinical information side colorbar
{
    idx <- match(colnames(m6a_rt), s_info$helab_id)
    s_info <- s_info[idx, ]
    
    N = length(idx)
    sex_c <- rep("pink", N)
    idx_t <- !is.na(match(s_info$sex, "M"))
    sex_c[idx_t] <- "blue"
    sex_c <- c(rep("grey", 10), sex_c)
    
    smoking_c <- rep("green", N)
    idx_t <- !is.na(match(s_info$smoking, "Current"))
    smoking_c[idx_t] <- "red"
    idx_t <- !is.na(match(s_info$smoking, "Ex-Smoker"))
    smoking_c[idx_t] <- "orange"
    smoking_c <- c(rep("grey", 10), smoking_c)
    
    tumor_c <- c(rep("purple", 10), rep("black", 51))
    col_s <- cbind(tumor_c, sex_c, smoking_c)
}

## globally heatmap
{
    hmcols<-colorRampPalette(c("blue","white","red"))(256)
    
    # for all selected genes
    png(file = "log2_m6a_level_heatmap_withcolbar.png", width = 4, height = 5, units = "in", res = 300)
    heatmap.plus(log2(m6a_rm), Rowv = T, Colv = T, scale = "row", col=hmcols, labCol = NA, ColSideColors=col_s, na.rm = T, margins = c(0.5, 0.5))
    dev.off()
    
    # for all selected protein coding genes
    idx_p <- match(rownames(m6a_rm), pc_gene$V2) 
    idx_p <- !is.na(idx_p)
    m6a_rm_p <- m6a_rm[idx_p, ]
    
    png(file = "log2_m6a_level_heatmap_withcolbar_protein.png", width = 4, height = 5, units = "in", res = 300)
    heatmap.plus(log2(m6a_rm_p), Rowv = T, Colv = T, scale = "row", col=hmcols, labCol = NA, ColSideColors=col_s, na.rm = T, margins = c(0.5, 0.5))
    dev.off()
    
    ## tumor-specific immune related genes 
    g_imm <- read.table("/Users/Yong/Yong/m6a_profiling/2_peak_calling/1_peaks/merge/bed30_61/2_micro_envi/Top2_MF_terms_unique_genes.txt", as.is = T)
    idx_imm <- match(rownames(m6a_rm), g_imm$V1) 
    m6a_rm_imm <- m6a_rm[!is.na(idx_imm), ]
    
    m6a_rm_imm_log2 <- log2(m6a_rm_imm)
    m6a_rm_imm_log2[is.infinite(m6a_rm_imm_log2)] <- NA
    
    hmcols<-colorRampPalette(c("blue","white","red"))(256)
    png(file = "log2_m6a_level_heatmap_tumor_specific_immune_top2_genes.png", width = 8, height = 8, units = "in", res = 300)
    heatmap.2(m6a_rm_imm_log2, Colv = F, scale = "row", density.info = "none", trace = "none", col = hmcols)
    dev.off()
    
    
}
}

#######################################
##  differentially m6A mehtylated genes 
#######################################
{
    library("Hmisc")         # for rcorr
    library(heatmap.plus)
    
    ##########################################
    ## load and pre-process data and funcitons 
    {
    load("/Users/Yong/Yong/m6a_profiling/3_m6a_level/m6a_level_keep_dup/gene_with_peak/m6A_level_4gene_passed_filtering.Rdata")  
    pc_gene <- read.table("/Users/Yong/Yong/m6a_profiling/1_mapping_quant/5_RNA-seq/pc_gene_list.txt")
    linc_gene <- read.table("/Users/Yong/Yong/m6a_profiling/1_mapping_quant/5_RNA-seq/lincRNA_gene_list.txt")
    s_info <- read.table("/Users/Yong/Yong/m6a_profiling/2_peak_calling/1_peaks/merge/sample_info.txt", header = T) 
    
    normal_m6a <- m6a_rn
    tumor_m6a <- m6a_rt
    m6a_rs <- cbind(m6a_rn, m6a_rt)
    
    ## m6a level and mRNA level fold change 
    m6a_fd_o <- log2(rowMeans(m6a_rt, na.rm = T) / rowMeans(m6a_rn, na.rm = T))        ## m6a ratio between normal and tumor
    m6a_mean_o <- log2((rowMeans(m6a_rt, na.rm = T) + rowMeans(m6a_rn, na.rm = T))/2)   ## m6a mean of normal and tumor
    
    mrna_fd_o <- log2(rowMeans(input_rpkm_t, na.rm = T) / rowMeans(input_rpkm_n, na.rm = T))
    names(m6a_fd_o) <- rownames(m6a_rt)
    names(mrna_fd_o) <- rownames(m6a_rt)
    save(m6a_fd_o, mrna_fd_o, file = "log2_m6a_fd_mrna_fd.RData")
    
    ### compared with peakcalling results
    g_match_cnt <- function(g_s, g_r)
    {
        idx <- match(g_s, g_r)
        idx_1 <- !is.na(idx)
        g_m <- g_s[idx_1]
        #print(length(g_m))
        return <- length(g_m)
    }
    
    }
    
    #################################
    ##  gene m6A  differenctial test  
    {
        ## for differential test: based on the IP/Input directly 
        {
            N <- nrow(tumor_m6a)
            mean_n <- vector()
            mean_t <- vector()
            pval <- vector()
            fd <- vector()
            pval_nn <- vector()
            pval_nt <- vector()
            
            for (i in 1:N) 
            {
                ## normality test
                tmp <- shapiro.test(normal_m6a[i, ])
                pval_nn[i] <- tmp$p.value
                
                tmp <- shapiro.test(tumor_m6a[i, ])
                pval_nt[i] <- tmp$p.value
                
                ## difference test 
                tmp <- wilcox.test(tumor_m6a[i, ], normal_m6a[i, ])
                #tmp <- t.test(tumor_m6a[i, ], normal_m6a[i, ])
                pval[i] <- tmp$p.value
                fd[i] <- mean(tumor_m6a[i, ], na.rm = T) / mean(normal_m6a[i, ], na.rm = T)        # fold-change  Tumor/normal !!!! 
                mean_t[i] <- mean(tumor_m6a[i, ], na.rm = T)
                mean_n[i] <- mean(normal_m6a[i, ], na.rm = T)
            }
            pval_a <- p.adjust(pval, method = "BH")         # multiple test adjustment 
        }
        
        m6a_de_out <- cbind(log2(mean_n), log2(mean_t), log2(fd), pval, pval_a) 
        rownames(m6a_de_out) <- rownames(normal_m6a)
        colnames(m6a_de_out) <- c("log2(mean_m6a_level_normal)", "log2(mean_m6a_level_tumor)", "log2(m6a_level_tumor/normla)", "p-value", "p-adj_BH")
        
        saveRDS(m6a_de_out, file = "m6a_de_test_all_8030_fc_pvalues.RDS")
        write.csv(m6a_de_out, file = "m6a_de_test_all_8030_fc_pvalues.csv")
        
        ## normality checking plots
        png(file = "normality_shapiro_wilk_test_for_each_gene_in_tumor.png", width = 4, height = 3, units = "in", res = 300)
        plot(density(pval_nt))
        abline(v = 0.05, col = "red", lty = 2)
        dev.off()
        
        png(file = "normality_shapiro_wilk_test_for_each_gene_in_normal.png", width = 4, height = 3, units = "in", res = 300)
        plot(density(pval_nn))
        abline(v = 0.05, col = "red", lty = 2)
        dev.off()
    }
    
    #################################
    ## significant gene with m6A diff
    ## choose fd = 1.5
    {
        m6a_cut = 3/2
        idx_1 <- pval_a < 0.05
        fd1 <- m6a_cut 
        fd2 <- 1/(m6a_cut)
        idx_2 <-  fd > fd1 | fd < fd2              # loose constrains
        idx_3 <- idx_1 & idx_2
        
        # non-diff and diff 
        m6a_nd <- m6a_rs[!idx_3, ]
        m6a_d <- m6a_rs[idx_3, ]
        fd_d <- fd[idx_3]
        pval_a_d <- pval_a[idx_3]
        
        # sorted based of m6A level fold change
        idx_s <- order(fd_d)
        fd_d_s <- fd_d[idx_s]
        pval_a_ds <- pval_a_d[idx_s]
        m6a_ds <- m6a_d[idx_s, ]
        
        m6a_gd <- cbind(fd_d_s, pval_a_ds)
        rownames(m6a_gd) <- rownames(m6a_ds)
        write.csv(m6a_gd, file = "all_gene_with_diff_m6a_level.csv", quote = F)
        
        idx_u <- fd_d_s > 1
        idx_d <- fd_d_s < 1
        g_u <- rownames(m6a_ds)[idx_u]
        g_d <- rownames(m6a_ds)[idx_d]
        write.table(g_u, file = "all_gene_with_diff_m6a_level_upregulated.txt", quote = F, col.names = F, row.names = F)
        write.table(g_d, file = "all_gene_with_diff_m6a_level_downregulated.txt", quote = F, col.names = F, row.names = F)
        
        ### diff-m6A gene compared to the m6a peak calling resutls 
        g_t <- read.table("/Users/Yong/Yong/m6a_profiling/2_peak_calling/1_peaks/merge/bed30_61/gene_with_tumor_only_peak.txt", as.is  = T)
        g_n <-  read.table("/Users/Yong/Yong/m6a_profiling/2_peak_calling/1_peaks/merge/bed30_61/gene_with_normal_only_peak.txt", as.is = T)
        g_c <- read.table("/Users/Yong/Yong/m6a_profiling/2_peak_calling/1_peaks/merge/bed30_61/gene_with_common_only_peak.txt", as.is = T)
        
        g_peak <- list(g_n$V1, g_c$V1, g_t$V1)
        L_peak <- length(g_peak)
        
        g_diff_cnt <- matrix(0, 2, L_peak + 1)
        colnames(g_diff_cnt) <- c("Normal_only", "Common_only", "Tumor_only", "others")
        rownames(g_diff_cnt) <- c("Hyper", "Hypo")
        
        for(i in 1: L_peak)
        {
            g_diff_cnt[1, i] <- g_match_cnt(g_u, g_peak[[i]])
            g_diff_cnt[2, i] <- g_match_cnt(g_d, g_peak[[i]])
        }
        
        g_diff_cnt[1, L_peak + 1] <- length(g_u) - sum(g_diff_cnt[1, 1:L_peak])
        g_diff_cnt[2, L_peak + 1] <- length(g_d) - sum(g_diff_cnt[2, 1:L_peak])
        
        write.table(g_diff_cnt, file = "all_gene_with_diff_m6a_level_cmp2_peakcalling_res.txt", quote = F, col.names = T, row.names = T, sep = "\t" )
        
        idx_uf <- is.na(match(g_u, g_n$V1))            ## rm 21
        g_uf <- g_u[idx_uf]
        idx_df <- is.na(match(g_d, g_t$V1))
        g_df <- g_d[idx_df]                            ## rm 43
        
        g_amb <- c(g_u[!idx_uf], g_d[!idx_df])
        g_amb_t <- c(rep("g_u_with_normal_only_peak", length(g_u[!idx_uf])), rep("g_d_with_tumor_only_peak", length(g_d[!idx_df])))
        g_amb_m <- data.frame(g_amb, g_amb_t) 
        
        ## exclude genes show opposite direction with peak calling results
        write.table(g_uf, file = "all_gene_with_diff_m6a_level_upregulated_corrected.txt", quote = F, col.names = F, row.names = F)
        write.table(g_df, file = "all_gene_with_diff_m6a_level_downregulated_corrected.txt", quote = F, col.names = F, row.names = F)
        write.table(g_amb_m, file = "all_gene_with_diff_m6a_level_corrected_amb.txt", quote = F, col.names = T, row.names = F, sep = "\t" )
        
        #################### without corretion #######
        # selected protein coding genes
        idx_p <- !is.na(match(rownames(m6a_ds), pc_gene$V2)) & is.na(match(rownames(m6a_ds), g_amb))
        m6a_dsp <- m6a_ds[idx_p, ]     
        fd_d_sp <- fd_d_s[idx_p]
        pval_a_dsp <- pval_a_ds[idx_p]
        m6a_gdp <- cbind(fd_d_sp, pval_a_dsp)
        rownames(m6a_gdp) <- rownames(m6a_dsp)
        write.csv(m6a_gdp, file = "pc_gene_with_diff_m6a_level.csv", quote = F)
        
        # selected lincRNA genes 
        idx_l <- !is.na(match(rownames(m6a_ds), linc_gene$V2)) & is.na(match(rownames(m6a_ds), g_amb))
        m6a_dsl <- m6a_ds[idx_l, ]     
        fd_d_sl <- fd_d_s[idx_l]
        pval_a_dsl <- pval_a_ds[idx_l]
        m6a_gdl <- cbind(fd_d_sl, pval_a_dsl)
        rownames(m6a_gdl) <- rownames(m6a_dsl)
        write.csv(m6a_gdl, file = "lincRNA_gene_with_diff_m6a_level.csv", quote = F)
        
    } 
    
    ######################################################
    ## heatmap for gene with differential m6A level (log2) 
    ## after filtering
    {
        library(RColorBrewer)
        col_l <- brewer.pal(n = 10, name = "RdBu")      
        
        sample_c <- c(rep("cadetblue3", 10), rep("coral2", 51))
        gene_c <- c(rep(col_l[8], length(g_df)), rep(col_l[2], length(g_uf)))
        idx_ff <- is.na(match(rownames(m6a_ds), g_amb))
        
        m6a_ds_log2 <- log2(m6a_ds[idx_ff, ])
        idx_inf <- is.infinite(m6a_ds_log2)
        m6a_ds_log2[idx_inf] <- NA
        hmcols<-colorRampPalette(c("blue","white","red"))(50)

        png(file = "m6a_level_DE_withcolbar.png", width = 5, height = 6, units = "in", res = 600)
        heatmap.2(m6a_ds_log2, Colv= T, Rowv= F, dendrogram = "column" , scale="row", na.color = "grey",
                  trace ="none", density.info="none", col=hmcols, RowSideColors= gene_c, ColSideColors= sample_c, colsep = 51, rowsep = length(g_df),
                  lhei = c(1,5),  lwid = c(1.3, 4), margins=c(0.5, 0.5), labRow = FALSE, labCol = F,  offsetCol = 0, key.title = NA)
        dev.off()
        
        sum(fd_d_s[idx_ff] > fd1)
        sum(fd_d_s[idx_ff] < fd2)
        
        # all m6A diff protein coding genes with filtering
        idx_pp <- !is.na(match(rownames(m6a_ds_log2), pc_gene$V2))
        png(file = "m6a_level_DE_withcolbar_protein.png", width = 5, height = 6, units = "in", res = 600)
        heatmap.2(m6a_ds_log2[idx_pp, ], Colv= T, Rowv= F, dendrogram = "column" , scale="row", na.color = "grey",
                  trace ="none", density.info="none", col=hmcols, RowSideColors= gene_c[idx_pp], ColSideColors= sample_c, colsep = 51, rowsep = sum(gene_c[idx_pp] == col_l[8]),
                  lhei = c(1,5),  lwid = c(1.3, 4), margins=c(0.5, 0.5), labRow = FALSE, labCol = F,  offsetCol = 0, key.title = NA)
        dev.off()
        
        # all m6A diff lincRNA genes
        idx_ll <- !is.na(match(rownames(m6a_ds_log2), linc_gene$V2))
        png(file = "m6a_level_DE_withcolbar_lincRNA.png", width = 6, height = 6, units = "in", res = 600)
        heatmap.2(m6a_ds_log2[idx_ll, ], Colv= T, Rowv= F, dendrogram = "column" , scale="row", na.color = "grey",
                  trace ="none", density.info="none", col=hmcols, RowSideColors= gene_c[idx_ll], ColSideColors= sample_c, rowsep = sum(gene_c[idx_ll] == col_l[8]),
                  lhei = c(1,5),  lwid = c(1.3, 4), margins=c(0.5, 8), labRow = rownames(m6a_ds_log2[idx_ll, ]), labCol = F,  offsetCol = 0, key.title = NA)
        dev.off()
        
        # all m6A non-diff genes (log2(m6a level))  :: without m6a-diff and m6a-amb
        png(file = "m6a_level_non_DE_withcolbar.png", width = 4, height = 5, units = "in", res = 300)
        hmcols<-colorRampPalette(c("blue","white","red"))(256)
        heatmap.plus(log2(m6a_nd), Rowv = T, Colv = T, scale = "row", col=hmcols, labCol = NA, ColSideColors=col_s, na.rm = T, margins = c(0.5, 0.5))
        #legend(0, 3, legend = c("Male", "Female", "NA"), fill = c("blue", "pink", "gray"), border = F, bty="n", y.intersp = 0.7, cex=0.7)
        dev.off()
    }
}


########################################
## gene m6A level diff vs RNA level diff 
########################################
{

########################################################################
##  means of m6A level and RNA level comparison between normal and tumor
{
    ## after filtering
    m6a_ds <-  m6a_ds[idx_ff, ]     
    input_rpkm <- cbind(input_rpkm_n, input_rpkm_t)
    idx_1 <- match(rownames(m6a_ds), rownames(input_rpkm))
    idx_2 <- match(rownames(m6a_nd), rownames(input_rpkm))
    mrna_ds <- input_rpkm[idx_1, ]                           # mrna_sorted
    mrna_nd <- input_rpkm[idx_2, ]
    
    m6a_s <- log2(rbind(m6a_nd, m6a_ds))           # log2 transformed
    m6a_sg <- c(rep("normal", 10), rep("tumor", 51))
    mrna_s <- log2(rbind(mrna_nd, mrna_ds) + 1)     # log2 transformed
    
    idx_u <- match(g_uf, rownames(m6a_ds))
    idx_d <- match(g_df, rownames(m6a_ds))
    col_m <- rep("red",nrow(m6a_ds))
    col_m[idx_d]  <- "blue"
    
    # mean m6A level scatterplot
    m6a_me_t <-  log2(rowMeans(m6a_ds[, 11:61], na.rm = T))
    m6a_me_n <-  log2(rowMeans(m6a_ds[, 1:10], na.rm = T))
    
    png(file = "Mean_m6A_level_in_Normal_vs_Tumor.png",  width = 4, height = 4, units = "in", res = 600)
    par(mar = c(3, 3, 0.5, 0.5), mgp = c(1.5, 0.5, 0))
    plot(m6a_me_n, m6a_me_t, col = col_m, 
         xlim = c(-5, 5), ylim = c(-5, 5), xlab = "Mean m6A level in Normal", ylab = "Mean m6A level in Tumor")
    abline(a = 0, b = 1, lty = 2, col = "gray")
    abline(v = c(0, 2), h = c(0, 2), lty = 2, col = "gray", lwd = 0.5)
    dev.off()
    
    # boxplot
    pdf(file = "Mean_m6A_level_in_Normal_vs_Tumor_boxplot_flipped.pdf",  width = 3, height = 3)
    par(mar = c(3, 3, 0.5, 0.5), mgp = c(2, 0.75, 0))
    boxplot(list(m6a_me_n[idx_d], m6a_me_t[idx_d], m6a_me_n[idx_u], m6a_me_t[idx_u]), col = rep(c("cadetblue3", "coral2"), 2), ylab = "Mean m6A level")
    #abline( h = 0, lty = 2)
    abline(v = 2.5, lty = 2)
    dev.off()
    
    wilcox.test(m6a_me_n[idx_u], m6a_me_t[idx_u], alternative = "less")
    wilcox.test(m6a_me_n[idx_d], m6a_me_t[idx_d], alternative = "greater")
    wilcox.test(m6a_me_n[idx_u], m6a_me_n[idx_d], alternative = "less")
    wilcox.test(m6a_me_t[idx_u], m6a_me_t[idx_d], alternative = "less")

    ## mean mRNA level
    mrna_me_t <-  log2(rowMeans(mrna_ds[, 11:61], na.rm = T) + 1)
    mrna_me_n <-  log2(rowMeans(mrna_ds[, 1:10], na.rm = T) + 1)
    
    png(file = "Mean_mRNA_in_Normal_vs_Tumor.png",  width = 4, height = 4, units = "in", res = 600)
    par(mar = c(3, 3, 0.5, 0.5), mgp = c(1.5, 0.5, 0))
    plot(mrna_me_n, mrna_me_t , col = col_m, xlab = "Mean mRNA in Normal", ylab = "Mean mRNA in Tumor")
    abline(a = 0, b = 1, lty = 2, col = "gray")
    dev.off()
    
    # boxplot
    pdf(file = "Mean_mRNA_in_Normal_vs_Tumor_boxplot_flipped.pdf",  width = 3, height = 3)
    par(mar = c(3, 3, 0.5, 0.5), mgp = c(2, 0.75, 0))
    boxplot(list( mrna_me_n[idx_d], mrna_me_t[idx_d], mrna_me_n[idx_u], mrna_me_t[idx_u]), col = rep(c("cadetblue3", "coral2"), 2), ylab = "Mean mRNA level")
    abline( v = 2.5, lty = 2)
    dev.off()
    
    wilcox.test(mrna_me_n[idx_u], mrna_me_t[idx_u])
    wilcox.test(mrna_me_n[idx_d], mrna_me_t[idx_d])
    wilcox.test(mrna_me_n[idx_u], mrna_me_n[idx_d], "greater")
    wilcox.test(mrna_me_t[idx_u], mrna_me_t[idx_d], "greater")
    
}

##################################
##### m6A level FC vs mRNA FC only
{
         ##############################
         ## all genes m6A FC vs RNA FC
         {
            m6a_sr <- rbind(m6a_nd, m6a_ds)
            mrna_sr <- rbind(mrna_nd, mrna_ds)
            
            m6a_fd <- log2(rowMeans(m6a_sr[, 11:61], na.rm = T) / rowMeans(m6a_sr[, 1:10], na.rm = T))
            mrna_fd <- log2(rowMeans(mrna_sr[, 11:61], na.rm = T) / rowMeans(mrna_sr[, 1:10], na.rm = T))
            
            save(m6a_fd, mrna_fd, file = "log2_m6a_fd_mrna_fd_sorted_by_de_and_nd.RData")
            
            idx_u1 <- !is.na(match(rownames(m6a_sr), g_uf))
            idx_d1 <- !is.na(match(rownames(m6a_sr), g_df))
            idx_u2 <- abs(mrna_fd) < 1  & idx_u1
            idx_d2 <- abs(mrna_fd) < 1  & idx_d1
            
            col_mf <- col_m[idx_ff]
            col_mf[col_mf == "blue"]  = col_l[8]
            col_mf[col_mf == "red"]  = col_l[2]
            col_c <- c(rep("gray", nrow(m6a_nd)), col_mf)
            col_c[idx_u2] <- col_l[2]   #"indianred1"       #
            col_c[idx_d2] <- col_l[8]   # "royalblue1"
            
            a <- sum(mrna_fd >=1 & idx_u1 )
            b <- sum(mrna_fd <= -1 & idx_u1)
            c <- sum(mrna_fd >=1 & idx_d1)
            d <- sum(mrna_fd <= -1 & idx_d1)
            mm <- matrix(c(a, b, c, d), 2, 2)
            chisq.test(mm)
            
            e <- sum(idx_u2)
            f <- sum(idx_d2)
            
            png(file = "log2FC_m6a_level_vs_mrna_without_cnt.png", width = 3, height = 3, units = "in", res = 600)
            par(mar = c(3, 3, 0.5, 0.5), mgp =c(2, 0.75, 0))
            plot(mrna_fd, m6a_fd, col = col_c, xlab = "log2FC_mRNA", ylab = "log2FC_m6A", ylim = c(-3, 3))
            abline(h = c(log2(m6a_cut), -log2(m6a_cut)), v= c(1, -1), lty = 2, lwd = 0.75)
            #text(x = c(5, -5, 5, -5), y = c(2.5, 2.5, -2.5, -2.5), labels = c(a, b, c, d), col = "black")
            #text(x = c(0, 0), y = c(2.5, -2.5), labels = c(e, f), col = "orange")
            dev.off()
        }
        
        #################################################################
        ## mRNA no big differene but m6a level has significant difference
        ## sm: slected middle area genes 
        {   
            m6a_sm_u <- m6a_sr[idx_u2, ]
            mrna_sm_u <- mrna_sr[idx_u2, ]
            m6a_sm_d <- m6a_sr[idx_d2, ]
            mrna_sm_d <- mrna_sr[idx_d2, ]
            
            g_sm <- c(rownames(m6a_sm_u), rownames(m6a_sm_d))
            write.table(g_sm, file = "all_gene_with_diff_m6a_level_corrected_nodiff_mrna.txt", quote = F, col.names = F, row.names = F)
            
            ## mRNA fd for stringent 
            mrna_sm_u_fd <- log2(rowMeans(mrna_sm_u[, 11:61], na.rm = T) / rowMeans(mrna_sm_u[, 1:10], na.rm = T))
            mrna_sm_d_fd <- log2(rowMeans(mrna_sm_d[, 11:61], na.rm = T) / rowMeans(mrna_sm_d[, 1:10], na.rm = T))
            
            idx_sm_u_fd <- abs(mrna_sm_u_fd) < 0.05
            idx_sm_d_fd <- abs(mrna_sm_d_fd) < 0.05
            sm_lu <- sum(idx_sm_u_fd)
            sm_ld <- sum(idx_sm_d_fd)
            col_sm_u <- c(rep("coral2", sm_lu), rep("cadetblue3", sm_lu))
            col_sm_d <- c(rep("coral2", sm_ld), rep("cadetblue3", sm_ld))
            
            ## for hyper
            mrna_sm_u_mean <- log2(rowMeans(mrna_sm_u[idx_sm_u_fd, ], na.rm = T))    ## mean across normal and tumor 
            m6a_sm_u_t_mean <- log2(rowMeans(m6a_sm_u[idx_sm_u_fd, 11:61], na.rm = T))  ## for tumor 
            m6a_sm_u_n_mean <- log2(rowMeans(m6a_sm_u[idx_sm_u_fd, 1:10], na.rm = T))  ## for tumor 
            
            pdf(file = "hyper_genes_with_stringent_nodiff_mrna.pdf", width = 3, height = 3)
            par(mar = c(3, 3, 0.5, 0.5), mgp =c(2, 0.75, 0))
            plot(rep(mrna_sm_u_mean, 2), c(m6a_sm_u_t_mean, m6a_sm_u_n_mean), col = col_sm_u, pch = 16, xlim = c(0, 6.5))
            ## draw link lines
            for(i in 1: sm_lu)
            {
                lines(x = c(mrna_sm_u_mean[i], mrna_sm_u_mean[i]), y = c(m6a_sm_u_t_mean[i] - 0.06, m6a_sm_u_n_mean[i] + 0.06), lty =2 , col = col_l[2])
            }
            text(x = mrna_sm_u_mean, y = m6a_sm_u_t_mean, labels = rownames(mrna_sm_u)[idx_sm_u_fd], pos = 1, cex = 0.5)
            dev.off()
            
            ## for hypo
            mrna_sm_d_mean <- log2(rowMeans(mrna_sm_d[idx_sm_d_fd, ], na.rm = T))    ## mean across normal and tumor 
            m6a_sm_d_t_mean <- log2(rowMeans(m6a_sm_d[idx_sm_d_fd, 11:61], na.rm = T))  ## for tumor 
            m6a_sm_d_n_mean <- log2(rowMeans(m6a_sm_d[idx_sm_d_fd, 1:10], na.rm = T))  ## for tumor 
            
            pdf(file = "hypo_genes_with_stringent_nodiff_mrna.pdf", width = 3, height = 3)
            par(mar = c(3, 3, 0.5, 0.5), mgp =c(2, 0.75, 0))
            plot(rep(mrna_sm_d_mean, 2), c(m6a_sm_d_t_mean, m6a_sm_d_n_mean), col = col_sm_d, pch = 16, xlim = c(0, 5))
            ## draw link lines
            for(i in 1: sm_ld)
            {
                lines(x = c(mrna_sm_d_mean[i], mrna_sm_d_mean[i]), y = c(m6a_sm_d_t_mean[i] + 0.06, m6a_sm_d_n_mean[i] - 0.06), lty =2 , col = col_l[8])
            }
            text(x = mrna_sm_d_mean, y = m6a_sm_d_n_mean, labels = rownames(mrna_sm_d)[idx_sm_d_fd], pos = 1, cex = 0.5)
            dev.off()
            
            g_sm_s <- c(rownames(mrna_sm_u)[idx_sm_u_fd], rownames(mrna_sm_d)[idx_sm_d_fd]) 
        }
        
        ######################
        ### plot lincRNA only 
        {
            idx_li <- !is.na(match(names(mrna_fd), linc_gene$V2))
            idx_lic <- col_c == col_l[8] | col_c == col_l[2]
            idx_lis <- idx_li & idx_lic
            
            pdf(file = "log2FC_m6a_level_vs_mrna_lincRNA.pdf", width = 4, height = 4)
            par(mar = c(3, 3, 0.5, 0.5), mgp =c(2, 0.7, 0))
            plot(mrna_fd[idx_li], m6a_fd[idx_li], col = col_c[idx_li], xlab = "log2FC_mRNA", ylab = "log2FC_m6A", ylim = c(-2, 2), xlim = c(-6, 6))
            abline(h = c(log2(m6a_cut), -log2(m6a_cut)), v= c(1, -1), lty = 2, lwd = 0.75)
            text(x = mrna_fd[idx_lis], y = m6a_fd[idx_lis], labels = names(mrna_fd)[idx_lis], col = "black", pos = 3, cex = 0.5)
            dev.off()
            
            ## lincRNA no big difference at mRNA level and with significant difference at m6A level 
            linc_genes <- names(mrna_fd)[idx_lis]
            idx_ls <- match(linc_genes, names(m6a_me_n))
            
            linc_RNA_o <- cbind(mrna_me_n[idx_ls],  mrna_me_t[idx_ls], m6a_me_n[idx_ls],  m6a_me_t[idx_ls],  mrna_fd[idx_lis], m6a_fd[idx_lis])
            colnames(linc_RNA_o) <- c("log2(mean_mRNA_normal)", "log2(mean_mRNA_tumor)",
                                      "log2(mean_m6A_normal)", "log2(mean_m6A_tumor)",
                                      "log2(mean_mRNA_FC)", "log2(mean_m6A_FC))")
            write.csv(linc_RNA_o, file = "lincRNA_sig_m6a_diff_with_no_diff_mRNA.csv") 
        }
    
}
    
############################################
### m6A level FC vs mRNA FC with DEG p value 
{
    ###########################
    ## based own data
    {
        de_res <- read.csv("/Users/Yong/Yong/m6a_profiling/1_mapping_quant/5_RNA-seq/htseq/count/61_sample_Input_DESeq_DE.csv")
        
        g_ff <- c(g_uf, g_df)                      ##!!! m6a_ds and mrna_ds sorted by g_df--> g_uf; to draw g_uf first then g_df
        idx_a <- match(g_ff, de_res$X)
        
        library(RColorBrewer)
        col_l <- brewer.pal(n = 10, name = "RdBu")  
        
        p_ff <- -log10(de_res$padj[idx_a])
        fc_ff <- de_res$log2FoldChange[idx_a]
        
        idx1 <- abs(fc_ff) <= 1 | p_ff <= -log10(0.05)
        idx_gu <- as.logical(c(rep(1, length(g_uf)), rep(0, length(g_df))))
        idx_gd <- as.logical(c(rep(0, length(g_uf)), rep(1, length(g_df))))
        
        sum(idx1)
        dm_gene_nodiff_rna <- rownames(de_res[idx_a, ])[idx1]
        
        sum(idx_gu & p_ff > -log10(0.05) & fc_ff < -1)
        sum(idx_gd & p_ff > -log10(0.05) & fc_ff < -1)
        sum(idx_gu & p_ff > -log10(0.05) & fc_ff > 1)
        sum(idx_gd & p_ff > -log10(0.05) & fc_ff > 1)
        
        col_gff <- c(rep(col_l[2], length(g_uf)), rep(col_l[8], length(g_df)))
        col_gff[idx1&idx_gu] <- col_l[5]
        col_gff[idx1&idx_gd] <- col_l[6]
        
        pdf(file = "DMG_mRNA_vocanal_plot.pdf", width = 3, height = 3)
        par(mar = c(3, 3, 0.5, 0.5), mgp =c(2, 0.75, 0))
        plot(fc_ff, p_ff, col =  col_gff , xlab = "log2(mRNA_FC)", ylab = "-log10(FDR)", xlim = c(-7, 7), pch = 16, cex = 0.75)
        abline(h = -log10(0.05), lty = 2, col = "gray")
        abline(v = c(-1, 1), lty = 2, col = "gray")
        dev.off()
        
        ## out put diff_m6A genes m6a level FC , p and mRNA level, logfFC p 
        g_ffs <- c(g_df, g_uf)
        idx_o1 <- match(g_ffs, rownames(m6a_de_out))
        idx_o2 <- match(g_ffs, de_res$X)
        idx_o3 <- match(g_ffs, names(mrna_me_n))
        
        m6a_de_m <- m6a_de_out[idx_o1, ]
        mrna_de_m <- de_res[idx_a, ]
        mrna_de_out <-cbind(mrna_me_n[idx_o3], mrna_me_t[idx_o3], mrna_de_m$log2FoldChange, mrna_de_m$pvalue, mrna_de_m$padj)
        colnames(mrna_de_out) <- c("log2(mean_mrna_level_normal)", "log2(mean_mrna_level_tumor)", "log2(mrna_level_tumor/normla)", "p-value_deseq2", "p-adj_deseq2")
        
        ## merge m6a_de and correspodning 
        out <- cbind(m6a_de_m, mrna_de_out)
        saveRDS(out, file = "Corrected_m6A_DE_gene_and_coressponding_mRNA_DESeq2_res.RDS")
        write.csv(out, file = "Corrected_m6A_DE_gene_and_coressponding_mRNA_DESeq2_res.csv")
    }
    
    ###########################
    ### based on TCGA RNA data
    {
        de_res_tcga <- read.csv("/Users/Yong/Yong/m6a_profiling/1_mapping_quant/5_RNA-seq/htseq/count/tcga_sample_DESeq_DE.csv")
        
        ##### impact vs tcga 
        ## 
        idx_tcga <- match(rownames(de_res_tcga), rownames(de_res))
        library(epitrans)
        corr_scatter_hex(de_res_tcga$log2FoldChange, de_res$log2FoldChange[idx_tcga], "pearson", "TCGA_log2FC",
                         "IMPACT_log2FC", 50, "TCGA_log2FC_vs_IMPACT_log2FC_tumor_against_normal.pdf")
        
        g_ff <- c(g_uf, g_df)                      ##!!! m6a_ds and mrna_ds sorted by g_df--> g_uf; to draw g_uf first then g_df
        idx_a <- match(g_ff, de_res_tcga$X)
        
        library(RColorBrewer)
        col_l <- brewer.pal(n = 10, name = "RdBu")  
        
        p_ff <- -log10(de_res_tcga$padj[idx_a])
        fc_ff <- de_res_tcga$log2FoldChange[idx_a]
        
        idx1 <- abs(fc_ff) <= 1 | p_ff <= -log10(0.05)
        idx_gu <- as.logical(c(rep(1, length(g_uf)), rep(0, length(g_df))))
        idx_gd <- as.logical(c(rep(0, length(g_uf)), rep(1, length(g_df))))
        
        sum(idx1, na.rm = T)
        dm_gene_nodiff_rna_tcga <- (de_res_tcga[idx_a, 1])[idx1]
        
        sum(idx_gu & p_ff > -log10(0.05) & fc_ff < -1, na.rm = T)
        sum(idx_gd & p_ff > -log10(0.05) & fc_ff < -1, na.rm = T)
        sum(idx_gu & p_ff > -log10(0.05) & fc_ff > 1, na.rm = T)
        sum(idx_gd & p_ff > -log10(0.05) & fc_ff > 1, na.rm = T)
        
        col_gff <- c(rep(col_l[2], length(g_uf)), rep(col_l[8], length(g_df)))
        col_gff[idx1&idx_gu] <- col_l[5]
        col_gff[idx1&idx_gd] <- col_l[6]
        
        pdf(file = "DMG_mRNA_TCGA_vocanal_plot.pdf", width = 3, height = 3)
        par(mar = c(3, 3, 0.5, 0.5), mgp =c(2, 0.75, 0))
        plot(fc_ff, p_ff, col =  col_gff , xlab = "log2(mRNA_FC)_TCGA", ylab = "-log10(FDR)", xlim = c(-2.5, 2.5), pch = 16, cex = 0.75)
        abline(h = -log10(0.05), lty = 2, col = "gray")
        abline(v = c(-1, 1), lty = 2, col = "gray")
        dev.off()
        
    }
}
 
}