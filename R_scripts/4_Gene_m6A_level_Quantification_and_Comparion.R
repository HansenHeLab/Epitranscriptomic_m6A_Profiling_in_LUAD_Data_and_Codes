################################################################
## This script implements:
## 1. gene_m6A_level_and_filtering
## 2. gene m6A level vs H1ESC percentage m6A level
## 3. gene m6A level distribution and global diff 
## 4. differential m6A methylated genes 
## 5. differential m6A analysis sub sampling
## 6. gene m6A level diff vs RNA level diff 
## 7. EML4 related analyses
#################################################################

rm(list = ls())
setwd("./data")
               
library("matrixStats")
library(ggplot2)

###############################
## gene_m6A_level_and_filtering
###############################
{
#### load data ######################################
# normalized by ERCC and DESeq with all duplications
{
    load("./paired_63_count_ercc_deseq.Rdata") 
    
    ## remove outliers tumor 45 and tumor9
    idx_s <- match(c("tumor9", "tumor45"), colnames(ip_rpkm))
   
    ip_count_norm <- ip_count_norm[, -idx_s]
    input_count_norm <- input_count_norm[, -idx_s]
    
    ip_rpkm_norm <- ip_rpkm[, -idx_s]
    input_rpkm_norm <- input_rpkm[, -idx_s]
    
    gene_a <- read.table("./0_merged_peaks/gene_with_m6a_peak_in_tumor_or_normal.txt")    # all gene with m6a Peak  
    
    pc_gene <- read.table("pc_gene_list.txt")
    lincRNA <- read.table("lincRNA_gene_list.txt")
    
    #s_info <- read.table("/Users/Yong/Yong/m6a_profiling/2_peak_calling/1_peaks/merge/sample_info.txt", header = T)       # sample information 
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
    pdf(file = "../Results/genes_filtering_based_RPKM_and_peaks.pdf", width = 3.5, height = 3)
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
    gene.ty(gene_m6a, pc_gene, lincRNA, "../Results/Gene_with_m6A_RPKM>1")
   # gene.ty(gene_nm6a, pc_gene, lincRNA, "../Results/Gene_without_m6A_RPKM>1")
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

#save(m6a_rn_all, m6a_rt_all, input_count_n_all, input_count_t_all, ip_count_n_all, ip_count_t_all,
#     input_rpkm_n_all, input_rpkm_t_all, ip_rpkm_n_all, ip_rpkm_t_all, file = "m6A_level_4gene_all.Rdata")

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
## comparison with the RPKM based IP / Input
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

pdf(file = "../Results/mean_m6a_level_log2_boxplot.pdf", width = 4, height = 3)
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
load("H1sec_m6a_level.RData")

idx_s <- match(names(m6a_rff), names(tumor_r_m6a))
sum(!is.na(idx_s))
m6a_rffs <- m6a_rff[!is.na(idx_s)]
tumor_r_m6a_s <- tumor_r_m6a[idx_s[!is.na(idx_s)]]
normal_r_m6a_s <- normal_r_m6a[idx_s[!is.na(idx_s)]]

## scatter plot
source("../R_scripts/functions/cor.plot.R")
cor.plot(normal_r_m6a_s, m6a_rffs, "../Results/Normal_m6a_level_vs_H1sec_m6a_level", "Normal_log2(mean(IP/Input))",  "H1ESC_log2(IP/Input)", "png")
cor.plot(tumor_r_m6a_s, m6a_rffs, "../Results/Tumor_m6a_level_vs_H1sec_m6a_level", "Tumor_log2(mean(IP/Input))",  "H1ESC_log2(IP/Input)", "png")

## scatter plot function
dat <- data.frame(normal_r_m6a_s, m6a_rffs)
reg <- lm(m6a_rffs ~ normal_r_m6a_s)
g <- ggplot(dat, aes(normal_r_m6a_s, m6a_rffs))
g <- g + geom_hex(bins = 50) +  geom_abline(intercept = reg$coefficients[1], slope = reg$coefficients[2], color="red", linetype="dashed")
g <- g + theme_classic()
ggsave("../Results/Normal_m6a_level_vs_H1sec_m6a_level_Hexagonal.pdf", width = 4, height = 3)

## for tumor
dat <- data.frame(tumor_r_m6a_s, m6a_rffs)
reg <- lm(m6a_rffs ~ tumor_r_m6a_s)
g <- ggplot(dat, aes(tumor_r_m6a_s, m6a_rffs))
g <- g + geom_hex(bins = 50) +  geom_abline(intercept = reg$coefficients[1], slope = reg$coefficients[2], color="red", linetype="dashed")
g <- g + theme_classic()
ggsave("../Results/Tumor_m6a_level_vs_H1sec_m6a_level_Hexagonal.pdf", width = 4, height = 3)



#########################################################
## m6a level global comparison between tumor and normal
##  scatter plot for inter- and intra group comparion 
{
  
  ## overlapping wtih % m6A level genes
  IN <- read.csv(file = "H1sec_percentage_m6a_level.csv")
  idx <- as.character(IN[, 8]) == "TRUE"
  m6a_p <- IN[idx, 4]                          # H1-ESC rep 1 
  names(m6a_p) <- IN[idx, 1]
  
  # normal vs tumor
  {
    #cor.plot(normal_r_m6a, tumor_r_m6a , "Tumor_vs_Normal_m6A_level", "Mean m6A level (Normal)", "Mean m6A level (Tumor)")
    
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
    pdf(file = "../Results/m6a_percentage_level_for_3group_genes_in_both_tumor_and_normal.pdf", width = 3, height = 3)
    par(mar = c(3, 3, 0.5, 0.5), mgp = c(2, 0.75, 0))
    boxplot(list(gene_lm_p, gene_mm_p, gene_hm_p), ylab = "m6A Percentage(%)", xlab = "In both normal and tumor", xaxt = "n", ylim = c(0, 1.1),
            col = c("steelblue4", "steelblue3", "steelblue1"), outcol = c("steelblue4", "steelblue3", "steelblue1"))
    axis(1, at = 1:3, labels = t_l) 
    text(x = 1:3, y = 1.05, labels = g_l)
    dev.off()
  }
  
  
}




}

################################################
## m6A level comparison at gene and peak level 
################################################
{
library(Hmisc)

load("peak_m6a_level_per_gene.Rdata")
# gene_m6a <-  list(gene_pSum, gene_pmMax, peak_max)  ## sum of peak  or max peak per gene 
gene_m6a_gn <- rownames(gene_m6a[[1]])

gene.cor <- function(x, y, group,  file)
{
  # x : raw m6a ratio
  # y : mRNA level 
  # group: smaples group information 
  # file:  prefix of output files 
  N_g <- table(group)
  x <- log2(x)            #  log2(IP/Input)
  idx_inf <- is.infinite(x)
  x[idx_inf] <- NA
  y <- log2(y + 1)        #  log2(RPKM +1)
  
  ## per gene level 
  {
    L <- nrow(x)       
    cor_gn <- vector()
    cor_gn_p <- vector()
    cor_gt <- vector()
    cor_gt_p <- vector()
    
    cor_me <- "pearson"
    
    for (i in 1:L) 
    {
      tmp <- rcorr(x[i,  1:N_g[1]], y[i,  1:N_g[1]], type = cor_me)
      cor_gn[i] <- tmp$r[1, 2]
      cor_gn_p[i] <- tmp$P[1, 2]
      tmp <- rcorr(x[i, 11: ncol(x)], y[i, 11 : ncol(y)], type = cor_me)
      cor_gt[i] <- tmp$r[1, 2]
      cor_gt_p[i] <- tmp$P[1, 2]
    }
    
    names(cor_gn) <- rownames(x)
    names(cor_gt) <- rownames(x)
    
    cor_ga <- c(cor_gn, cor_gt)
    cor_gty <- c(rep("Normal", length(cor_gn)), rep("Tumor", length(cor_gt)))
    
    gene_name <- names(cor_ga)
    dat <- data.frame(gene_name, cor_ga, cor_gty)
    #write.csv(dat, paste(file, "_per_gene.csv", sep = ""))
    
    g <- ggplot(dat, aes(x = cor_gty, y = cor_ga, fill = cor_gty, show.legend=F)) + geom_boxplot(outlier.shape = NA) +  theme_bw() + theme(axis.line = element_line(colour = "black"),panel.grid.major = element_blank(),panel.grid.minor = element_blank(),panel.border =  element_blank(),  panel.background = element_blank()) 
    g <- g + scale_fill_manual( values = c("cadetblue3", "coral2"))
    g <- g + geom_hline(yintercept = 0, lty = 2, col = "red") + labs(x = "Sample Group", y = "Cor(peak_m6A, whole_gene_m6A)") + theme(legend.position="none")
    name_2 <- paste(file, "_per_gene.pdf", sep = "")
    ggsave(name_2,  width = 2, height = 3, units = "in")
    
    out <- list("cor.gn" = cor_gn, "cor.gn.pval" = cor_gn_p, "cor.gt" = cor_gt, "cor.gt.pval" = cor_gt_p)
    return(out)
  }   
}

## matched_genes

m6a_pp <- cbind(m6a_rn, m6a_rt)
s_type <- c(rep("Normal", 10), rep("Tumor", 51))
idx_gg <- match(gene_m6a_gn, rownames(m6a_pp))

cor_ntg_pmMax_gene <- gene.cor(gene_m6a[[2]][!is.na(idx_gg), ], m6a_pp[idx_gg[!is.na(idx_gg)], ], s_type, "../Results/cor_gene_pmMax_and_whole_gene_m6A_per_gene_metpeak")
cor_ntg_pSum_gene <- gene.cor(gene_m6a[[1]][!is.na(idx_gg), ], m6a_pp[idx_gg[!is.na(idx_gg)], ], s_type, "../Results/cor_gene_pSum_and_whole_gene_m6A_per_gene_metpeak")

}


################################################
##  gene m6A level distribution and global diff 
################################################
{
    library("Hmisc")                 # for correlaiton with NA
    library("ggplot2")
    library("matrixStats")

    ###############################
    ## mean and cv of m6A per gene 
    {
      
        ## mean m6a level per  genes 
        {
            normal_r_m6a <- log2(rowMeans(m6a_rn, na.rm = T))
            tumor_r_m6a<- log2(rowMeans(m6a_rt, na.rm = T))
            
            ## densidy plot for gene with m6A only
            pdf(file = "../Results/mean_m6a_level_log2_density_mean.pdf", width = 4, height = 4)
            par(mar = c(3, 3, 2, 1), mgp =c(2, 0.7, 0))
            plot(density(tumor_r_m6a, na.rm = T), col = "coral2", xlab = "log2(mean(IP/Input))", main = "Gene with m6A")
            lines(density(normal_r_m6a, na.rm = T), col = "cadetblue3")
            legend(2.5, 0.32, legend = c("Tumor", "Normal"), col = c("red", "blue"), lty = c(1, 1), cex = 0.8)
            abline(v = mean(tumor_r_m6a, na.rm = T), col = "coral2", lty = 2)
            abline(v = mean(normal_r_m6a, na.rm = T), col = "cadetblue", lty = 2)
            dev.off()
            
            library(BSDA)
            z.test(normal_r_m6a, tumor_r_m6a, sigma.x = sd(normal_r_m6a), sigma.y = sd(tumor_r_m6a), alternative = "greater")
            #t.test(normal_r_m6a, tumor_r_m6a, alternative = "greater")
            #wilcox.test(normal_r_m6a, tumor_r_m6a, alternative = "greater")
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
            ggsave("../Results/m6a_level_cv_cdf.pdf",  width = 3, height = 3)
        }
      
        
    }
    
    #########################################
    ## samples m6A levels based correlation 
    {
        library("Hmisc")         
        library("corrplot")
        library("s2dverification")       # for color bar separately
        
        m6a_rm  <- cbind(m6a_rn, m6a_rt)
        cor_m <- rcorr(m6a_rm)
        
        pdf("../Results/corrplot_based_on_tumor_normal_m6A.pdf", width = 6, height = 6)
        corrplot(cor_m$r,  order = "hclust", addrect = 20, method = "square")
        dev.off()
    }
    
   
}


#######################################
##  differential m6A methylated genes 
#######################################
{
    library("Hmisc")         # for rcorr
    library("gplots")
    ##########################################
    ## load and pre-process data and funcitons 
    {
    s_info <- read.table("sample_info.txt", header = T) 
    
    normal_m6a <- m6a_rn
    tumor_m6a <- m6a_rt
    m6a_rs <- cbind(m6a_rn, m6a_rt)
    
    ## m6a level and mRNA level fold change 
    m6a_fd_o <- log2(rowMeans(m6a_rt, na.rm = T) / rowMeans(m6a_rn, na.rm = T))        ## m6a ratio between normal and tumor
    m6a_mean_o <- log2((rowMeans(m6a_rt, na.rm = T) + rowMeans(m6a_rn, na.rm = T))/2)   ## m6a mean of normal and tumor
    
    mrna_fd_o <- log2(rowMeans(input_rpkm_t, na.rm = T) / rowMeans(input_rpkm_n, na.rm = T))
    names(m6a_fd_o) <- rownames(m6a_rt)
    names(mrna_fd_o) <- rownames(m6a_rt)
    #save(m6a_fd_o, mrna_fd_o, file = "log2_m6a_fd_mrna_fd.RData")
    
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
    ##  gene m6A  differential test  
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
        
        if(FALSE){
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
        
        #write.csv(m6a_gd, file = "all_gene_with_diff_m6a_level.csv", quote = F)
        
        idx_u <- fd_d_s > 1
        idx_d <- fd_d_s < 1
        g_u <- rownames(m6a_ds)[idx_u]
        g_d <- rownames(m6a_ds)[idx_d]
        
        #write.table(g_u, file = "all_gene_with_diff_m6a_level_upregulated.txt", quote = F, col.names = F, row.names = F)
        #write.table(g_d, file = "all_gene_with_diff_m6a_level_downregulated.txt", quote = F, col.names = F, row.names = F)
        
        ### diff-m6A gene compared to the m6a peak calling resutls 
        g_t <- read.table("./0_merged_peaks/gene_with_tumor_only_peak.txt", as.is  = T)
        g_n <-  read.table("./0_merged_peaks/gene_with_normal_only_peak.txt", as.is = T)
        g_c <- read.table("./0_merged_peaks/gene_with_common_only_peak.txt", as.is = T)
        
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
        
        #write.table(g_diff_cnt, file = "all_gene_with_diff_m6a_level_cmp2_peakcalling_res.txt", quote = F, col.names = T, row.names = T, sep = "\t" )
        
        idx_uf <- is.na(match(g_u, g_n$V1))            ## rm 21
        g_uf <- g_u[idx_uf]
        idx_df <- is.na(match(g_d, g_t$V1))
        g_df <- g_d[idx_df]                            ## rm 43
        
        g_amb <- c(g_u[!idx_uf], g_d[!idx_df])
        g_amb_t <- c(rep("g_u_with_normal_only_peak", length(g_u[!idx_uf])), rep("g_d_with_tumor_only_peak", length(g_d[!idx_df])))
        g_amb_m <- data.frame(g_amb, g_amb_t) 
        
        ## exclude genes show opposite direction with peak calling results
        write.table(g_uf, file = "all_gene_with_diff_m6a_level_upregulated.txt", quote = F, col.names = F, row.names = F)
        write.table(g_df, file = "all_gene_with_diff_m6a_level_downregulated.txt", quote = F, col.names = F, row.names = F)
        #write.table(g_amb_m, file = "all_gene_with_diff_m6a_level_corrected_amb.txt", quote = F, col.names = T, row.names = F, sep = "\t" )
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

        png(file = "../Results/m6a_level_DE_withcolbar.png", width = 5, height = 6, units = "in", res = 600)
        heatmap.2(m6a_ds_log2, Colv= T, Rowv= F, dendrogram = "column" , scale="row", na.color = "grey",
                  trace ="none", density.info="none", col=hmcols, RowSideColors= gene_c, ColSideColors= sample_c, colsep = 51, rowsep = length(g_df),
                  lhei = c(1,5),  lwid = c(1.3, 4), margins=c(0.5, 0.5), labRow = FALSE, labCol = F,  offsetCol = 0, key.title = NA)
        dev.off()
        
        sum(fd_d_s[idx_ff] > fd1)
        sum(fd_d_s[idx_ff] < fd2)
        
        # all m6A diff protein coding genes with filtering
        pc_gene <- read.table("pc_gene_list.txt")
        idx_pp <- !is.na(match(rownames(m6a_ds_log2), pc_gene$V2))
        png(file = "../Results/m6a_level_DE_withcolbar_protein.png", width = 5, height = 6, units = "in", res = 600)
        heatmap.2(m6a_ds_log2[idx_pp, ], Colv= T, Rowv= F, dendrogram = "column" , scale="row", na.color = "grey",
                  trace ="none", density.info="none", col=hmcols, RowSideColors= gene_c[idx_pp], ColSideColors= sample_c, colsep = 51, rowsep = sum(gene_c[idx_pp] == col_l[8]),
                  lhei = c(1,5),  lwid = c(1.3, 4), margins=c(0.5, 0.5), labRow = FALSE, labCol = F,  offsetCol = 0, key.title = NA)
        dev.off()
        
        # all m6A diff lincRNA genes
        linc_gene <- read.table("lincRNA_gene_list.txt")
        idx_ll <- !is.na(match(rownames(m6a_ds_log2), linc_gene$V2))
        
        hm_lt <- t(m6a_ds_log2[idx_ll, ])
        fd_lt <- colMeans(hm_lt[1:10, ], na.rm = T) - colMeans(hm_lt[11:61, ], na.rm = T)    ## log2 tranformed
        idx_fd_lt <- order(fd_lt, decreasing = T)
        
        png(file = "../Results/m6a_level_DE_withcolbar_lincRNA.png", width = 10, height = 5, units = "in", res = 600)
        heatmap.2(hm_lt[, idx_fd_lt], Colv= F, Rowv= F, dendrogram = "none" , scale="column", na.color = "grey", labRow = F,
                  trace ="none", density.info="none", col=hmcols, ColSideColors= gene_c[idx_ll], RowSideColors= sample_c,
                  margins=c(8, 0.5), offsetCol = 0, key.title = NA, cexCol = 0.75, rowsep =11, colsep = 25)
        
        dev.off()
        
    }
  
  
}

###########################################
##  differential m6A analysis sub sampling
## sample sampling 10 tumor samples 100 times for the :
###########################################
{
  ## without correction
  g_u_ref <- g_u
  g_d_ref <- g_d
  
  hyper_1k_cnt <- hypo_1k_cnt <- vector()                  ## number of genes hyper or hypo methylated 
  hyper_1k_ol_cnt <- hypo_1k_ol_cnt <- vector()            ## number of genes hyper or hypo methylated, overlapping wtih ref
  hyper_1k_ol_frac <- hypo_1k_ol_frac <- vector()            ## fraction of genes hyper or hypo methylated, overlapping wtih ref
  
  k = 1
  while (k <= 100){
    
    print(paste0(k , "th subampling 10 tumors ...."))
    
    set.seed(k)
    idx_st <- sample(1:51, 10)
    idx_st
    
    normal_m6a <- m6a_rn
    tumor_m6a  <- m6a_rt[, idx_st]
    
    ## for differential test: based on the IP/Input directly 
    N <- nrow(tumor_m6a)
    mean_n <- vector()
    mean_t <- vector()
    pval <- vector()
    fd <- vector()
    
    for (i in 1:N) 
    {
      ## difference test 
      tmp <- wilcox.test(tumor_m6a[i, ], normal_m6a[i, ])
      #tmp <- t.test(tumor_m6a[i, ], normal_m6a[i, ])
      pval[i] <- tmp$p.value
      fd[i] <- mean(tumor_m6a[i, ], na.rm = T) / mean(normal_m6a[i, ], na.rm = T)        # fold-change  Tumor/normal !!!! 
      #mean_t[i] <- mean(tumor_m6a[i, ], na.rm = T)
      #$mean_n[i] <- mean(normal_m6a[i, ], na.rm = T)
    }
    pval_a <- p.adjust(pval, method = "BH")         # multiple test adjustment 
    
    
    ## significant diff genes and intersection with ref diff gene  
    m6a_cut = 3/2
    
    idx_u <- fd  >  m6a_cut     &  pval_a < 0.05
    idx_d <- fd  < 1/(m6a_cut)  &  pval_a< 0.05
    idx_u[is.na(idx_u)] <- FALSE
    idx_d[is.na(idx_d)] <- FALSE
    
    g_u <- rownames(normal_m6a)[idx_u]
    g_d <- rownames(normal_m6a)[idx_d]
    
    hyper_1k_cnt[k] <- length(g_u) 
    hypo_1k_cnt[k]  <- length(g_d)
    
    ## overlapping
    hyper_1k_ol_cnt[k] <- length(intersect(g_u, g_u_ref))
    hypo_1k_ol_cnt[k] <- length(intersect(g_d, g_d_ref))
    
    ## fraction 
    hyper_1k_ol_frac[k] <- length(intersect(g_u, g_u_ref)) / length(g_u_ref)
    hypo_1k_ol_frac[k] <- length(intersect(g_d, g_d_ref)) /  length(g_d_ref)
    
    k = k + 1
  }
  
  
  ## visualization
  library(RColorBrewer)
  library(ggplot2)
  
  col_l <- brewer.pal(n = 10, name = "RdBu") 
  
  ol_frac <- c(hyper_1k_ol_frac, hypo_1k_ol_frac)
  ol_group <- rep(c("Hyper", "Hypo"), each = 100)
  
  dat <- data.frame(ol_frac, ol_group)
  
  g <- ggplot(dat, aes(x = ol_group, y = ol_frac, col = ol_group)) + geom_boxplot(outlier.shape = NA) 
  g <- g + geom_point(aes(col = ol_group), position = position_jitter(width = 0.3)) 
  g <- g + labs(y = "Overlapped farction against original test", x = "Subsampling 10 tumors 100 times") + theme_classic() 
  g <- g + scale_color_manual( values = c(col_l[2], col_l[8])) + theme(legend.position = "none")
  ggsave("../Results/Overlapped_farction_of_Hyper_Hypo_genes_by_subsampling.pdf", width = 4, height = 4, units = "in")
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
    
    # boxplot
    pdf(file = "../Results/Mean_m6A_level_in_Normal_vs_Tumor_boxplot_flipped.pdf",  width = 3, height = 3)
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
    pdf(file = "../Results/Mean_mRNA_in_Normal_vs_Tumor_boxplot_flipped.pdf",  width = 3, height = 3)
    par(mar = c(3, 3, 0.5, 0.5), mgp = c(2, 0.75, 0))
    boxplot(list( mrna_me_n[idx_d], mrna_me_t[idx_d], mrna_me_n[idx_u], mrna_me_t[idx_u]), col = rep(c("cadetblue3", "coral2"), 2), ylab = "Mean mRNA level")
    abline( v = 2.5, lty = 2)
    dev.off()
    
    wilcox.test(mrna_me_n[idx_u], mrna_me_t[idx_u])
    wilcox.test(mrna_me_n[idx_d], mrna_me_t[idx_d])
    wilcox.test(mrna_me_n[idx_u], mrna_me_n[idx_d], "greater")
    wilcox.test(mrna_me_t[idx_u], mrna_me_t[idx_d], "greater")
    
}

############################################
### m6A level FC vs mRNA FC with DEG p value 
{
  ###########################
  ## based own data
  {
        de_res <- read.csv("61_sample_Input_DESeq_DE.csv")
        
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
        dm_gene_nodiff_rna <- (de_res[idx_a, 1])[idx1]
        
        length(dm_gene_nodiff_rna)
        
        sum(idx_gu & p_ff > -log10(0.05) & fc_ff < -1)
        sum(idx_gd & p_ff > -log10(0.05) & fc_ff < -1)
        sum(idx_gu & p_ff > -log10(0.05) & fc_ff > 1)
        sum(idx_gd & p_ff > -log10(0.05) & fc_ff > 1)
        
        col_gff <- c(rep(col_l[2], length(g_uf)), rep(col_l[8], length(g_df)))
        col_gff[idx1&idx_gu] <- col_l[5]
        col_gff[idx1&idx_gd] <- col_l[6]
        
        pdf(file = "../Results/DMG_mRNA_vocanal_plot.pdf", width = 3, height = 3)
        par(mar = c(3, 3, 0.5, 0.5), mgp =c(2, 0.75, 0))
        plot(fc_ff, p_ff, col =  col_gff , xlab = "log2(mRNA_FC)", ylab = "-log10(FDR)", xlim = c(-7, 7), pch = 16, cex = 0.75)
        abline(h = -log10(0.05), lty = 2, col = "gray")
        abline(v = c(-1, 1), lty = 2, col = "gray")
        dev.off()
        
    }
    
    ## m6A levels in DEGs
  {
    
    ###mean m6A level per gene
    m6a_mt <- rowMeans(m6a_rt, na.rm = T )
    all_g <- c(g_uf, g_df) 
    idx_a <- match(all_g, names(m6a_mt))
    m6a_mt<- m6a_mt[idx_a]     ##  in tumor only 
    
    ## dm genes with nodiff_RNA
    idx_g <- match(all_g, dm_gene_nodiff_rna)
    sig <- m6a_mt[!is.na(idx_g)]
    non_sig <- m6a_mt[is.na(idx_g)]
    
    #tmp <- t.test(sig, non_sig)
    tmp <- wilcox.test(sig, non_sig)
    m_fc <- mean(sig, na.rm = T)/mean(non_sig, na.rm = T)
    tt <- paste("sig/non_sig=", formatC(m_fc, format = "e", digits = 2), "; p=", formatC(tmp$p.value, format = "e", digits = 2), sep = "");
    tt 
    
    m6a_m <- c(log2(non_sig), log2(sig))
    m6a_tag <- c(rep("non_sig", length(non_sig)), rep("sig", length(sig)))
    dat <- data.frame(m6a_m, m6a_tag)
    
    g <- ggplot(dat, aes(x = m6a_tag, y = m6a_m, fill = m6a_tag)) + geom_boxplot()
    g <- g + scale_fill_manual(values = c("gray", "brown4")) + ylim(-4, 4)
    g <- g + theme_classic() + theme(legend.position = "none") 
    ggsave("../Results/DMGs_mRNA_level_non-diff_vs_diff_m6A_levels_in_tumor.pdf", width = 2.7, height = 3)
    
  }
}
 
}


########################################
## EML4 related analysese
########################################
{
  gene_s <- "EML4"
  
  ##################
  ## m6A FC vs Pval 
  {
  library(ggrepel)
  ## combine m6a and gene diff res 
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
  de <- cbind(m6a_de_m, mrna_de_out)
  
  ## non-diff RNA
  idx_no <- de[, 10] >= 0.05 | abs(de[, 8]) <= 1    ## 447
  idx_k <- idx_no
  
  dat <- data.frame(de[idx_k, ])
  
  ## add gene color 
  library(RColorBrewer)
  col_l <- brewer.pal(n = 10, name = "RdBu")      
  idx_hyper <- dat$log2.m6a_level_tumor.normla. > 0
  
  gene_dm <- rep("Hypo", nrow(dat))
  gene_dm[idx_hyper] <- "Hyper"
  dat <- cbind(dat, gene_dm)
  
  
  idx_g <- match("EML4", rownames(dat))
  g_lab <- rep("", nrow(dat))
  g_lab[idx_g] <- "EML4"
  
  ## without RNA level
  g <- ggplot(dat, aes(x = log2.m6a_level_tumor.normla., y = -log10(p.value), alpha = 0.9,
                       col = gene_dm,  label = g_lab)) 
  g <- g  + geom_point()   + geom_text_repel(color = "red") 
  g <- g + scale_color_manual(values = c(col_l[2], col_l[8]))
  g <- g + labs(x = "log2FC_m6A", y = "-log10(pval_m6A)")  + theme_classic() 
  ggsave("../Results/EML4_m6A_log2FC_vs_pvals_all_447_DMG_nodiff_RNA_nGene_labels.pdf", width = 3.75, height = 3)
  }
  
  
  ###########################
  ## EML4 peak m6a level
  {
    load("EML4_peaks_m6a_level.Rdata")
    
    idx_sn <- match(colnames(m6a_rn), colnames(peak_m6a))
    idx_st <- match(colnames(m6a_rt), colnames(peak_m6a))
    
    ## using EML4_peak_region_last_exon
    peak_m6a_n <- peak_m6a[2, idx_sn]
    peak_m6a_t <- peak_m6a[2, idx_st]
    
    peak_s <- c(log2(peak_m6a_n),  log2(peak_m6a_t))
    group_s <- c(rep("Normal", 10), rep("Tumor", 51))
    
    
    wilcox.test(peak_m6a_n, peak_m6a_t, alternative = "less")
    
    dat <- data.frame(peak_s, group_s)

    g <- ggplot(dat, aes(x = group_s, y = peak_s, col = group_s)) + geom_boxplot(outlier.shape = NA) 
    g <- g + geom_point(aes(col = group_s), position = position_jitter(width = 0.3)) 
    #g <- g + geom_text_repel(label = names(peak_s))
    g <- g + theme_bw() + theme(axis.line = element_line(colour = "black"),panel.grid.major = element_blank(),panel.grid.minor = element_blank(),panel.border =  element_blank(),  panel.background = element_blank()) 
    g <- g + scale_color_manual( values = c("cadetblue3", "coral2")) + theme(legend.position = "none")
    ggsave(paste("../Results/", gene_s, "_peak_m6A_level_boxplot_with_samples.pdf", sep = ""), width = 2.5, height = 3)

    
  }
  
  
  ###########################
  ## EML4 mRNA levels in TCGA
  {
    #load("/Users/yong/OneDrive - UHN/Projects/Epitranscriptomics/LUAD_m6A/1_mapping_quant/5_RNA-seq/TCGA/GDC_LUAD/tcga_htseq_fpkm_log2.Rdata")
    #idx_g <- match(gene_s, rownames(tcga))
    #tcga <- tcga[idx_g, ]
    
    # save(tcga, tcga_type, file = "TCGA_LUAD_EML4_expr.Rdata")
    load("TCGA_LUAD_EML4_expr.Rdata")
    
    idx_n <- tcga_type == "normal"
    idx_t <- tcga_type == "tumor"
    
    #tmp <- t.test(tcga[idx_g, idx_n], tcga[idx_g, idx_t])
    
    tmp <- wilcox.test(tcga[idx_n], tcga[idx_t], alternative = "less")
    
    m_fc <- mean(2^(tcga[idx_t] - 1)) / mean(2^(tcga[idx_n] - 1))
    tt <- paste("T/N_FC=", formatC(m_fc, format = "e", digits = 2), "; p=", formatC(tmp$p.value, format = "e", digits = 2), sep = "")
    
    name=paste("../Results/", gene_s, "_mRNA_level_boxplot_TCGA.pdf", sep = "")
    #png(file = name, width = 3, height = 3, units = "in", res = 300)
    pdf(file = name, width = 2.5, height = 3)
    par(mar = c(3, 3, 2, 0.5), mgp =c(2, 0.75, 0), bty = "l", las = 1)
    boxplot(list(tcga[idx_n], tcga[idx_t]),  main = tt, cex.main = 0.75, ylab = "TCGA mRNA level", xlab = gene_s, col =  c("cadetblue3", "coral2"))
    axis(1, at = seq(1, 2, 1), labels = c("Normal", "Tumor") )     # for x axis
    dev.off()
  
}
  
  ##########################################################
  ## EML4 with matched rna and protein in CPTAC LUAD corhort
  {
    
    ## load RNA and protein data for matched NAT and Tumor
    ## log2 transformed with pseudo count 
    rna_nm <- read.csv("CPTAC_LUAD_RNA_expr_4matched_NAT.csv")
    rna_tm <- read.csv("CPTAC_LUAD_RNA_expr_4matched_Tumor.csv")
    gene_lm <- rna_nm$X
    rna_nm <- as.matrix(rna_nm[, -1])
    rna_tm <- as.matrix(rna_tm[, -1])
    
    ## protein 
    protein_nm <- read.csv("CPTAC_LUAD_Protein_expr_4matched_NAT.csv")
    protein_tm <- read.csv("CPTAC_LUAD_Protein_expr_4matched_Tumor.csv")
    gene_lp <- protein_nm$X
    protein_nm <- as.matrix(protein_nm[, -1])
    protein_tm <- as.matrix(protein_tm[, -1])
    
    ###################################
    ## RNA and protein levels for EML4
    ## paired test
    ##################################
    {
      gene_s = "EML4"
      
      ## RNA
      {
        idx_g <- match(gene_s, gene_lm)
        
        tmp <- wilcox.test(rna_nm[idx_g, ], rna_tm[idx_g, ], paired = TRUE)
        #m_fc <- mean(2^gene_e$EML4[idx_t] ) / mean(2^gene_e$EML4[idx_n] - 1)
        m_fc <- mean((2^rna_tm[idx_g, ] ) / (2^rna_nm[idx_g, ]))
        tt <- paste("T/N_FCm=", formatC(m_fc, format = "e", digits = 2), "; wc_p=", formatC(tmp$p.value, format = "e", digits = 2), sep = "")
        
        ## data frame
        type <-  as.factor(rep(c("Normal","Tumor"), each = ncol(rna_nm)))
        group <- rep(c(1,3), each = ncol(rna_nm))
        score <- c(rna_nm[idx_g, ], rna_tm[idx_g, ])
        dat <- data.frame(type, group, score)
        
        g <- ggplot(dat, aes(x = group, y = score)) 
        g <- g + geom_point(aes(col = type))  + xlim(0, 3.5) + ylim(0, 7)
        g <- g + scale_color_manual(values = c("blue","red")) + theme_classic() 
        g <- g + labs(title = tt, y = "RNA levels")
        
        ## segment   
        L <- length(rna_nm[idx_g, ])
        for(i in 1:L)
        {
          g <- g + geom_segment(x = 1 + 0.2, y = rna_nm[idx_g, i], xend = 3 - 0.2, yend = rna_tm[idx_g, i], 
                                size = 0.25, color = "gray90", alpha = 0.7)
        }
        ggsave("../Results/EML4_mRNA_level_boxplot_in_CPTAC_LUAD_paired.pdf", width = 3.5, height = 3)
      }
      
      ## protein
      {
        idx_g <- match(gene_s, gene_lp)
        
        tmp <- wilcox.test(protein_nm[idx_g, ], protein_tm[idx_g, ], paired = TRUE)
        m_fc <- mean((2^protein_tm[idx_g, ]) / (2^protein_nm[idx_g, ]))
        tt <- paste("T/N_FCm=", formatC(m_fc, format = "e", digits = 2), "; wc_p=", formatC(tmp$p.value, format = "e", digits = 2), sep = "")
        
        ## data frame
        type <-  as.factor(rep(c("Normal","Tumor"), each = ncol(protein_nm)))
        group <- rep(c(1,3), each = ncol(protein_nm))
        score <- c(protein_nm[idx_g, ], protein_tm[idx_g, ])
        dat <- data.frame(type, group, score)
        
        g <- ggplot(dat, aes(x = group, y = score)) 
        g <- g + geom_point(aes(col = type))  + xlim(0, 3.5) + ylim(-2, 3)
        g <- g + scale_color_manual(values = c("blue","red")) + theme_classic() 
        g <- g + labs(title = tt, y = "Protein levels")
        
        ## segment   
        L <- length(protein_nm[idx_g, ])
        for(i in 1:L)
        {
          g <- g + geom_segment(x = 1 + 0.2, y = protein_nm[idx_g, i], xend = 3 - 0.2, yend = protein_tm[idx_g, i], 
                                size = 0.25, color = "gray90", alpha = 0.7)
        }
        ggsave("../Results/EML4_Protein_level_boxplot_in_CPTAC_LUAD_paired.pdf", width = 3.5, height = 3)
      }
    }
    
     
  }
  
}
