#######################################################
##   1_cor_m6a_mRNA.R
##   1: at sample level
##   2: at gene level
#######################################################

rm(list = ls())
setwd("/Users/Yong/Yong/m6a_profiling/5_m6a_mrna/1_cor_m6a_mRNA_per_sample/")
library("Hmisc")

# load m6a level and match mrna level

load("/Users/Yong/Yong/m6a_profiling/3_m6a_level/m6a_level_keep_dup/gene_with_peak/m6A_level_4gene_passed_filtering.Rdata")  
pc_gene <- read.table("/Users/Yong/Yong/m6a_profiling/1_mapping_quant/5_RNA-seq/pc_gene_list.txt")


##################################################################################
#### correlation between the m6A level (IP/Input) correlation and mRNA  ##########

###############
## sample level
{
 # log2 transformed
 # function for drawing the boxplot cor(m6a_level, mRNA) at sample and gene level
   
     sample.cor <- function(x, y, group,  file)
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
        
        ## per sample level   
        {    
            N <- ncol(x)       # sample size
            cor_a <- vector()
            cor_p <- vector()
            cor_me <- "pearson"
            
            for (i in 1:N) 
            {
                tmp <- rcorr(x[ , i], y[, i], type = cor_me)
                cor_a[i] <- tmp$r[1, 2]
                cor_p[i] <- tmp$P[1, 2]
            }
            
            names(cor_a) <- colnames(x)
            out <- list("r" = cor_a, "group" = group, "p.value" = cor_p)
           
            
            dat <- data.frame(cor_a, group)
            g <- ggplot(dat, aes(x = group, y = cor_a)) + geom_boxplot(outlier.shape=NA)  +geom_point(aes( col = group), position = position_jitter(width = 0.3)) +  theme_bw() + theme(axis.line = element_line(colour = "black"),panel.grid.major = element_blank(),panel.grid.minor = element_blank(),panel.border =  element_blank(),  panel.background = element_blank()) 
            g <- g + scale_color_manual( values = c("cadetblue3", "coral2")) + labs(x = "Sample type", y = "Cor(m6A, mRNA)") + theme(legend.position="none")
            name_1 <- paste(file, "_per_sample.pdf", sep = "")
            ggsave(name_1,  width = 2, height = 3)
            
            
            ## examples for max and min
            
            source("/Users/Yong/Yong/R/functions/cor.plot.R")
            ii_n <- which(cor_a == min(cor_a[1:N_g[1]]))
            ii_t <- which(cor_a == min(cor_a[(N_g[2] + 1) : sum(N_g)]))
            cor.plot(y[, ii_n], x[, ii_n], "Normal Example (per sample)_min", "mRNA level", "log2(IP/Input)", "PNG")
            cor.plot(y[, ii_t], x[, ii_t], "Tumor Example (per sample)_min", "mRNA level", "log2(IP/Input)", "PNG")
            
            # for max
            ii_n <- which(cor_a == max(cor_a[1:N_g[1]]))
            ii_t <- which(cor_a == max(cor_a[(N_g[2] + 1) : sum(N_g)]))
            cor.plot(y[, ii_n], x[, ii_n], "Normal Example (per sample)_max", "mRNA level", "log2(IP/Input)", "PNG")
            cor.plot(y[, ii_t], x[, ii_t], "Tumor Example (per sample)_max", "mRNA level", "log2(IP/Input)", "PNG")
        }
        
    return(out)
    }
 
# cor(m6a, mrna)
    s_type <- c(rep("Normal", 10), rep("Tumor", 51))
    
    # all gene
    m6a_rs <- cbind(m6a_rn, m6a_rt)                     # merger normal and tumor
    mrna_s <- cbind(input_rpkm_n, input_rpkm_t)
    
    # protein coding gene only
    idx_p <- !is.na(match(rownames(mrna_s), pc_gene$V2))
    mrna_sp <- mrna_s[idx_p, ]
    m6a_rsp <- m6a_rs[idx_p, ]
      
    cor_s <- sample.cor(m6a_rs, mrna_s, s_type, "cor_m6A_level_and_mRNA_RPKM_all_gene")
    
 #   cor_s_pc <- sample.cor(m6a_rsp, mrna_sp, s_type, "cor_m6A_level_and_mRNA_RPKM_PC_gene")
}

## cor test 
shapiro.test(cor_s[[1]][1:10])
shapiro.test(cor_s[[1]][11:61])

r <- cor_s$r
group <- cor_s$group
dat <- data.frame(r, group)
leveneTest(r, group, data = dat)

#t.test(cor_s[[1]][1:10], cor_s[[1]][11:61], alternative = "less")
wilcox.test(cor_s[[1]][1:10], cor_s[[1]][11:61], alternative = "less")



if(FALSE){
###############################################
##############################################
global.cor <- function(x, y, group,  file)
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
    
    
    source("/Users/Yong/Yong/R/functions/cor.plot.R")
    
    ## per sample level   
    {    
        N <- ncol(x)       # sample size
        cor_a <- vector()
        cor_me <- "pearson"
        
        for (i in 1:N) 
        {
            tmp <- rcorr(x[ , i], y[, i], type = cor_me)
            cor_a[i] <- tmp$r[1, 2]
        }
        
        dat <- data.frame(cor_a, group)
        g <- ggplot(dat, aes(x = group, y = cor_a)) + geom_boxplot(outlier.shape=NA)  +geom_point(aes( col = group), position = position_jitter(width = 0.3)) +  theme_bw() + theme(axis.line = element_line(colour = "black"),panel.grid.major = element_blank(),panel.grid.minor = element_blank(),panel.border =  element_blank(),  panel.background = element_blank()) 
        g <- g + scale_color_manual( values = c("cadetblue3", "coral2")) + labs(x = "Sample Group", y = "Cor(m6A, mRNA)") + theme(legend.position="none")
        name_1 <- paste(file, "_per_sample.png", sep = "")
        ggsave(name_1,  width = 3, height = 3, units = "in", dpi = 300)
        
        ## examples
        ii_n <- which(cor_a == min(cor_a[1:N_g[1]]))
        ii_t <- which(cor_a == min(cor_a[(N_g[2] + 1) : sum(N_g)]))
        cor.plot(y[, ii_n], x[, ii_n], "Normal Example (per sample)", "mRNA level", "log2(IP/Input)")
        cor.plot(y[, ii_t], x[, ii_t], "Tumor Example (per sample)", "mRNA level", "log2(IP/Input)")
    }
    
    ## per gene level 
    {
        L <- nrow(x)       
        cor_gn <- vector()
        cor_gt <- vector()
        
        for (i in 1:L) 
        {
            tmp <- rcorr(x[i,  1:N_g[1]], y[i,  1:N_g[1]], type = cor_me)
            cor_gn[i] <- tmp$r[1, 2]
            tmp <- rcorr(x[i, 11: 61], y[i, 11 : 61], type = cor_me)
            cor_gt[i] <- tmp$r[1, 2]
        }
        
        names(cor_gn) <- rownames(x)
        names(cor_gt) <- rownames(x)
        
        save(cor_gn, cor_gt, file = paste(file, ".Rdata", sep = ""))     # save correaltion between m6A level and mRNA
        
        cor_ga <- c(cor_gn, cor_gt)
        cor_gty <- c(rep("Normal", length(cor_gn)), rep("Tumor", length(cor_gt)))
        
        dat <- data.frame(cor_ga, cor_gty)
        g <- ggplot(dat, aes(x = cor_gty, y = cor_ga, fill = cor_gty, show.legend=F)) + geom_boxplot(outlier.shape = NA) +  theme_bw() + theme(axis.line = element_line(colour = "black"),panel.grid.major = element_blank(),panel.grid.minor = element_blank(),panel.border =  element_blank(),  panel.background = element_blank()) 
        g <- g + scale_fill_manual( values = c("cadetblue3", "coral2"))
        g <- g + geom_hline(yintercept = 0, lty = 2, col = "red") + labs(x = "Sample Group", y = "Cor(m6A, mRNA)") + theme(legend.position="none")
        name_2 <- paste(file, "_per_gene.png", sep = "")
        ggsave(name_2,  width = 3, height = 3, units = "in", dpi = 300)
        
        ## examples
        if(FALSE){
            cor_tmp <- cor_gt - cor_gn
            ii_ma <- which(cor_tmp == max(cor_tmp))
            ii_mi <- which(cor_tmp == min(cor_tmp))
            
            eg_nn <- paste("Normal_", rownames(y)[ii_ma], sep = "")
            eg_nt <- paste("Tumor_", rownames(y)[ii_ma], sep = "")
            cor.plot(y[ii_ma, 1:10], x[ii_ma, 1:10], eg_nn, "mRNA level", "log2(IP/Input)")
            cor.plot(y[ii_ma, 11:61], x[ii_ma, 11:61], eg_nt, "mRNA level", "log2(IP/Input)")
            
            eg_nn <- paste("Normal_", rownames(y)[ii_mi], sep = "")
            eg_nt <- paste("Tumor_", rownames(y)[ii_mi], sep = "")
            cor.plot(y[ii_mi, 1:10], x[ii_mi, 1:10], eg_nn, "mRNA level", "log2(IP/Input)")
            cor.plot(y[ii_mi, 11:61], x[ii_mi, 11:61], eg_nt, "mRNA level", "log2(IP/Input)")
        }
    }   
    
}

##############################################
## comparsion of cor between tumor and normal
{
## out put corrlation table
{
load("cor_m6A_level_and_mRNA_RPKM_all_gene.Rdata")

m6a_nm <- log2(rowMeans(m6a_rs[, 1:10], na.rm = T))
mrna_nm <- log2(rowMeans(mrna_s[, 1:10]) + 1)
m6a_tm <- log2(rowMeans(m6a_rs[, 11:61], na.rm = T))
mrna_tm <- log2(rowMeans(mrna_s[, 11:61]) + 1)

cor_out <- cbind(rownames(m6a_rs), mrna_nm, mrna_tm, m6a_nm, m6a_tm, cor_gn, cor_gt)
colnames(cor_out) <- c("Gene_name", "log2(mRNA_mean)_normal", "log2(mRNA_mean)_tumor", "log2(m6A_mean)_normal", "log2(m6A_mean)_tumor", "Cor(m6A, mRNA)_normal", "Cor(m6A, mRNA)_tumor")
write.table(cor_out, file = "Cor_m6A_mRNA_per_gene_for_normal_tumor.tsv", col.names = T, row.names = F, quote = F, sep = "\t")
}

## scater plot
{
L <- length(cor_gn)
cor_c <- 0.3                # corelation cutoff
col_r <- rep("gray80", L)
idx1 <- cor_gn < -cor_c & cor_gt < -cor_c
idx2 <- cor_gn > cor_c & cor_gt > cor_c
idx3 <- cor_gn < -cor_c & cor_gt > cor_c
idx4 <- cor_gn > cor_c & cor_gt < -cor_c

col_r[idx1] <- "blue"
col_r[idx2] <- "red"
col_r[idx3] <- "orange"
col_r[idx4] <- "purple"

png(file = "cor_m6a_mrna_normal_vs_tumor.png", width = 4, height = 4, units = "in", res = 300)
par(mar = c(3, 3, 0.5, 0.5), mgp =c(2, 0.7, 0))
plot(cor_gn, cor_gt, col = col_r, xlab = "cor(m6A, mRNA)_Normal", ylab = "cor(m6A, mRNA)_Tumor", xlim = c(-0.8, 0.8), ylim = c(-1, 1))
abline(h = c(cor_c, -cor_c), v= c(cor_c, -cor_c), lty = 2, lwd = 0.75)
text(x = c(-0.8, 0.8, -0.8, 0.8), y = c(-0.7, 0.7, 0.7, -0.7), labels = c(sum(idx1), sum(idx2), sum(idx3), sum(idx4)))
dev.off()

g1 <- rownames(m6a_rs)[idx1]
g2 <- rownames(m6a_rs)[idx2]
g3 <- rownames(m6a_rs)[idx3]
g4 <- rownames(m6a_rs)[idx4]
ga <- rownames(m6a_rs)

write.table(g1, file = "cor_both_neg_normal_vs_tumor_genes.txt", row.names = F, col.names = F, quote = F)
write.table(g2, file = "cor_both_pos_normal_vs_tumor_genes.txt", row.names = F, col.names = F, quote = F)
write.table(g3, file = "cor_normal_neg_tumor_pos_genes.txt", row.names = F, col.names = F, quote = F)
write.table(g4, file = "cor_normal_pos_tumor_neg_genes.txt", row.names = F, col.names = F, quote = F)
write.table(ga, file = "all_gene_with_m6a_after_filtering.txt", row.names = F, col.names = F, quote = F)
}

}
}