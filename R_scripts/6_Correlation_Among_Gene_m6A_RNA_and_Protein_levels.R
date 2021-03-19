################################################################
## This script implements:
## 1. correlation between the m6A level and RNA at sample level
## 2. correlation between the m6A level and RNA at gene level
## 3. compare correlation coefficent between normal and tumor per gene 
## 4. to be added ...
#################################################################

rm(list = ls())
setwd("/work/dir/")


############################################################
## correlation between the m6A level and RNA at sample level
############################################################
{
    load("/Users/Yong/Yong/m6a_profiling/3_m6a_level/m6a_level_keep_dup/gene_with_peak/m6A_level_4gene_passed_filtering.Rdata")  
    
    ################################################################################
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
    
    ##########################
    ## correlation comparison
    {
    s_type <- c(rep("Normal", 10), rep("Tumor", 51))
    
    # all gene
    m6a_rs <- cbind(m6a_rn, m6a_rt)                     # merger normal and tumor
    mrna_s <- cbind(input_rpkm_n, input_rpkm_t)
    
    # protein coding gene only
    pc_gene <- read.table("/Users/Yong/Yong/m6a_profiling/1_mapping_quant/5_RNA-seq/pc_gene_list.txt")
    idx_p <- !is.na(match(rownames(mrna_s), pc_gene$V2))
    mrna_sp <- mrna_s[idx_p, ]
    m6a_rsp <- m6a_rs[idx_p, ]
    
    cor_s <- sample.cor(m6a_rs, mrna_s, s_type, "cor_m6A_level_and_mRNA_RPKM_all_gene")
    
    ## cor test 
    shapiro.test(cor_s[[1]][1:10])
    shapiro.test(cor_s[[1]][11:61])
    
    r <- cor_s$r
    group <- cor_s$group
    dat <- data.frame(r, group)
    leveneTest(r, group, data = dat)
    
    #t.test(cor_s[[1]][1:10], cor_s[[1]][11:61], alternative = "less")
    wilcox.test(cor_s[[1]][1:10], cor_s[[1]][11:61], alternative = "less")
    }
    
    
    ######################################################
    ## for tumor only and subsampling to 10 for 1000 times
    {
    cor_me <- "pearson" 
    ip_count <- cbind(ip_count_n, ip_count_t)
    input_count <- cbind(input_count_n, input_count_t)
    
    m6a_log2 <- log2(ip_count / input_count)
    m6a_log2[is.infinite(m6a_log2 )] <- NA
    
    input_rpkm <- cbind(input_rpkm_n, input_count_t)
    input_rpkm_log2 <- log2(input_rpkm + 1)
    
    ### permuated IP 1000 for cor(input, ip.input)
    {
        N <- nrow(ip_count)
        M <- ncol(ip_count)
        S <- 1000
        
        ## cor(inpt, ip/input) by permutated gene IP 
        cor_mat_r <- matrix(0, M, S)
        cor_mat_p <- matrix(0, M, S)
        
        for (i in 1:S)
        {
            idx_s <- sample(N)
            m6a_s <- ip_count[idx_s, ] / input_count
            m6a_s_log2 <- log2(m6a_s)
            m6a_s_log2[is.infinite(m6a_s_log2 )] <- NA
            
            for(j in 1:M)
            {
                tmp <- rcorr(input_rpkm_log2[, j],  m6a_s_log2[, j], type = cor_me)
                cor_mat_r[j, i] <- tmp$r[1, 2]
                cor_mat_p[j, i] <- tmp$P[1, 2]
            }
            
        }
        
        
        cor_mat_r[is.na(cor_mat_r)] <- 0       ## assign no cor for NA
        cor_mat_p[is.na(cor_mat_p)] <- 1     
        
        save(cor_mat_r, cor_mat_p, file = "cor_rna_m6a_by_permutated_IP_1000_times_per_sample.Rdata")
        
        pdf("cor_rna_m6a_egNL_by_permutated_IP_1000_times_per_sample.pdf", width = 4, height =  3)
        par(mar = c(4, 3, 2, 0.5), mgp = c(2, 0.75, 0))
        plot(density(cor_mat_r), xlab = "cor(IP, IP/INPUT)",
             main = "cor_rna_m6a_by_permutated_IP_1000_times_per_sample", cex.main = 0.75, cex.lab = 0.75)
        abline(h = 0, lty = 2, col = "gray")
        dev.off()
    }
    }
}


############################################################
## correlation between the m6A level and RNA at gene level
############################################################
{
    library("Hmisc")         # for rcorr
    library("VennDiagram")
    
    ######################################
    ## load m6a level and match mrna level
    {
    load("/Users/Yong/Yong/m6a_profiling/3_m6a_level/m6a_level_keep_dup/gene_with_peak/m6A_level_4gene_passed_filtering.Rdata")  
    pc_gene <- read.table("/Users/Yong/Yong/m6a_profiling/1_mapping_quant/5_RNA-seq/pc_gene_list.txt")
    
    m6a_rt_log2 <- log2(m6a_rt)
    m6a_rn_log2 <- log2(m6a_rn)
    input_rpkm_t_log2 <- log2(input_rpkm_t + 1)
    input_rpkm_n_log2 <- log2(input_rpkm_n + 1)
    
    ## merge load 
    # cor(m6a, mrna) per sample
    s_type <- c(rep("Normal", 10), rep("Tumor", 51))
    m6a_rs <- cbind(m6a_rn, m6a_rt)                     # merger normal and tumor
    mrna_s <- cbind(input_rpkm_n, input_rpkm_t)
    }
    
    ###################################
    ## subsampling to same sample size
    {
        ## subsampling tumor sample size 
        {
            ## subsampling tumor samples to 10 for 100 times
            sub_t <- 100               # subsampling times
            sub_s <- seq(5, 50, 5)     # subsampling sample size
            L <- length(sub_s)
            cor_sub <- matrix(0, sub_t, L)
            
            ## for single gene 
            {
                #for single gene 
                k = 1;
                for (i in sub_s)
                {
                    for (j in 1:sub_t)
                    {
                        
                        g_id <- 200                         # random gene  
                        
                        # mRNA
                        idx_s <- sample(1:51, i)
                        en_s <- input_rpkm_t_log2[g_id, idx_s]                  
                        
                        # for 1 gene subsampling sub_t tiems
                        tmp = cor.test(m6a_rt_log2[g_id, idx_s], en_s)      # test gene 
                        cor_sub[j, k] = tmp$estimate
                    }
                    k = k + 1;
                    
                }
                
                png(file = "Cor_m6a_MRNA_single_gene_subsampl100_withdiff_sample_size.png", width = 5, height = 4, units = "in", res = 300)
                boxplot(cor_sub, xaxt = "n", xlab = "Sample Size", ylab = "Pearson correlation", main = "Single gene subsampling 100 times")
                axis(1, at = 1:L, labels = sub_s)     # for x axis
                dev.off()
            }
            
            ## for all genes sample down sample size once 
            {
                #for all genes
                
                k = 1;
                sub_s <- seq(5, 50, 5)
                L_n <- nrow(m6a_rt_log2)
                cor_sub_a <- matrix(0, L_n, L)
                for (i in sub_s)
                {
                    
                    for (j in 1:L_n)
                        
                    {
                        # mRNA
                        idx_s <- sample(1:51, i)
                        en_s <- input_rpkm_t_log2[j, idx_s]
                        
                        # for 1 gene subsampling sub_t tiems
                        if (sum(!is.na(m6a_rt_log2[j, idx_s])) > 2)
                        {
                            tmp = cor.test(m6a_rt_log2[j, idx_s], en_s)      # test gene 
                            cor_sub_a[j, k] = tmp$estimate
                        } else
                        {cor_sub_a[j, k] = NA }
                    }
                    k = k + 1;
                    
                }
                
                png(file = "Cor_m6a_MRNA_all_gene_subsample_1_withdiff_sample_size.png", width = 5, height = 4, units = "in", res = 300)
                boxplot(cor_sub_a, xaxt = "n", xlab = "Sample Size", ylab = "Pearson correlation", main = "all gene subsampling once")
                axis(1, at = 1:L, labels = sub_s)     # for x axis
                dev.off()
            }
            
        }
        
        ## subsampling normal sample size 
        {
            
            ## subsampling tumor samples to 10 for 100 times
            sub_t <- 100               # subsampling times
            sub_s <- seq(3, 9, 3)     # subsampling sample size
            L <- length(sub_s)
            cor_sub <- matrix(0, sub_t, L)
            
            ## for single gene 
            {
                #for single gene 
                k = 1;
                for (i in sub_s)
                {
                    
                    for (j in 1:sub_t)
                        
                    {
                        
                        g_id <- 200                         # random gene  
                        
                        # mRNA
                        idx_s <- sample(1:10, i)
                        en_s <- input_rpkm_n_log2[g_id, idx_s]                  
                        
                        # for 1 gene subsampling sub_t tiems
                        tmp = cor.test(m6a_rn_log2[g_id, idx_s], en_s)      # test gene 
                        cor_sub[j, k] = tmp$estimate
                    }
                    k = k + 1;
                    
                }
                
                png(file = "Cor_m6a_MRNA_single_gene_subsampl100_withdiff_sample_size_normal.png", width = 5, height = 4, units = "in", res = 300)
                boxplot(cor_sub, xaxt = "n", xlab = "Sample Size", ylab = "Pearson correlation", main = "Single gene subsampling 100 times")
                axis(1, at = 1:L, labels = sub_s)     # for x axis
                dev.off()
            }
            
            ## for all genes sample down sample size once 
            {
                #for all genes
                
                k = 1;
                sub_s <- seq(3, 9, 3)
                L_n <- nrow(m6a_rn_log2)
                cor_sub_a <- matrix(0, L_n, L)
                for (i in sub_s)
                {
                    
                    for (j in 1:L_n)
                        
                    {
                        # mRNA
                        idx_s <- sample(1:10, i)
                        en_s <- input_rpkm_n_log2[j, idx_s]
                        
                        # for 1 gene subsampling sub_t tiems
                        if (sum(!is.na(m6a_rn_log2[j, idx_s])) > 2)
                        {
                            tmp = cor.test(m6a_rn_log2[j, idx_s], en_s)      # test gene 
                            cor_sub_a[j, k] = tmp$estimate
                        } else
                        {cor_sub_a[j, k] = NA }
                    }
                    k = k + 1;
                    
                }
                
                png(file = "Cor_m6a_MRNA_all_gene_subsample_1_withdiff_sample_size_normal.png", width = 5, height = 4, units = "in", res = 300)
                boxplot(cor_sub_a, xaxt = "n", xlab = "Sample Size", ylab = "Pearson correlation", main = "all gene subsampling once")
                axis(1, at = 1:L, labels = sub_s)     # for x axis
                dev.off()
            }
            
        }
        
    }
    
    ##########################
    #### correlation per gene
    {
    ## cor funciton 
    gene.cor <- function(x, y, group,  file)
    {
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
            
            dat <- data.frame(cor_ga, cor_gty)
            g <- ggplot(dat, aes(x = cor_gty, y = cor_ga, fill = cor_gty, show.legend=F)) + geom_boxplot(outlier.shape = NA) +  theme_bw() + theme(axis.line = element_line(colour = "black"),panel.grid.major = element_blank(),panel.grid.minor = element_blank(),panel.border =  element_blank(),  panel.background = element_blank()) 
            g <- g + scale_fill_manual( values = c("cadetblue3", "coral2"))
            g <- g + geom_hline(yintercept = 0, lty = 2, col = "red") + labs(x = "Sample Group", y = "Cor(m6A, mRNA)") + theme(legend.position="none")
            name_2 <- paste(file, "_per_gene.pdf", sep = "")
            ggsave(name_2,  width = 3, height = 3, units = "in")
            
            out <- list("cor.gn" = cor_gn, "cor.gn.pval" = cor_gn_p, "cor.gt" = cor_gt, "cor.gt.pval" = cor_gt_p)
            return(out)
        }   
    }
    
    # for all tumor samples 
    cor_ntg <- gene.cor(m6a_rs, mrna_s, s_type, "cor_m6A_level_and_mRNA_RPKM_all_gene")
    save(cor_ntg, file = "cor_m6A_level_and_mRNA_per_gene_all_samples.Rdata")
    
    # for randomly selected 10 tumor samples
    idx_rs <- sample(1:51, 10)
    m6a_rs_20 <- cbind(m6a_rn, m6a_rt[, idx_rs])                     # merger normal and tumor
    mrna_s_20 <- cbind(input_rpkm_n, input_rpkm_t[, idx_rs])
    s_type_20 <- c(rep("Normal", 10), rep("Tumor", 10))
    gene.cor(m6a_rs_20, mrna_s_20, s_type_20, "cor_m6A_level_and_mRNA_RPKM_all_gene_tumor10")
    }
    
    ##############################
    # cor_n and cor_t separately 
    {
    r_cut <- 0.3
    p_cut <- 0.05     
    
    # normal 
    cor_gn <- cor_ntg$cor.gn
    cor_gn_p <- cor_ntg$cor.gn.pval
    cor_gn_p_adj <- p.adjust(cor_gn_p)      
    idx_n <- abs(cor_gn) > r_cut & cor_gn_p < p_cut
    
    idx_np <- cor_gn > r_cut & cor_gn_p < p_cut     # positive
    idx_nn <- cor_gn < -r_cut & cor_gn_p < p_cut    # negative
    gene_np <- names(cor_gn)[idx_np]
    gene_nn <- names(cor_gn)[idx_nn]
    write.table(gene_np, file = "cor_m6a_mRNA_sig_pos_in_normal.txt", quote = F, col.names = F, row.names = F)
    write.table(gene_nn, file = "cor_m6a_mRNA_sig_neg_in_normal.txt", quote = F, col.names = F, row.names = F)
    
    # tumor 
    cor_gt <- cor_ntg$cor.gt
    cor_gt_p <- cor_ntg$cor.gt.pval
    cor_gt_p_adj <- p.adjust(cor_gt_p)
    idx_t <- abs(cor_gt) > r_cut & cor_gt_p < p_cut
    
    idx_tp <- cor_gt > r_cut & cor_gt_p < p_cut     # positive
    idx_tn <- cor_gt < -r_cut & cor_gt_p < p_cut    # negative
    gene_tp <- names(cor_gt)[idx_tp]
    gene_tn <- names(cor_gt)[idx_tn]
    write.table(gene_tp, file = "cor_m6a_mRNA_sig_pos_in_tumor.txt", quote = F, col.names = F, row.names = F)
    write.table(gene_tn, file = "cor_m6a_mRNA_sig_neg_in_tumor.txt", quote = F, col.names = F, row.names = F)
    }
    
    #########################
    ## volcano plot for cor
    {
        # plot function
        cor.vol.plot <- function(cor_r, cor_p, r_c, p_c, ytxt, nytxt, name)
        {
            L <- length(cor_r)
            idx_u <- cor_r > r_c & cor_p < p_c
            idx_d  <- cor_r < -r_c & cor_p < p_c
            idx_nu <- (cor_r >= 0 & cor_r <= r_c ) | (cor_r > r_c & cor_p >= p_c)
            idx_nd  <- (cor_r <= 0 & cor_r >= -r_c) | (cor_r < -r_c & cor_p >= p_c)
            
            ## chisquare test
            t <- matrix(c(sum(idx_u), sum(idx_d), sum(idx_nu), sum(idx_nd)), 2, 2)
            tc <- chisq.test(t)
            pval <- as.numeric(tc$p.value)
            
            if(pval > 2.2e-16)
            {
                pval <- formatC(pval, format = "e", digits = 2)
                sub_t <- paste("Chisq.test: p_value = ", pval, sep = "")        # correlation coefficient
                
            } else{
                
                sub_t <- paste("Chisq.test: p_value < 2.2e-16") 
            }
            
            ## color
            col_c <- rep("gray", L)
            col_c[idx_u] <- "#FFC000"    #viridis(2)[2]  # "red"
            col_c[idx_d] <- "#702FA0"    #viridis(2)[1]  #"blue"
            
            library("viridis")  
            
            ## 
            #idx_p <- -log10(cor_p) > 6
            #cor_p[idx_p] <- 10^(-6)
            
            pdf(file = name, width = 3, height = 3)
            par(mar = c(3, 3, 1, 1), mgp = c(2, 0.75, 0))
            plot(cor_r, -log10(cor_p), col = col_c, xlab = "Correlation", ylab = "-log10(p)", 
                 xlim = c(-1, 1), ylim = c(0, max(-log10(cor_p))), main = sub_t, cex.main = 0.6, font.main = 3)
            #abline(v = c(-r_c, r_c), h = -log10(p_c), lty = 2)
            abline(h = -log10(p_c), lty = 2)
            lines(x = c(-r_c, -r_c), y = c(-log10(p_c), max(-log10(cor_p))), lty = 2)                   #  max(-log10(cor_p))
            lines(x = c(r_c, r_c), y = c(-log10(p_c),  max(-log10(cor_p))), lty = 2)
            lines(x = c(0, 0), y = c(0, -log10(p_c)), lty = 2)
            text(x = c(-0.6, 0.6), y = c(ytxt, ytxt), labels = c(sum(idx_d), sum(idx_u)), col = viridis(2))                 ## c("blue", "red") 
            text(x = c(-0.6, 0.6), y = c(nytxt, nytxt), labels = c(sum(idx_nd), sum(idx_nu)), col = c("gray", "gray"))
            
            dev.off()
            
        }
        
        cor.vol.plot(cor_gn, cor_gn_p, 0.3, 0.05, 4.5, 0.5,  "Cor_sig_per_gene_in_normal.pdf")
        cor.vol.plot(cor_gt, cor_gt_p, 0.3, 0.05, 4.5, 0.5, "Cor_sig_per_gene_in_tumor.pdf")
        
        ## venn plot for shared genes
        pdf(file = "cor_per_gene_sig_cmp_vennplot.pdf", width = 2, height = 2)
        par(mar = c(0, 0, 0, 0))
        draw.pairwise.venn(sum(idx_n), sum(idx_t), sum(idx_n & idx_t), fill = c("cadetblue3", "coral2"), alpha = rep(0.9, 2))
        dev.off()
    }
}


####################################################################
## compare correlation coefficent between normal and tumor per gene 
#####################################################################
{
    library("Hmisc")         # for rcorr
    library("matrixStats")
    library("gplots")
    library("psych")           # r.test  function 
    
    #####
    load("cor_m6A_level_and_mRNA_per_gene_all_samples.Rdata")
    load("cor_m6a_mRNA_tumor_subsample10_1000.Rdata")
    
    ################################
    # compare correlation coefficent
    ## per gene 
    cor_a <- cbind(cor_ntg$cor.gn, cor_ntg$cor.gt, cor_10mat_mean, rowMedians(cor_10mat_sort))
    colnames(cor_a) <- c("nomal_10", "tumor_51", "tumor_10_mean", "tumor_10_median")
    png("cor_per_gene_tumor_vs_normal_with_subsampling.png", width = 7, height = 4, units = "in", res = 300)
    par(mar = c(3, 3, 0.5, 0.5), mgp = c(2.3, 1, 0))
    boxplot(cor_a, ylab = "correlation")
    abline(h = 0, lty = 2)
    dev.off()
    
    ## permutation test based on the mean/median of all genes
    pdf("cor_m6A_mRNA_per_gene_permutaiton_subsampling_1000times_10tumor_Mean.pdf", width = 3, height = 3)
    par(mar = c(4, 3, 1, 0.5), mgp = c(2, 0.75, 0))
    plot(density(colMeans(cor_10mat, na.rm = T)), main = "subsmapling 10 tumors 1000 times", 
         xlab = "Mean of cor(m6A, mRNA)_g per sampling", y = "Density", xlim = c(-0.3, 0.1), cex.main = 0.75, cex.lab = 0.75)
    abline(v = c(mean(cor_ntg$cor.gn), mean(cor_ntg$cor.gt)), col = c("red", "blue"), lty = 2)
    abline(h = 0, col = "gray", lty = 2)
    text(x = mean(cor_ntg$cor.gn), y = 15, labels = "corG_NL", pos = 4, col = "red", cex = 0.5)
    text(x = mean(cor_ntg$cor.gt), y = 15, labels = "corG_Tumor", pos = 4, col = "blue", cex = 0.5)
    dev.off()
    
    pdf("cor_m6A_mRNA_per_gene_permutaiton_subsampling_1000times_10tumor_Median.pdf", width = 3, height = 3)
    par(mar = c(4, 3, 1, 0.5), mgp = c(2, 0.75, 0))
    plot(density(colMedians(cor_10mat, na.rm = T)), main = "subsmapling 10 tumors 1000 times", 
         xlab = "Median of cor(m6A, mRNA)_g per sampling", y = "Density", xlim = c(-0.3, 0.1), cex.main = 0.75, cex.lab = 0.75)
    abline(v = c(median(cor_ntg$cor.gn), median(cor_ntg$cor.gt)), col = c("red", "blue"), lty = 2)
    abline(h = 0, col = "gray", lty = 2)
    text(x = median(cor_ntg$cor.gn), y = 13, labels = "corG_NL", pos = 4, col = "red", cex = 0.5)
    text(x = median(cor_ntg$cor.gt), y = 13, labels = "corG_Tumor", pos = 4, col = "blue", cex = 0.5)
    dev.off()
    
    ######################################
    ## fisher r to z test: single tail test
    ## cutoff 0.05 
    
    cor_nr <- cor_ntg$cor.gn
    cor_tr <- cor_ntg$cor.gt
    
    ## r to z transform
    cor_nr_z <- fisherz(cor_nr)
    cor_tr_z <- fisherz(cor_tr)
    
    t.test(cor_nr_z, cor_tr_z)
    wilcox.test(cor_nr_z, cor_tr_z)
    
    cor_ga <- c(cor_nr_z, cor_tr_z)
    cor_gty <- c(rep("Normal", length(cor_nr_z)), rep("Tumor", length(cor_tr_z)))
    
    library(BSDA)
    z.test(cor_nr_z, cor_tr_z, sigma.x = sd(cor_nr_z), sigma.y = sd(cor_tr_z), alternative = "less")
    
    outlier.shape = NA
    
    dat <- data.frame(cor_ga, cor_gty)
    g <- ggplot(dat, aes(x = cor_gty, y = cor_ga, fill = cor_gty, show.legend=F)) + geom_boxplot(outlier.shape = NA) +  theme_bw() + theme(axis.line = element_line(colour = "black"),panel.grid.major = element_blank(),panel.grid.minor = element_blank(),panel.border =  element_blank(),  panel.background = element_blank()) 
    g <- g + scale_fill_manual( values = c("cadetblue3", "coral2")) + ylim(-1.5, 1)
    g <- g + geom_hline(yintercept = 0, lty = 2, col = "red") + labs(x = "Sample Group", y = "Cor(m6A, mRNA)_Z") + theme(legend.position="none")
    name_2 <- "cor_m6A_level_and_mRNA_RPKM_all_gene_r2z.pdf"
    ggsave(name_2,  width = 2.5, height = 3)
    
    
    ## fisher r to z test 
    N <- nrow(cor_10mat_sort)
    M_n <- 10
    M_t <- 51
    r_z <- vector()
    r_zp <- vector()
    for(i in 1:N)
    {
        tmp <- r.test(M_t, cor_tr[i], cor_nr[i], n2 = M_n)
        r_z[i] <- tmp$z
        r_zp[i] <- tmp$p
    }
    
    r_z_s <- r_z*(cor_tr_z - cor_nr_z)
    plot(r_z, -log10(r_zp))
    
    #r_zp_adj <- p.adjust(r_zp)
    ## considering significance separately
    idx_ns <- abs(cor_nr) > 0.3 & cor_ntg$cor.gn.pval < 0.05
    idx_ts <- abs(cor_tr) > 0.3 & cor_ntg$cor.gt.pval < 0.05
    idx_sig <- idx_ns | idx_ts       # cor are significant in at leat one group
    idx_s <- r_zp < 0.05 & idx_sig
    
    cor_cut = 0   ##  or 0            # consider fisher r to z test only
    col_c <- rep("gray", N)
    idx_1 <- cor_nr_z - cor_tr_z > cor_cut
    idx_2 <- cor_nr_z - cor_tr_z < -cor_cut
    idx_s1 <- idx_1 & idx_s
    idx_s2 <- idx_2 & idx_s
    col_c[idx_s1] <- "red"
    col_c[idx_s2] <- "blue"
    
    g_r <- names(cor_nr)[idx_s1]
    g_b <- names(cor_nr)[idx_s2]
    save(g_r, g_b, file = "cor_diff_genes_fisher_r2z_test_filter.Rdata")
    
    table(col_c)
    
    x <- c(cor_nr_z[!idx_s], cor_nr_z[idx_s1], cor_nr_z[idx_s2])
    y <- c(cor_tr_z[!idx_s], cor_tr_z[idx_s1], cor_tr_z[idx_s2])
    col_cc <- c(col_c[!idx_s], col_c[idx_s1], col_c[idx_s2])
    
    z_mm <- cor_tr_z - cor_nr_z
    
    png(paste0("cor_per_gene_tumor_vs_normal_transformed_r_filter_", cor_cut, ".png"), width = 4, height = 4, units = "in", res = 300)
    par(mar = c(3, 3, 0.5, 0.5), mgp = c(2, 0.5, 0))
    plot(x, y, col = col_cc, xlim = c(-2.3, 2.3), ylim = c(-1.3, 1.3),  xlab = "Z transformed r (NL)", ylab = "Z transformed r (Tumor)")
    abline( a = 0, b = 1,  lty = 2)
    abline( a = max(z_mm[idx_s1]), b = 1,  lty = 3)
    abline(  a = min(z_mm[idx_s2]), b = 1,  lty = 3)
    text(x = -1.5, y = 1, labels = paste0("N = ", table(col_c)[1]), col = "blue")
    text(x = 1.5, y = 1, labels = paste0("N = ", table(col_c)[2]), col = "gray")
    text(x = 1.5, y = -1, labels = paste0("N = ", table(col_c)[3]), col = "red")
    dev.off()
    
    
    cor_g1 <- cbind(cor_nr_z[idx_s1], cor_tr_z[idx_s1])
    g1 <- names(cor_nr_z)[idx_s1]
    colnames(cor_g1) <- c("Normal", "Tumor")
    png(paste0("cor_per_gene_tumor_vs_normal_transformed_group1_boxplot_filter_", cor_cut, ".png"), width = 2.3, height = 2, units = "in", res = 300)
    par(mar = c(2, 3, 0.5, 0.5), mgp = c(1.5, 0.5, 0))
    boxplot(cor_g1, col = c("cadetblue3", "coral2"), ylab = "transformed R")
    abline( h = 0,  lty = 2)
    dev.off()
    
    cor_g2 <- cbind(cor_nr_z[idx_s2], cor_tr_z[idx_s2])
    g2 <- names(cor_nr_z)[idx_s2]
    colnames(cor_g2) <- c("Normal", "Tumor")
    png(paste0("cor_per_gene_tumor_vs_normal_transformed_group2_boxplot_filter_", cor_cut, ".png"), width = 2.3, height = 2, units = "in", res = 300)
    par(mar = c(2, 3, 0.5, 0.5), mgp = c(1.5, 0.5, 0))
    boxplot(cor_g2, col = c("cadetblue3", "coral2"), ylab = "transformed R")
    abline( h = 0,  lty = 2)
    dev.off()
    
    cor_g3 <- cbind(cor_nr_z[!idx_s], cor_tr_z[!idx_s])
    g3 <- names(cor_nr_z)[!idx_s]
    colnames(cor_g3) <- c("Normal", "Tumor")
    png(paste0("cor_per_gene_tumor_vs_normal_transformed_group3_boxplot_filter_", cor_cut, ".png"), width = 2.3, height = 2, units = "in", res = 300)
    par(mar = c(2, 3, 0.5, 0.5), mgp = c(1.5, 0.5, 0))
    boxplot(cor_g3, col = c("cadetblue3", "coral2"), ylab = "transformed R")
    abline( h = 0,  lty = 2)
    dev.off()
}