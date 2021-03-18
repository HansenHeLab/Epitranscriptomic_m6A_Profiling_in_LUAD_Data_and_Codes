################################################################
## This script implements:
## 1. m6A related enzymes mRNA level comparison
## 2. correlation between m6a level and enzyme per each sample
## 3. correlation between m6a level and enzyme per each gene
#################################################################

rm(list = ls())
setwd("/work/dir/")

############################################
## m6A related enzymes mRNA level comparison
############################################
{
    m6a_en <- read.table("/Users/Yong/Yong/m6a_profiling/4_enzymes/m6a_enzyme.txt", header = T)      # enzyme list  
    load("/Users/Yong/Yong/m6a_profiling/1_mapping_quant/5_RNA-seq/htseq/count/paired_61_rpkm.Rdata")        # rpkm for IP and input 
    load("/Users/Yong/Yong/m6a_profiling/1_mapping_quant/5_RNA-seq/TCGA/GDC_LUAD/tcga_htseq_fpkm_log2.Rdata")        # FPKM from TCGA
    
    N_n <- 10     # number of normal sample
    N_t <- 51     # number of tumor sample
    s_type <- c(rep("normal", N_n), rep("tumor", N_t))
    
    idx_en <- match(m6a_en$gene_name, rownames(input_rpkm))
    en_cmb <- as.matrix(input_rpkm[idx_en, ])          
    idx <- m6a_en$function. == "writer" | m6a_en$function. == "eraser"          # writer or eraser
    en_cmb_we <- log2(en_cmb[idx, ] + 1)                         # writer and eraser
    en_cmb_re <- log2(en_cmb[!idx, ] + 1)                        # Reader
    en_cmb_all <- log2(en_cmb + 1)
    mrna_all  <- log2(input_rpkm + 1)                           # all 
    
    save(en_cmb_we, en_cmb_re, en_cmb_all, file = "en_cmb_we_re_log2_rpkm.Rdata")
    
    #### GDC TCGA LUAD data  : Normal: 59
    tcga_g <- rownames(tcga)
    idx_ki <- match("VIRMA", tcga_g)         # KIAA1429 uses name VIRMA
    tcga_g[idx_ki] <- "KIAA1429"
    idx_we <- match(rownames(en_cmb_we), tcga_g)
    idx_re <- match(rownames(en_cmb_re), tcga_g)
    en_cmb_we_t <- tcga[idx_we, ]
    en_cmb_re_t <- tcga[idx_re, ] 
    idx_1 <- tcga_type == "normal"
    idx_2 <- tcga_type == "tumor"
    
    ## boxplot for both our own data and TCGA data
    en.boxplot <- function(en_expr, sample_type, en_name, ww, hh)
    {
        L <- nrow(en_expr)
        idx_n <- sample_type == "normal"
        idx_t <- sample_type == "tumor"
        
        en_s <- list()
        pval <- vector()
        for (i in 1:L)
        {
            en_s[[2*i - 1]] <- en_expr[i, idx_n]
            en_s[[2*i]] <- en_expr[i, idx_t]
            
            ## test 
            tmp <- wilcox.test(en_expr[i, idx_n], en_expr[i, idx_t])
            pval[i] <- tmp$p.value
            
        }
        
        pdf(file = en_name, width = ww, height = hh)
        par(mar = c(5, 5, 0.5, 0.5))
        boxplot(en_s, at = 1:(2*L) + c(0.05, -0.1), col = c("cadetblue3", "coral2"), outline = F, xaxt = "n", ylab = "log2(RPKM + 1)")
        xl <- rownames(en_expr)
        axis(1, at = seq(1, 2*L, 2) + 0.5, labels = xl , las = 2)     # for x axis
        dev.off()
        
        return(pval)
    }
    
    en.boxplot(en_cmb_we, s_type, "m6a_writer_eraser_mRNA_boxplot.pdf", 6, 4)
    en.boxplot(en_cmb_we_t, tcga_type, "m6a_writer_eraser_mRNA_TCGA_boxplot.pdf", 6, 4)
    en.boxplot(en_cmb_re, s_type, "m6a_reader_mRNA_boxplot.pdf", 7, 4)
    en.boxplot(en_cmb_re_t, tcga_type, "m6a_reader_mRNA_TCGA_boxplot.pdf", 7, 4)
   
    ## only includes 8we
     en_cmb_8we <- en_cmb_we[-c(3, 5, 6), ]
     en_cmb_8we_t <- en_cmb_we_t[-c(3, 5, 6), ]
     
     en.boxplot(en_cmb_8we, s_type, "m6a_8writer_eraser_mRNA_boxplot.pdf", 6, 3)
     en.boxplot(en_cmb_8we_t, tcga_type, "m6a_8writer_eraser_mRNA_TCGA_boxplot.pdf", 6, 3)
    
     ##########################################################
     ## boxplot with shade 
     ## 8 writers and 2 erasers, not inclued METTL16 for snRNA
     {
     ## for TCGA
     idx_n <- tcga_type == "normal"
     tcga_type_s <- c(tcga_type[idx_n], tcga_type[!idx_n])
     col_c <- c( "coral2", "cadetblue3")
     en_cmb_8we_ts <- cbind(en_cmb_8we_t[, idx_n], en_cmb_8we_t[, !idx_n])
     
     ## for both 
     expr_we_list <- list(en_cmb_8we_ts, en_cmb_8we)
     outname_list <- c("m6a_8writer_eraser_mRNA_TCGA_boxplot_h.pdf", "m6a_8writer_eraser_mRNA_boxplot_h.pdf")
     
     for(i in 1:2)
     {
     ## transfer to dataform
     expr_we <- expr_we_list[[i]]
     x <- as.vector(t(expr_we)) 
     y1 <- rep(rownames(expr_we), each = ncol(expr_we))
     y2 <- paste(y1, sample_type, sep = "_")
     y <- factor(y2, levels = rev(unique(y2)))
     dat <- data.frame(x, y)

     ## plot
     library("BoutrosLab.plotting.general")
     pdf(file = outname[i] , width = 4, height = 6)
     create.boxplot(formula = y ~ x,
                    data = dat,
                    col = col_c,
                    xaxis.cex = 1,
                    xlimits = c(1, 6),
                    yaxis.cex = 1,
                    outliers = FALSE,   
                    
                    ## for shade
                    # draw rectangles
                    add.rectangle = TRUE,
                    # coordinates of rectangles given by four parameters
                    xleft.rectangle = 0,
                    xright.rectangle = 9,
                    ybottom.rectangle = seq(0.5, 18.5, 4),
                    ytop.rectangle = seq(2.5, 19.5, 4),
                    # set rectangle colour
                    col.rectangle = "grey",
                    # set rectangle alpha (transparency)
                    alpha.rectangle = 0.25)
    dev.off()
     }
     
    }
     
    ##################### 
    ## enzyme correlation
    {
    corr.plot <- function(expr, name, groups)
    {
        library(corrplot)
        M <- cor(expr)
        png(file = name, res = 600, width = 5, height = 5, units = "in")
        par(mar = c(2, 5, 5, 2 ))
        corrplot(M, order = "hclust", addrect = groups)
        dev.off()
    }
    corr.plot(t(en_cmb_we[, 1:10]), "WE_normal_corrplot.png", 3)
    corr.plot(t(en_cmb_we[, 11:61]), "WE_tumor_corrplot.png", 3)
    corr.plot(t(en_cmb_we_t[, idx_1]), "WE_normal_TCGA_corrplot.png", 3)
    corr.plot(t(en_cmb_we_t[, idx_2]), "WE_tumor_TCGA_corrplot.png", 3)
    }
     
    ################################
    #### enzyme correlatioin network
    {
    en_8we_tumor <- t(en_cmb_we[-c(3, 5, 6), 11:61])
    
    library(tidyverse)  
    library(corrr)
    
    cor_res <- correlate( en_8we_tumor)
    pdf("cor_between_8we_network_plot.pdf", width = 3, height = 3)
    network_plot(cor_res, min_cor = 0.0,  colours = c("blue", "gray", "red"),  curved = T)
    dev.off()
    
    ## for TCGA 
    en_8we_tumor_t <- t(en_cmb_we_t[-c(3, 5, 6), idx_2])
    cor_res <- correlate( en_8we_tumor_t)
    
    pdf("cor_between_8we_network_plot_TCGA.pdf", width = 3, height = 3)
    network_plot(cor_res, min_cor = 0.0,  colours = c("blue", "gray", "red"),  curved = T)
    dev.off()
    }
    
}

###########################################################
## correlation between m6a level and enzyme per each sample
###########################################################
{
    #############################
    ## loading data and functions
    {
    library("Hmisc")
    library("ggplot2")
    library("s2dverification")       # for color bar separately
    library("corrplot")
    
    source("/Users/Yong/Yong/R/functions/cor.plot.R")
    
    ## total IP count based corrlation 
    load("/Users/Yong/Yong/m6a_profiling/4_enzymes/1_m6a_enzyme_mRNA_protein/en_cmb_we_re_log2_rpkm.Rdata")     # log 2 transformed
    uniq_reads <- read.csv("/Users/Yong/Yong/m6a_profiling/1_mapping_quant/61samples_input_ip_uniquely_mapped_reads.csv", row.names = 1)
    uniq_reads_log <- log2(uniq_reads)
    
    ## cor return with p value 
    pval.cor <- function(m6a, enzyme)
    {
        N <- nrow(enzyme)
        cor_p <- matrix(0, N, 2)
        rownames(cor_p) <- rownames(enzyme)
        colnames(cor_p) <- c("cor_R", "cor_pval")
        
        for (i in 1 : N )
        {
            tmp <- cor.test(m6a, enzyme[i, ])
            cor_p[i, 1] <- tmp$estimate
            cor_p[i, 2] <- tmp$p.value
        }
        
        return(cor_p)
    }
    
    }
    
    #############################################################################################
    ## Global m6A level meausred by the uniquely mapped IP read, normalized by totoal ek12 reads
    {
        cor_n <- cor(uniq_reads[1:10, ])
        cor_t <- cor(uniq_reads[11:61, ])
        
        png(file = "cor_uniqely_mapped_read_in_normal.png", width = 4, height = 4,  unit = "in", res = 600)
        corrplot.mixed(cor_n)
        dev.off()
    
        cor.plot(uniq_reads$IP_ek12[1:10], uniq_reads$IP_human[1:10], "Normal_IP_read_ek12_vs_human", "ek12", "human")
        cor.plot(uniq_reads_log$IP_ek12[1:10], uniq_reads_log$IP_human[1:10], "Normal_IP_read_log_ek12_vs_human", "ek12", "human")
        
        png(file = "cor_uniqely_mapped_read_in_tumor.png", width = 4, height = 4,  unit = "in", res = 600)
        corrplot.mixed(cor_t)
        dev.off()
        
        cor.plot(uniq_reads$IP_ek12[11:61], uniq_reads$IP_human[11:61], "Tumor_IP_read_ek12_vs_human", "ek12", "human")
        cor.plot(uniq_reads_log$IP_ek12[11:61], uniq_reads_log$IP_human[11:61], "Tumor_IP_read_log_ek12_vs_human", "ek12", "human")
        
        ## ip reads normalized by the total ek12 reads
        ip_ek12_norm <- uniq_reads$IP_human * (uniq_reads$IP_ek12[1] /  uniq_reads$IP_ek12)
        names(ip_ek12_norm) <- rownames(uniq_reads)
        
        png("global_m6a_all_IP_reads_normalized_by_total_ek12.png", height = 3, width = 3, unit = "in", res = 600)
        par(mar = c(4, 4, 0.5, 0.5))
        boxplot(list(ip_ek12_norm[1:10], ip_ek12_norm[11:61]), xlab = ("Nromal      Tumor"),  ylab = "Normalized reads")
        dev.off()
        
        png("cor_we_global_m6a_all_IP_reads_normalized_by_total_ek12_normal.png", height = 3, width = 3, unit = "in", res = 600)
        par(mar = c(5, 3, 0.5, 0.5))
        barplot(cor(log2(ip_ek12_norm[1:10]), t(en_cmb_we[, 1:10])), las = 2)
        dev.off()
        
        png("cor_we_global_m6a_all_IP_reads_normalized_by_total_ek12_tumor.png", height = 3, width = 3, unit = "in", res = 600)
        par(mar = c(5, 3, 0.5, 0.5))
        barplot(cor(log2(ip_ek12_norm[11:61]), t(en_cmb_we[, 11:61])), las = 2)
        dev.off()
        
        
        
    }
    
    #############################################################################################
    ## Global m6A level meausred by the uniquely mapped IP read, normalized by ek12 reads regressioin
    {
        
        load("/Users/Yong/Yong/m6a_profiling/1_mapping_quant/5_RNA-seq/htseq/count/uniquely_mapped_ek12_reads.Rdata")
        ip_all_count_log2_ercc <- uniq_reads_log$Input_human
        names(ip_all_count_log2_ercc) <- rownames(uniq_reads_log)
        N3 <- length(ip_all_count_log2_ercc)
        
        ip_ek12_count_log2_61 <- ip_ek12_count_log2[, -c(11, 49, 64)]
        name_l <- strsplit(colnames(ip_ek12_count_log2_61), "_")
        name_v <- vector()
        for(i in 1:N3)
        {
            name_v[i] <- name_l[[i]][1]
        }
        
        idx_n <- match(names(ip_all_count_log2_ercc), name_v)
        ip_ek12_count_log2_61_s <- ip_ek12_count_log2_61[, idx_n]
        
        for(i in 2:N3)
        {
            
            idx1 <- rowSums(ip_ek12_count_log2_61_s[, c(1, i)] > 3) == 2              # selected cut off 3
            tmp1 <- lm(ip_ek12_count_log2_61_s[idx1, i]~ip_ek12_count_log2_61_s[idx1, 1])
            ip_all_count_log2_ercc[i] <- (uniq_reads_log$Input_human[i] - tmp1$coefficients[1]) / tmp1$coefficients[2]
        }
        
        ip_all_count_ercc <- 2^ip_all_count_log2_ercc
        
        png("global_m6a_all_IP_reads_normalized_by_ek12_regression.png", height = 3, width = 3, unit = "in", res = 600)
        par(mar = c(4, 4, 0.5, 0.5))
        boxplot(list(ip_all_count_ercc[1:10], ip_all_count_ercc[11:61]), xlab = ("Nromal      Tumor"),  ylab = "Normalized reads")
        dev.off()
        
        sample_type <- c(rep("normal", 10), rep("tumor", 51))
        dat <- data.frame(ip_all_count_log2_ercc, sample_type )
        g <- ggplot(dat, aes(x = sample_type, y = ip_all_count_log2_ercc, col = sample_type)) + geom_boxplot(outlier.shape = NA) 
        g <- g + geom_point(aes(col = sample_type), position = position_jitter(width = 0.3)) 
        g <- g + theme_bw() + theme(axis.line = element_line(colour = "black"),panel.grid.major = element_blank(),panel.grid.minor = element_blank(),panel.border =  element_blank(),  panel.background = element_blank()) 
        g <- g + scale_color_manual( values = c("cadetblue3", "coral2")) + theme(legend.position = "none")
        ggsave("Global_m6A_level_by_IP_reads_regression_ek12_boxplot_61samples.pdf", width = 2.5, height = 3, units = "in")
        
        png("cor_we_global_m6a_all_IP_reads_normalized_by_ek12_regression_normal.png", height = 3, width = 3, unit = "in", res = 600)
        par(mar = c(5, 3, 0.5, 0.5))
        barplot(cor(log2(ip_all_count_ercc[1:10]), t(en_cmb_we[, 1:10])), las = 2)
        dev.off()
        
        cor_pval_n <- pval.cor(log2(ip_all_count_ercc[1:10]), en_cmb_we[ ,1:10])
        cor_pval_n
        
        png("cor_we_global_m6a_all_IP_reads_normalized_by_ek12_regression_tumor.png", height = 3, width = 3, unit = "in", res = 600)
        par(mar = c(5, 3, 0.5, 0.5))
        barplot(cor(ip_all_count_ercc[11:61], t(en_cmb_we[, 11:61])), las = 2)
        dev.off()
        
        cor_pval_t <- pval.cor(log2(ip_all_count_ercc[11:61]), en_cmb_we[ ,11:61])
        cor_pval_t
        
        png("cor_we_global_m6a_all_IP_reads_normalized_by_ek12_regression_tumor_normal.png", height = 3, width = 3, unit = "in", res = 600)
        par(mar = c(5, 3, 0.5, 0.5))
        barplot(cor(ip_all_count_ercc, t(en_cmb_we)), las = 2)
        dev.off()
        
        
    }
    
    #############################################################################################
    ## Global m6A level measured by totoal IP read count (from htseq-count, based on annotation)
    {
        load("/Users/Yong/Yong/m6a_profiling/1_mapping_quant/5_RNA-seq/htseq/count/paired_61_count_ercc_deseq.Rdata")      # for all genes ~50,000  
        ip_sum <- colSums(ip_count_norm)
        input_sum <- colSums(input_count_norm)
        m6a_g  <- ip_sum          ## ip reads sum for global m6A levle 
        # m6a_g <- ip_ek12_norm     ## all IP uniquely mapped reads normalized by ek12 reads
        
        sample_type <- c(rep("normal", 10), rep("tumor", 51))
        dat <- data.frame(m6a_g, sample_type)
        
        g <- ggplot(dat, aes(x = sample_type, y = m6a_g, col = sample_type)) + geom_boxplot(outlier.shape = NA) 
        g <- g + geom_point(aes(col = sample_type), position = position_jitter(width = 0.3)) 
        g <- g + theme_classic() +  scale_color_manual( values = c("cadetblue3", "coral2")) + theme(legend.position = "none")
        ggsave("Global_m6A_level_by_IP_reads_sum_boxplot_61samples.pdf", width = 2.5, height = 3, units = "in")
        #ggsave("Global_m6A_level_by_IP_reads_sum_ek12_boxplot_61samples.pdf", width = 2.5, height = 3, units = "in")
    }
    
    #############################################
    ## multiple linear regression: based on m6a_g  
    {
        library(relaimpo)
        m6a_g <- ip_all_count_ercc
        
        x_n <- data.frame(t(en_cmb_we[-c(3,5, 6), 1:10]))
        x_t <- data.frame(t(en_cmb_we[-c(3,5, 6), 11:61]))
        y_m6a <- log2(m6a_g)           # 
        dat_n <- data.frame(y = y_m6a[1:10], x_n)
        dat_t <- data.frame(y = y_m6a[11:61], x_t)
        
        fit_n  <- lm(y ~ ., data = dat_n)            # remove three varabiles as small sample size
        fit_t <-  lm(y ~ ., data = dat_t)
        
        if (FALSE){
            summary(fit)
            layout(matrix(c(1,2,3,4),2,2))
            plot(fit)
        }
        
        png(file = "Multiple_LM_8WE_m6A_per_sample_normal.png", width = 8, height = 4, units = "in", res = 600)
        par(mfrow = c(1, 2), mar = c(5,4,0.5,0.5))
        barplot(fit_n$coefficients, las = 2, ylab = "Coefficient")
        tmp <- calc.relimp(fit_n,type= "lmg", rela=TRUE)
        barplot(tmp$lmg, las = 2, ylab = "% of R2")
        dev.off()
        
        png(file = "Multiple_LM_8WE_m6A_per_sample_tumor.png", width = 8, height = 4, units = "in", res = 600)
        par(mfrow = c(1, 2), mar = c(5,4,0.5,0.5))
        barplot(fit_t$coefficients, las = 2, ylab = "Coefficient")
        tmp <- calc.relimp(fit_t,type= "lmg", rela=TRUE)
        barplot(tmp$lmg, las = 2, ylab = "% of R2")
        dev.off()
    
        # Bootstrap Measures of Relative Importance (100 samples) 
        boot <- boot.relimp(fit_t, b = 100, type = "lmg", rank = TRUE,  diff = TRUE, rela = TRUE)
        booteval.relimp(boot)       # print result
        
        png(file = "Multiple_LM_8WE_m6A_per_sample_tumor_R2.png", width = 4, height = 5, units = "in", res = 600)
        par(mar = c(5,4,0.5,0.5))
        plot(booteval.relimp(boot,sort=F), las = 2)   # plot result
        dev.off()
    
        # Bootstrap Measures of Relative Importance (100 samples) 
        boot <- boot.relimp(fit_n, b = 100, type = "lmg", rank = TRUE,  diff = TRUE, rela = TRUE)
        booteval.relimp(boot)       # print result
        
        png(file = "Multiple_LM_8WE_m6A_per_sample_normal_R2.png", width = 4, height = 5, units = "in", res = 600)
        par(mar = c(5,4,0.5,0.5))
        plot(booteval.relimp(boot,sort=F), las = 2)   # plot result
        dev.off()
        
    }
    
}


###########################################################
## correlation between m6a level and enzyme per each gene
###########################################################
{
    #############################
    ## loading data and functions
    {
    library("gplots")       # for heatmap.2
    library("ggplot2")
    load("/Users/Yong/Yong/m6a_profiling/4_enzymes/1_m6a_enzyme_mRNA_protein/en_cmb_we_re_log2_rpkm.Rdata")     # log 2 transformed
    load("/Users/Yong/Yong/m6a_profiling/3_m6a_level/m6a_level_keep_dup/gene_with_peak/m6A_level_4gene_passed_filtering.Rdata")
    
    write.csv(en_cmb_we, file = "IMPACT_LUAD_m6A_writers_erasers_log2_RPKM.csv")
    m6a_s <- cbind(m6a_rn, m6a_rt)
    write.csv(m6a_s, file = "IMPACT_LUAD_methylated_gene_m6A_level.csv")
    }
    
    #########################################
    ## correaltion per gene: cor(m6a, enzyme)
    {
        ## cor(m6A level, enzyme) per gene
        m6a_rn_log2 <- log2(m6a_rn)
        idx <- is.infinite(m6a_rn_log2)
        m6a_rn_log2[idx] <- NA                  # for correlation 
        m6a_rt_log2 <- log2(m6a_rt)
        idx <- is.infinite(m6a_rt_log2)
        m6a_rt_log2[idx] <- NA    
        m6a_mm_log2 <- cbind(m6a_rn_log2, m6a_rt_log2)
        
        N1 <- nrow(m6a_rn)
        N2 <- nrow(en_cmb_all)
        cor_mn <- matrix(0, N1, N2)
        cor_mn_p <- matrix(0, N1, N2)
        cor_mt <- matrix(0, N1, N2)
        cor_mt_p <- matrix(0, N1, N2)
        cor_mm <- matrix(0, N1, N2)
        cor_mm_p <- matrix(0, N1, N2)
        cor_mm_60 <- matrix(0, N1, N2)
        cor_mm_60_p <- matrix(0, N1, N2)
        
        cor_method  = "pearson"
        # cor_method  = "spearman"
        
        for (i in 1: N1)
        {
            for (j in 1:N2)
            {
                # cor for normal 
                tmp = cor.test(m6a_rn_log2[i, ], en_cmb_all[j, 1:10], method = cor_method)
                cor_mn[i, j] = tmp$estimate
                cor_mn_p[i, j] = tmp$p.value
                
                # cor for tumor 
                tmp = cor.test(m6a_rt_log2[i, ], en_cmb_all[j, 11:61], method = cor_method)
                cor_mt[i, j] = tmp$estimate
                cor_mt_p[i, j] = tmp$p.value
                
                # cor for merge 
                tmp = cor.test(m6a_mm_log2[i, ], en_cmb_all[j, ], method = cor_method)
                cor_mm[i, j] = tmp$estimate
                cor_mm_p[i, j] = tmp$p.value
                
                ## cor or merge without tumor18 (FTO outlier)
                tmp = cor.test(m6a_mm_log2[i, -19], en_cmb_all[j, -19], method = cor_method)
                cor_mm_60[i, j] = tmp$estimate
                cor_mm_60_p[i, j] = tmp$p.value
                
            }
        }
        
        rownames(cor_mn) <- rownames(m6a_rn)
        colnames(cor_mn) <- rownames(en_cmb_all)
        rownames(cor_mt) <- rownames(m6a_rn)
        colnames(cor_mt) <- rownames(en_cmb_all)
        rownames(cor_mm) <- rownames(m6a_rn)
        colnames(cor_mm) <- rownames(en_cmb_all)
        rownames(cor_mm_60) <- rownames(m6a_rn)
        colnames(cor_mm_60) <- rownames(en_cmb_all)
        
        save(cor_mn, cor_mn_p, cor_mt, cor_mt_p, file = paste(cor_method, "_cor_m6a_and_enzymes_mrna.Rdata", sep = "")) 
        save(cor_mm, cor_mm_p, cor_mm_60, cor_mm_60_p, file = paste(cor_method, "_cor_m6a_and_enzymes_mrna_merged_normal_tumor.Rdata", sep = "")) 
    }
    
    #################################################################
    ##  overall patterns for correatlion between the m6A and enzymes
    {
    ## heatmap for correlation with row group color bar
    {
        hm_cor_bar <- function(cor_t, name, colv, w_in, h_in, kcex)
        {
            png(file = name, width = w_in, height = h_in, units = "in", res = 300)
            hmcols<-colorRampPalette(c("blue","white","red"))(256)
            
            # for rowg roup sidebar
            hr <- hclust(dist(cor_t, method="euclidean"), method="complete")     # heatmap.2 clustring based on euclidean distance by default 
            mycl <- cutree(hr, k=5)
            clusterCols <- c("dodgerblue2", "yellow2", "springgreen2", "red", "magenta2") 
            #clusterCols <- rainbow(length(unique(mycl)))
            myClusterSideBar <- clusterCols[mycl]  
            
            # for column group sidebar
            hr_col <- hclust(dist(t(cor_t), method="euclidean"), method="complete")
            mycl_col <- cutree(hr_col, k=2)
            clusterCols_col <- c("orange", "gray")
            #clusterCols_col <- rainbow(length(unique(mycl_col)))
            
            library(plotrix)    # for the finction of color.id  from the hex code 
            clusterCols_col_name <- sapply(clusterCols_col,color.id)
            myClusterSideBar_col <- clusterCols_col[mycl_col]  
            
            heatmap.2(cor_t, Colv= as.dendrogram(hr_col), Rowv=as.dendrogram(hr), dendrogram = "both" , scale="none", na.color = "grey",
                      trace ="none", density.info="none", col=hmcols, RowSideColors= myClusterSideBar, ColSideColors= myClusterSideBar_col,
                      lhei = c(1,5),  lwid = c(1, 4), margins=c(5, 1), labRow = FALSE, labCol = colnames(cor_t), cexCol = 1, offsetCol = 0, 
                      srtCol = 75, key.title = NA, key.xlab = "Correlation", key.par =list(cex = kcex))
            
            dev.off()
            
        }
        
        ## for 8 writer complex and erasers enzymes without filtering: not include RBM15, RBM15B, METTL16, SETD2
        cor_mn_8 <- cor_mn[, c(1:2, 4, 7:9, 11:12)]
        cor_mn_p_8 <- cor_mn_p[, c(1:2, 4, 7:9, 11:12)]
        cor_mt_8 <- cor_mt[, c(1:2, 4, 7:9, 11:12)]
        cor_mt_p_8 <- cor_mt_p[, c(1:2, 4, 7:9, 11:12)]
        
        hm_cor_bar(cor_mn_8, "hm_cor_m6a_8we_enzyme_normal_with_color_bar.png", T, 5, 4, 0.5)
        hm_cor_bar(cor_mt_8, "hm_cor_m6a_8we_enzyme_tumor_with_color_bar.png", T, 5, 4, 0.5)
        
    }
    
    ## re order heatmap for correlation with row group color bar
    {
            hmcols<-colorRampPalette(c("blue","white", "red"))(10)
            cor_t <- cor_mt_8
            idx <- abs(cor_t) > 0.5    
            sum(idx)
            cor_t[idx] <- 0.5
            
            # for rowg roup sidebar
            hr <- hclust(dist(cor_t, method="euclidean"), method="complete")     # heatmap.2 clustring based on euclidean distance by default 
            mycl <- cutree(hr, k=3)
            clusterCols <- c("darkseagreen2", "darkseagreen4", "darkslategray")
            myClusterSideBar <- clusterCols[mycl]  
            
            # for column group sidebar
            hr_col <- hclust(dist(t(cor_t), method="euclidean"), method="complete")
            mycl_col <- cutree(hr_col, k=2)
            clusterCols_col <- c("darkorchid3", "orange2")
            myClusterSideBar_col <- clusterCols_col[mycl_col]  
            
            pdf(file = "hm_cor_m6a_8we_enzyme_tumor_with_color_bar_new.pdf", width = 4, height = 5.5)
            heatmap.2(cor_t, Colv= as.dendrogram(hr_col), Rowv=as.dendrogram(hr), dendrogram = "both" , scale="none", na.color = "grey",
                      trace ="none", density.info="none", col=hmcols, RowSideColors= myClusterSideBar, ColSideColors= myClusterSideBar_col,
                      lhei = c(1,5),  lwid = c(1, 4), margins=c(5, 1), labRow = FALSE, labCol = colnames(cor_t), cexCol = 1, offsetCol = 0, 
                      srtCol = 75, key.title = NA, key.xlab = "Correlation", key.par =list(cex = 0.5))
            dev.off()
            
            ###########
            ## reorder 
            
            cor_t_col_s <- cor_t[, c(1, 2, 5, 4, 6, 3, 8, 7)]
            table(mycl)
            tmp <- table(myClusterSideBar)
            tmp_m  <- cor_t_col_s[hr$order, ]
            
            cor_t_s <- rbind( tmp_m[1:2776, ], tmp_m[4750:8030, ],   tmp_m[2777:4749,  ])      ## original order : from bottom to top [1, 3, 2]
            
            library(heatmap.plus) 
            col_c <- cbind(as.matrix(myClusterSideBar_col), as.matrix(myClusterSideBar_col))
            col_rv <- rev(c(rep("darkslategray", tmp["darkslategray"]), rep("darkseagreen4", tmp["darkseagreen4"]), rep("darkseagreen2", tmp["darkseagreen2"])))
            col_r <- as.matrix(cbind(col_rv, col_rv))
            
            #pdf(file = "hm_cor_m6a_8we_enzyme_tumor_with_color_bar_resort.pdf", width = 3, height = 6)
            png(file = "hm_cor_m6a_8we_enzyme_tumor_with_color_bar_resort.png", width = 3, height = 6, units = "in", res = 600)
            heatmap.plus(cor_t_s, scale = "none", labRow = NA,  Rowv = NA, Colv = NA, col = hmcols, ColSideColors = col_c, RowSideColors = col_r) 
            dev.off()
            
        }
    
    ## for tumor row groups
    {
            table(mycl)
            group1 <- names(mycl)[mycl == 3]
            group2 <- names(mycl)[mycl == 2]
            group3 <- names(mycl)[mycl == 1]
            write.table(group1, file = "cor_enzyme_m6A_tumor_gene_group1.txt", row.names = F, col.names = F, quote = F )
            write.table(group2, file = "cor_enzyme_m6A_tumor_gene_group2.txt", row.names = F, col.names = F, quote = F )
            write.table(group3, file = "cor_enzyme_m6A_tumor_gene_group3.txt", row.names = F, col.names = F, quote = F )
            
            ## run cluster profiler 
            library(clusterProfiler)
            group1_ID <- bitr(group1, 'SYMBOL', "ENTREZID", "org.Hs.eg.db") [, "ENTREZID"]
            group1_KEGG <- enrichKEGG(group1_ID, keyType = "kegg", pvalueCutoff = 0.05, qvalueCutoff = 0.05, pAdjustMethod = "BH",
                                      minGSSize = 5, maxGSSize = 500, organism = "hsa", use_internal_data = FALSE )
            
            write.csv(summary(group1_KEGG), "Group1_KEGG_clusterprofiler.csv", row.names = F)
            
            pdf("Group1_KEGG_clusterprofiler_top5.pdf", width = 5, height = 5)
            dotplot(group1_KEGG, showCategory = 5)
            dev.off()
            
            
            
            group2_ID <- bitr(group2, 'SYMBOL', "ENTREZID", "org.Hs.eg.db") [, "ENTREZID"]
            group2_KEGG <- enrichKEGG(group2_ID, keyType = "kegg", pvalueCutoff = 0.05, qvalueCutoff = 0.05, pAdjustMethod = "BH",
                                      minGSSize = 5, maxGSSize = 500, organism = "hsa", use_internal_data = FALSE )
            write.csv(summary(group2_KEGG), "Group2_KEGG_clusterprofiler.csv", row.names = F)
            
            group3_ID <- bitr(group3, 'SYMBOL', "ENTREZID", "org.Hs.eg.db") [, "ENTREZID"]
            group3_KEGG <- enrichKEGG(group3_ID, keyType = "kegg", pvalueCutoff = 0.05, qvalueCutoff = 0.05, pAdjustMethod = "BH",
                                      minGSSize = 5, maxGSSize = 500, organism = "hsa", use_internal_data = FALSE )
            
            write.csv(summary(group3_KEGG), "Group3_KEGG_clusterprofiler.csv", row.names = F)
            
            
            kegg_top5_3 <- rbind(head(summary(group1_KEGG), 5), head(summary(group2_KEGG), 5), head(summary(group3_KEGG), 5))
            kegg_top2_3 <- rbind(head(summary(group1_KEGG), 2), head(summary(group2_KEGG), 2), head(summary(group3_KEGG), 2))
            
            
            ###  dot plot for top5 KEGG for G1 ,G2 and G3
            {
                library("BoutrosLab.plotting.general")
                
                ## top5 for each group    
                {
                    spot_color <- rep( rev(c("darkseagreen2", "darkseagreen4", "darkslategray")), each = 5)
                    spot_color <- spot_color[-6]        ## only 4 for group2
                    
                    fdr_log <- -log10(kegg_top5_3$p.adjust)
                    gene_cnt <- kegg_top5_3$Count
                    names(gene_cnt) <- kegg_top5_3$Description
                    
                    pdf(file = "Tumor_3Groups_Top5_KEGG.pdf", heigh = 6, width = 2)
                    {
                        spot.size.function <- function(x) { 0.1 + (0.02 * abs(x)); }
                        
                        create.dotmap(
                            
                            gene_cnt,  yaxis.cex = 0.5,  xaxis.cex = 0.5, xlab.cex = 0.5, ylab.cex = 0.7,
                            
                            spot.size.function = spot.size.function,
                            
                            key = list(space = "right", 
                                       
                                       points = list(cex = spot.size.function(seq(0, 150, 30)),
                                                     
                                                     pch = 19),
                                       
                                       # dot labels
                                       text = list(lab = c( "0", "30","60", "90", "120", "150"),
                                                   cex = 0.5, adj = 1.0, fontface = "normal")
                            ),
                            
                            # add borders to points
                            pch = 21,
                            pch.border.col = "white",
                            
                            # add the background
                            bg.data = fdr_log,
                            
                            # add a colourkey
                            colourkey = TRUE,
                            
                            # set colour scheme for background data
                            colour.scheme = c("gray", "black"),
                            
                            #### back groud alpha        !!!!!
                            bg.alpha = 0.5,
                            
                            # make bg colour scheme a discrete colour scheme, with breaks at these places
                            at = seq(1, 5),
                            xaxis.rot = 90
                        )
                    }
                    dev.off()
                }
                
                ## top2 for each group   
                {
                    spot_color <- rep( rev(c("darkseagreen2", "darkseagreen4", "darkslategray")), each = 2)
                    
                    fdr_log <- -log10(kegg_top2_3$p.adjust)
                    gene_cnt <- kegg_top2_3$Count
                    names(gene_cnt) <- kegg_top2_3$Description
                    
                    pdf(file = "Tumor_3Groups_Top2_KEGG.pdf", heigh = 6, width = 2)
                    {
                        spot.size.function <- function(x) { 0.1 + (0.02 * abs(x)); }
                        
                        create.dotmap(
                            
                            gene_cnt,  yaxis.cex = 0.5,  xaxis.cex = 0.5, xlab.cex = 0.5, ylab.cex = 0.7,
                            
                            spot.size.function = spot.size.function,
                            spot.colour.function = spot_color,
                            
                            key = list(space = "right", 
                                       
                                       points = list(cex = spot.size.function(seq(0, 150, 50)),
                                                     pch = 19),
                                       
                                       # dot labels
                                       text = list(lab = c( "0", "50", "100", "150"),
                                                   cex = 0.5, adj = 1.0, fontface = "normal")
                            ),
                            
                            # add borders to points
                            pch = 21,
                            pch.border.col = "white",
                            
                            # add the background
                            bg.data = fdr_log,
                            
                            # add a colourkey
                            colourkey = TRUE,
                            
                            # set colour scheme for background data
                            colour.scheme = c("gray", "black"),
                            
                            #### back groud alpha        !!!!!
                            bg.alpha = 0.5,
                            
                            # make bg colour scheme a discrete colour scheme, with breaks at these places
                            at = seq(1, 5),
                            xaxis.rot = 90
                        )
                    }
                    dev.off()
                }
                
                
                
                
            } 
            
            ###########################
            g_l <- c(group1, group2, group3)
            L1 <- length(group1)
            L2 <- length(group2)
            L3 <- length(group3)
            L <- length(g_l)
            g_lab <- c(rep("Group1", L1), rep("Group2", L2), rep("Group3", L3)) 
            
            
            ## box plot for mRNA
            idx_m <-   match(g_l, rownames(input_rpkm_t))
            mrna_t <- log2(rowMeans(input_rpkm_t[idx_m, ]) + 1) 
            dat <- data.frame(mrna_t, g_lab)
            
            wilcox.test(mrna_t[g_lab == "Group1"], mrna_t[g_lab == "Group2"], alternative = "greater")
            wilcox.test(mrna_t[g_lab == "Group1"], mrna_t[g_lab == "Group3"], alternative = "greater")
            wilcox.test(mrna_t[g_lab == "Group2"], mrna_t[g_lab == "Group3"], alternative = "greater")
            
            g <- ggplot(dat, aes(x = g_lab, y = mrna_t,  fill = g_lab)) + geom_boxplot() + theme_classic() + scale_y_continuous(limits = c(0, 15), breaks = seq(0, 15, 3))
            g <- g + scale_fill_manual(values = c("darkslategray", "darkseagreen4", "darkseagreen2")) + theme(legend.position = "none")
            ggsave("Tumor_only_3group_genes_mRNA.pdf", width = 2, height = 3, units = "in")
            
            
            ### boxplot for m6A levle 
            idx_m <-   match(g_l, rownames(m6a_rt))
            m6a_t <- log2(rowMeans(m6a_rt[idx_m, ], na.rm = T)) 
            dat <- data.frame(m6a_t, g_lab)
            
            wilcox.test(m6a_t[g_lab == "Group1"], m6a_t[g_lab == "Group2"], alternative = "less")
            wilcox.test(m6a_t[g_lab == "Group1"], m6a_t[g_lab == "Group3"], alternative = "less")
            wilcox.test(m6a_t[g_lab == "Group2"], m6a_t[g_lab == "Group3"], alternative = "less")
            
            g <- ggplot(dat, aes(x = g_lab, y = m6a_t,  fill = g_lab)) + geom_boxplot() + theme_classic() 
            g <- g + scale_fill_manual(values = c("darkslategray", "darkseagreen4", "darkseagreen2")) + theme(legend.position = "none")
            ggsave("Tumor_only_3group_genes_m6A.pdf", width = 2, height = 3, units = "in")
            
        } 
           
    }
   
     ################################
    ## pusblished PAR-CLIP targets
    {
        M3 <- read.table("/Users/Yong/Yong/m6a_profiling/4_enzymes/Published_data/METTL3_14_WTAP/METTL3_PAR_CLIP_targets.txt")
        M14 <- read.table("/Users/Yong/Yong/m6a_profiling/4_enzymes/Published_data/METTL3_14_WTAP/METTL14_PAR_CLIP_targets.txt")
        WTAP <- read.table("/Users/Yong/Yong/m6a_profiling/4_enzymes/Published_data/METTL3_14_WTAP/WTAP_PAR_CLIP_targets.txt")
        A5 <- read.table("/Users/Yong/Yong/m6a_profiling/4_enzymes/Published_data/ALKBH5/ALKBH5_PAR_CLIP_targets.txt")
        FTO <- read.table("/Users/Yong/Yong/m6a_profiling/4_enzymes/Published_data/FTO/FTO_PAR_CLIP_targets.txt")
        
        ## gene match function 
        g_match_cnt <- function(g_s, g_r)
        {
            idx <- match(g_s, g_r)
            idx_1 <- !is.na(idx)
            g_m <- g_s[idx_1]
            return <- length(g_m)
        }
        
        que_list <- list(group1, group2, group3)
        #que_list <- list(gene_nn, gene_np, gene_tn, gene_tp, gene_nt_n, gene_nt_p)
        ref_list <- list(M3$V1, M14$V1, WTAP$V1, A5$V1, FTO$V1)
        
        ## table 
        L1 <- length(que_list)
        L2 <- L1*length(ref_list)
        target <- matrix(0, 3, L2)     # targets vs non targets
        
        colnames(target) <- rep(c("Group 1", "Group 2", "Group 3"), L2/L1)
        rownames(target) <- c("target", "non_target", "target_pe")
        for(i in seq(1, L2, L1))
        {
            for(k in 1:L1)
            {
                target[1, (i-1) + k] <- g_match_cnt(que_list[[k]], ref_list[[(i + L1 - 1)/L1]])
                target[2, (i-1) + k] <- length(que_list[[k]]) - target[1, (i-1) + k]
                target[3, (i-1) + k] <- target[1, (i-1) + k] / length(que_list[[k]])
            }
            
        }
        
        write.csv(target, file = "cor_enzyme_m6A_tumor_gene_3group_WE_targets.csv")
        
        # test red group and the rest groups
        p_val <- vector()
        p_val_p <- vector()
        k = 1
        for(i in seq(1, L2, L1))
        {
            for (j in 0:2)
            {
                ## chi-square test 
                tmp <- chisq.test(target[1:2, c( i, i+j)])
                p_val[k] <- tmp$p.value
                
                
                ## propotion test
                tmp <- prop.test(target[1, c(i, i +j)], colSums(target[1:2, c(i, i +j)]))
                p_val_p[k] <- tmp$p.value
                
                k = k + 1
            }
        }
        
        p_val 
        ## p_val and p_val_p are the same for the 2 by 2 contingency table
        
        ### grouped barplot
        pen_m <- matrix(target[3, ], 3, 5)
        rownames(pen_m) <- c("Group 1", "Group 2", "Group 3")
        colnames(pen_m) <- c("METTL3", "METTL14", "WTAP", "ALKBH5", "FTO")
        col_bar <-  c("darkslategray", "darkseagreen4", "darkseagreen2")
        
        pdf(file = "cor_enzyme_m6A_tumor_gene_3group_WE_targets.pdf", width = 4, height = 3)
        par(mar = c(3, 3, 0.5, 0.5), mgp = c(2, 0.75, 0))
        barplot(pen_m, beside = T, col = col_bar, ylim = c(0, 0.6),  ylab = "PAR-CLIP targets (%)")
        dev.off()
        
    }
}
