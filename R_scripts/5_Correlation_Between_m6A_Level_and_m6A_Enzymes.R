################################################################
## This script implements:
## 1. m6A related enzymes mRNA level comparison
## 2. correlation between m6A levels and related enzymes per gene
## 3. characteristics of correlation based subgroups: mRNA, m6A and regulators targets
## 4. characteristics of correlation based subgroups: peak density and gene features 
#################################################################

rm(list = ls())
setwd("./data")

############################################
## m6A related enzymes mRNA level comparison
############################################
{
    ###############################
    ## loading and extracting data 
    {
    m6a_en <- read.table("./1_enzymes/m6a_enzyme.txt", header = T)      # enzyme list  
    load("./paired_63_count_ercc_deseq.Rdata") 
    
    ## remove outliers tumor 45 and tumor9
    idx_s <- match(c("tumor9", "tumor45"), colnames(ip_rpkm))

    ip_rpkm <- ip_rpkm[, -idx_s]
    input_rpkm <- input_rpkm[, -idx_s]
    
    
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
    
    #save(en_cmb_we, en_cmb_re, en_cmb_all, file = "en_cmb_we_re_log2_rpkm.Rdata")
    
    #### GDC TCGA LUAD data  : Normal: 59
    if(FALSE){
    load("/Users/yong/OneDrive - UHN/Projects/Epitranscriptomics/LUAD_m6A/1_mapping_quant/5_RNA-seq/TCGA/GDC_LUAD/tcga_htseq_fpkm_log2.Rdata")        # FPKM from TCGA
    tcga_g <- rownames(tcga)
    idx_ki <- match("VIRMA", tcga_g)         # KIAA1429 uses name VIRMA
    tcga_g[idx_ki] <- "KIAA1429"
    idx_we <- match(rownames(en_cmb_we), tcga_g)
    idx_re <- match(rownames(en_cmb_re), tcga_g)
    tcga <- tcga[c(idx_we, idx_re), ]
    save(tcga, tcga_type, file = "TCGA_LUAD_ezymes_expr.Rdata")
    }
    
    load("TCGA_LUAD_ezymes_expr.Rdata")
    
    tcga_g <- rownames(tcga)
    idx_ki <- match("VIRMA", tcga_g)         # KIAA1429 uses name VIRMA
    tcga_g[idx_ki] <- "KIAA1429"
    idx_we <- match(rownames(en_cmb_we), tcga_g)
    idx_re <- match(rownames(en_cmb_re), tcga_g)

    en_cmb_we_t <- tcga[idx_we, ]
    en_cmb_re_t <- tcga[idx_re, ] 
    
    idx_1 <- tcga_type == "normal"
    idx_2 <- tcga_type == "tumor"
    
    ## only includes 8we
     en_cmb_8we <- en_cmb_we[-c(3, 5, 6), ]
     en_cmb_8we_t <- en_cmb_we_t[-c(3, 5, 6), ]
     }
     
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
     outname <- c("../Results/m6a_8writer_eraser_mRNA_TCGA_boxplot_h.pdf", "../Results/m6a_8writer_eraser_mRNA_boxplot_h.pdf")
     sample_type  <- list(tcga_type_s, s_type)
     
     for(i in 1:2)
     {
     ## transfer to dataform
     expr_we <- expr_we_list[[i]]
     x <- as.vector(t(expr_we)) 
     y1 <- rep(rownames(expr_we), each = ncol(expr_we))
     y2 <- paste(y1, sample_type[[i]], sep = "_")
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
     
     
     ## test 
     {
       enzymes.test <- function(en_expr, sample_type)
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
         
         return(pval)
       }
       
       enzymes.test(en_cmb_we, s_type)
       enzymes.test(en_cmb_we_t, tcga_type)
     }
     
    }
}     
 
##############################################################
## correlation between m6A levels and related enzymes per gene 
##############################################################
{
  #################################################
  ## load m6A levels and calcuated the correlations
  {
  load("m6A_level_4gene_passed_filtering.Rdata")
  
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
    }
  }
  
  rownames(cor_mn) <- rownames(m6a_rn)
  colnames(cor_mn) <- rownames(en_cmb_all)
  rownames(cor_mt) <- rownames(m6a_rn)
  colnames(cor_mt) <- rownames(en_cmb_all)
  
  save(cor_mn, cor_mn_p, cor_mt, cor_mt_p, file = "cor_m6a_and_enzymes_mrna.Rdata") 
  }
  
  #################################################
  ## correaltion patterns
  {
    library(gplots)
    ## for 8 writer complex and erasers enzymes without filtering: not include RBM15, RBM15B, METTL16, SETD2 in tumor
    cor_mt_8 <- cor_mt[, c(1:2, 4, 7:9, 11:12)]
    cor_mt_p_8 <- cor_mt_p[, c(1:2, 4, 7:9, 11:12)]
    
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

    pdf(file = "../Results/hm_cor_m6a_8we_enzyme_tumor_with_color_bar_new.pdf", width = 4, height = 5.5)
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
    
    pdf(file = "hm_cor_m6a_8we_enzyme_tumor_with_color_bar_resort.pdf", width = 3, height = 6)
    #png(file = "hm_cor_m6a_8we_enzyme_tumor_with_color_bar_resort.png", width = 3, height = 6, units = "in", res = 600)
    heatmap.plus(cor_t_s, scale = "none", labRow = NA,  Rowv = NA, Colv = NA, col = hmcols, ColSideColors = col_c, RowSideColors = col_r) 
    dev.off()
    
  }
  
  #################################
  ## subgroups based on correlation 
  {
    table(mycl)
    group1 <- names(mycl)[mycl == 3]
    group2 <- names(mycl)[mycl == 2]
    group3 <- names(mycl)[mycl == 1]
    
    write.table(group1, file = "../Results/cor_enzyme_m6A_tumor_gene_group1.txt", row.names = F, col.names = F, quote = F )
    write.table(group2, file = "../Results/cor_enzyme_m6A_tumor_gene_group2.txt", row.names = F, col.names = F, quote = F )
    write.table(group3, file = "../Results/cor_enzyme_m6A_tumor_gene_group3.txt", row.names = F, col.names = F, quote = F )
    
    ## run cluster profiler
    ## result might be slight different due the clusterprofiling enrichment anaysis  !!!
    if(FALSE){
    library(clusterProfiler)
    group1_ID <- bitr(group1, 'SYMBOL', "ENTREZID", "org.Hs.eg.db") [, "ENTREZID"]
    group1_KEGG <- enrichKEGG(group1_ID, keyType = "kegg", pvalueCutoff = 0.05, qvalueCutoff = 0.05, pAdjustMethod = "BH",
                              minGSSize = 5, maxGSSize = 500, organism = "hsa", use_internal_data = FALSE )
    
   
    group2_ID <- bitr(group2, 'SYMBOL', "ENTREZID", "org.Hs.eg.db") [, "ENTREZID"]
    group2_KEGG <- enrichKEGG(group2_ID, keyType = "kegg", pvalueCutoff = 0.05, qvalueCutoff = 0.05, pAdjustMethod = "BH",
                              minGSSize = 5, maxGSSize = 500, organism = "hsa", use_internal_data = FALSE )

    group3_ID <- bitr(group3, 'SYMBOL', "ENTREZID", "org.Hs.eg.db") [, "ENTREZID"]
    group3_KEGG <- enrichKEGG(group3_ID, keyType = "kegg", pvalueCutoff = 0.05, qvalueCutoff = 0.05, pAdjustMethod = "BH",
                              minGSSize = 5, maxGSSize = 500, organism = "hsa", use_internal_data = FALSE )
    
    #kegg_top5_3 <- rbind(head(summary(group1_KEGG), 5), head(summary(group2_KEGG), 5), head(summary(group3_KEGG), 5))
    kegg_top2_3 <- rbind(head(summary(group1_KEGG), 2), head(summary(group2_KEGG), 2), head(summary(group3_KEGG), 2))
    
    ###  dot plot for top5 KEGG for G1 ,G2 and G3
    {
      library("BoutrosLab.plotting.general")
      
      ## top2 for each group   
      {
        spot_color <- rep( rev(c("darkseagreen2", "darkseagreen4", "darkslategray")), each = 2)
        
        fdr_log <- -log10(kegg_top2_3$p.adjust)
        gene_cnt <- kegg_top2_3$Count
        names(gene_cnt) <- kegg_top2_3$Description
        
        pdf(file = "../Results/Tumor_3Groups_Top2_KEGG.pdf", heigh = 6, width = 2)
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
    }
  }
}

##############################################################
## characteristics of correlation based subgroups: mRNA, m6A and regulators targets
##############################################################
{
  library(ggplot2)
  ##############################
  ## RNA and m6A levels for G1-3
  {
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
    ggsave("../Results/Tumor_only_3group_genes_mRNA.pdf", width = 2, height = 3, units = "in")

    
    ### boxplot for m6A levle 
    idx_m <-   match(g_l, rownames(m6a_rt))
    m6a_t <- log2(rowMeans(m6a_rt[idx_m, ], na.rm = T)) 
    dat <- data.frame(m6a_t, g_lab)
    
    wilcox.test(m6a_t[g_lab == "Group1"], m6a_t[g_lab == "Group2"], alternative = "less")
    wilcox.test(m6a_t[g_lab == "Group1"], m6a_t[g_lab == "Group3"], alternative = "less")
    wilcox.test(m6a_t[g_lab == "Group2"], m6a_t[g_lab == "Group3"], alternative = "less")
    
    g <- ggplot(dat, aes(x = g_lab, y = m6a_t,  fill = g_lab)) + geom_boxplot() + theme_classic() 
    g <- g + scale_fill_manual(values = c("darkslategray", "darkseagreen4", "darkseagreen2")) + theme(legend.position = "none")
    ggsave("../Results/Tumor_only_3group_genes_m6A.pdf", width = 2, height = 3, units = "in")
    
  }
  
  ##############################
  ## write and eraser targets 
  {
    M3 <- read.table("./1_enzymes/METTL3_PAR_CLIP_targets.txt")
    M14 <- read.table("./1_enzymes/METTL14_PAR_CLIP_targets.txt")
    WTAP <- read.table("./1_enzymes/WTAP_PAR_CLIP_targets.txt")
    A5 <- read.table("./1_enzymes/ALKBH5_PAR_CLIP_targets.txt")
    FTO <- read.table("./1_enzymes/FTO_PAR_CLIP_targets.txt")
    
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
    
    #write.csv(target, file = "cor_enzyme_m6A_tumor_gene_3group_WE_targets.csv")
    
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
    
    pdf(file = "../Results/cor_enzyme_m6A_tumor_gene_3group_WE_targets.pdf", width = 6, height = 3)
    par(mar = c(3, 3, 0.5, 0.5), mgp = c(2, 0.75, 0))
    barplot(pen_m, beside = T, col = col_bar, ylim = c(0, 0.6),  ylab = "PAR-CLIP targets (%)")
    dev.off()
  }
  
  ##############################
  ## reader targets
  {
    ## YTH
    YTHDF1 <- read.table("./1_enzymes/YTHDF1_targets.txt", as.is = T)           # promote tranlation
    YTHDF2 <- read.table("./1_enzymes/YTHDF2_targets.txt", as.is = T)           # RNA decay 
    YTHDF3 <- read.table("./1_enzymes/YTHDF3_targets.txt", as.is = T)           # promote tranlation
    YTHDC2 <- read.table("./1_enzymes/YTHDC2_target.txt", as.is = T)
    
    ## GF2BP1-3
    IGF2BP1 <- read.table("./1_enzymes/IGF2BP1_targets.txt", as.is = T)
    IGF2BP2 <- read.table("./1_enzymes/IGF2BP2_targets.txt", as.is = T)
    IGF2BP3 <- read.table("./1_enzymes/IGF2BP3_targets.txt", as.is = T)
    
    que_list <- list(group1, group2, group3)
    ref_list <- list(YTHDF1$V1, YTHDF2$V1, YTHDF3$V1, YTHDC2$V1, IGF2BP1$V1, IGF2BP2$V1, IGF2BP3$V1)
    
    ref_len <- vector()
    for (i in 1:7){
      ref_len[i] <- length(ref_list[[i]])
    }
    
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
    
    #write.csv(target, file = "cor_enzyme_m6A_tumor_gene_3group_Readers_targets.csv")
    
    # test red group and the rest groups
    ## G1 vs G2 or G3
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
    if (FALSE){
      ## G2 vs G1, G3
      ## G3 vs G1, G2
      
      p_val <- vector()
      p_val_p <- vector()
      k = 1
      for(i in seq(3, L2, L1))
      {
        for (j in 0:2)
        {
          ## chi-square test 
          tmp <- chisq.test(target[1:2, c( i, i+j-2)])
          p_val[k] <- tmp$p.value
          
          
          ## propotion test
          tmp <- prop.test(target[1, c(i, i + j - 2)], colSums(target[1:2, c(i, i +j - 2)]))
          p_val_p[k] <- tmp$p.value
          
          k = k + 1
        }
      }
      
      p_val 
      
      
    }
    
    
    ### grouped barplot
    pen_m <- matrix(target[3, ], 3, 7)
    rownames(pen_m) <- c("Group 1", "Group 2", "Group 3")
    colnames(pen_m) <- c("YTHDF1", "YTHDF2", "YTHDF3", "YTHDC2", "IGF2BP1", "IGF2BP2", "IGF2BP3")
    col_bar <-  c("darkslategray", "darkseagreen4", "darkseagreen2")
    
    pdf(file = "../Results/cor_enzyme_m6A_tumor_gene_3group_Reader_targets.pdf", width = 7, height = 4)
    par(mar = c(3, 3, 0.5, 0.5), mgp = c(2, 0.75, 0))
    barplot(pen_m, beside = T, col = col_bar, ylim = c(0, 0.4),  ylab = "PAR-CLIP targets (%)")
    dev.off()
    p_val 
  }
  
   
  
}

##############################################################
## characteristics of correlation based subgroups: peak density and gene features 
##############################################################
{
  ############################
  ## peaks meta distributions
  {
    comLength = c(0.136, 0.459, 0.405)
    weight <- comLength/sum(comLength)
    sep <- cumsum(weight)
    
    ## function
    metaPeakDistributions <-  function(ct, comLength, color, name) 
    {
      library("ggplot2")
      source("../R_scripts/functions/summarySE.R")     # barplot with error bar
      
      ## for lincRNA
      ct$weight <- ct$count
      ct1 <- ct[ct$category == "mRNA", ]
      ct2 <- ct[ct$category == "lncRNA", ]
      pos = Feature = weight = NULL
      id1 <- which(match(ct1$comp, c("Front", "Back")) > 0)
      ct1 <- ct1[-id1, ]
      id2 <- which(match(ct2$comp, c("Front", "Back")) > 0)
      ct2 <- ct2[-id2, ]
      featureSet <- as.character(unique(ct$Feature))
      for (i in 1:length(featureSet)) {
        id <- (ct1$Feature == featureSet[i])
        ct1$weight[id] <- ct1$weight[id]/sum(ct1$weight[id])
        id <- (ct2$Feature == featureSet[i])
        ct2$weight[id] <- ct2$weight[id]/sum(ct2$weight[id])
      }
      
      p2 <- ggplot(ct2, aes(x = pos, group = Feature, weight = weight)) + geom_density()
      p2_out <- ggplot_build(p2)
      
      lnc <- p2_out$data[[1]]
      lnc_se <- summarySE(data = lnc, measurevar = "y", groupvars="x")
      
      ## plot mean with shadow for sd
      g2 <-  ggplot(lnc_se, aes(x = x, y = y, fill = color)) + geom_ribbon(aes( ymin = y - sd, ymax = y + sd), alpha = 1) + geom_line(col = "black", lty = 2)
      g2 <- g2 + theme_bw() + theme(axis.line = element_line(colour = "black"),panel.grid.major = element_blank(),panel.grid.minor = element_blank(),panel.border =  element_blank(),  panel.background = element_blank()) 
      g2 <- g2 + annotate("rect", xmin = 0, xmax = 1, ymin = -0.05, ymax = -0.00, alpha = 0.99, colour = "black") 
      g2 <-  g2 + scale_fill_manual(values = color)   + theme(legend.position = "none")
      
      ggsave(paste(name, "_peaks_distribution_in_lncRNA.pdf", sep = ""), width = 3, height = 2, units = "in")
      
      
      ## for mRNA
      adjust = 1
      weight <- comLength/sum(comLength)
      names(weight) <- c("5'UTR", "CDS", "3'UTR")
      cds_id <- which(ct1$comp == "CDS")
      utr3_id <- which(ct1$comp == "UTR3")
      utr5_id <- which(ct1$comp == "UTR5")
      ct1$count[utr5_id] <- ct1$count[utr5_id] * weight["5'UTR"]
      ct1$count[cds_id] <- ct1$count[cds_id] * weight["CDS"]
      ct1$count[utr3_id] <- ct1$count[utr3_id] * weight["3'UTR"]
      featureSet <- as.character(unique(ct$Feature))
      for (i in 1:length(featureSet)) {
        id <- (ct1$Feature == featureSet[i])
        ct1$weight[id] <- ct1$count[id]/sum(ct1$count[id])
      }
      x <- cumsum(weight)
      ct1$pos[utr5_id] <- ct1$pos[utr5_id] * weight["5'UTR"] +  0
      ct1$pos[cds_id] <- ct1$pos[cds_id] * weight["CDS"] + x[1]
      ct1$pos[utr3_id] <- ct1$pos[utr3_id] * weight["3'UTR"] + x[2]
      
      p1 <- ggplot(ct1, aes(x = pos, group = Feature, weight = weight)) + geom_density()
      p1_out <- ggplot_build(p1)
      
      pc <- p1_out$data[[1]]
      pc_se <- summarySE(data = pc, measurevar = "y", groupvars="x")
      
      g1 <-  ggplot(pc_se, aes(x = x, y = y, fill = color)) + geom_ribbon(aes( ymin = y - sd, ymax = y + sd), alpha = 1) + geom_line(col = "black", lty = 2)
      g1 <- g1 + theme_bw() + theme(axis.line = element_line(colour = "black"),panel.grid.major = element_blank(),panel.grid.minor = element_blank(),panel.border =  element_blank(),  panel.background = element_blank()) 
      g1 <- g1 +  geom_vline(xintercept = x[2], linetype = "dotted")
      g1 <- g1 + annotate("rect", xmin = 0, xmax = x[1], ymin = -0.12, ymax = -0.08, alpha = 0.99, colour = "black") + annotate("rect", xmin = x[2], xmax = 1, ymin = -0.12, ymax = -0.08, alpha = 0.99, colour = "black") + annotate("rect", xmin = x[1], xmax = x[2], ymin = -0.16, ymax = -0.04, alpha = 0.2, colour = "black")
      g1 <- g1 + scale_fill_manual(values = color)   + theme(legend.position = "none")
      ggsave(paste(name, "_peaks_distribution_in_PC.pdf", sep = ""), width = 3, height = 2, units = "in")
      
      se_out <- list(pc_se, lnc_se)
      return(se_out)
    }
  
    ## loading processed data
    G1 <- readRDS("G1_peaks_distribution.RDS")
    G2 <- readRDS("G2_peaks_distribution.RDS")
    G3 <- readRDS("G3_peaks_distribution.RDS")
    
    col_G1 = "darkslategray"
    col_G2 = "darkseagreen4"
    col_G3 = "darkseagreen2"
    
    G1_se <- metaPeakDistributions(G1, comLength, col_G1, "Tumor_G1")
    G2_se <- metaPeakDistributions(G2, comLength, col_G2, "Tumor_G2")
    G3_se <- metaPeakDistributions(G3, comLength, col_G3, "Tumor_G3")
    
    ## merge plot for PC
    g <-  ggplot() + geom_line(data = G1_se[[1]], aes(x = x, y = y), col = col_G1)
    g <- g + geom_line(data = G2_se[[1]], aes(x = x, y = y), col = col_G2)
    g <- g + geom_line(data = G3_se[[1]], aes(x = x, y = y), col = col_G3)
    g <- g + theme_bw() + theme(axis.line = element_line(colour = "black"),panel.grid.major = element_blank(),panel.grid.minor = element_blank(),panel.border =  element_blank(),  panel.background = element_blank()) 
    g <- g +  geom_vline( xintercept = sep[2], linetype = "dotted")
    g <- g + annotate("rect", xmin = 0, xmax = sep[1], ymin = -0.12, ymax = -0.08, alpha = 0.99, colour = "black") + annotate("rect", xmin = sep[2], xmax = 1, ymin = -0.12, ymax = -0.08, alpha = 0.99, colour = "black") + annotate("rect", xmin = sep[1], xmax = sep[2], ymin = -0.16, ymax = -0.04, alpha = 0.2, colour = "black")
    ggsave(paste("../Results/Tumor_3Groups_peaks_distribution_in_PC.pdf", sep = ""), width = 3, height = 2, units = "in")
    
    # KS test
    ks.test(G1_se[[1]]$y, G2_se[[1]]$y)
    ks.test(G1_se[[1]]$y, G3_se[[1]]$y)
    ks.test(G2_se[[1]]$y, G3_se[[1]]$y)
    
    ## merge plot for LncRNA
    g <-  ggplot() + geom_line(data = G1_se[[2]], aes(x = x, y = y), col = col_G1)
    g <- g + geom_line(data = G2_se[[2]], aes(x = x, y = y), col = col_G2)
    g <- g + geom_line(data = G3_se[[2]], aes(x = x, y = y), col = col_G3)
    g <- g + theme_bw() + theme(axis.line = element_line(colour = "black"),panel.grid.major = element_blank(),panel.grid.minor = element_blank(),panel.border =  element_blank(),  panel.background = element_blank()) 
    g <- g +  annotate("rect", xmin = 0, xmax = 1, ymin = -0.05, ymax = -0.00, alpha = 0.99, colour = "black") 
    #ggsave(paste("Tumor_3Groups_peaks_distribution_in_lincRNA.pdf", sep = ""), width = 3, height = 2, units = "in")
    
    # KS test
    ks.test(G1_se[[2]]$y, G2_se[[2]]$y)
    ks.test(G1_se[[2]]$y, G3_se[[2]]$y)
    ks.test(G2_se[[2]]$y, G3_se[[2]]$y)
    
    
    ## adding ZNF gens
    ZNF <- readRDS("ZNF_peaks_distribution.RDS")
    col_ZNF = "red"
    
    ZNF_se <- metaPeakDistributions(ZNF, comLength, col_ZNF, "Tumor_ZNF")
    
    
    ## merge plot for PC
    g <-  ggplot() + geom_line(data = G1_se[[1]], aes(x = x, y = y), col = col_G1)
    g <- g + geom_line(data = G2_se[[1]], aes(x = x, y = y), col = col_G2)
    g <- g + geom_line(data = G3_se[[1]], aes(x = x, y = y), col = col_G3)
    g <- g + geom_line(data = ZNF_se[[1]], aes(x = x, y = y), col = col_ZNF)
    g <- g + theme_bw() + theme(axis.line = element_line(colour = "black"),panel.grid.major = element_blank(),panel.grid.minor = element_blank(),panel.border =  element_blank(),  panel.background = element_blank()) 
    g <- g +  geom_vline( xintercept = sep[2], linetype = "dotted")
    g <- g + annotate("rect", xmin = 0, xmax = sep[1], ymin = -0.12, ymax = -0.08, alpha = 0.99, colour = "black") + annotate("rect", xmin = sep[2], xmax = 1, ymin = -0.12, ymax = -0.08, alpha = 0.99, colour = "black") + annotate("rect", xmin = sep[1], xmax = sep[2], ymin = -0.16, ymax = -0.04, alpha = 0.2, colour = "black")
    ggsave(paste("../Results/Tumor_3Groups_peaks_distribution_withZNF_in_PC.pdf", sep = ""), width = 3, height = 2, units = "in")
    
    
    
  }
  
  #############################
  ## genomic features for G1-G3
  {
    gene_f <- read.table("GENCODE_V25_12_Transcript_with_most_exon_4gene_PC_features.txt", header = T, as.is = T)
    
    idx_g1 <- match(rownames(gene_f), group1)
    idx_g2 <- match(rownames(gene_f), group2)
    idx_g3 <- match(rownames(gene_f), group3)

    g_info <- c(rep("G1", sum(!is.na(idx_g1))), rep("G2", sum(!is.na(idx_g2))), rep("G3", sum(!is.na(idx_g3))) )
    
    gene_f3g <- rbind(gene_f[!is.na(idx_g1), ], gene_f[!is.na(idx_g2), ], gene_f[!is.na(idx_g3), ]) 
    gene_f3g <- data.frame(cbind(gene_f3g, g_info))
    
    ## exon numbers
    g <- ggplot(gene_f3g, aes(x = g_info, y = log2(exon_num),  fill = g_info)) + geom_boxplot() + theme_classic()    ##+ scale_y_continuous(limits = c(0, 15), breaks = seq(0, 15, 3))
    g <- g + scale_fill_manual(values = c("darkslategray", "darkseagreen4", "darkseagreen2")) + theme(legend.position = "none")
    ggsave("../Results/Tumor_only_3group_genes_exon_numbers.pdf", width = 2, height = 3, units = "in")
    wilcox.test(gene_f3g$exon_num[idx1], gene_f3g$exon_num[idx2], alternative = "greater")
    wilcox.test(gene_f3g$exon_num[idx2], gene_f3g$exon_num[idx3], alternative = "greater")
    wilcox.test(gene_f3g$exon_num[idx1], gene_f3g$exon_num[idx3], alternative = "greater")
    
    ## last CDS
    g <- ggplot(gene_f3g, aes(x = g_info, y = log2(last_CDS_len),  fill = g_info)) + geom_boxplot() + theme_classic()    ##+ scale_y_continuous(limits = c(0, 15), breaks = seq(0, 15, 3))
    g <- g + scale_fill_manual(values = c("darkslategray", "darkseagreen4", "darkseagreen2")) + theme(legend.position = "none")
    ggsave("../Results/Tumor_only_3group_genes_last_CDS_len.pdf", width = 2, height = 3, units = "in")
    wilcox.test(gene_f3g$last_CDS_len[idx1], gene_f3g$last_CDS_len[idx2], alternative = "less")
    wilcox.test(gene_f3g$last_CDS_len[idx2], gene_f3g$last_CDS_len[idx3], alternative = "less")
    wilcox.test(gene_f3g$last_CDS_len[idx1], gene_f3g$last_CDS_len[idx3], alternative = "less")
    
    
    ## 5 UTR CDS and #UTR
    group <- c(as.character(gene_f3g$g_info), as.character(gene_f3g$g_info), as.character(gene_f3g$g_info))
    length <- c(log2(gene_f3g$X5UTR_len), log2(gene_f3g$CDS_len),  log2(gene_f3g$X3UTR_len))
    features <- rep(c("A_5UTR","B_CDS","C_3UTR"), each = nrow(gene_f3g))
    dat <- data.frame(group, length, features)
    
    g <- ggplot(dat, aes(x = group, y = length, col = group)) + geom_boxplot() 
    g <- g + scale_color_manual( values = c("darkslategray", "darkseagreen4", "darkseagreen2")) + theme_bw()  + theme(legend.position = "top") 
    g <- g + facet_grid(cols = vars(features)) + theme(legend.position = "none")
    ggsave("../Results/5UTR_CDS_3UTR_len_for_3groups_genes_facet_plot.pdf", width = 5, height = 3, units = "in")

    ## comparisons
    idx1 <- gene_f3g$g_info == "G1"
    idx2 <- gene_f3g$g_info == "G2"
    idx3 <- gene_f3g$g_info == "G3"
    
    wilcox.test(gene_f3g$X5UTR_len[idx1], gene_f3g$X5UTR_len[idx2])
    wilcox.test(gene_f3g$X5UTR_len[idx2], gene_f3g$X5UTR_len[idx3])
    wilcox.test(gene_f3g$X5UTR_len[idx1], gene_f3g$X5UTR_len[idx3])
    
    wilcox.test(gene_f3g$CDS_len[idx1], gene_f3g$CDS_len[idx2])
    wilcox.test(gene_f3g$CDS_len[idx2], gene_f3g$CDS_len[idx3])
    wilcox.test(gene_f3g$CDS_len[idx1], gene_f3g$CDS_len[idx3])
    
    wilcox.test(gene_f3g$X3UTR_len[idx1], gene_f3g$X3UTR_len[idx2])
    wilcox.test(gene_f3g$X3UTR_len[idx2], gene_f3g$X3UTR_len[idx3])
    wilcox.test(gene_f3g$X3UTR_len[idx1], gene_f3g$X3UTR_len[idx3])
    
}
  
  ####################################################
  ## diffrene grouop of genes' features with ZNF_G2G3
  ## g2 and g3 without ZNF
  {
    znf <- read.table("ZNF_genes_G2G3.txt", as.is = T)
    
    idx1 <- match(group2, znf$V1)
    g2_f <- group2[is.na(idx1)] 
    idx2 <- match(group3, znf$V1)
    g3_f <- group3[is.na(idx2)] 
    
    idx_g1 <- match(rownames(gene_f), group1)
    idx_g2 <- match(rownames(gene_f), g2_f)
    idx_g3 <- match(rownames(gene_f), g3_f)
    idx_znf <- match(rownames(gene_f), znf$V1)
    
    g_info <- c( rep("G2_no_ZNF", sum(!is.na(idx_g2))), rep("G3_no_ZNF", sum(!is.na(idx_g3))), rep("ZNF_G23", sum(!is.na(idx_znf))))
    gene_f3g <- rbind( gene_f[!is.na(idx_g2), ], gene_f[!is.na(idx_g3), ], gene_f[!is.na(idx_znf), ]) 
    
    gene_f3g <- data.frame(cbind(gene_f3g, g_info))
    #write.csv(gene_f3g, "Tumor_only_3group_genes_CDS_len_with_ZNFG23_noG1.csv")
    
    idx2 <- gene_f3g$g_info == "G2_no_ZNF"
    idx3 <- gene_f3g$g_info == "G3_no_ZNF"
    idx4 <- gene_f3g$g_info == "ZNF_G23"
    
    wilcox.test(gene_f3g$CDS_len[idx4], gene_f3g$CDS_len[idx3])
    wilcox.test(gene_f3g$CDS_len[idx4], gene_f3g$CDS_len[idx2])
    
    g <- ggplot(gene_f3g, aes(x = g_info, y = log2(CDS_len),  fill = g_info)) + geom_boxplot() + theme_classic()    ##+ scale_y_continuous(limits = c(0, 15), breaks = seq(0, 15, 3))
    g <- g + scale_fill_manual(values = c("darkseagreen4", "darkseagreen2", "darkseagreen1")) + theme(legend.position = "none")
    ggsave("../Results/Tumor_only_3group_genes_CDS_len_with_ZNFG23_noG1.pdf", width = 2.5, height = 2, units = "in")
    
    
  }
  
  
}
