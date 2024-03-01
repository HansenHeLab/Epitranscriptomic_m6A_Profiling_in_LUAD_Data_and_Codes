################################################################
## This script implements:
## 1. Association analysis among m6a, mRNA, protein clinic info at gene level
## 2. Survival analysis and KM plots based on m6a, mRNA, protein
## 3. BLVRA related analyses
## 4. BLVRA related analyses in validation cohort
#################################################################

rm(list = ls())
setwd("./data")

###########################################################################
##  Association analysis among m6a, mRNA, protein clinic info at gene level 
###########################################################################
{
 
  library("s2dverification")       # for color bar separately
  library("ggplot2")
  library("dplyr")
  library("ggpubr")
  
  ####################
  # load m6a level 
  if(False){
  load ("/Users/Yong/Yong/m6a_profiling/6_mrna_m6a_protein/1_all_overlap/log2_mrna_m6a_protein_shared_all.Rdata")
  s_info <- read.table("/Users/Yong/Yong/m6a_profiling/2_peak_calling/1_peaks/merge/sample_info.txt", header = T) 
  
  ### slected gene with m6A and protein across all samples 
  idx_1 <- rowSums(is.na(m6a_ms)) == 0                ## all samples with m6A
  idx_2 <- rowSums(is.na(protein_ms)) < 11            ## at least 40 samples with protein
  idx_k <- idx_1 & idx_2
  sum(idx_k)
  
  m6a_ms <- m6a_ms[idx_k, ]
  mrna_ms <- mrna_ms[idx_k, ]
  protein_ms  <- protein_ms[idx_k, ]
  
  dim(mrna_ms)
  
  save(m6a_ms, mrna_ms, protein_ms, file = "log2_mrna_m6a_protein_shared_more_than_40samples.Rdata")
  }
  
  load("log2_mrna_m6a_protein_shared_more_than_40samples.Rdata")
  s_info <- read.table("sample_info.txt", header = T) 

  ####################################################
  ## pre-defined function and run association analysis
  {
  cor.clinical4gene <- function(expr, info, prefix)
  {
    library("dplyr")
    library("ggpubr")
    library("RColorBrewer") 
    
    N <- nrow(expr)
    M <- ncol(info)-1                                    # clinical features, first column is the tumorIDs
    gene_p <- matrix(0, N, M)                            # p values for the clinical and m6a level
    colnames(gene_p) <- colnames(info)[2:(M + 1)] 
    rownames(gene_p) <- rownames(expr)
    
    ## anova test 
    for(i in 1:N)
    {
      expr_ge <- expr[i, ]
      dat <- data.frame(info, expr_ge)
      
      for(j in 2:(M + 1))
      {
        tmp <- aov(expr_ge ~ dat[, j], data = dat)
        gene_p[i, j-1] <- summary(tmp)[[1]][["Pr(>F)"]][1]
      }
    }
    
    ## multiple test 
    gene_fdr <- gene_p
    L <- ncol(gene_fdr)
    for (i in 1:L )
    {
      gene_fdr[, i] <- p.adjust(gene_fdr[, i], method = "BH")
    }
    
    ## merge pvale and fdr to vector for plotting
    p_merge <- -log10(as.vector(gene_p))
    fdr_merge <- -log10(as.vector(gene_fdr))
    
    ## -log10 pvalue > 5  <- 5
    p_merge[p_merge >= 5] <- 5
    
    p_group <- rep(colnames(gene_p), each = nrow(gene_p))
    dat_p <- data.frame(p_merge, p_group)
    dat_fdr <- data.frame(fdr_merge, p_group)
    p_max <- max(p_merge) + 1
    fdr_max <- max(fdr_merge) + 1
    
    ##  significantly associated with clinical infomation based on P value
    g <- ggplot(dat_p, aes(x = p_group, y = p_merge))  +geom_point(aes( col = p_group),stat="identity", position = position_jitter(width = 0.3), show.legend = FALSE) + theme_bw() + theme(axis.line = element_line(colour = "black"),panel.grid.major = element_blank(),panel.grid.minor = element_blank(),panel.border =  element_blank(),  panel.background = element_blank()) 
    # g <- g + geom_hline(yintercept=-log10(0.01), linetype="dashed", color = "red", size=0.5)
    g <- g + geom_hline(yintercept=-log10(0.05), linetype="dashed", color = "red", size=0.5)
    g <- g + scale_y_continuous(limits = c(0, p_max), expand = c(0, 0)) + labs(x = "Clinical Features", y = "-log10(p_value)") 
    g <- g + scale_color_manual(values = brewer.pal(n = 8, name = "Accent")[5:8])
    name <- paste(prefix, "_p_value.pdf", sep = "")
    ggsave(name,  width = 4, height = 3)
    #ggsave(name,  width = 4, height = 3, units = "in", dpi = 600)
    
    
    ## significantly associated with clinical infomation based on FDR
    g <- ggplot(dat_fdr, aes(x = p_group, y = fdr_merge))  +geom_point(aes( col = p_group),stat="identity", position = position_jitter(width = 0.3), show.legend = FALSE) + theme_bw() + theme(axis.line = element_line(colour = "black"),panel.grid.major = element_blank(),panel.grid.minor = element_blank(),panel.border =  element_blank(),  panel.background = element_blank()) 
    g <- g + geom_hline(yintercept=-log10(0.05), linetype="dashed", color = "red", size=0.5)
    g <- g + scale_y_continuous(limits = c(0, fdr_max), expand = c(0, 0)) + labs(x = "Clinical Features", y = "-log10(FDR)")
    g <- g + scale_color_manual(values = brewer.pal(n = 8, name = "Accent")[5:8])
    name <- paste(prefix, "_FDR.pdf", sep = "")
    ggsave(name,  width = 4, height = 3)
    #ggsave(name,  width = 4, height = 3, units = "in", dpi = 600)
    
    p_return <- list(gene_p, gene_fdr)
    return(p_return)
    
  }
  
  idx_si <- match(colnames(m6a_ms), s_info$helab_id)
  s_info_ss <- s_info[idx_si, ]
  s_info_ss <- s_info_ss[, c(1, 2, 4, 8)]              # 3 featuress
  
  ## somking levels order 
  s_info_ss$smoking <- factor(s_info_ss$smoking, 
                              levels =c("Current", "Ex-Smoker", "Never"))
  
  
  m6a_cli <- cor.clinical4gene(m6a_ms, s_info_ss, "../Results/gene_m6A_assoc_clinical_info")
  mrna_cli <- cor.clinical4gene(mrna_ms, s_info_ss, "../Results/gene_mRNA_assoc_clinical_info")
  protein_cli <- cor.clinical4gene(protein_ms, s_info_ss, "../Results/gene_protein_assoc_clinical_info")
  }
  
  ####################################################
  ## plot the top candidates 
  {
    
    library(RColorBrewer) 
    
    ## for sex
    col_set <- list(brewer.pal(n = 8, name = "Set2")[4:3], brewer.pal(n = 11, name = "RdYlGn")[c(2, 5, 9)],
                    brewer.pal(n = 11, name = "BrBG")[c(4, 3, 1)])
    
  
    gene_p <- m6a_cli[[1]]
    gene_mp <- mrna_cli[[1]]
    gene_pp <- protein_cli[[1]] 
    
    for(i in 1: 3)
    {
      idx <- which(gene_p[, i] == min(gene_p[, i]) )
      p_min <- formatC(gene_p[idx, i], digits = 2)
      
      g_n <- rownames(gene_p)[idx]
      idx_g <- match(g_n, rownames(m6a_ms))
      dat <- data.frame(m6a_ms[idx_g, ], s_info_ss[, i + 1])
      colnames(dat) <- c("m6a_ge", colnames(s_info_ss)[i + 1])
      
      tt <- paste(g_n, ": p_value = ", p_min, sep = "")
      g <- ggline(dat, x = colnames(s_info_ss)[i + 1] , y = "m6a_ge",  col = colnames(dat)[2], 
                  add = c("mean_se", "jitter"), add.params = list(shape = 1),
                  ylab = "m6A level", xlab = colnames(s_info_ss)[i + 1])
      g <- g + scale_color_manual(values =  col_set[[i]]) + theme(legend.position = "none") 
      g <- g + ggtitle(tt) + theme(plot.title = element_text(size = 8))
      
      name <- paste("../Results/m6A_Most_sig_gene_under_", colnames(s_info_ss)[i + 1], ".pdf", sep = "")
      ggsave(name,  width = 2.5, height = 3, units = "in")
      
      
      ## for match mRNA levl
      p_mrna <- formatC(gene_mp[idx, i], digits = 2)
      
      idx_m <- match(g_n, rownames(mrna_ms))
      dat <- data.frame(mrna_ms[idx_m, ], s_info_ss[, i + 1])
      colnames(dat) <- c("mrna", colnames(s_info_ss)[i + 1])
      tt <- paste(g_n, ": p_value = ", p_mrna, sep = "")
      
      g <- ggline(dat, x = colnames(s_info_ss)[i + 1] , y = "mrna", col = colnames(dat)[2],
                  add = c("mean_se", "jitter"), add.params = list(shape = 1), 
                  ylab = "mRNA level", xlab = colnames(s_info_ss)[i + 1])
      
      g <- g + scale_color_manual(values =  col_set[[i]]) + theme(legend.position = "none") 
      g <- g + ggtitle(tt) + theme(plot.title = element_text(size = 8))
      name <- paste("../Results/m6A_Most_sig_gene_mRNA_under_", colnames(s_info_ss)[i + 1], ".pdf", sep = "")
      ggsave(name,  width = 2.5, height = 3, units = "in")
      
      ## for match protein levle 
      p_protein <- formatC(gene_pp[idx, i], digits = 2)
      
      idx_m <- match(g_n, rownames(protein_ms))
      dat <- data.frame(protein_ms[idx_m, ], s_info_ss[, i + 1])
      colnames(dat) <- c("Protein", colnames(s_info_ss)[i + 1])
      tt <- paste(g_n, ": p_value = ", p_protein, sep = "")
      
      g <- ggline(dat, x = colnames(s_info_ss)[i + 1] , y = "Protein", col = colnames(dat)[2], 
                  add = c("mean_se", "jitter"),  add.params = list(shape = 1), 
                  ylab = "Protein level", xlab = colnames(s_info_ss)[i + 1])
      
      g <- g + scale_color_manual(values =  col_set[[i]]) + theme(legend.position = "none") 
      g <- g + ggtitle(tt) + theme(plot.title = element_text(size = 8))
      name <- paste("../Results/m6A_Most_sig_gene_Protein_under_", colnames(s_info_ss)[i + 1], ".pdf", sep = "")
      ggsave(name,  width = 2.5, height = 3, units = "in")
    }
      
      
    
  
  }
  
}


##################################################
##  Survival analysis based on m6a, mRNA, protein
##################################################
{
  ### Load required packages and process survial info
  {
  library(survival)
  library(survminer)
  library(dplyr)
  library(ggplot2)
  library(VennDiagram)
  library(RColorBrewer)  
  
  idx <- !is.na(s_info$surv_time)
  s_infos <- s_info[idx, ]
  
  ## for death 
  N <- nrow(s_infos)
  censor_d <- rep(0, N)
  idx_d <- s_infos$status == "D"
  censor_d[idx_d] = 1
  
  ## for reccurence
  censor_rec <- rep(0, N)
  idx_rec <- s_infos$rec == "Y"
  censor_rec[idx_rec] = 1
  s_infos <- cbind(s_infos, censor_d, censor_rec)
  
  ######### matching 
  idx <- match(s_infos$helab_id, colnames(m6a_ms))
  s_info_ss <- s_infos[!is.na(idx), ]  
  
  ## resorted colnames to match with s_info_ss
  m6a_ms_f <- m6a_ms[, idx[!is.na(idx)]]
  mrna_ms_f <- mrna_ms[, idx[!is.na(idx)]]
  protein_ms_f <- protein_ms[, idx[!is.na(idx)]]
  
  ### all samples with protein 
  idx_ff <- rowSums(is.na(protein_ms_f)) == 0
  m6a_ms_f <- m6a_ms_f[idx_ff, ]
  mrna_ms_f <- mrna_ms_f[idx_ff, ]
  protein_ms_f <- protein_ms_f[idx_ff, ]
  }
  
  ## run survival analysis with top and bottom 22 smaples 
  expr_list <- list(m6a_ms_f, mrna_ms_f, protein_ms_f)
  sur_p_list <- list()
  
  ##########################
  ## all availabe samples 45
  ## for top 10, 15, 22
  for(sub_s in c(22))
  {    
    ############################
    ## survial pvalue 
    {   
      for(k in 1: length(expr_list))
      {
        expr <- expr_list[[k]] 
        info <- s_info_ss
        
        ##
        {
          library(survival)
          library(survminer)
          library(dplyr)
          
          N <- nrow(expr)
          M <- ncol(expr)
          gene_p <- vector()                            
          
          ## survival test 
          for(i in 1:N)
          {
            expr_g <- expr[i, ] 
            
            LL <- length(expr_g)
            idx_gs <- order(expr_g) 
            
            ## low 1/3 and high 1/3
            expr_gs <- c(expr_g[idx_gs[1:sub_s]], expr_g[idx_gs[(LL-sub_s+1):LL]])
            relative <- rep("High", length(expr_gs))
            relative[1:sub_s] <- "Low"
            relative <- as.factor(relative)
            
            ## corresponding samples 
            idx_ss <- match(names(expr_gs), info$helab_id)
            info_t <- info[idx_ss, ] 
            dat_t <- data.frame(info_t, relative)
            surv_d <- Surv(time = info_t$surv_time, event = info_t$censor_d)
            fit_d  <- survfit(surv_d ~ relative , data = dat_t)
            gene_p[i] <- surv_pvalue(fit_d)[2][1, 1]
            
          }
          
          names(gene_p) <- rownames(expr)
        }
        
        sur_p_list[[k]] <- gene_p
      }
      #save(m6a_ms_f, protein_ms_f, mrna_ms_f, s_info_ss, sur_p_list, file = paste("top_bottom_", sub_s, "_samples_based_survial_analysis_related.RData", sep = ""))
    }  
  }
  
  ################################
  ## pval | FDR  for each platform    
  {
    p_merge <- c(sur_p_list[[1]], sur_p_list[[2]], sur_p_list[[3]])
    p_merge_fdr <- c(p.adjust(sur_p_list[[1]], method = "BH"), p.adjust(sur_p_list[[2]], method = "BH"), p.adjust(sur_p_list[[3]], method = "BH"))
    p_merge <- -log10(p_merge)
    p_merge_fdr <-  -log10(p_merge_fdr)
    
    p_group <- c(rep("m6A", length(sur_p_list[[1]])), rep("mRNA", length(sur_p_list[[2]])), rep("Protein", length(sur_p_list[[3]])))
    dat_p <- data.frame(p_merge, p_group)
    p_max <- max(p_merge, na.rm = T) + 1
    
    ## based on p value 
    ## genes' m6A , mRNA, Protein level are significantly associated with survial information based on P value
    g <- ggplot(dat_p, aes(x = p_group, y = p_merge))  +geom_point(aes( col = p_group),stat="identity", position = position_jitter(width = 0.3), show.legend = FALSE) + theme_bw() + theme(axis.line = element_line(colour = "black"),panel.grid.major = element_blank(),panel.grid.minor = element_blank(),panel.border =  element_blank(),  panel.background = element_blank()) 
    #g <- g + geom_hline(yintercept=-log10(0.01), linetype="dashed", color = "red", size=0.5)
    g <- g + geom_hline(yintercept=-log10(0.05), linetype="dashed", color = "red", size= 0.5)
    g <- g + scale_y_continuous(limits = c(0, p_max), expand = c(0, 0)) + labs(x = "Features", y = "-log10(p_value)") 
    g <- g + scale_color_manual(values = brewer.pal(n = 8, name = "Set1")[1:3])
    ggsave("../Results/gene_m6A_mRNA_protein_assoc_survail_info.pdf" ,  width = 3, height = 3)
    
    ## based on FDR 
    dat_p <- data.frame(p_merge_fdr, p_group)
    p_max <- max(p_merge_fdr, na.rm = T) + 0.5
    g <- ggplot(dat_p, aes(x = p_group, y = p_merge_fdr))  +geom_point(aes( col = p_group),stat="identity", position = position_jitter(width = 0.3), show.legend = FALSE) + theme_bw() + theme(axis.line = element_line(colour = "black"),panel.grid.major = element_blank(),panel.grid.minor = element_blank(),panel.border =  element_blank(),  panel.background = element_blank()) 
    #g <- g + geom_hline(yintercept=-log10(0.01), linetype="dashed", color = "red", size=0.5)
    g <- g + geom_hline(yintercept=-log10(0.05), linetype="dashed", color = "red", size=0.5)
    g <- g + scale_y_continuous(limits = c(0, p_max), expand = c(0, 0)) + labs(x = "Features", y = "-log10(FDR)") 
    g <- g + scale_color_manual(values = brewer.pal(n = 8, name = "Set1")[1:3])
    ggsave(paste("top_bottom_", sub_s, "_samples_", "gene_m6A_mRNA_protein_assoc_survail_info_FDR.pdf", sep = "") ,  width = 3, height = 3)
  }  
  
 
}

##################################################
## 3 BLVRA related analyses
##################################################
{

  #################################
  ## KM plot for BLVRA
  {
    m6a_p <- sur_p_list[[1]]
    mrna_p <- sur_p_list[[2]]
    protein_p <- sur_p_list[[3]] 
    
    #gene_l <- c(names(m6a_p)[order(m6a_p)[1:5]])
    gene_l <"BLVRA"
    expr_l <- list(m6a_ms_f, mrna_ms_f, protein_ms_f)
    expr_type <- c("m6A", "mRNA", "Protein")
    
    ## loop structure
    {
      L <- length(gene_l)
      M <- length(expr_l)
      for (i in 1: L)
      {
        for (j in 1: M)
        {
          gene_n <- gene_l[i];
          expr_t <- expr_l[[j]]
          type <- expr_type[j]
          info <- s_info_ss
          idx <- match(gene_n, rownames(expr_t))
          
          ### with records
          if (!is.na(idx))
          {
            
            expr_g <- expr_t[idx, ] 
            LL <- length(expr_g)
            idx_gs <- order(expr_g) 
            
            ## top/bottom sub_s
            expr_gs <- c(expr_g[idx_gs[1:sub_s]], expr_g[idx_gs[(LL-sub_s+1):LL]])
            relative <- rep("High", length(expr_gs))
            relative[1:sub_s] <- "Low"
            relative <- as.factor(relative)
            
            ## corresponding samples 
            idx_ss <- match(names(expr_gs), info$helab_id)
            info_t <- info[idx_ss, ] 
            
            dat <- data.frame(info_t, relative)
            surv <- Surv(time = info_t$surv_time, event = info_t$censor_d)
            fit  <- survfit(surv ~ relative , data = dat) 
            
            #g <- ggsurvplot(fit, data = dat, pval = TRUE, title = paste(gene_n, type, sep = "_"), risk.table = TRUE)
            ggsurvplot(fit, data = dat, pval = TRUE, palette = brewer.pal(n = 12, name = "Paired")[10:9], title = paste(gene_n, type, sep = "_"), font.title = c(10)) 
            ggsave(paste("../Results/", gene_n, "_", type, "_based_survival_KM_plot.pdf", sep = ""), width = 3, height = 3)
            
          }
        }
      }
    }
  }
  
  
  #############################
  ## KM plot for BLVRA in TCGA
  {
    gene_s <- "BLVRA"
    s_info_tcga <- read.csv("TCGA_LUAD_surv.csv")
    
    idx_rm <- s_info_tcga$OS_MONTHS == 0 | s_info_tcga$OS_MONTHS == "[Not Available]"
    s_info_tcga <- s_info_tcga[!idx_rm, ]
    event <- rep(0, nrow(s_info_tcga))
    event[s_info_tcga$OS_STATUS == "DECEASED"] = 1    
    s_info_tcga  <-  cbind(s_info_tcga , event)
    
    tcga_expr_sn <- substr(colnames(tcga_t), 0 , 12)
    
    idx_ss <- match(s_info_tcga$PATIENT_ID, tcga_expr_sn)
    
    s_info_tcga_s <- s_info_tcga[!is.na(idx_ss), ]
    
    tcga_t_s <- tcga_t[, idx_ss[!is.na(idx_ss)]]
    sur_time <-  as.numeric(as.character(s_info_tcga_s$OS_MONTHS))
    s_info_tcga_s  <-  cbind(s_info_tcga_s , sur_time)
    
    ## for BLVRA 
    idx_gs <- match(gene_s, rownames(tcga_t_s))
    
    idx_ss <- match(colnames(tcga_t_s), colnames(tcga_t))
    tcga_expr_gene_s <- tcga_t[idx_gs, idx_ss]
    
    
    expr_l <- list(tcga_expr_gene_s)
    expr_type <- paste(gene_s, "_TCGA_mRNA_All", sep = "")
    
    ##########
    ## Km Plot
    {
      M <- length(expr_l)
      for (j in 1: M)
      {
        expr_t <- expr_l[[j]]
        type <- expr_type[j]
        info <- s_info_tcga_s
        
        ### with records
        if (TRUE)
        {
          idx_na <- !is.na(expr_t)
          expr_g <- expr_t[idx_na] 
          info_tt <- info[idx_na, ] 
          idx_h <- expr_g >=  median(expr_g)
          
          
          relative <- rep("Low", length(expr_g))
          relative[idx_h] <- "High"
          relative <- as.factor(relative)
          
          dat <- data.frame(info_tt, relative)
          surv <- Surv(time = info_tt$sur_time, event = info_tt$event)
          fit  <- survfit(surv ~ relative , data = dat) 
          
          #g <- ggsurvplot(fit, data = dat, pval = TRUE, title = paste(gene_n, type, sep = "_"), risk.table = TRUE)
          ggsurvplot(fit, data = dat, pval = TRUE, palette = brewer.pal(n = 12, name = "Paired")[10:9], title = type, font.title = c(8)) 
          ggsave(paste("../Results/", type, "_based_survival_KM_plot_TCGA.pdf", sep = ""), width = 3, height = 3)
          
        }
      }
    }
    
  }
  
  #################################
  ## BLVRA gene and peak m6a levels
  {
    load("BLVRA_peaks_level_m6A_intensity.RData")
    
    gene_s <- "BLVRA"
    idx_g <- match(gene_s, rownames(m6a_ms_f))
    
    idx_s <- match(colnames(m6a_ms_f), rownames(peak2))   
    
    peak2_f <- peak2[idx_s, ]    ## BLVRA narrow peak peak sum ip / input
    
    #cor.plot(m6a_ms_f[idx_g, ], log2(peak2_f$peak_sum_fd), paste0("../Results/", gene_s, "_peak_vs_gene_m6A_level_added"), "gene_m6a_level", "peak_m6a_level_sum", "PDF")
    
    file_name <- paste("../Results/", gene_s, "_peak_vs_gene_m6A_level.pdf", sep = "");
    m6a_g <- m6a_ms_f[idx_g, ]
    m6a_p <- log2(peak2_f$peak_sum_fd)
   
    dat <- data.frame(m6a_g, m6a_p)
    g <- ggplot(dat, aes(x=m6a_g, y=m6a_p)) + geom_point(color='#2980B9', size = 2)  
    g <- g + geom_smooth(method=lm, color='#2C3E50') + theme_classic()
    ggsave(file_name, width = 3.5, height = 3)
    
  }
  
  #######################################
  ## BLVRA peak m6a levels based KM plots
  {
    peak_m6a <- peak2_f$peak_sum_fd
    names(peak_m6a) <- rownames(peak2_f)
    
    for(sub_s in c(22))
      
    {    
      ############################
      ## for surival KMPlot
      {
        
        expr_l <- list(peak_m6a)
        expr_type <- c("../Results/BLVRA_peak_level_m6A")
        
        ## loop structure
        {
          
          M <- length(expr_l)
          for (j in 1: M)
          {
            expr_t <- expr_l[[j]]
            type <- expr_type[j]
            info <- s_info_ss
            
            ### with records
            if (TRUE)
            {
              
              expr_g <- expr_t 
              LL <- length(expr_g)
              idx_gs <- order(expr_g) 
              
              ## top/bottom sub_s
              expr_gs <- c(expr_g[idx_gs[1:sub_s]], expr_g[idx_gs[(LL-sub_s+1):LL]])
              relative <- rep("High", length(expr_gs))
              relative[1:sub_s] <- "Low"
              relative <- as.factor(relative)
              
              ## corresponding samples 
              idx_ss <- match(names(expr_gs), info$helab_id)
              info_t <- info[idx_ss, ] 
              
              dat <- data.frame(info_t, relative)
              surv <- Surv(time = info_t$surv_time, event = info_t$censor_d)
              fit  <- survfit(surv ~ relative , data = dat) 
              
              #g <- ggsurvplot(fit, data = dat, pval = TRUE, title = paste(gene_n, type, sep = "_"), risk.table = TRUE)
              ggsurvplot(fit, data = dat, pval = TRUE, palette = brewer.pal(n = 12, name = "Paired")[10:9], title = paste(type, sep = "_"), font.title = c(10)) 
              ggsave(paste(type, "_based_survival_KM_plot.pdf", sep = ""), width = 3, height = 3)
              
            }
          }
        }
        
        
      }
      
    }
  }
  
  ##################################################
  ## correlation across BLVRA m6a, mRNA and protein
  {
    ### ggplot with regression confidence intervals 
    
      library(ggplot2)
      
      gene_s <- "BLVRA"
      idx_g <- match(gene_s, rownames(m6a_ms))
      
      name2 = paste0("../Results/", gene_s, "_cor m6A_mRNA in Tumor")
      name3 = paste0("../Results/", gene_s, "_cor m6a_Protein in Tumor")
      name4 = paste0("../Results/", gene_s, "_cor mRNA_Protein in Tumor")
      
      ## m6A, mRNA
      file_name <- paste(name2, "_with_interval.pdf", sep = "");
      m6a <- m6a_ms[idx_g, ]
      mrna <- mrna_ms[idx_g, ]
      dat <- data.frame(m6a, mrna)
      g <- ggplot(dat, aes(x=m6a, y=mrna)) + geom_point(color='#2980B9', size = 2)  
      g <- g + geom_smooth(method=lm, color='#2C3E50') + theme_classic()
      ggsave(file_name, width = 3.5, height = 3)
    
      ## m6A, protein
      
      file_name <- paste(name3, "_with_interval.pdf", sep = "");
      m6a <- m6a_ms[idx_g, ]
      protein <- protein_ms[idx_g, ]
      dat <- data.frame(m6a, protein)
      g <- ggplot(dat, aes(x=m6a, y=protein)) + geom_point(color='#2980B9', size = 2)  
      g <- g + geom_smooth(method=lm, color='#2C3E50') + theme_classic()
      ggsave(file_name, width = 3.5, height = 3)
      
      ## mRNA and protein
      file_name <- paste(name4, "_with_interval.pdf", sep = "");
      mrna <- mrna_ms[idx_g, ] 
      protein <- protein_ms[idx_g, ]
      dat <- data.frame(mrna, protein)
      g <- ggplot(dat, aes(x=mrna, y=protein)) + geom_point(color='#2980B9', size = 2)  
      g <- g + geom_smooth(method=lm, color='#2C3E50') + theme_classic()
      ggsave(file_name, width = 3.5, height = 3)
    }
}

##################################################
##  BLVRA related analyses in validation cohort
##################################################
{
  
  ####### loading and matching data
  {
  sample_80 <- read.csv("MB_80_samples_for_validation.csv") 
  m6a_80 <-read.csv("MB_80_samples_BLVRA_m6A_levels.csv")
  
  ##### 5 probes 
  gene_s_expr <- read.table("MB_BLVRA_array_expr.txt", as.is = T, header = T, row.names = 1) 
  
  idx_s <- match(sample_80$sample, colnames(gene_s_expr))
  expr <- t(gene_s_expr[, idx_s])
  
  ## BLVRA mean value across probes
  BLVRA_mean <- rowMeans(expr)
  
  ## selected probe for the best probe, which is with maximum mean value 
  BLVRA_max <- expr[, which(colMeans(expr) == max(colMeans(expr)))]
  
  ## BLVRA best probes based on jetset score
  library(jetset)
  p_score <- jscores("hgu133plus2", symbol = gene_s)
  p_best <- rownames(p_score)[which(p_score$overall == max(p_score$overall))]
  idx_best <- match(p_best, colnames(expr))
  BLVRA_best <- expr[, idx_best]
  
  expr_80 <- cbind(expr, BLVRA_mean, BLVRA_max, BLVRA_best)
  
  ############
  ## m6A, mRNA
  file_name <- "../Results/BLVRA_cor_log2_m6A_vs_RNA_in_MB_cohort.pdf"
  m6a <- log2(m6a_80$Corrected_relative_methylated_level)
  mrna <- expr_80[, 8]   
  
  dat <- data.frame(m6a, mrna)
  g <- ggplot(dat, aes(x=m6a, y=mrna)) + geom_point(color='#2980B9', size = 2)  
  g <- g + geom_smooth(method=lm, color='#2C3E50') + theme_classic()
  ggsave(file_name, width = 3.5, height = 3)
  
  cor.test(m6a, mrna)
  }
  
  
  ######################################################################
  ## based on top/bottom m6a level 40 samples survival analysis
  {
    ## different subsample size FOR 2 sub groups
    
    expr_all <- cbind(expr_80[, 8], m6a_80$Corrected_relative_methylated_level)      ## check best probe for each gene only
   
    colnames(expr_all) <- c("BLVRA_mRNA", "BLVRA_m6A")
    
    BLVRA_sur <- cbind(sample_80, expr_all)
    
    #save(BLVRA_sur, file = "MB_BLVRA_sur_analysis_data.Rdata")
    
    ## for log-rank test
    for(sub_s in c(40))
    {    
    
      ############################
      ## for surival KMPlot
      {
        
        gene_l <- c(colnames(expr_all))
        gene_l
        expr_l <- list(expr_all)
        expr_type <- c("")
        
        ## loop structure
        {
          L <- length(gene_l)
          M <- length(expr_l)
          for (i in 1: L)
          {
            for (j in 1: M)
            {
              gene_n <- gene_l[i];
              expr_t <- t(expr_l[[j]])     ##
              type <- expr_type[j]
              info <- sample_80
              idx <- match(gene_n, rownames(expr_t))
              
              ### with records
              if (!is.na(idx))
              {
                expr_g <- expr_t[idx, ] 
                LL <- length(expr_g)
                idx_gs <- order(expr_g) 
                
                ## top/bottom sub_s
                expr_gs <- c(expr_g[idx_gs[1:sub_s]], expr_g[idx_gs[(LL-sub_s+1):LL]])
                relative <- rep("High", length(expr_gs))
                relative[1:sub_s] <- "Low"
                relative <- as.factor(relative)
                
                ## corresponding samples 
                idx_ss <- match(names(expr_gs), info$sample)
                info_t <- info[idx_ss, ] 
                
                dat <- data.frame(info_t, relative)
                surv <- Surv(time = info_t$surv_time, event = info_t$censor_d)
                fit  <- survfit(surv ~ relative , data = dat) 
                
                #g <- ggsurvplot(fit, data = dat, pval = TRUE, title = paste(gene_n, type, sep = "_"), risk.table = TRUE)
                ggsurvplot(fit, data = dat, pval = TRUE, palette = brewer.pal(n = 12, name = "Paired")[10:9], title = paste(gene_n, type, sub_s, sep = "_"), font.title = c(10)) 
                ggsave(paste("../Results/", gene_n, "_based_survival_KM_plot_in_MB_cohort.pdf", sep = ""), width = 3, height = 3)
                
              }
            }
          }
        }
        
        
      }
    }
    
  }
  
}
  