################################################################
## This script implements:
## 1. clustering based on m6A, RNA and protein 
## 2. m6A based clusters heatmap 
## 3. m6A based clusters characteristics
## 4. m6A based clusters survival analyses
## 5. p2 vs P1,3, 4, 5 DEGs 
#################################################################

rm(list = ls())
setwd("./data")

###############################################
##  clustering based on m6A, RNA and protein 
###############################################
{

  ################################
  ## load m6a, protein and rna
  ## select the most variable ones
  {
  library(ConsensusClusterPlus)
  library(matrixStats)
  ## m6a 
  load("m6A_level_4gene_passed_filtering.Rdata")
  
  ## rna
  load("./paired_63_count_ercc_deseq.Rdata") 

  ## remove outliers tumor 45 and tumor9
  idx_s <- match(c("tumor9", "tumor45"), colnames(input_rpkm))
  input_rpkm <- input_rpkm[, -idx_s]
  
  ## protein
  protein_ori <- readRDS("protein_expr_51.rds")
  
  ## focusing on Tumor 
  idx_k = rowSums(is.na(m6a_rt)) == 0 
  m6a <- m6a_rt[idx_k, ]                   ## ip/input, without NA
  
  ## rna for match samples 
  idx_s <- match(colnames(m6a), colnames(input_rpkm))
  rna <- input_rpkm[, idx_s]
  idx_k <- rowMedians(rna) > 1
  rna <- rna[idx_k, ]
  
  ## protein
  idx_k <- rowSums(is.na(protein_ori)) == 0 
  idx_s <- match(colnames(m6a), colnames(protein_ori))
  protein <- 2^(protein_ori[idx_k, idx_s])                 ## protein without NA
  
  #filter based on variance (top 20%)   
  {
    
    prob_cutoff = 0.80
    
    rna_iqr<- apply(rna, 1,  IQR)
    cutoff <- quantile(rna_iqr, probs = prob_cutoff)
    rna_iqr <- rna[rna_iqr > cutoff, ]   
    rna_iqr_z <- t(scale(t(rna_iqr)))       ## z transform
    
    m6a_iqr<- apply(m6a, 1,  IQR)
    cutoff <- quantile(m6a_iqr, probs = prob_cutoff)
    m6a_iqr <- m6a[m6a_iqr > cutoff, ]  
    m6a_iqr_z <- t(scale(t(m6a_iqr)))       ## z transform
    
    protein_iqr<- apply(protein, 1,  IQR)
    cutoff <- quantile(protein_iqr, probs = prob_cutoff)
    protein_iqr <- protein[protein_iqr > cutoff, ] 
    protein_iqr_z <- t(scale(t(protein_iqr)))       ## z transform
    
    #save(rna_iqr, m6a_iqr, protein_iqr, rna_iqr_z, m6a_iqr_z, protein_iqr_z, 
    #     file = "m6A_RNA_Protein_IQR_top_20perc.rda")
    
  }
  
}

  #########################
  ## sample wise clustering
  ## m6A, RNA, Protein
  {
  expr <- list(m6a_iqr_z, rna_iqr_z, protein_iqr_z)
  type <- c("m6A", "RNA", "Protein")
  
  for (i in 1:3)
  {
    epxr_s <- as.matrix(expr[[i]])
    
    sample_cluster <- ConsensusClusterPlus(
      epxr_s,
      maxK = 10,
      reps = 1000,
      pItem = 0.8,
      pFeature = 0.8,
      clusterAlg = 'hc',
      distance = 'euclidean',
      innerLinkage = 'ward.D',
      finalLinkage = 'ward.D',
      seed = 17,
      title = paste0(type[i], '_based_sample_ConsensusClusterPlus'),
      writeTable = TRUE,
      corUse = 'complete.obs',
      verbose = TRUE,
      plot = 'pdf'
    )
    
    saveRDS(sample_cluster, file = paste0(type[i], "_based_sample_ConsensusClusterPlus.RDS"))
    
  }
  }

  #########################
  ## gene wise clustering
  ## m6A only
  {
  gene_cluster <- ConsensusClusterPlus(
    as.matrix(t(m6a_iqr_z)),
    maxK = 10,
    reps = 1000,
    pItem = 0.8,
    pFeature = 0.8,
    clusterAlg = 'hc',
    distance = 'euclidean',
    innerLinkage = 'ward.D',
    finalLinkage = 'ward.D',
    seed = 17,
    title = paste0('m6A_gene_ConsensusClusterPlus'),
    writeTable = TRUE,
    corUse = 'complete.obs',
    verbose = TRUE,
    plot = 'pdf'
  )
  saveRDS(gene_cluster, file = paste0("m6A_based_gene_ConsensusClusterPlus.RDS"))
  
}
}

#############################
## m6A based clusters heatmap
#############################
{
  library(BoutrosLab.plotting.general)
  library(getopt)
  
  load("m6A_RNA_Protein_IQR_top_20perc.rda")
  sample_cluster_m6a <- readRDS("m6A_based_sample_ConsensusClusterPlus.RDS")
  gene_cluster_m6a  <-  readRDS("m6A_based_gene_ConsensusClusterPlus.RDS")
  
  ## clinical features 
  clin_f <- read.csv("sample_info_53.csv")
  idx_s <- match(colnames(m6a_iqr_z), clin_f$helab_id)
  clin_fs <- clin_f[idx_s, c('helab_id','age', 'sex','smoking','pathStage','mutation')] 
  
  nClu_s <- 5    ## number of selected clusters for samples
  nClu_g <- 5     ## number of selected clusters for genes
  
  s_idx <- sample_cluster_m6a[[nClu_s]]$consensusTree$order
  sample_ordered <- (sample_cluster_m6a[[nClu_s]]$consensusClass[s_idx])
  
  saveRDS(sample_ordered, file = paste0("m6A_based_sample_in_", nClu_s, "_clusters.RDS"))
  
  g_idx <- gene_cluster_m6a[[nClu_g]]$consensusTree$order
  gene_ordered <- (gene_cluster_m6a[[nClu_g]]$consensusClass[g_idx])
  #saveRDS(gene_ordered, file = paste0("m6A_based_gene_in_", nClu_g, "_clusters.RDS"))

  
  ## resort scaled m6A 
  idx_s <- match(names(sample_ordered), colnames(m6a_iqr_z))
  idx_g <- match(names(gene_ordered), rownames(m6a_iqr_z))
  data <- as.data.frame(m6a_iqr_z[idx_g, idx_s])
  
  #################################################
  ### heatmap without spaces between the clusters
  ## convert >=4 to 4 for data (heatmap only)
  {
    range(data)
    data[data >= 4 ] = 4
    
    key.min <- -4  
    key.max <-  4
    key.colour.interval.num <- 50;
    key.scale <- c(
      seq(key.min, 0, -key.min / key.colour.interval.num),
      seq(0, key.max, key.max / key.colour.interval.num)
    );
    key.scale <- unique(key.scale);
    
    sample_col <- default.colours(nClu_s)[sample_ordered];
    gene_col <- default.colours(nClu_g)[gene_ordered];
    
    sample_covariate <- list(
      rect = list(
        col = 'transparent',
        fill = sample_col,
        lwd = 1.5
      )
    )
    
    gene_covariate <- list(
      rect = list(
        col = 'transparent',
        fill = rev(gene_col),
        lwd = 1.5
      )
    )
    
    
    create.heatmap(
      x = data,
      file = paste0("../Results/m6A_clustering_heatmap.png"),
      clustering.method = 'none',
      at = key.scale,
      colour.scheme = c('darkorchid4', 'white', 'darkgreen'),
      print.colour.key = TRUE,
      colourkey.cex = 1,
      same.as.matrix = TRUE,
      covariates.top = sample_covariate,
      covariates = gene_covariate,
      resolution = 500,
      scale.data = FALSE
    )
  }
  
  #################################################
  ### heatmap  spaces between the clusters
  {
    ###############################################
    # split samples and genes into separate matrices and save in list
    gene.split <- split(data, as.factor(gene_ordered));
    #head(gene.split)
    sample.and.gene <- lapply(gene.split, function(x) split(as.data.frame(t(x)), as.factor(sample_ordered)))
    
    #head(sample.and.gene)
    rownames(clin_fs) <- clin_fs$helab_id
    clin_fs$"helab_id" <- NULL
    clin_fs[clin_fs$sex == "F",]$sex <- "female"
    clin_fs[clin_fs$sex == "M",]$sex <- "male"
    clin_fs[clin_fs$pathStage == "1A",]$pathStage <- "I"
    clin_fs[clin_fs$pathStage == "1B",]$pathStage <- "I"
    clin_fs[clin_fs$pathStage == "2A",]$pathStage <- "II"
    clin_fs[clin_fs$pathStage == "2B",]$pathStage <- "II"
    clin_fs[clin_fs$pathStage == "3A",]$pathStage <- "III"
    clin_fs[clin_fs$pathStage == "3B",]$pathStage <- "III"
    clin_fs[clin_fs$pathStage == "4",]$pathStage <- "IV" 
    clin_fs$mutation[is.na(clin_fs$mutation)] <- "No Data"
    clin_fs$mutation[clin_fs$mutation == "KRAS "] <- "KRAS"
    
    # create m6A and covariate heatmaps 
    source('../R_scripts/functions/heatmap_colouring.R')
    
    plot.objects <- list()
    counter <- 1
    for (i in 1:nClu_s) {
      for (j in 1:nClu_g) {
        plot.objects[[counter]] <- create.subset.heatmap(sample.and.gene[[j]][[i]])
        counter <- counter + 1;
      }
      annot <- clin_fs[match(rownames(sample.and.gene[[1]][[i]]), rownames(clin_fs)), ]
      plot.objects[[counter]] <- create.covariate.heatmap(annot);
      counter <- counter + 1;
    }
    
    
    cov.length <- 500
    gene.length <- nrow(m6a_iqr_z) + cov.length
    sample.length <- ncol(m6a_iqr_z)
    
    # need to add white space after final sample.k etc is picked
    ylab.label <- rev(paste0('P', 1:nClu_s))
    xlab.label <- paste0('M', 1:nClu_g)
    
    # want each row / col to be the same height so scale based on how many samples or gene are in each cluster
    panel.heights <- rep(NA, nClu_s);
    for (i in 1:nClu_s) {
      panel.heights[i] <- nrow(sample.and.gene[[1]][[i]]) / sample.length
    }
    
    panel.widths <- rep(NA, nClu_g);
    for (i in 1:nClu_g) {
      panel.widths[i] <- ncol(sample.and.gene[[i]][[1]]) / gene.length
    }
    
    
    # rev panel heights because boxes are plotted from the bottom
    create.multiplot(
      plot.objects = plot.objects,
      file = paste0("../Results/m6A_clustering_heatmap_with_clinical_features.png"),
      plot.layout = c(nClu_g + 1, nClu_s),
      panel.heights = rev(panel.heights),
      panel.widths = c(panel.widths, cov.length / gene.length),
      y.spacing = -0.6,
      x.spacing = 0.2,
      x.relation = 'free',
      y.relation = 'free',
      ylab.label = ylab.label,
      xlab.label = xlab.label,
      xlab.padding = 0,
      xlab.to.xaxis.padding = -1.5,
      ylab.padding = 1,
      xlab.cex = 1,
      ylab.cex = 1,
      yaxis.tck = 0,
      xaxis.tck = 0,
      print.new.legend = TRUE,
      right.padding = 14,
      use.legacy.settings = TRUE,
      legend = list(
        inside = list(
          x = 1.05,
          y = 0.99,
          fun = legend.grob(
            covariates.legend,
            between.row = 0.75,
            size = 2.25,
            label.cex = 0.80,
            title.cex = 0.80,
            title.just = 'left'
          )
        )
      ),
      resolution = 500
    )
    
  }

}

######################################
## m6A based clusters characteristics
#####################################
{
  
  ##############################################
  ## m6A, RNA and protein based clusters overlap
  {
  library(ggplot2)
  library(ggpubr)
  library(dplyr)
  library(RColorBrewer)
  library(pheatmap)
  
  ## all 8030 gene m6a level Z 
  m6a_rt_z <- t(scale(t(m6a_rt)))  
  sample.order_5 <- sample_ordered
  
  clin_f <- read.csv("sample_survival_info_53.csv")
  idx_s <- match(colnames(m6a_iqr_z), clin_f$helab_id)
  clin_fs <- clin_f[idx_s, ] 
   
  clin_fs$clusters_5 <- sample.order_5[match(clin_fs$helab_id, names(sample.order_5))]
  
  ## protein and RNA based clusters
  cluster_rna <- readRDS("RNA_based_sample_ConsensusClusterPlus.RDS")
  cluster_protein <- readRDS("Protein_based_sample_ConsensusClusterPlus.RDS")
  
  nClu_s <- 5
  sample.order_m6a <- sort(sample.order_5) 
  
  s_idx <- cluster_rna[[nClu_s]]$consensusTree$order
  sample.order_rna <- (cluster_rna[[nClu_s]]$consensusClass[s_idx])
  idx <- match(names(sample.order_m6a), names(sample.order_rna))
  sample.order_rna  <- sample.order_rna[idx]
  
  s_idx <- cluster_protein[[nClu_s]]$consensusTree$order
  sample.order_protein <- (cluster_protein[[nClu_s]]$consensusClass[s_idx])
  idx <- match(names(sample.order_m6a), names(sample.order_protein))
  sample.order_protein  <- sample.order_protein[idx]
  
  ## m6a_rna_protein clusters
  clus <- data.frame(sample.order_m6a, sample.order_rna, sample.order_protein)
  for(i in 1:3){clus[, i] <- paste0("P", clus[, i])}
  
  #write.csv(clus, file = "m6a_rna_protein_clusters_assigment.csv")
  m6a_col <- brewer.pal(n = 8, name = "Dark2")[1:5][sample.order_m6a]
  mrna_col <- brewer.pal(n = 12, name = "Set3")[1:5][sample.order_rna]
  protein_col <-  brewer.pal(n = 12, name = "Set3")[5:9][sample.order_protein]
  
  col_v <- c(m6a_col, mrna_col, protein_col)
  
  ##########################
  ## check grops overlapping
  
  cluster_intersect <- function(a, b, name_a, name_b)
  {
    N <- length(unique(a))
    M <- length(unique(b))
    mm <- matrix(0, N, M)
    for(i in 1:N)
    {
      s1 <- names(a)[a == i]
      for (j in 1:M)
      {
        s2 <- names(b)[b == j]
        mm[i, j] = length(intersect(s1, s2))
      }
    }
    
    rownames(mm) <- paste0(name_a, "_P",1:N)
    colnames(mm) <- paste0(name_b, "_P",1:M)
    
    ## fisher's exact test
    ft <- fisher.test(mm, simulate.p.value = TRUE)
    tt <- paste0("Fisher's test: Pval = ", round(ft$p.value, 4))
    
    ## plot 
    library(pheatmap)
    pdf(paste0(name_a, "_", name_b, "_based_clusters_overlapping.pdf"), width = 3, height = 3)
    pheatmap(mm, display_numbers = T, number_format = "%.0f",
             color = colorRampPalette(c('white','coral2'))(20),
             cluster_rows = F, cluster_cols = F, fontsize_number = 15,
             fontsize_row = 5, fontsize_col = 5,angle_col = 0, 
             legend = FALSE, main = tt, cex.main = 0.5)
    dev.off()
    
    ## reture matrix 
    return(mm)
    
  }
  
  m6m <- cluster_intersect(sample.order_m6a, sample.order_rna, "../Results/m6A", "RNA")
  m6p <- cluster_intersect(sample.order_m6a, sample.order_protein, "../Results/m6A", "Protein")
  mp  <- cluster_intersect(sample.order_rna, sample.order_protein, "../Results/RNA", "Protein")
  }
  
  
  ##################
  ## mean m6A level
  ## for selected variable genes
  {
    m6a_rt_z <- t(scale(t(m6a_rt)))  
    
    mean_m6a <- list()
    mean_m6a_iqr <- list()
    sample.order <- sample.order_5
    
    for (i in 1:5)
    {
      idx <- sample.order == i
      
      ## all 8030 genes
      idx_s <- match(names(sample.order)[idx], colnames(m6a_rt_z))
      tmp <- m6a_rt_z[, idx_s]
      mean_m6a[[i]] <- rowMeans(tmp)
      
      ## IQR genes
      idx_s <- match(names(sample.order)[idx], colnames(m6a_iqr_z))
      tmp <- m6a_iqr_z[, idx_s]
      mean_m6a_iqr[[i]] <- rowMeans(tmp)
    }
    
    
    ## with comparison 
    cmp <- list( c("P2", "P1"), c("P2", "P3"), c("P2", "P4"), c("P2", "P5"))
    mean_m6a <- unlist(mean_m6a)
    group <- rep(c("P1", "P2", "P3", "P4", "P5"), each = nrow(m6a_rt_z))
    dat <- data.frame(mean_m6a, group)
    
    g <- ggboxplot(dat, x = "group", y = "mean_m6a",color = "group")
    g <- g + stat_compare_means(comparisons = cmp) 
    g <- g + scale_color_manual(values = brewer.pal(n = 12, name = "PRGn")[c(7,2,4,9,10)])
    ggsave("../Results/Mean_m6A_levels_for_5_clusters_8030_genes_with_Pval.pdf", width = 5, height = 4.5)            
    
    mean_m6a_iqr <- unlist(mean_m6a_iqr)
    group_iqr <- rep(c("P1", "P2", "P3", "P4", "P5"), each = nrow(m6a_iqr_z)) 
    dat <- data.frame(mean_m6a_iqr, group_iqr)
    
    g <- ggboxplot(dat, x = "group_iqr", y = "mean_m6a_iqr",color = "group_iqr")
    g <- g + stat_compare_means(comparisons = cmp) 
    g <- g + scale_color_manual(values = brewer.pal(n = 12, name = "PRGn")[c(7,2,4,9,10)])
    ggsave("../Results/Mean_m6A_levels_for_5_clusters_IQR_genes_with_Pval.pdf", width = 5, height = 4.5)            
    
  }
  
  
  ################
  ## enzymes mRNA
  {

    m6a_en <- read.table("./1_enzymes/m6a_enzyme.txt", header = T)      # enzyme list  
    idx_g <- match(m6a_en$gene_name, rownames(input_rpkm))
    idx_s <- match(colnames(m6a_rt), colnames(input_rpkm))
    
    m6a_en_expr <- input_rpkm[idx_g, idx_s]
    m6a_en_expr_z <- t(scale(t(m6a_en_expr)))
    
    ## m6a enzymes mean expr (Z) per groups
    m6a_en_pz <- matrix(0, nrow(m6a_en), 5)
    for(i in 1:nrow(m6a_en))
    {
      for(j in 1:5)
      {
        idx <- sample.order == j
        idx_ss <- match(names(sample.order)[idx], colnames(m6a_en_expr_z))
        m6a_en_pz[i, j] <- mean(m6a_en_expr_z[i, idx_ss], na.rm = T) 
      }
    }
    
    rownames(m6a_en_pz) <- rownames(m6a_en_expr)
    colnames(m6a_en_pz) <- c("P1", "P2", "P3", "P4", "P5")
    
    library("gplots")
    
    hmcols<-colorRampPalette(c("blue","gray95","red"))(25)
    
    ## pdf for selected enzymes to be consistent with 
    m6a_en_s <- m6a_en_pz[c(1, 2, 4, 7, 8, 9, 11, 12), ]
    pdf("../Results/m6A_8_writer_eraser_mean_RNA_in_clusters.pdf", height = 4 , width = 8)
    heatmap.2(t(m6a_en_s), Colv=F, Rowv=F, dendrogram = "none" , scale='none', na.color = "grey", 
              trace ="none",col=hmcols, density.info="none",  
              cexCol = 1, cexRow = 0.75, margins=c(3, 5), 
              lhei = c(1.5, 3), lwid = c(2, 6), 
              keysize = 1, key.title = NULL, key.xlab = "Mean Z intensity")
    dev.off()
    
  }
  
  
  ##################
  ## enzymes protein
  {
    p_expr <- readRDS("protein_expr_51.rds")
    
    idx_g <- match(m6a_en$gene_name, rownames(p_expr), nomatch = 0)
    idx_s <- match(colnames(m6a_rt), colnames(p_expr))
    p_en_expr <- p_expr[idx_g, idx_s]  
    
    p_en_expr_z <- t(scale(t(p_en_expr)))
    
    p_en_pz <- matrix(0, nrow(p_en_expr_z), 5)
    for(i in 1:nrow(p_en_expr_z))
    {
      for(j in 1:5)
      {
        idx <- sample.order == j
        idx_ss <- match(names(sample.order)[idx], colnames(p_en_expr_z))
        p_en_pz[i, j] <- mean(p_en_expr_z[i, idx_ss], na.rm = T) 
      }
    }
    
    rownames(p_en_pz) <- rownames(p_en_expr)
    colnames(p_en_pz) <- c("P1", "P2", "P3", "P4", "P5")
    
    hmcols<-colorRampPalette(c("blue","gray95","red"))(25)
    
    p_en_s <- p_en_pz[c(1, 2, 4, 5), ]
    pdf("../Results/m6A_4_writer_eraser_mean_Protein_in_clusters.pdf", height = 4 , width = 8)
    heatmap.2(t(p_en_s), Colv=F, Rowv=F, dendrogram = "none" , scale='none', na.color = "grey", 
              trace ="none",col=hmcols, density.info="none",  
              cexCol = 1, cexRow = 0.75, margins=c(3, 5), 
              lhei = c(1.5, 3), lwid = c(2, 4), 
              keysize = 1, key.title = NULL, key.xlab = "Mean Z intensity")
    dev.off()
    
  }
  
  
  #########################
  ## with clinical features
  ## p2 vs other groups
  {
    idx_p2 <- clin_fs$clusters_5 == 2
    
    ## age
    age_c <- clin_fs$age > 65
    s1 = table(age_c[idx_p2])
    s2 = table(age_c[!idx_p2])
    p_age = fisher.test(rbind(s1, s2))$p.value
    
    ## sex 
    s1 = table(clin_fs$sex[idx_p2])
    s2 = table(clin_fs$sex[!idx_p2])
    p_sex = fisher.test(rbind(s1, s2))$p.value
    
    ## smoking
    s1 = table(clin_fs$smoking[idx_p2])
    s2 = table(clin_fs$smoking[!idx_p2])
    p_smoking = fisher.test(rbind(s1, s2))$p.value
    
    ## path stages
    idx1 <- clin_fs$pathStage == "1A" | clin_fs$pathStage == "1B"
    idx2 <- clin_fs$pathStage == "1A" | clin_fs$pathStage == "2B"
    stage <- rep("III_IV", nrow(clin_fs))
    stage[idx1] <- "I"
    stage[idx2] <- "II"
    
    s1 = table(stage[idx_p2])
    s2 = table(stage[!idx_p2])
    p_stage = fisher.test(rbind(s1, s2))$p.value
    
    ## mutation
    idx_na = is.na(clin_fs$mutation)
    mut <- clin_fs$mutation[!idx_na]
    clu <- clin_fs$clusters_5[!idx_na]
    idx_p2f <- clu == 2
    
    ## yes vs no
    mut_1 <- mut
    idx_no <- mut == "No"
    mut_1[!idx_no] = "Yes"
    s1 = table(mut_1[idx_p2f])
    s2 = table(mut_1[!idx_p2f])
    p_mut_yn = fisher.test(rbind(s1, s2))$p.value
    
    ## KRAS 
    mut_1 <- mut
    idx_no <- mut == "KRAS" 
    mut_1[!idx_no] = "No"
    s1 = table(mut_1[idx_p2f])
    s2 = table(mut_1[!idx_p2f])
    p_mut_kras = fisher.test(rbind(s1, s2))$p.value
    
    ### bar plot 
    pval <- -log10(c(p_age, p_sex, p_smoking, p_stage,  p_mut_yn, p_mut_kras))
    names(pval) <- c("Age","Sex","Smoking Staus","pathStage", "No Mutation", "KRAS mutated")
    
    pdf("../Results/m6A_P2_vs_others_assoc_with_clinical_features.pdf", width = 3, height = 4)
    par(mar = c(5, 7, 1, 1))
    barplot(pval, horiz = T, las = 1, xlab = "-log10(P value)")
    dev.off()
    
  }
  
  
}

#######################################
## m6A based clusters survival analyses
######################################
{
  library(survival)
  library(survminer)
  
  {
    survival_data <- data.frame(cluster_cat_5 = clin_fs$clusters_5, 
                                survival_time = clin_fs$surv_time,status = clin_fs$status)
    
    survival_data$survival_time[survival_data$survival_time == "No Data"] <- NA
    survival_data <- survival_data[complete.cases(survival_data),]
    survival_data$status[survival_data$status == "A"] <- 0
    survival_data$status[survival_data$status == "D"] <- 1
    survival_data$survival_time <- as.numeric(survival_data$survival_time)
    survival_data$status <- as.numeric(survival_data$status)
    survival_data$cluster_cat <- as.numeric(survival_data$cluster_cat)
    survival_data
    
    #km analysis for all 5 clusters separately:
    {
      kmsurvival <- survfit(Surv(survival_data$survival_time,survival_data$status) ~ survival_data$cluster_cat_5)
      summary(kmsurvival)
      
      png(paste0("../Results/5cluster_KM_plots.png"))
      par(mar=c(5,6,4,1)+.1)
      plot(kmsurvival, col = 1:5, lwd=4,mark.time = TRUE, xlab = "Time (years)", ylab = "Survival Probability",cex.axis=1.8,
           cex.lab = 2, main = "m6A Patient Clusters", cex.main = 2)
      legend(
        "bottomleft",
        legend=c("P1","P2","P3","P4","P5"),
        col=1:5,lty = 1,lwd = 4,
        horiz=FALSE,
        bty="n",cex = 1.8)
      dev.off()
      
      #P2 cluster vs rest, combine 1,3,4,5. 
      survival_data$cluster_cat_5[survival_data$cluster_cat_5 != 2] <- 1
      survival_data
      kmsurvival <- survfit(Surv(survival_data$survival_time,survival_data$status) ~ survival_data$cluster_cat_5)
      summary(kmsurvival)
      
      png(paste0("../Results/5_2cluster_KM_plots.png"))
      par(mar=c(5,6,4,1)+.1)
      plot(kmsurvival, col = 1:2, lwd=4,mark.time = TRUE, xlab = "Time (years)", ylab = "Survival Probability",cex.axis=1.8,
           cex.lab = 2,main = "m6A Patient Clusters",cex.main=2)
      legend(
        "bottomleft",
        legend=c("P1,3,4,5","P2"),
        col=1:2,lty = 1,lwd = 4,
        horiz=FALSE,
        bty="n",cex = 1.8)
      p_km <- survdiff(Surv(survival_data$survival_time, survival_data$status) ~ survival_data$cluster_cat_5)
      p_km
      dev.off()
    }
    

  }
}

#######################################
##  p2 vs P1,3, 4, 5 DEGs analyses
######################################
{
  ######################################################
  ### DE test beased m6a_leve high low group comparison
  ### rwa count with ERCC correction 
  {
  library(DESeq2)
  
  load("paired_63_count_ercc.Rdata")
  
  ## p2 vs others 
  idx_2 <- clin_fs$clusters_5 == 2
  
  idx_s1 <- match(clin_fs$helab_id[idx_2], colnames(input_count_norm))
  idx_s2 <- match(clin_fs$helab_id[!idx_2], colnames(input_count_norm))
  
  expr <- round(cbind(input_count_norm[, idx_s1], input_count_norm[, idx_s2]))
  type <- factor(c(rep("P2", length(idx_s1)), rep("Others", length(idx_s2))))
  
  
  dds <- DESeqDataSetFromMatrix(expr, DataFrame(type), ~type)
  
  ## load exonBygene_f
  load("exonBygene_f.RData")
  
  rowRanges(dds) <- exonBygene_f
  dds <- DESeq(dds)
  res <- results(dds, contrast=c("type", "P2","Others"))
  
  #write.csv(res, "P2_vs_Others_DESEQ2_all.csv")
  fc <- res$log2FoldChange
  pval <- res$pvalue
  fdr <- res$padj
  
  ## for pre-ranked GSEA analysis 
  stat <- fc * -log10(pval)
  names(stat) <- rownames(res)
  stat_s <- sort(stat, decreasing =  T) 
  #write.table(stat_s, paste("m6A_level_based_P2_vs_others_samples_GSEA_ranked.rnk", sep = ""), row.names = T, col.names = F, quote = F,  sep = "\t") 
  }
  
  ##########################
  ## for significantly DEGs
  {
   
    ## DEG volcano plot based on FDR
    idx_up <- fc > log2(3/2) & fdr < 0.05 & !is.na(fdr)
    idx_do <- fc < -log2(3/2) & fdr < 0.05 & !is.na(fdr)
    
    col_p <- brewer.pal(n = 10, name = "RdBu")  
    col_deg <- rep("gray", length(fc))
    col_deg[idx_up] <- col_p[2]
    col_deg[idx_do] <- col_p[8]
    sub_title = paste0("Down: N= ", sum(idx_do, na.rm = T), ";  Up: N=", sum(idx_up, na.rm = T))
    
    pdf(file = paste0("../Results/P2_vs_Others_Volcano_plot_DEGs_FDR.pdf"), width = 4, height = 4)
    par(mar = c(4, 3, 2, 0.5), mgp =c(2, 0.75, 0))
    plot(fc, -log10(fdr), col =  col_deg , pch = 16, xlab = "log2(RNA_FC)", 
         ylab = "-log10(FDR)", main = "P2 vs Others", sub = sub_title, cex = 0.75, xlim = c(-8, 8))
    abline(h = -log10(0.05), lty = 2, col = "gray")
    abline(v = c(-log2(1.5), log2(1.5)), lty = 2, col = "gray")
    dev.off()
    
    
    
    #################################
    ## functional enrichment for DEGs
    if(TRUE){
      
      library(clusterProfiler)
      library(enrichplot)
      
      ## for KEGG and BP enrichment 
      KEGG_gmt <- read.gmt("./2_GMT/c2.cp.kegg.v7.1.symbols.gmt")
      BP_gmt <- read.gmt("./2_GMT/c5.bp.v7.1.symbols.gmt")
      CC_gmt <- read.gmt("./2_GMT/c5.cc.v7.1.symbols.gmt")
      MF_gmt <- read.gmt("./2_GMT/c5.mf.v7.1.symbols.gmt")
      
      idx_up <- fc > log2(3/2) & fdr < 0.05 & !is.na(fdr)
      idx_do <- fc < -log2(3/2) & fdr < 0.05 & !is.na(fdr)
      
      gene_up <- rownames(res)[idx_up]
      gene_down <- rownames(res)[idx_do]
      
      
      ############################################
      ## DEGs KEGG and BP functional enrichment analysis 
      ## for upreguated genes 
      ## for reference 
      if(FALSE){
        
        kegg_enrich <-  enricher(gene_up, TERM2GENE = KEGG_gmt) 
        if(length(kegg_enrich$geneID)){
          g <- clusterProfiler::dotplot(kegg_enrich) + ggtitle(paste0("P2_vs_others_pval_based_DEGs_Up_KEGG"))
          ggsave(paste0("P2_vs_others_pval_based_DEGs_Up_KEGG.pdf"), width = 7, height = 4)
        }
        
        bp_enrich <-  enricher(gene_up, TERM2GENE = BP_gmt) 
        if(length(bp_enrich$geneID)){
          g <- clusterProfiler::dotplot(bp_enrich) + ggtitle(paste0("P2_vs_others_pval_based_DEGs_Up_GO_BP"))
          ggsave(paste0("P2_vs_others_pval_based_DEGs_Up_GO_BP.pdf"), width = 7, height = 4)
        }
        
        cc_enrich <-  enricher(gene_up, TERM2GENE = CC_gmt) 
        if(length(cc_enrich$geneID)){
          g <- clusterProfiler::dotplot(cc_enrich) + ggtitle(paste0("P2_vs_others_pval_based_DEGs_Up_GO_CC"))
          ggsave(paste0("P2_vs_others_pval_based_DEGs_Up_GO_CC.pdf"), width = 7, height = 4)
        }
        
        mf_enrich <-  enricher(gene_up, TERM2GENE = MF_gmt) 
        if(length(mf_enrich$geneID)){
          g <- clusterProfiler::dotplot(mf_enrich) + ggtitle(paste0("P2_vs_others_pval_based_DEGs_Up_GO_MF"))
          ggsave(paste0("P2_vs_others_pval_based_DEGs_Up_GO_MF.pdf"), width = 7, height = 4)
        }
        
        
        ## for downregulated genes
        kegg_enrich <-  enricher(gene_down, TERM2GENE = KEGG_gmt) 
        if(length(kegg_enrich$geneID)){
          g <- clusterProfiler::dotplot(kegg_enrich) + ggtitle(paste0("P2_vs_others_pval_based_DEGs_Down_KEGG"))
          ggsave(paste0("P2_vs_others_pval_based_DEGs_Down_KEGG.pdf"), width = 7, height = 4)
        }
        
        bp_enrich <-  enricher(gene_down, TERM2GENE = BP_gmt) 
        if(length(bp_enrich$geneID)){
          g <- clusterProfiler::dotplot(bp_enrich) + ggtitle(paste0("P2_vs_others_pval_based_DEGs_Down_GO_BP"))
          ggsave(paste0("P2_vs_others_pval_based_DEGs_Down_GO_BP.pdf"), width = 7, height = 4)
        }
        
        cc_enrich <-  enricher(gene_down, TERM2GENE = CC_gmt) 
        if(length(cc_enrich$geneID)){
          g <- clusterProfiler::dotplot(cc_enrich) + ggtitle(paste0("P2_vs_others_pval_based_DEGs_Down_GO_CC"))
          ggsave(paste0("P2_vs_others_pval_based_DEGs_Down_GO_CC.pdf"), width = 7, height = 4)
        }
        
        mf_enrich <-  enricher(gene_down, TERM2GENE = MF_gmt) 
        if(length(mf_enrich$geneID)){
          g <- clusterProfiler::dotplot(mf_enrich) + ggtitle(paste0("P2_vs_others_pval_based_DEGs_Down_GO_MF"))
          ggsave(paste0("P2_vs_others_pval_based_DEGs_Down_GO_MF.pdf"), width = 7, height = 4)
        }
        
      }
      
    }
    
    
  
  
}
}
