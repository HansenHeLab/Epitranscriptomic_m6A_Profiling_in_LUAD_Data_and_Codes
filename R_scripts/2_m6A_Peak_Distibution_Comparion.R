##########################################################################
## This script implements:
## 1. m6A peak counts, overlapping comparison bwteen normal and tumor
## 2. peaks distribution in protein and lincRNA & peak length distribution
## 3. gene with common peaks and unique peaks
## 4. reoccurence rate for m6a residuals in tumor and normal 
## 5. gene expression comparion among gene with common and specific peak
## 6. association between peak counts and clinical information
##########################################################################

rm(list=ls())
setwd("/work/dir/")

library("heatmap.plus")
library("gplots")
library("ggplot2")
library("matrixStats")
library("car")    ## for leveneTest

####################################################################
# m6A peak counts and overlapping comparison bwteen normal and tumor
####################################################################
{

############################    
# m6A peaks count per sample
{
peak_cnt <- read.table("peak_cnt_per_sample.txt", header = T)     # tumor45 and tumor9 have been removed
leveneTest(count ~ sample_type, data = peak_cnt)
shapiro.test(peak_cnt$count[1:10])
shapiro.test(peak_cnt$count[11:61])

## failed to meet the assumptions
wilcox.test(peak_cnt$count[1:10], peak_cnt$count[11:61], alternative = "greater")
median(peak_cnt$count[1:10])
median(peak_cnt$count[11:61])

g <- ggplot(peak_cnt, aes(x = sample_type, y = count, col = sample_type)) + geom_boxplot(outlier.shape = NA) 
g <- g + geom_point(aes(col = sample_type), position = position_jitter(width = 0.3)) 
g <- g + theme_bw() + theme(axis.line = element_line(colour = "black"),panel.grid.major = element_blank(),panel.grid.minor = element_blank(),panel.border =  element_blank(),  panel.background = element_blank()) 
g <- g + scale_color_manual( values = c("cadetblue3", "coral2")) + theme(legend.position = "none")
ggsave("peak_count_boxplot_61samples.pdf", width = 2.5, height = 3, units = "in")
}

#######################################
# for pairwise peaks overlap percentage 
{
normal_oc <- read.table("normal_pairwise_overlap_cnt.txt", header = T)
normal_perc <- normal_oc$overlaped / normal_oc$all
tumor_oc <- read.table("tumor_pairwise_overlap_cnt.txt", header = T)
tumor_perc <- tumor_oc$overlaped / tumor_oc$all

length(normal_perc)
length(tumor_perc)
overlap_perc<- c(normal_perc, tumor_perc)
group <- c(rep("Normal", length(normal_perc)), rep("Tumor", length(tumor_perc)))
dat <- data.frame(overlap_perc, group)    

shapiro.test(normal_perc)
shapiro.test(tumor_perc)
leveneTest(overlap_perc ~ group, data = dat)
wilcox.test(normal_perc, tumor_perc, alternative = "greater")

g <- ggplot(dat, aes(x = group, y = overlap_perc, col = group, alpha = 1)) + 
    geom_violin(aes( fill=group), trim=FALSE) + geom_boxplot(width=0.2) + theme_classic()
g <- g + scale_color_manual( values = c("cadetblue3", "coral2"))  + scale_fill_manual( values = c("cadetblue3", "coral2"))
g <- g + theme(legend.position = "none")
ggsave("peak_pairwise_overlap_percentage_violin_plot.pdf", width = 2.5, height = 3, units = "in")

## boxplot
perc <- list(normal_perc, tumor_perc)
col_f <- c("cadetblue3", "coral2")
pdf(file = "peak_pairwise_overlap_percentage_boxplot.pdf", height = 4, width = 3)
boxplot(perc, col = col_f)
dev.off()
}

#########################################################
# for subsmapling tumor sample to the same size of normal 
{
    peak_cmp <- read.table("subsample_10tumor_peak_count_cmp.txt", header = T)
    idx <- peak_cmp$type == "total"
    peak_cmp <- peak_cmp[!idx, ]
    
    g <- ggplot(peak_cmp, aes(x = type, y = peak_count, fill = type)) + geom_boxplot(outlier.shape = NA)  +  theme_bw() + theme(axis.line = element_line(colour = "black"),panel.grid.major = element_blank(),panel.grid.minor = element_blank(),panel.border =  element_blank(),  panel.background = element_blank()) 
    g <- g + scale_fill_manual( values = c("cadetblue3", "grey", "coral2"))
    ggsave("peak_cmp_boxplot_subsampling1000.pdf", width = 10, height = 8, units = "cm")
    
    peak_cmp_se <- summarySE(data = peak_cmp, measurevar = "peak_count", groupvars = "type" )
    
    g <- ggplot(peak_cmp_se, aes(x = type, y = peak_count, fill = type)) + geom_bar(position = position_dodge(), stat = "identity") 
    g <- g +  geom_errorbar(aes(ymin = peak_count - sd, ymax = peak_count + sd), width = .2, position=position_dodge(.9))
    g <- g + theme_bw() + theme(axis.line = element_line(colour = "black"),panel.grid.major = element_blank(),panel.grid.minor = element_blank(),panel.border =  element_blank(),  panel.background = element_blank()) 
    g <- g + scale_fill_manual( values = c("cadetblue3", "grey", "coral2"))  + theme(legend.position = "none")
    
    ggsave("peak_cmp_barplot_subsampling1000.pdf", width = 3, height = 1.6, units = "in")
}

}


######################################################################
# peaks distribution in protein and lincRNA & peak length distribution
######################################################################
{
#########################################################
## look at the pc gene, lincRNA and other genes with m6A
{
gene_dis <- read.table("cat_peaks_gene_type_summary_barplot.txt", header = T, sep = "\t")
source("/Users/Yong/Yong/R/functions/summarySE.R")     # barplot with error bar

mean(gene_dis$gene_count[1:10])
mean(gene_dis$gene_count[11:61])
wilcox.test(gene_dis$gene_count[gene_dis$class == 1],gene_dis$gene_count[gene_dis$class == 2])
wilcox.test(gene_dis$gene_count[gene_dis$class == 3],gene_dis$gene_count[gene_dis$class == 4])
wilcox.test(gene_dis$gene_count[gene_dis$class == 5],gene_dis$gene_count[gene_dis$class == 6])

## summarySE function 
summarySE <- function(data=NULL, measurevar, groupvars=NULL, na.rm=FALSE,
                      conf.interval=.95, .drop=TRUE) {
    library(plyr)
    
    # New version of length which can handle NA's: if na.rm==T, don't count them
    length2 <- function (x, na.rm=FALSE) {
        if (na.rm) sum(!is.na(x))
        else       length(x)
    }
    
    # This does the summary. For each group's data frame, return a vector with
    # N, mean, and sd
    datac <- ddply(data, groupvars, .drop=.drop,
                   .fun = function(xx, col) {
                       c(N    = length2(xx[[col]], na.rm=na.rm),
                         mean = mean   (xx[[col]], na.rm=na.rm),
                         sd   = sd     (xx[[col]], na.rm=na.rm)
                       )
                   },
                   measurevar
    )
    
    # Rename the "mean" column    
    datac <- rename(datac, c("mean" = measurevar))
    
    datac$se <- datac$sd / sqrt(datac$N)  # Calculate standard error of the mean
    
    # Confidence interval multiplier for standard error
    # Calculate t-statistic for confidence interval: 
    # e.g., if conf.interval is .95, use .975 (above/below), and use df=N-1
    ciMult <- qt(conf.interval/2 + .5, datac$N-1)
    datac$ci <- datac$se * ciMult
    
    return(datac)
}
gene_dis_se <- summarySE(data = gene_dis, measurevar = "gene_count", groupvars=c("class"))
gene_dis_se$class <- as.character(gene_dis_se$class)

write.csv(gene_dis_se, file = "cat_peaks_gene_type_summary_stat.csv")

g <- ggplot(gene_dis_se, aes(x = class, y = gene_count, fill = class)) + geom_bar(position = position_dodge(), stat = "identity") 
g <- g +  geom_errorbar(aes(ymin = gene_count - sd, ymax = gene_count + sd), width = .2, position=position_dodge(.9))
g <- g + theme_bw() + theme(axis.line = element_line(colour = "black"),panel.grid.major = element_blank(),panel.grid.minor = element_blank(),panel.border =  element_blank(),  panel.background = element_blank()) 
g <- g + scale_fill_manual(values = rep(c("cadetblue3", "coral2"), 3)) + ylim(0, 6000)
ggsave("cat_peaks_gene_type_summary_barplot.pdf", width = 4, height = 3, units = "in")
}

############################
## peak length distribution
{
tumor_a  <- read.table("/Users/Yong/Yong/m6a_profiling/2_peak_calling/1_peaks/merge/bed30_61/0_merged_peaks/tumor_merge.bed")
normal_a <- read.table("/Users/Yong/Yong/m6a_profiling/2_peak_calling/1_peaks/merge/bed30_61/0_merged_peaks/normal_merge.bed")
tumor_o  <- read.table("/Users/Yong/Yong/m6a_profiling/2_peak_calling/1_peaks/merge/bed30_61/0_merged_peaks/tumor_only.bed")
normal_o <- read.table("/Users/Yong/Yong/m6a_profiling/2_peak_calling/1_peaks/merge/bed30_61/0_merged_peaks/normal_only.bed")
tumor_s  <- read.table("/Users/Yong/Yong/m6a_profiling/2_peak_calling/1_peaks/merge/bed30_61/0_merged_peaks/tumor_shared.bed")
normal_s <- read.table("/Users/Yong/Yong/m6a_profiling/2_peak_calling/1_peaks/merge/bed30_61/0_merged_peaks/normal_shared.bed")

tumor_ap  <- read.table("/Users/Yong/Yong/m6a_profiling/2_peak_calling/1_peaks/merge/bed30_61/0_merged_peaks/tumor_merge_protein.bed")
normal_ap <- read.table("/Users/Yong/Yong/m6a_profiling/2_peak_calling/1_peaks/merge/bed30_61/0_merged_peaks/normal_merge_protein.bed")
tumor_al  <- read.table("/Users/Yong/Yong/m6a_profiling/2_peak_calling/1_peaks/merge/bed30_61/0_merged_peaks/tumor_merge_lincRNA.bed")
normal_al <- read.table("/Users/Yong/Yong/m6a_profiling/2_peak_calling/1_peaks/merge/bed30_61/0_merged_peaks/normal_merge_lincRNA.bed")

ag <- c(length(unique(normal_a$V4)), length(unique(tumor_a$V4)))
apg <- c(length(unique(normal_ap$V4)), length(unique(tumor_ap$V4)))
alg <- c(length(unique(normal_al$V4)), length(unique(tumor_al$V4)))
all <- rbind(ag, apg, alg)

tumor_op  <- read.table("/Users/Yong/Yong/m6a_profiling/2_peak_calling/1_peaks/merge/bed30_61/0_merged_peaks/tumor_only_protein.bed")
normal_op <- read.table("/Users/Yong/Yong/m6a_profiling/2_peak_calling/1_peaks/merge/bed30_61/0_merged_peaks/normal_only_protein.bed")
tumor_ol  <- read.table("/Users/Yong/Yong/m6a_profiling/2_peak_calling/1_peaks/merge/bed30_61/0_merged_peaks/tumor_only_lincRNA.bed")
normal_ol <- read.table("/Users/Yong/Yong/m6a_profiling/2_peak_calling/1_peaks/merge/bed30_61/0_merged_peaks/normal_only_lincRNA.bed")
og <- c(length(unique(normal_o$V4)), length(unique(tumor_o$V4)))
opg <- c(length(unique(normal_op$V4)), length(unique(tumor_op$V4)))
olg <- c(length(unique(normal_ol$V4)), length(unique(tumor_ol$V4)))
only <- rbind(og, opg, olg)


## merged peaks length distribution
normal_len <- normal_a$V3-normal_a$V2
hist(log10(normal_len))
sum(normal_len<500)/length(normal_len)

tumor_len <- tumor_a$V3-tumor_a$V2
hist(log10(tumor_len))
sum(tumor_len<500)/length(tumor_len)
}

######################################################################
# peak summit distribution in protein coding 5 non-overlapping region  
{
        file <- list.files(path = ".", pattern = "stat");
        
        count <- vector()
        dens <- vector()
        class <- vector()
        group <- vector()
        sample <- vector()
        N <- length(file)
        
        ref_den <- vector()
        
        for (i in 1: N)
        {
            tmp <-read.table(file[i], skip=1)
            
            tmp[3,2:3] = tmp[3, 2:3] - tmp[4, 2:3]
            tmp[3,4] = tmp[3, 2]/tmp[3, 3]
            
            tmp[5,2:3] = colSums(tmp[4:6,2:3])
            tmp[5, 4] = tmp[5, 2]/tmp[5, 3]
            
            tmp[7,2:3] = tmp[7, 2:3] - tmp[6, 2:3]
            tmp[7,4] = tmp[7, 2]/tmp[7, 3]
            
            tmp_1 = tmp[-c(4,6), ];
            
            count <-  c(count, tmp_1$V2)
            dens <- c(dens, tmp_1$V4)
            class <- c(class, as.character(tmp_1$V1))
            group <- c(group, rep(substr(file[i], 1, 12), 5))
            sample[i] <- substr(file[i], 1, 12)
            
            ref_den[i] <- sum(tmp$V2)/sum(tmp$V3)
            
        }
        
        count_m <- matrix(count, 5, N)
        dens_lm <- -1/log10(matrix(dens, 5, N))
        enrich <- matrix(dens/(3*ref_den), 5, N)
        
        rownames(count_m) <- class[1:5]
        rownames(dens_lm) <- class[1:5]
        rownames(enrich) <- class[1:5]
        colnames(count_m) <- sample
        colnames(dens_lm) <- sample
        colnames(enrich) <- sample
        
        #write.csv(count_m, file="count_dis.csv")
        #write.csv(dens_lm, file="log10_dens.csv")
        write.csv(enrich, file="enrich.csv")
    }    

###################
# merged Guitarplot
{
## function to draw merged Guitarplot
metaPeakDistributions <-  function(ct, comLength, color, name) 
{
        library("ggplot2")
        source("/Users/Yong/Yong/R/functions/summarySE.R")     # barplot with error bar
        
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
        
}

## for normal
load("~/Yong/m6a_profiling/2_peak_calling/6_metagene/peak_distribution/Normal_Metagene_analysis.Rdata")

ct = count
col_n = "cadetblue3"
col_t = "coral2"
comLength = c(0.136, 0.459, 0.405)
metaPeakDistributions(count, comLength, col_n, "Normal")

## for tumor
load("~/Yong/m6a_profiling/2_peak_calling/6_metagene/peak_distribution/Tumor_25S_Metagene_analysis.Rdata")
t1 <- count
load("~/Yong/m6a_profiling/2_peak_calling/6_metagene/peak_distribution/Tumor_26S_Metagene_analysis.Rdata")
t2 <- count
count_t <- rbind(t1, t2)
metaPeakDistributions(count_t, comLength, col_t, "Tumor")

}
       
}


##########################################
# gene with common peaks and unique peaks
##########################################
{
# genelist with merged gene
# because peak merge will lead to multiple gene assign 

gene.split <- function(gene){
    N <- length(gene)
    gene_l <- gene[1]
    
    for (i in 2:N)
    {
        gene_l <- paste(gene_l, gene[i], sep = ",")
    }
    
    gene_list <- strsplit(gene_l, split = ",")
    gene_v <- unique(as.character(gene_list[[1]]))
    return(gene_v)
}

common_g <- read.table("/Users/Yong/Yong/m6a_profiling/2_peak_calling/1_peaks/merge/bed30_61/0_merged_peaks/tumor_shared.bed")
common_gl <- gene.split(unique(common_g$V4))
tumor_gl <- gene.split(unique(tumor_o$V4))
normal_gl <- gene.split(unique(normal_o$V4))

# venn plot for gene with common or unique peaks
a <- length(common_gl)
b <- length(tumor_gl)
c <- length(normal_gl)
ab <- sum(!is.na(match(tumor_gl, common_gl)))
ab_g <- as.character(tumor_gl[!is.na(match(tumor_gl, common_gl))])
ac <- sum(!is.na(match(normal_gl, common_gl)))
ac_g <- as.character(normal_gl[!is.na(match(normal_gl, common_gl))])
bc <- sum(!is.na(match(normal_gl, tumor_gl)))
bc_g <- as.character(normal_gl[!is.na(match(normal_gl, tumor_gl))])
abc <- sum(!is.na(match(bc_g, common_gl)))
abc_g <- as.character(bc_g[!is.na(match(bc_g, common_gl))])

library(VennDiagram)
grid.newpage()
venn.plot <- draw.triple.venn(a, b, c, ab, bc, ac, abc, fill = c( "grey", "coral2", "cadetblue3"), euler.d = T, scaled = T)

# all gene with m6A peak
m6a_g <- unique(c(as.character(common_gl), as.character(tumor_gl), as.character(normal_gl)))

# gene only with tumor only peak 
over_t <- unique(c(ab_g, bc_g))
idx_l <- match(over_t, tumor_gl)
tumor_og <- tumor_gl[-idx_l]

# gene only with normal only peak 
over_t <- unique(c(ac_g, bc_g))
idx_l <- match(over_t, normal_gl)
normal_og <- normal_gl[-idx_l]

# gene only with common only peak 
over_t <- unique(c(ac_g, ab_g))
idx_l <- match(over_t, common_gl)
common_og <- common_gl[-idx_l]

# for Gprofiler
## all genes
write.table(tumor_gl, quote = F, row.names = F, col.names = F, file = "all_gene_with_tumor_only_peak.txt")
write.table(normal_gl, quote = F, row.names = F, col.names = F, file = "all_gene_with_normal_only_peak.txt")
write.table(common_gl, quote = F, row.names = F, col.names = F, file = "all_gene_with_common_peak.txt")
write.table(m6a_g, quote = F, row.names = F, col.names = F, file = "gene_with_m6a_peak_in_tumor_or_normal.txt")

## gen 
write.table(tumor_og, quote = F, row.names = F, col.names = F, file = "gene_with_tumor_only_peak.txt")
write.table(normal_og, quote = F, row.names = F, col.names = F, file = "gene_with_normal_only_peak.txt")
write.table(common_og, quote = F, row.names = F, col.names = F, file = "gene_with_common_only_peak.txt")
}


#########################################################
## reoccurence rate for m6a residuals in tumor and normal 
#########################################################
{
    ## for common peaks tumor / nomal separately
    {
    samples <- strsplit(as.character(tumor_s$V5), ",")
    N <- length(samples)
    tumor_peak_ot <- vector()
    for (i in 1: N)
    {
        tumor_peak_ot[i] <- length(unique(samples[[i]]))
    }
    tumor_peak_of <- tumor_peak_ot/51    # 
    
    samples <- strsplit(as.character(normal_s$V5), ",")
    M <- length(samples)
    normal_peak_ot <- vector()
    for (i in 1: M)
    {
        normal_peak_ot[i] <- length(unique(samples[[i]]))
    }
    normal_peak_of <- normal_peak_ot/10
    
    peak_of <- c(tumor_peak_of, normal_peak_of)
    peak_g <- c(rep("tumor", N), rep("normal",M)) 
    dat <- data.frame(peak_of, peak_g)
    
    g <- ggplot(dat, aes(peak_of,fill = peak_g)) + geom_density(alpha = 0.9)  + theme(legend.direction = 'horizontal', legend.position = 'top')+  theme_bw() + theme(axis.line = element_line(colour = "black"),panel.grid.major = element_blank(),panel.grid.minor = element_blank(),panel.border =  element_blank(),panel.background = element_blank()) 
    g <- g + scale_fill_manual( values = c("cadetblue3", "coral2")) 
    g <- g + scale_x_continuous(expand = c(0, 0), name = "m6A peak reoccurence rate") + scale_y_continuous(expand = c(0, 0)) + theme(legend.position = "none")
    ggsave(filename="m6A peak reoccurence rate_shared_new.pdf" , width = 3, height = 3)
    }
    
    
    ## for shared common peaks
    shared_com  <- read.table("/Users/Yong/Yong/m6a_profiling/2_peak_calling/1_peaks/merge/bed30_61/0_merged_peaks/shared_merge.bed")
    samples <- strsplit(as.character(shared_com$V5), ",")
    N <- length(samples)
    shared_com_ot <- vector()
    for (i in 1: N)
    {
        shared_com_ot[i] <- length(unique(samples[[i]]))
    }
    idx_r <- shared_com_ot == 2 | shared_com_ot == 3
    shared_com_ot <- shared_com_ot[!idx_r]
    shared_com_of <- shared_com_ot/61     ## merge normal and tumor 
    
    
    # for unique peak
    {
    samples <- strsplit(as.character(tumor_o$V5), ",")
    N <- length(samples)
    tumor_peak_ot <- vector()
    for (i in 1: N)
    {
        tumor_peak_ot[i] <- length(unique(samples[[i]]))
    }
    idx_r <- tumor_peak_ot == 1
    tumor_peak_ot <- tumor_peak_ot[!idx_r]
    tumor_peak_of <- tumor_peak_ot/51    # 
    
    samples <- strsplit(as.character(normal_o$V5), ",")
    M <- length(samples)
    normal_peak_ot <- vector()
    for (i in 1: M)
    {
        normal_peak_ot[i] <- length(unique(samples[[i]]))
    }
    idx_r <- normal_peak_ot == 1
    normal_peak_ot <- normal_peak_ot[!idx_r]
    normal_peak_of <- normal_peak_ot/10
    
    peak_of <- c(tumor_peak_of, normal_peak_of)
    peak_g <- c(rep("tumor", length(tumor_peak_of)), rep("normal",length(normal_peak_of))) 
    dat <- data.frame(peak_of, peak_g)
    
    g <- ggplot(dat, aes(peak_of, fill = peak_g)) + geom_density(alpha = 0.9)  + theme(legend.direction = 'horizontal', legend.position = 'top')+  theme_bw() + theme(axis.line = element_line(colour = "black"),panel.grid.major = element_blank(),panel.grid.minor = element_blank(),panel.border =  element_blank(),                                                                                                                                                                 panel.background = element_blank()) 
    g <- g + scale_fill_manual( values = c("cadetblue3", "coral2"))
    g <- g + scale_x_continuous(expand = c(0, 0), name = "m6A peak reoccurence rate") + scale_y_continuous(expand = c(0, 0)) + theme(legend.position = "none")
    
    ggsave(filename="m6A peak reoccurence rate_unique.pdf", width = 3, height = 3)
    
    
    ### percentage
    tumor_peak_ot_p <- table(tumor_peak_ot) / length(tumor_peak_of) 
    normal_peak_ot_p <-  table(normal_peak_ot) / length(normal_peak_of) 
    shared_com_of_p <- table(shared_com_ot) / length(shared_com_of) 
        
    p_cdf <- function(p_v)
    {
        L <- length(p_v)
        cdf <- vector()
        for(i in 1:L )
        {
            cdf[i] <- sum(p_v[1:i])
        }
        return(cdf)
    }
     
    tumor_peak_ot_p_cdf <- p_cdf(tumor_peak_ot_p)
    normal_peak_ot_p_cdf <- p_cdf(normal_peak_ot_p)
    shared_com_of_cdf <- p_cdf(shared_com_of_p)
    
    }
    
    ## merge uniuqe peaks and common peaks
    peak_of <- c(tumor_peak_of, normal_peak_of, shared_com_of)
    peak_g <- c(rep("tumor", length(tumor_peak_of)), rep("normal",length(normal_peak_of)), rep("common",length(shared_com_of))) 
    dat <- data.frame(peak_of, peak_g)
    
    g <- ggplot(dat, aes(peak_of, fill = peak_g)) + geom_density(alpha = 0.9)  + theme(legend.direction = 'horizontal', legend.position = 'top')+  theme_bw() + theme(axis.line = element_line(colour = "black"),panel.grid.major = element_blank(),panel.grid.minor = element_blank(),panel.border =  element_blank(),                                                                                                                                                                      panel.background = element_blank()) 
    g <- g + scale_fill_manual( values = c("#C79C8D", "cadetblue3", "coral2"))
    g <- g + scale_x_continuous(expand = c(0, 0), name = "m6A peak reoccurence rate") + scale_y_continuous(expand = c(0, 0)) + theme(legend.position = "none")
    
    ggsave(filename="m6A_peak_reoccurence_rate_for_3types_genes.pdf", width = 3, height = 3)
    
}


#######################################################################
#### gene expression comparion among gene with common and specific peak 
#######################################################################
{

####################
## all related genes
{
load("/Users/Yong/Yong/m6a_profiling/1_mapping_quant/5_RNA-seq/htseq/count/paired_61_rpkm.Rdata")

og <- c(normal_og, tumor_og, common_og)
l <- length(og)
l1 <- length(normal_og)
l2 <- length(tumor_og)
l3 <- length(common_og)

idx_og <- match(og, rownames(input_rpkm))
og_expr <- input_rpkm[idx_og, ]
og_nm <- log2(rowMeans(og_expr[, 1:10]) + 1)
og_tm <- log2(rowMeans(og_expr[, 11:61]) +1)

og_ntm <- c(og_nm, og_tm)
og_ntm_g <- as.character(c(rep(3, l1), rep(5, l2), rep(1, l3), rep(4, l1), rep(6, l2), rep(2, l3)))

dat <- data.frame(og_ntm, og_ntm_g)
g <- ggplot(dat, aes(x = og_ntm_g, y = og_ntm, fill = og_ntm_g, col = og_ntm_g, alpha = 1)) + geom_boxplot(width=0.2, outlier.shape=NA) +  geom_violin(trim=FALSE) 
g <- g + theme_bw()  + theme(axis.line = element_line(colour = "black"),panel.grid.major = element_blank(),panel.grid.minor = element_blank(),panel.border =  element_blank(),  panel.background = element_blank()) 
g <- g + scale_fill_manual( values = rep(c("cadetblue3", "coral2"), 3)) + scale_color_manual( values =  rep(c("cadetblue3", "coral2"), 3)) + ylim(0, 8)
g <- g + theme(legend.position = "none")
ggsave("gene_expr_tumor_normal_unique_or_common_violin.pdf",  width = 5, height = 3, units = "in")

g <- ggplot(dat, aes(x = og_ntm_g, y = og_ntm, fill = og_ntm_g)) + geom_boxplot(outlier.shape=NA)
g <- g + theme_bw()  + theme(axis.line = element_line(colour = "black"),panel.grid.major = element_blank(),panel.grid.minor = element_blank(),panel.border =  element_blank(),  panel.background = element_blank()) 
g <- g + scale_fill_manual( values = rep(c("cadetblue3", "coral2"), 3)) + scale_color_manual( values =  rep(c("cadetblue3", "coral2"), 3)) + ylim(0, 8)
g <- g + theme(legend.position = "none")
ggsave("gene_expr_tumor_normal_unique_or_common_boxplot.pdf",  width = 4, height = 3, units = "in")

## wilcox.test 
wilcox.test(og_nm[1:l1], og_tm[1:l1], alternative = "greater")     # normal only
wilcox.test(og_nm[(l1+1) : (l-l3)], og_tm[(l1+1) : (l-l3)], alternative = "less")    # tumor only
wilcox.test(og_nm[(l1+l2+1) : l], og_tm[(l1+l2+1) : l])        # shared
}

########################################################################
## by DESeq2 results  : overall gene expression between normal and tumor 
{
de_res <- read.csv("/Users/Yong/Yong/m6a_profiling/1_mapping_quant/5_RNA-seq/htseq/count/61_sample_Input_DESeq_DE.csv")
idx_de <- match(og, de_res$X)       ## matched gene list

## by own test
og_pad <- de_res$padj [idx_de]
og_fd <- de_res$log2FoldChange[idx_de]

sum(is.na(og_fd))           # do not include "NA" and "Infinite"
idx_p <- og_pad < 0.05       #with FDR < 0.05 
idx_fd <- abs(og_fd) > 1         # 
idx_de <- idx_p & idx_fd

de <- idx_de
n_p <- table(de[1:l1])      # for normal
t_p <- table(de[(l1+1) : (l-l3)])
c_p <- table(de[(l1+l2+1) : l])
cmb_p <- rbind(c_p, n_p, t_p)

## foldchange dis
fd_list <- list(og_fd[1:l1], og_fd[(l1+1) : (l-l3)], og_fd[(l1+l2+1) : l])
names(fd_list) <- c("normal_only", "tumor_only", "Common")

boxplot(fd_list, outline = F)
abline(h = c(-1, 1), col = "red", lty = 2)

# non-de de genes
normal_og_nde <- normal_og[!idx_de[1:l1]]
normal_og_nde <- normal_og_nde[!is.na(normal_og_nde)]
tumor_og_nde <- tumor_og[!idx_de[(l1+1):(l-l3)]]
tumor_og_nde <- tumor_og_nde[!is.na(tumor_og_nde)]
common_og_nde <- common_og[!idx_de[(l1+l2+1) : l]]
common_og_nde <- common_og_nde[!is.na(common_og_nde)]

write.table(tumor_og_nde, quote = F, row.names = F, col.names = F, file = "gene_with_tumor_only_peak_nde.txt")
write.table(normal_og_nde, quote = F, row.names = F, col.names = F, file = "gene_with_normal_only_peak_nde.txt")
write.table(common_og_nde, quote = F, row.names = F, col.names = F, file = "gene_with_common_only_peak_nde.txt")
}

###################################################################################
## split each group of genes into group with called peaks and without called peaks
{
        ## for tumor only 
        {
            gene_to <- read.table("/Users/Yong/Yong/m6a_profiling/2_peak_calling/1_peaks/merge/bed30_61/0_merged_peaks/tumor_only.bed", as.is = T)
            eg_g <- tumor_og
            
            idx_g <- match(eg_g, gene_to$V4)
            idx_gs <- idx_g[!is.na(idx_g)]
            eg_gs <- eg_g[!is.na(idx_g)]
            gene_tos <- gene_to[idx_gs, ]
            gene_tos_sl <- strsplit(gene_tos$V5, ",")       ## samples burdern tumor only called peaks 
            
            ano_diff_p <- vector()
            
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
                res.aov <- aov(expr ~ group, data = dat)
                # Summary of the analysis
                
                ano_diff_p[i] <-  summary(res.aov)[[1]]["Pr(>F)"][[1]][1]
                
            }
            ano_diff_padj <- p.adjust(ano_diff_p, method = "BH")
            sum(ano_diff_padj < 0.05, na.rm = T)
            sum(ano_diff_padj >= 0.05, na.rm = T)
        }
        
        ## for normal only 
        {
            gene_to <- read.table("/Users/Yong/Yong/m6a_profiling/2_peak_calling/1_peaks/merge/bed30_61/0_merged_peaks/normal_only.bed", as.is = T)
            eg_g <- normal_og
            
            idx_g <- match(eg_g, gene_to$V4)
            idx_gs <- idx_g[!is.na(idx_g)]
            eg_gs <- eg_g[!is.na(idx_g)]
            gene_tos <- gene_to[idx_gs, ]
            gene_tos_sl <- strsplit(gene_tos$V5, ",")       ## samples burdern normal only called peaks 
            
            ano_diff_p <- vector()
            
            L <- length(eg_gs)
            expr_no_mean <- matrix(0, L, 3)
            
            for(i in 1:L)
            {
                idx_gg <- match(eg_gs[i], rownames(input_rpkm))
                
                expr_t <- input_rpkm[idx_gg, 11:61]       
                
                idx_np <- match(gene_tos_sl[[i]], colnames(input_rpkm))          ## tumor samples with tumor only methylated genes
                expr_np <- input_rpkm[idx_gg, idx_np]            
                
                expr_nnp <- input_rpkm[idx_gg, -c(11:61, idx_np)]
                
                expr_no_mean[i, ] <- c(log2(mean(expr_nnp + 1)), log2(mean(expr_np + 1)), log2(mean(expr_t + 1)))
                
                
                expr <- c(expr_np, expr_nnp, expr_t)
                group <- c(rep(1, length(expr_np)), rep(2, length(expr_nnp)), rep(3, length(expr_t)))
                
                dat <- data.frame(expr, group)
                res.aov <- aov(expr ~ group, data = dat)
                # Summary of the analysis
                
                ano_diff_p[i] <-  summary(res.aov)[[1]]["Pr(>F)"][[1]][1]
                
            }
            ano_diff_padj <- p.adjust(ano_diff_p, method = "BH")
            sum(ano_diff_padj < 0.05, na.rm = T)
            sum(ano_diff_padj >= 0.05, na.rm = T)
        }
        
        ## for common peaks
        {
            gene_to <- read.table("/Users/Yong/Yong/m6a_profiling/2_peak_calling/1_peaks/merge/bed30_61/0_merged_peaks/shared_merge.bed", as.is = T)
            eg_g <- common_og
            
            idx_g <- match(eg_g, gene_to$V4)
            idx_gs <- idx_g[!is.na(idx_g)]
            eg_gs <- eg_g[!is.na(idx_g)]
            gene_tos <- gene_to[idx_gs, ]
            gene_tos_sl <- strsplit(gene_tos$V5, ",")       ## samples burdern tumor only called peaks 
            
            ano_diff_p <- vector()
            
            L <- length(eg_gs)
            
            expr_com_mean <- matrix(0, L, 4)
            
            for(i in 1:L)
            {
                idx_gg <- match(eg_gs[i], rownames(input_rpkm))
                
                ## normal
                idx_np <- match(gene_tos_sl[[i]], colnames(input_rpkm)[1:10], nomatch = 0)          ## tumor samples with tumor only methylated genes
                expr_np <- input_rpkm[idx_gg, idx_np]            
                expr_nnp <- input_rpkm[idx_gg, -c(11:61, idx_np)]
                
                ## tumor
                idx_tp <- match(gene_tos_sl[[i]], colnames(input_rpkm)[11:61], nomatch = 0)          ## tumor samples with tumor only methylated genes
                expr_tp <- input_rpkm[idx_gg, 11:61][idx_tp]            
                
                expr_tnp <- input_rpkm[idx_gg, -c(1:10)][-idx_tp]
                
                expr_com_mean[i, ] <- c(log2(mean(expr_nnp + 1)), log2(mean(expr_np + 1)), log2(mean(expr_tnp + 1)), log2(mean(expr_tp + 1)))
                
                
                expr <- c(expr_nnp, expr_np, expr_tnp, expr_tp)
                group <- c(rep(1, length(expr_nnp)), rep(2, length(expr_np)) ,  rep(3, length(expr_tnp)), rep(4, length(expr_tp)))
                
                dat <- data.frame(expr, group)
                res.aov <- aov(expr ~ group, data = dat)
                # Summary of the analysis
                
                ano_diff_p[i] <-  summary(res.aov)[[1]]["Pr(>F)"][[1]][1]
                
            }
            ano_diff_padj <- p.adjust(ano_diff_p, method = "BH")
            sum(ano_diff_padj < 0.05, na.rm = T)
            sum(ano_diff_padj >= 0.05, na.rm = T)
        }
        
        ##PLot with splited subgroups
        {
            og_ntm <- c(c(expr_com_mean), c(expr_no_mean), c(expr_to_mean))
            og_ntm_g <- as.character(c(rep(1:4, each = nrow(expr_com_mean)), rep(5:7, each = nrow(expr_no_mean)), rep(8:10, each = nrow(expr_to_mean))))
            col_f <- c("cadetblue3", "cadetblue3", "coral2", "coral2", "cadetblue3", "cadetblue3", "coral2", "cadetblue3", "coral2", "coral2")
            lty_b <- c(2, 1, 2, 1, 2, 1, 2, 2, 2, 1 ) 
            
            dat <- data.frame(og_ntm, og_ntm_g)
            
            g <- ggplot(dat, aes(x = og_ntm_g, y = og_ntm, fill = og_ntm_g)) + geom_boxplot(outlier.shape=NA, lty = lty_b)
            g <- g + theme_bw()  + theme(axis.line = element_line(colour = "black"),panel.grid.major = element_blank(),panel.grid.minor = element_blank(),panel.border =  element_blank(),  panel.background = element_blank()) 
            g <- g + scale_fill_manual( values = col_f) + scale_color_manual( values =  col_f) + ylim(0, 8)
            g <- g + theme(legend.position = "none")
            ggsave("gene_expr_tumor_normal_unique_or_common_boxplot_splitgroups.pdf",  width = 4, height = 3, units = "in")
        }
    }

}
    

###########################################################
## association between peak counts and clinical information
###########################################################
{
peak_cnt_t <- peak_cnt[11:61, ]
idx <- order(peak_cnt_t$count)
peak_cnt_ts <- peak_cnt_t[idx, ]

s_info <- read.table("/Users/Yong/Yong/m6a_profiling/2_peak_calling/1_peaks/merge/sample_info.txt", header = T)  
idx_s <- match(peak_cnt_ts$sample, s_info$helab_id)
s_infos <- s_info[idx_s, ]

count_s <- peak_cnt_ts$count
names(count_s) <- s_infos$pathstage
barplot(count_s, las = 2, cex.names = 0.8)

library("s2dverification")       # for color bar

## sex
N = length(idx_s)
lims <- 1:52
sex_c <- rep("pink", N)                           # for Female
idx_t <- !is.na(match(s_infos$sex, "M"))          # for Male
sex_c[idx_t] <- "blue"
ColorBar(lims, sex_c, vertical = F, draw_ticks = F, draw_separators = T)
h <- table(sex_c[27:51])
l <- table(sex_c[1:25])
sex_tt <- cbind(h, l)
chisq.test(sex_tt)

## smoking history 
smoking_c <- rep("green", N)
idx_t <- !is.na(match(s_infos$smoking, "Current"))
smoking_c[idx_t] <- "red"
idx_t <- !is.na(match(s_infos$smoking, "Ex-Smoker"))
smoking_c[idx_t] <- "orange"
ColorBar(lims, smoking_c, vertical = F, draw_ticks = F, draw_separators = T)
h <- table(smoking_c[27:51])
l <- table(smoking_c[1:25])
smoking_tt <- cbind(h, l)
chisq.test(smoking_tt)

## xenograph outcome
xeno_c <- rep("gray", N)
idx_t <- !is.na(match(s_infos$xeno, "Y"))
xeno_c[idx_t] <- "red"
idx_t <- !is.na(match(s_infos$xeno, "N"))
xeno_c[idx_t] <- "blue"
ColorBar(lims, xeno_c, vertical = F, draw_ticks = F, draw_separators = T)
h <- table(xeno_c[27:51])
l <- table(xeno_c[1:25])
xeno_tt <- cbind(h, l)
chisq.test(xeno_tt)
}