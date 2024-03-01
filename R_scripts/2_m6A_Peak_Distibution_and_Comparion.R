##########################################################################
## This script implements:
## 1. m6A peak counts, overlapping comparison between normal and tumor
## 2. peaks distribution in protein and lincRNA & peak length distribution
## 3. gene with common peaks and unique peaks
## 4. recurrence rate for m6a residuals in tumor and normal 
##########################################################################

rm(list=ls())
setwd("./data")

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
ggsave("../Results/peak_count_boxplot.pdf", width = 2.5, height = 3, units = "in")
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
ggsave("../Results/peak_pairwise_overlap_percentage_violin_plot.pdf", width = 2.5, height = 3, units = "in")

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
source("../R_scripts/functions/summarySE.R")     # barplot with error bar

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
ggsave("../Results/cat_peaks_gene_type_summary_barplot.pdf", width = 4, height = 3, units = "in")
}

############################
## peak length distribution
{
tumor_a  <- read.table("./0_merged_peaks/tumor_merge.bed")
normal_a <- read.table("./0_merged_peaks/normal_merge.bed")
tumor_o  <- read.table("./0_merged_peaks/tumor_only.bed")
normal_o <- read.table("./0_merged_peaks/normal_only.bed")
tumor_s  <- read.table("./0_merged_peaks/tumor_shared.bed")
normal_s <- read.table("./0_merged_peaks/normal_shared.bed")

tumor_ap  <- read.table("./0_merged_peaks/tumor_merge_protein.bed")
normal_ap <- read.table("./0_merged_peaks/normal_merge_protein.bed")
tumor_al  <- read.table("./0_merged_peaks/tumor_merge_lincRNA.bed")
normal_al <- read.table("./0_merged_peaks/normal_merge_lincRNA.bed")

ag <- c(length(unique(normal_a$V4)), length(unique(tumor_a$V4)))
apg <- c(length(unique(normal_ap$V4)), length(unique(tumor_ap$V4)))
alg <- c(length(unique(normal_al$V4)), length(unique(tumor_al$V4)))
all <- rbind(ag, apg, alg)

tumor_op  <- read.table("./0_merged_peaks/tumor_only_protein.bed")
normal_op <- read.table("./0_merged_peaks/normal_only_protein.bed")
tumor_ol  <- read.table("./0_merged_peaks/tumor_only_lincRNA.bed")
normal_ol <- read.table("./0_merged_peaks/normal_only_lincRNA.bed")
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


###################
# merged Guitarplot
{
## function to draw merged Guitarplot
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
        
}

## for normal
load("./0_merged_peaks/Normal_Metagene_analysis.Rdata")

ct = count
col_n = "cadetblue3"
col_t = "coral2"
comLength = c(0.136, 0.459, 0.405)
metaPeakDistributions(count, comLength, col_n, "../Results/Normal")

## for tumor
load("./0_merged_peaks/Tumor_Metagene_analysis.Rdata")

ct = count_t
metaPeakDistributions(count_t, comLength, col_t, "../Results/Tumor")

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

common_g <- read.table("./0_merged_peaks/tumor_shared.bed")
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
if(FALSE){
write.table(tumor_gl, quote = F, row.names = F, col.names = F, file = "all_gene_with_tumor_only_peak.txt")
write.table(normal_gl, quote = F, row.names = F, col.names = F, file = "all_gene_with_normal_only_peak.txt")
write.table(common_gl, quote = F, row.names = F, col.names = F, file = "all_gene_with_common_peak.txt")
write.table(m6a_g, quote = F, row.names = F, col.names = F, file = "gene_with_m6a_peak_in_tumor_or_normal.txt")

## gen 
write.table(tumor_og, quote = F, row.names = F, col.names = F, file = "gene_with_tumor_only_peak.txt")
write.table(normal_og, quote = F, row.names = F, col.names = F, file = "gene_with_normal_only_peak.txt")
write.table(common_og, quote = F, row.names = F, col.names = F, file = "gene_with_common_only_peak.txt")
}
}


#########################################################
## recurrence rate for m6a residuals in tumor and normal 
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
    ggsave(filename="../Results/m6A peak reoccurence rate_shared_new.pdf" , width = 3, height = 3)
    }
    
    ## for shared common peaks
    shared_com  <- read.table("./0_merged_peaks/shared_merge.bed")
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
    
    ggsave(filename="../Results/m6A peak reoccurence rate_unique.pdf", width = 3, height = 3)
    
    
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
    
}
