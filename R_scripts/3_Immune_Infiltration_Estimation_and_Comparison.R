########################################################################################
## This script implements:
## 1. immume infiltration estimation comparison in UHN PLCME
## 2. immume infiltration estimation comparison in TCGA LUAD
########################################################################################

rm(list = ls())
setwd("./data")

library(ggplot2)
library(gplots)

############################################
## immume infiltration estimation comparison
############################################
{
###################################################
## read-in immune celss estimated fraction by TIMER
{
IN <- read.csv(check.names = F, "immuneEstimation_PLCME.csv")
rownames(IN) <- IN$sampleID
frac<- as.matrix(IN[, 2: 7])

}

###############################################
## test difference between the tumor and normal
{
p_val <- vector()
L <-  ncol(frac)
for (i in 1:L)
{
    tmp <- wilcox.test(frac[1:10, i], frac[11:61, i])
    p_val[i] <- tmp$p.value
}
idx_s <- p_val  < 1  #< 0.05              ### pvalue cutoff 0.01

## draw  boxplot
frac_v <- as.vector(frac)
g_1 <- rep(c(rep("normal",  10), rep("tumor", 51)), 6)
g_2 <- rep(colnames(frac), each = 61)

frac_go <- paste(g_2, g_1, sep = "_")  
frac_gou <- unique(frac_go)
# remian sort order 
frac_gf <- as.factor(frac_go)
frac_g <- factor(frac_gf, level = frac_gou)

## significant ones with color 
idx_ss = which(idx_s == TRUE)
col_s <- rep("gray", 12)
col_s[idx_ss*2] <- "coral2"        ## Tumor
col_s[idx_ss*2-1] <- "cadetblue3"  ## normal 

dat_g <- data.frame(frac_v, frac_g)

## for violin plots
g <- ggplot(dat_g, aes(x = frac_g, y = frac_v, col = frac_g, alpha = 1))  
g <- g + geom_violin(aes( fill=frac_g), trim=FALSE, scale = "width", lwd = 0.2) + geom_boxplot(width=0.5, lwd = 0.3) + theme_classic()
g <- g + scale_color_manual( values = col_s)  + scale_fill_manual( values = col_s) + theme(legend.position = "none")
ggsave("../Results/6_types_immune_cells_fraction_TIMER_violin_plot_PLCME.pdf",  width = 4, height = 3)

p_val

}

}


#########################################################
## immune infiltration estimation comparison in TCGA LUAD
#########################################################
{
###################################################
## read-in immune celss estimated fraction by TIMER
{
luad <- read.csv("immuneEstimation_TCGA.csv", row.names = 1)   ## all TCGA samples

## load LUAD normal  and tumor smaple list
load("./TCGA_LUAD_normal_tumor_sample_ID.Rdata")
idx_n <- match(tcga_normal_sample, rownames(luad), nomatch = 0)
sum(idx_n != 0)
frac_n <- luad[idx_n, ]

idx_t <- match(tcga_tumor_sample, rownames(luad), nomatch = 0)
sum(idx_t != 0)
frac_t <- luad[idx_t, ]
frac <- as.matrix(rbind(frac_n, frac_t))
}

#################################
## draw  boxplot and violin plots
{
p_val <- vector()
L <-  ncol(frac)
for (i in 1:L)
{
        tmp <- wilcox.test(frac[1:nrow(frac_n), i], frac[59:nrow(frac), i])
        p_val[i] <- tmp$p.value
    }
idx_s <- p_val < 1             ### pvalue cutoff 0.01    
    
frac_v <- as.vector(frac)
g_1 <- rep(c(rep("normal",  nrow(frac_n)), rep("tumor", nrow(frac_t))), 6)
g_2 <- rep(colnames(frac), each = nrow(frac))

frac_go <- paste(g_2, g_1, sep = "_")  
frac_gou <- unique(frac_go)
# remian sort order 
frac_gf <- as.factor(frac_go)
frac_g <- factor(frac_gf, level = frac_gou)

## significant ones with color 
idx_ss = which(idx_s == TRUE)
col_s <- rep("gray", 12)
col_s[idx_ss*2] <- "coral2"        ## Tumor
col_s[idx_ss*2-1] <- "cadetblue3"  ## normal 

## for violin plots
g <- ggplot(dat_g, aes(x = frac_g, y = frac_v, col = frac_g, alpha = 1))  
g <- g + geom_violin(aes( fill=frac_g), trim=FALSE, scale = "width", lwd = 0.2) + geom_boxplot(width=0.5, lwd = 0.3) + theme_classic()
g <- g + scale_color_manual( values = col_s)  + scale_fill_manual( values = col_s) + theme(legend.position = "none")
ggsave("../Results/6_types_immune_cells_fraction_TIMER_violin_plot_TCGA_LUAD.pdf",  width = 4, height = 3)

p_val
}

}
