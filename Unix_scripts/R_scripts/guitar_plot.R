
args <- commandArgs(trailingOnly = TRUE)
sample_summit <- args[1];                     # for one sample's called peakss
sample_name <-args[2];
ref <- args[3]  			## ref txdb 

load(ref)

name_out <- paste(sample_name, "_summit_Guitar_plot", sep = "")

library(Guitar)
library(regioneR)
library(ggplot2)

print(sample_summit)
tmp_bed <- toGRanges(sample_summit)
sum_GR <- split(tmp_bed, f=tmp_bed$score)	 
names(sum_GR) <- sample_name
head(tmp_bed)
## plot for onesample 

GuitarPlot(gfeatures = sum_GR, returnCount = TRUE,  GuitarCoordsFromTxDb = guitar_txdb_hg38, rescaleComponent=T, saveToPDFprefix = name_out)

#######################################
## out put density for mRNA and lincRNA
######################################


## for lincRNA
    adjust = 1
    ct$weight <- ct$count
    
    ct1 <- ct[ct$category == "mRNA", ]
    ct2 <- ct[ct$category == "lncRNA", ]
    pos = Feature = weight = NULL
    id1 <- which(match(ct1$comp, c("Front", "Back")) > 0)
    ct1 <- ct1[-id1, ]
    id2 <- which(match(ct2$comp, c("Front", "Back")) > 0)
    ct2 <- ct2[-id2, ]
    
    ct1$weight <- ct1$weight/sum(ct1$weight)
    ct2$weight <- ct2$weight/sum(ct2$weight)

    p2 <- ggplot(ct2, aes(x = pos, group = Feature, weight = weight)) + geom_density(adjust = adjust)
    p2_out <- ggplot_build(p2)
    
    lnc <- cbind(p2_out$data[[1]]$x, p2_out$data[[1]]$y)
    colnames(lnc) <- c("pos", "density") 
   
    write.table(lnc, file = paste(name_out, "_lncRNA__density.txt", sep = ""))


## for mRNA 
if(FALSE){
    comLength = c(0.136, 0.459, 0.405)     ## relative length the 5'UTR, CDS, and 3'UTR, based on Guitar
    weight <- comLength/sum(comLength)
    names(weight) <- c("5'UTR", "CDS", "3'UTR")
    cds_id <- which(ct1$comp == "CDS")
    utr3_id <- which(ct1$comp == "UTR3")
    utr5_id <- which(ct1$comp == "UTR5")
    ct1$count[utr5_id] <- ct1$count[utr5_id] * weight["5'UTR"]
    ct1$count[cds_id] <- ct1$count[cds_id] * weight["CDS"]
    ct1$count[utr3_id] <- ct1$count[utr3_id] * weight["3'UTR"]
    ct1$weight <- ct1$count/sum(ct1$count)
    
    # scale position
    x <- cumsum(weight)
    ct1$pos[utr5_id] <- ct1$pos[utr5_id] * weight["5'UTR"] + 0
    ct1$pos[cds_id] <- ct1$pos[cds_id] * weight["CDS"] + x[1]
    ct1$pos[utr3_id] <- ct1$pos[utr3_id] * weight["3'UTR"] +  x[2]
    
    p1 <- ggplot(ct1, aes(x = pos, group = Feature, weight = weight))  + geom_density(adjust = adjust)  
    p1_out <- ggplot_build(p1)

    pc <- cbind(p1_out$data[[1]]$x, p1_out$data[[1]]$y)
    colnames(pc) <- c("pos", "density")
    write.table(pc, file = paste(name_out, "_PC__density.txt", sep = ""))   
    


}







 
