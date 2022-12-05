library(edgeR)
library(Rsamtools)
library(GenomicAlignments)
library(regioneR)
library(limma)
library(ggplot2)
library(ggrepel)
###-----Cell line-----###
##1. RNA seq data processing
RQ10902A07.gene <- read.table("/mnt/nas/wtlien/All_CCA/1.Cell_line/2.RNA_analysis/RQ10902A07_count", row.names = 1)
RQ10902A08.gene <- read.table("/mnt/nas/wtlien/All_CCA/1.Cell_line/2.RNA_analysis/RQ10902A08_count", row.names = 1)
# RQ10902A09.gene <- read.table("/mnt/nas/wtlien/All_CCA/1.Cell_line/2.RNA_analysis/RQ10902A09_count", row.names = 1)
# RQ10902A10.gene <- read.table("/mnt/nas/wtlien/All_CCA/1.Cell_line/2.RNA_analysis/RQ10902A10_count", row.names = 1)
# resistance.total.gene <- data.frame(RQ10902A07.gene,RQ10902A08.gene,RQ10902A09.gene,RQ10902A10.gene)
# colnames(resistance.total.gene) <- c("SSP.25", "SSP.25_GR","SNU.1196","SNU.1196_G")
resistance.total.gene <- data.frame(RQ10902A07.gene,RQ10902A08.gene)
colnames(resistance.total.gene) <- c("SSP.25", "SSP.25_GR")
dge_resistance.gene <- DGEList(resistance.total.gene)
dge_resistance.gene.list <- filterByExpr(dge_resistance.gene)
dge_resistance.gene <- dge_resistance.gene[dge_resistance.gene.list,keep.lib.sizes=FALSE]
# dge_resistance.gene <- dge_resistance.gene[dge_resistance.gene,keep.lib.sizes=FALSE]
dge_resistance.gene <- calcNormFactors(dge_resistance.gene)
resistance.gene.table <- data.frame(resistance.total.gene)
                                    
####################################################################
resistance.gene.table$Gene_ID <- rownames(resistance.gene.table)
resistance.gene.table$Gene_ID<- gsub("\\..*","",resistance.gene.table$Gene_ID)
Clinical.only1.5FC.RNA.ATAC$Gene_ID<- gsub("\\..*","",Clinical.only1.5FC.RNA.ATAC$Gene_ID)
resistance.gene.table <- merge(resistance.gene.table, Clinical.only1.5FC.RNA.ATAC , by="Gene_ID")
resistance.gene.table <- merge(resistance.gene.table, gene.feature , by="Gene_ID")

resistance.gene.table$Description <- c("NA")
resistance.gene.table <- data.frame(resistance.gene.table$Symbol,resistance.gene.table$Description,resistance.gene.table$SSP.25,resistance.gene.table$SSP.25_GR)
colnames(resistance.gene.table) <- c("Name","Description","SSP25","SSP25_GR")

write.table(resistance.gene.table,"/mnt/nas/yh/9.CCA/2.ATAC_seq/1.final/GSEA.txt",
 sep = "\t",
 col.names = TRUE,
 row.names = FALSE,
 quote = FALSE)    
