library(edgeR)
library(Rsamtools)
library(GenomicAlignments)
library(regioneR)
library(limma)
setwd("/mnt/nas/yh/1.RNA_seq/20220616_FADU/")
FADU_time.list <- read.table("FADU_condition.txt", header = TRUE, sep="\t")
FADU_time.list <- data.frame(FADU_time.list$sample_ID,FADU_time.list$sample_name)
colnames(FADU_time.list) <- c("sample_ID","condition")
hg19_gene.feature <- read.table(
  "/data2/reference/annotation/hg19/v25/gene_feature.table", header = FALSE) 
  hg19_gene.feature <- data.frame(hg19_gene.feature$V1,
                                hg19_gene.feature$V2,
                                hg19_gene.feature$V4)
  colnames(hg19_gene.feature) <- c("Gene_ID","type","symbol")
  hg19_gene.feature$Gene_ID<- gsub("\\..*","",hg19_gene.feature$Gene_ID)
  
####################################################################################
RQ10907A57.count <- read.table("RQ10907A57.count", header = FALSE, sep="\t")
table <- data.frame(RQ10907A57.count$V1)
colnames(table ) <- c("Gene_ID")
sample_ID <- FADU_time.list$sample_ID
condition.list <- FADU_time.list$condition

for(k in 1:nrow(FADU_time.list)){
  sample <- sample_ID[k]
  condition <- condition.list[k]
  count_table <- read.table(paste0(sample,".count"),header = F,sep = "\t")
  colnames(count_table) <- c("Gene_ID",paste0(condition))
  table <- merge(count_table,table,by="Gene_ID")
}
rownames(table) <- table$Gene_ID
table <- table[,-c(1)]
#########################################################################################
##normalize counts
group_FADU_HandN_RNA <-as.factor(FADU_time.list$condition)
design.edge_FADU_HandN_RNA <- model.matrix(~0+group_FADU_HandN_RNA)
colnames(design.edge_FADU_HandN_RNA) <- levels(group_FADU_HandN_RNA)

dge_HandN_RNA <- DGEList(table, group = FADU_time.list$condition)
dge_HandN_RNA_filter <- filterByExpr(dge_HandN_RNA)
dge_HandN_RNA <- dge_HandN_RNA[dge_HandN_RNA_filter,keep.lib.sizes=FALSE]
dge_HandN_RNA = calcNormFactors(dge_HandN_RNA)
dge_HandN_RNA = cpm(dge_HandN_RNA, normalized= T)
FADU_HandN_RNA.table <- data.frame(dge_HandN_RNA)
FADU_HandN_RNA.table$Gene_ID <- rownames(FADU_HandN_RNA.table)

FADU_HandN_Log2FC <- data.frame(  
  log2(FADU_HandN_RNA.table$FADU_4hr_H/FADU_HandN_RNA.table$FADU_4hr_N),  
  log2(FADU_HandN_RNA.table$FADU_8hr_H/FADU_HandN_RNA.table$FADU_8hr_N),  
  log2(FADU_HandN_RNA.table$FADU_12hr_H/FADU_HandN_RNA.table$FADU_12hr_N),  
  log2(FADU_HandN_RNA.table$FADU_16hr_H/FADU_HandN_RNA.table$FADU_16hr_N),  
  log2(FADU_HandN_RNA.table$FADU_20hr_H/FADU_HandN_RNA.table$FADU_20hr_N),  
  log2(FADU_HandN_RNA.table$FADU_24hr_H/FADU_HandN_RNA.table$FADU_24hr_N),
  log2(FADU_HandN_RNA.table$FADU_28hr_H/FADU_HandN_RNA.table$FADU_28hr_N),
  log2(FADU_HandN_RNA.table$FADU_32hr_H/FADU_HandN_RNA.table$FADU_32hr_N),
  log2(FADU_HandN_RNA.table$FADU_36hr_H/FADU_HandN_RNA.table$FADU_36hr_N),
  log2(FADU_HandN_RNA.table$FADU_40hr_H/FADU_HandN_RNA.table$FADU_40hr_N),
  log2(FADU_HandN_RNA.table$FADU_44hr_H/FADU_HandN_RNA.table$FADU_44hr_N),
  log2(FADU_HandN_RNA.table$FADU_48hr_H/FADU_HandN_RNA.table$FADU_48hr_N),
  row.names(FADU_HandN_RNA.table))

colnames(FADU_HandN_Log2FC) <- c("FADU_4hr",  
                                "FADU_8hr",  
                                "FADU_12hr",  
                                "FADU_16hr",  
                                "FADU_20hr",  
                                "FADU_24hr",
                                "FADU_28hr",
                                "FADU_32hr",
                                "FADU_36hr",
                                "FADU_40hr",
                                "FADU_44hr",
                                "FADU_48hr",
                                "Gene_ID")

FADU_HandN_Log2FC <- merge(FADU_HandN_Log2FC,hg19_gene.feature,by="Gene_ID")

###############################################################################
write.table(FADU_HandN_Log2FC,"/mnt/nas/yh/1.RNA_seq/20220616_FADU/FADU_HandN_Log2FC.xls",
             sep = "\t",
             col.names = TRUE,
             row.names = FALSE,
             quote = FALSE)
###############################################################################

###############################################################################
RQ10907A57.count <- read.table("RQ10907A57.count", header = FALSE, sep="\t")
table <- data.frame(RQ10907A57.count$V1)
colnames(table ) <- c("Gene_ID")
sample_ID <- FADU_time.list$sample_ID
condition.list <- FADU_time.list$condition
for(k in 1:nrow(FADU_time.list)){
  sample <- sample_ID[k]
  condition <- condition.list[k]
  count_table <- read.table(paste0(sample,".count"),header = F,sep = "\t")
  colnames(count_table) <- c("Gene_ID",paste0(condition))
  table <- merge(count_table,table,by="Gene_ID")
}
rownames(table) <- table$Gene_ID
table <- table[,-c(1)]

normoxia.table <-data.frame(table[,1:12])
##normalize counts
condition.table <-data.frame(sample_ID = c(colnames(table)),
                  condition = c("a","a","a","a","a","a","a","a","a","a","a","a",
                                "A","B","C","D","E","F","G","H","I","J","K","L"))

group_FADU_normoxia_RNA <-as.factor(condition.table$condition)
design.edge_FADU_HandN_RNA <- model.matrix(~0+group_FADU_normoxia_RNA)
colnames(design.edge_FADU_HandN_RNA) <- levels(group_FADU_normoxia_RNA)

dge_HandN_RNA <- DGEList(table, group = group_FADU_normoxia_RNA)
dge_HandN_RNA_filter <- filterByExpr(dge_HandN_RNA)
dge_HandN_RNA <- dge_HandN_RNA[dge_HandN_RNA_filter,keep.lib.sizes=FALSE]
dge_HandN_RNA = calcNormFactors(dge_HandN_RNA)
#dge_HandN_RNA = cpm(dge_HandN_RNA, normalized= T)
dge_HandN_RNA <- estimateDisp(dge_HandN_RNA, design.edge_FADU_HandN_RNA , robust = TRUE)


x <- c("A-a")
contrast <- makeContrasts(contrasts=x, levels = group_FADU_normoxia_RNA)
dge_HandN_RNA <- glmQLFit(dge_HandN_RNA, robust = TRUE)
res <- glmQLFTest(dge_HandN_RNA, contrast = contrast)

FADU_HandN_RNA.table <- data.frame(res)
colnames(FADU_HandN_RNA.table)<- c("FADU_48hr_HvsN_FC","FADU_44hr_HvsN_FC","FADU_36hr_HvsN_FC",
                                  "FADU_20hr_HvsN_FC","FADU_16hr_HvsN_FC","FADU_40hr_HvsN_FC",
                                  "FADU_32hr_HvsN_FC","FADU_28hr_HvsN_FC","FADU_24hr_HvsN_FC",
                                  "FADU_12hr_HvsN_FC","FADU_8hr_HvsN_FC","FADU_4hr_HvsN_FC","logPM","F","p-value")


FADU_HandN_RNA.table$Gene_ID <- rownames(FADU_HandN_RNA.table)
FADU_HandN_RNA.table <- merge(FADU_HandN_RNA.table,hg19_gene.feature,by="Gene_ID")
