rm(list = ls())
library(Seurat)
library(stringr)
library(tidyverse)
load("SAVE/pbmc_NSFK_AN.Rdata")
ident_df <- data.frame(cell=names(Idents(pbmc_NSFK)), cluster=Idents(pbmc_NSFK))
pbmc_NSFK2 <- subset(pbmc_NSFK, cells = as.vector(ident_df[!ident_df$cluster=="Mononeuclear phagocytes",1]))
Mphs <- readRDS("SAVE/pbmc_NSFK_Mononeuclear_phagocytes.Rds")

test_counts_1 <- as.data.frame(pbmc_NSFK2@assays$RNA@data)
test_counts_2 <- as.data.frame(Mphs@assays$RNA@data)
test_counts <- cbind(test_counts_1,test_counts_2)
#整合meta数据
table(pbmc_NSFK2@active.ident)
test_meta1 <- data.frame(Cell = rownames(pbmc_NSFK2@meta.data), 
                          cell_type = pbmc_NSFK2@active.ident ) 

table(Mphs@active.ident)
test_meta2 <- data.frame(Cell = rownames(Mphs@meta.data),cell_type = Mphs@active.ident) 
test_meta2$cell_type <- ifelse(str_detect(test_meta2$cell_type, "DC"),"DC",as.character(test_meta2$cell_type))
test_meta <- rbind(test_meta1,test_meta2)

#整理导出数据
head(test_meta)
test_meta$cell_type=gsub(' ','_',test_meta$cell_type)
test_meta$cell_type=gsub('\\+','',test_meta$cell_type) 
table(test_meta$cell_type)
length(unique(test_meta$Cell))

identical(colnames(test_counts),test_meta$Cell) 
test_counts=cbind(rownames(test_counts),test_counts)
colnames(test_counts)[1]='Gene'
#转换成ENSEMBL 并去重
library(clusterProfiler)
library(org.Hs.eg.db)
ids <- bitr(test_counts$Gene,'SYMBOL','ENSEMBL','org.Hs.eg.db')
test_counts <- merge(test_counts,ids,by.x='Gene',by.y='SYMBOL')
test_counts <- test_counts[!duplicated(test_counts$ENSEMBL),]
rownames(test_counts)<- test_counts$ENSEMBL
test_counts$Gene <- rownames(test_counts)
library(dplyr)
test_counts <- dplyr::select(test_counts,-ENSEMBL)

write.table(test_counts, "int3/test_counts.txt",  row.names=F, sep='\t',quote = F)
write.table(test_meta, "int3/test_meta.txt", row.names=F, sep='\t',quote = F)

