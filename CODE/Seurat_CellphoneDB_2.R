rm(list = ls())
library(Seurat)
Fib <- readRDS("SAVE/pbmc_NSFK_Fibroblasts.Rds")
Mac <- readRDS("SAVE/NSFK_Macrophage.Rds")

####整合count数据
test_counts_1 <- as.data.frame(Fib@assays$RNA@data)
test_counts_2 <- as.data.frame(Mac@assays$RNA@data)
test_counts <- cbind(test_counts_1,test_counts_2)
#整合meta数据
table(Fib@active.ident)
test_meta1 <-  data.frame(Cell = rownames(Fib@meta.data), 
                          cell_type = Fib@active.ident ) 

table(Mac@active.ident)
test_meta2 <-  data.frame(Cell = rownames(Mac@meta.data), 
                          cell_type = Mac@active.ident ) 

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

write.table(test_counts, "int1/test_counts.txt",  row.names=F, sep='\t',quote = F)
write.table(test_meta, "int1/test_meta.txt", row.names=F, sep='\t',quote = F)

#conda 代码
#conda activate cellphonedb
#cellphonedb method statistical_analysis test_meta.txt test_counts.txt  
#--counts-data hgnc_symbol --threads 4 用基因名的时候用这个

rm(list = ls())
Endo <- readRDS("SAVE/pbmc_NSFK_Endothelial_cells.Rds")
Mac <- readRDS("SAVE/NSFK_Macrophage.Rds")
####整合count数据
test_counts_1 <- as.data.frame(Endo@assays$RNA@data)
test_counts_2 <- as.data.frame(Mac@assays$RNA@data)
test_counts <- cbind(test_counts_1,test_counts_2)
#整合meta数据
table(Endo@active.ident)
test_meta1 <-  data.frame(Cell = rownames(Endo@meta.data), 
                          cell_type = Endo@active.ident ) 

table(Mac@active.ident)
test_meta2 <-  data.frame(Cell = rownames(Mac@meta.data), 
                          cell_type = Mac@active.ident ) 

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
test_counts <- dplyr::select(test_counts,-ENSEMBL)

write.table(test_counts, "int2/test_counts.txt",  row.names=F, sep='\t',quote = F)
write.table(test_meta, "int2/test_meta.txt", row.names=F, sep='\t',quote = F)

