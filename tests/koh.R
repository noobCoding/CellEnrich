library(DuoClustering2018)
library(SingleCellExperiment)
library(dplyr)
library(tidyr)
dat <- sce_full_Koh()
CountData <- dat@assays$data$normcounts
colnames(CountData) <- CellEnrich::cellNames(dat)

BiocManager::install('EnsDb.Hsapiens.v86')
library(EnsDb.Hsapiens.v86)
CountData[1:5,1:5]

# 1. Convert from ensembl.gene to gene.symbol

sce <- sce_filteredExpr10_Koh()

CountData <- counts(sce)

genes <- as.character(rownames(CountData))
genes <- sapply(genes, function(i){ strsplit(i, '\\.')[[1]][1] }, USE.NAMES = F)
genes <- ensembldb::select(EnsDb.Hsapiens.v86, keys= genes, keytype = "GENEID", columns = c("SYMBOL","GENEID"))
newGenes <- genes

genes <- as.character(rownames(CountData))
genes <- sapply(genes, function(i){ strsplit(i, '\\.')[[1]][1] }, USE.NAMES = F)

rem <- sapply(setdiff(genes, newGenes[,2]), function(i){which(genes == i)})
CountData <- CountData[-rem,]
rownames(CountData) <- newGenes[,1]
colnames(CountData) <- CellEnrich::cellNames(sce)
CountData[1:5,1:5]

groupinfo <- colnames(CountData)
CellNames <- sapply(groupinfo, function(i){strsplit(i, '-')[[1]][1]}, USE.NAMES = F)
colnames(CountData) <- CellNames
GroupInfo <- sapply(groupinfo, function(i){strsplit(i, '-')[[1]][2]}, USE.NAMES = F)

# minimal value is 1 :

library(Matrix)

load("koh.RData")
CountData <- CountData - 1
load("kohinfo.RData")
library(CellEnrich)
CellEnrich(CountData, GroupInfo)

# duplicated row remove

# load("koh.RData")
dup <- sort(rownames(CountData)[which(duplicated(rownames(CountData)))])
udup <- unique(dup)

for(i in 1:length(udup)){
  if(i%%10==0) print(i)
  idx <- which(rownames(CountData)==udup[i])
  if(length(idx)==1){next}
  adds <- colMeans(CountData[idx,])
  CountData <- rbind(CountData, adds)
  rownames(CountData)[nrow(CountData)] = udup[i]
  CountData <- CountData[-idx,]
}

save(CountData, file='koh_dgc.RData')

clres <- clustering_summary_filteredExpr10_Koh_v1()
clres <- clres %>% filter(method == 'PCAHC')
clres$cluster = as.character(as.numeric(clres$cluster) + 1)


kohTRUE <- clres %>% group_by(cell) %>% arrange(run) %>% filter(row_number()==1)
kohTRUE <- kohTRUE %>% select(cell, cluster, trueclass)


kohTRUE2 <- kohTRUE %>% inner_join(ctr2)
mclust::adjustedRandIndex(kohTRUE2$trueclass, kohTRUE2$cl)
mclust::adjustedRandIndex(kohTRUE2$trueclass, kohTRUE2$cluster)
