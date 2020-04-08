library(biomaRt)

listMarts()

ensembl <- useMart("ENSEMBL_MART_ENSEMBL")
dataset <- listDatasets(ensembl)
head(dataset)
dataset[grep("Mouse", dataset$description), ]
ensembl <- useDataset("mmusculus_gene_ensembl", mart = ensembl)
attributes <- listAttributes(ensembl)

attributes[grep("GO", attributes$description), ]

dataGO <- getBM(attributes = c("external_gene_name", "name_1006"), mart = ensembl)

dataGO <- dataGO[order(dataGO$name_1006), ]
dataGO <- dataGO[which(dataGO$name_1006 != ""), ]
dataGO$external_gene_name <- toupper(dataGO$external_gene_name)
dataGO$name_1006 <- toupper(dataGO$name_1006)

genesets <- list()
terms <- unique(dataGO$name_1006)
for (i in 1:length(terms)) {
  genesets[[i]] <- dataGO$external_gene_name[which(dataGO$name_1006 == terms[i])]
}
names(genesets) <- terms
lgs <- sapply(1:length(genesets), function(i) {
  length(genesets[[i]])
})
genesets <- genesets[intersect(which(lgs >= 15), which(lgs <= 500))]
save(genesets, file = "mouseGO.RData")

# ------ don't use

BiocManager::install("pathview")

datakegg <- CHRONOS::downloadPathways("mmu", "All")

library(CHRONOS)
datakegg <- downloadPathways("mmu", "All")
KEGGS <- downloadPathways("mmu", datakegg)
pathview::download.kegg(datakegg, "mmu", "xml")


# ------ KEGG

BiocManager::install("KEGGREST")
library(KEGGREST)

k <- names(keggList("pathway", "mmu"))

kv <- sapply(unname(keggList("pathway", "mmu")), function(i) {
  strsplit(i, " - ")[[1]][1]
}, USE.NAMES = FALSE)

k <- lapply(k, function(i) {
  g <- keggGet(i)[[1]]$GENE
  g <- g[which(1:length(g) %% 2 == 0)]
  g <- sapply(g, function(j) {
    strsplit(j, ";")[[1]][1]
  }, USE.NAMES = FALSE)
  return(unlist(g))
})

names(k) <- kv
genesets <- k

save(genesets, file= 'mouseKEGG.RData')
# ------ GO
