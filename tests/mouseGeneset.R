library(biomaRt)

listMarts()

ensembl <- useMart("ENSEMBL_MART_ENSEMBL")
dataset <- listDatasets(ensembl)
head(dataset)
dataset[grep("Mouse", dataset$description), ]
ensembl <- useDataset("mmusculus_gene_ensembl", mart = ensembl)

dataset[grep("Human", dataset$description), ]
ensembl <- useDataset("hsapiens_gene_ensembl", mart = ensembl)

attributes <- listAttributes(ensembl)

attributes[grep("GO", attributes$description), ]

dataGO <- getBM(attributes = c("namespace_1003", "external_gene_name", "name_1006"), mart = ensembl)

dataGO <- dataGO[which(dataGO$name_1006 != ""), ]

library(dplyr)
dataGOMF <- dataGO %>% filter(namespace_1003 =='molecular_function')
dataGOBP <- dataGO %>% filter(namespace_1003 =='biological_process')
dataGOCC <- dataGO %>% filter(namespace_1003 =='cellular_component')

buildGeneset = function(dataGO, filename){
  dataGO <- dataGO %>% select(-namespace_1003)

  dataGO <- dataGO[order(dataGO$name_1006), ]

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
  save(genesets, file = filename)
}

buildGeneset(dataGOMF, 'mouseGOMF.RData')
buildGeneset(dataGOBP, 'mouseGOBP.RData')
buildGeneset(dataGOCC, 'mouseGOCC.RData')

buildGeneset(dataGOMF, 'humanGOMF.RData')
buildGeneset(dataGOBP, 'humanGOBP.RData')
buildGeneset(dataGOCC, 'humanGOCC.RData')

buildGeneset(dataGO, 'humanGO.RData')

# ------ don't use

# BiocManager::install("pathview")

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
