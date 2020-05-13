library(ggplot2)
library(dplyr)

# hyper

getHyperPvalue <- function(genes, genesets) {
  A <- length(unique(unlist(genesets)))
  pv <- rep(0, length(genesets))

  for (i in 1:length(genesets)) {
    q <- length(intersect(genesets[[i]], genes))
    m <- length(genesets[[i]])
    n <- A - m
    k <- length(genes)
    pv[i] <- 1 - phyper(q - 1, m, n, k)
  }
  names(pv) <- names(genesets)

  res <- data.frame(name = names(pv), qv = p.adjust(pv, 'fdr'))
  if(nrow(res %>% filter(qv <= 0.05)) > 10){
    res <- res %>% filter(qv <= 0.05) %>% arrange(qv)
    write.table(res, file='tmp.csv', sep=',', append = FALSE, row.names = FALSE, col.names = FALSE)
    return()

  }
  else{
    res <- res %>% arrange(qv)
    write.table(res[1:10,], file='tmp.csv', sep = ',', append = FALSE, row.names = FALSE, col.names = FALSE)
  }

}

########################### pbmc

# B
genes <- c('TYROBP','S100A9','S100A8','CST3','LYZ',
           'FCER1G','FTL','TYMP','FTH1','HLA-DRA',
           'GSTP1','SAT1','CTSS','GPX1','CFD',
           'CD14','SERPINA1','CDA','VIM','APOBEC3A',
           'FOLR3','NCF1','RPL29','RPL10','CD99',
           'RPS24','RPS4Y1')
# CD14+MONO
genes <- c('CCL5','GZMA','IL32','CD3D','CTSW',
           'CD3E','RPS27A','HCST','LCK','PRF1',
           'B2M','HLA-C','CD8A','RPL21','HLA-A',
           'CCL4','CD7','GZMB','RPS3A','CD99',
           'RPS27','RPS12','NCR3','MATK','RPL13',
           'RPS25','RPS3','RPL3','RPL32','RPS2',
           'RPL14','RPL7','RPS24')

# CD8 T
genes <- c('FCGR3A','CD68','SPI1','CFD','FCER1G',
           'IFI30','SERPINA1','TYROBP','TYMP','CST3',
           'HLA-DPA1','CTSS','NPC2','SAT1','TIMP1',
           'HES4','CEBPB','FTH1','FTL','CDKN1C',
           'TMSB4X','RPS19','CD79B','ADA','B2M',
           'VIM','RPL21','RPL10','RPL14','RPS24',
           'RPL7','RPS4Y1')

# DC

genes <- c('GZMB','PRF1','GZMA','FCGR3A','CTSW',
           'CCL4','TYROBP','FCER1G','CD247','CD7',
           'XCL2','KLRD1','AKR1C3','CCL5','CD99',
           'RPS27A','RPS3A','RPS24')

# FCGR3A+ Mono

genes <- c('CD3D','RPS27A','RPS27','LTB','RPL9',
           'RPL21','RPS3A','IL32','LDHB','RPS25',
           'CD3E','RPL31','RPL3','RPL30','RPS12',
           'RPS3','RPS23','CD7','IL7R','RPLP2',
           'RPL32','RPL13','RPS18','LCK','RPL18A',
           'RPL11','CCR7','RPS14','RPS28','TMSB4X',
           'CD27','RPL10','CD3G','LEF1','SELL',
           'CD8B','FHIT','ACTN1','PFN2','NOG')

# Memory CD4 T

genes <- c('CD79A','HLA-DRB1','HLA-DQA1','HLA-DRA','HLA-DPB1',
           'HLA-DPA1','HLA-DQB1','CD79B','MS4A1','CD74',
           'LTB','TCL1A','CD37','CXCR4','FCER2',
           'RPS27A','RPL21','RPL9','RPS3A','RPL30',
           'BIRC3','RASGRP2')

# Naive CD4 T

genes <- c('IL32','LTB','CD3D','RPS27A','RPS27',
           'CD3E','LDHB','RPSA','RPL9','RPL21',
           'RPS25','RPS3A','IL7R','TMSB4X','RPS3',
           'CD2','RPL3','RPL31','RPL30','RPS12',
           'B2M','LCK','VIM','RPLP2','HLA-A',
           'RPL32','CD7','RPL13','RPL11','RPL10',
           'RPLP1','HLA-C','ACTB','CD3G','FTH1',
           'RPS4Y1')


# NK

genes <- c('HLA-DQA1','IFI30','HLA-DQA2','HLA-DQB1','HLA-DRB5',
           'FCER1A','CD33','SPI1','HLA-DMA','SERPINF1',
           'HLA-DMB','CTSH','PLD4','CD1C','ALDH2',
           'KLF4','FCER1G','CFP','TYROBP','CST3',
           'SAMHD1','HLA-DPA1','HLA-DPB1','PPP1R14A','HLA-DRA',
           'BLNK','CD74','CD72','VIM','RPL5',
           'RPL10','RPL7A','RPL22','RPS24')

# Platelet

genes <- c('GP9','GNG11','PF4','ITGA2B','PTCRA',
           'PPBP','CA2','CLU','SNCA','MYL9',
           'CD9','GP1BA','F13A1','TUBB1','RUFY1',
           'CCL5')


# compareplot

# cellenrich
# fisher
# intersect

# pbmc 0.3 mean

int <- data.frame(
  val = c(
    25, 5, 3,
    28, 15, 6,
    12, 4, 1,
    37, 10, 5,
    35, 7, 1,
    7, 25, 1,
    1, 10, 1,
    20, 25, 7,
    2, 5, 0
  ),
  method = rep(c('CellEnrich', 'Fisher','Overlap')),
  group = rep(c('B','CD14+ Mono','CD8 T','DC','FCGR3A+ Mono','Memory CD4 T','Naive CD4 T','NK','Platelet'), each = 3)
)

ggplot(int, aes(x = group, y= val, fill = method)) +
  geom_bar(stat = 'identity', position = 'dodge') +
  xlab('Group') +
  ylab('Number of Pathway') +
  ggtitle(label = 'pbmc OddRatio 0.3 mean')

# pbmc 0.1 mean

int <- data.frame(
  val = c(
    26, 5, 3,
    42, 15, 8,
    26, 4, 1,
    48, 7, 4,
    49, 7, 1,
    8, 25, 1,
    1, 10, 1,
    31, 25, 7,
    15, 5, 1
  ),
  method = rep(c('CellEnrich', 'Fisher','Overlap')),
  group = rep(c('B','CD14+ Mono','CD8 T','DC','FCGR3A+ Mono','Memory CD4 T','Naive CD4 T','NK','Platelet'), each = 3)
)

ggplot(int, aes(x = group, y= val, fill = method)) +
  geom_bar(stat = 'identity', position = 'dodge') +
  xlab('Group') +
  ylab('Number of Pathway') +
  ggtitle(label = 'pbmc OddRatio 0.1 mean')

# pbmc 0.3 median

int <- data.frame(
  val = c(
    25, 5, 3,
    28, 15, 6,
    12, 4, 1,
    37, 7, 4,
    35, 7, 1,
    8, 25, 1,
    1, 10, 1,
    20, 25, 7,
    11, 5, 1
  ),
  method = rep(c('CellEnrich', 'Fisher','Overlap')),
  group = rep(c('B','CD14+ Mono','CD8 T','DC','FCGR3A+ Mono','Memory CD4 T','Naive CD4 T','NK','Platelet'), each = 3)
)

ggplot(int, aes(x = group, y= val, fill = method)) +
  geom_bar(stat = 'identity', position = 'dodge') +
  xlab('Group') +
  ylab('Number of Pathway') +
  ggtitle(label = 'pbmc OddRatio 0.3 median')

# pbmc 0.1 median

int <- data.frame(
  val = c(
    26, 5, 3,
    41, 15, 7,
    22, 4, 1,
    48, 7, 4,
    48, 7, 1,
    9, 25, 1,
    1, 10, 1,
    31, 25, 7,
    25, 25, 1
  ),
  method = rep(c('CellEnrich', 'Fisher','Overlap')),
  group = rep(c('B','CD14+ Mono','CD8 T','DC','FCGR3A+ Mono','Memory CD4 T','Naive CD4 T','NK','Platelet'), each = 3)
)

ggplot(int, aes(x = group, y= val, fill = method)) +
  geom_bar(stat = 'identity', position = 'dodge') +
  xlab('Group') +
  ylab('Number of Pathway') +
  ggtitle(label = 'pbmc OddRatio 0.1 median')


# pancreas

# acinar
genes <- c('MRC2','PDGFRB','GNG11','COL18A1','ITGA11',
           'LTBP1','FMOD','COL12A1','SNAI2','SFRP2',
           'LAMA2','COL6A2','COL5A1','THBS2','COL5A3',
           'EDNRA','ZEB2','LUM','COL4A1','TPM2',
           'VCAN','FBN1','COL5A2','COL4A2','COL15A1',
           'ANTXR1','VIM','COL3A1','COL6A3','COL6A1',
           'COL1A2','FN1','CREB3L1','UNC5B','LAMB1',
           'COL1A1','CALD1','MSN','PSAT1')

# activated_stellate
genes <- c('ELOVL2','COL22A1','NPY1R','BHMT','BMP4',
           'MYO1A','PDE2A','CALCR','CALY','VTN',
           'AKR1C4','COL2A1','GHRL','ASGR1','OPRK1',
           'FGF14','HMGCS2','HRH2','SUCNR1','NKX2-2',
           'UGT2B4','L1CAM','GLDC','F10','CPLX2',
           'DEFB1','APOH','ADRB1','PPP2R2B','KIF5C',
           'SNAP25','HNF4A','CALD1','ENPP2','CDKN2A',
           'FOXA2','SERPINA1','RASD1','ANK3','ACSL1',
           'CLDN7','PAPSS2','SOX4')
# alpha
genes <- c('GCG','PPP1R1A','TTR','CPE','GNG4',
           'MAFB','PAX6','PTPRN','NEUROD1','CNTN1',
           'ENPP2','ISL1','KCNMB2','PDE3B','PTPRN2',
           'GRIA3','PGR','F10','SSTR2','FEV',
           'ALDH1A1','MYO10','SH3GL2','PLCE1','CAMK2B',
           'CD99L2','FABP5','GLS','SSX2IP','DPP4',
           'PAPSS2','SDC2','TUSC3','LDHA','PFN2',
           'VIM','CLU','LAPTM4B','CKB','FOXA2',
           'CLDN7')

# beta

genes <- c('HHEX','GAD2','ISL1','NEUROD1','SST',
           'PARVB','PTPRN','RFX6','CAMK2B','ABCC8',
           'SNAP25','CASR','UCHL1','KCNJ8','LEPR',
           'CPE','GABRB3','SSTR1','CADM1','GHSR',
           'HADH','EDN3','UNC5B','PPP1R1A','MAFB',
           'CD9','PDX1','TJP2','PCLO','BAMBI',
           'SEC11C','BHLHE41','SCD5','RGS2','SERPINA1',
           'SLC38A1','AQP3','CALD1','PFN2','FOXA2',
           'MARCKS')

# delta

genes <- c('ABCC8','INS','UCHL1','ADCYAP1','GAD2',
           'HADH','MAFA','PDX1','IAPP','PAX6',
           'PTPRN','NEUROD1','G6PC2','CASR','MAFB',
           'ATP2A3','GNG4','SCD5','CPE','KIF5C',
           'SH3GL2','PPP1R1A','PFKFB2','RASD1','ROBO2',
           'GSN','RASGRF1','PFN2','VEGFA','SIPA1L2',
           'ALCAM','RRAGD','FOXA2','CLDN7')

# ductal

genes <- c('CFTR','KRT19','SDC4','SERPING1','ZFP36L1',
           'CLDN1','SLC4A4','DEFB1','LITAF','SOX4',
           'CD9','CFB','COL18A1','SPP1','CAV2',
           'SERPINA5','CTSH','CD24','LCN2','MMP7',
           'AQP1','ALDH1A3','BACE2','VTCN1','CD44',
           'KRT18','CXADR','DAB2','ENAH','TMSB4X',
           'CLDN4','ANXA2','LAMC1','ADCY5','JUP',
           'CLDN7','VAMP8','CDH1','ATP1B1','PERP'
           )

# endothelial
genes <- c('KDR','ERG','S1PR1','GNG11','CDH5',
           'ITGA5','MEF2C','RHOJ','CLEC2B','ESAM',
           'SNAI1','RAMP2','ACVRL1','DLL4','A2M',
           'CXCR4','PDE2A','COL4A1','PODXL','ARHGDIB',
           'SHANK3','CALCRL','COL4A2','PTPRB','RAPGEF4',
           'FLT1','VIM','EFNA1','ID3','MSN',
           'WWTR1','PLXNA2')
# epsilon
genes <- c('ITGB3','GLI2','EDNRB','PDGFRB','SNAI1',
           'GYPC','FABP4','CAV1','TRPV2','IFI16',
           'EDNRA','A2M','MRC2','HEYL','TBX3',
           'ENPEP','BST2','COLEC11','CD4','ITGA4',
           'COL6A2','ZEB2','AVPR1A','COL5A1','LAMA4',
           'WNT5A','COL4A1','ESAM','NOTCH2','COL5A2',
           'ARHGDIB','TRPC6','COL4A2','COX4I2','COL15A1',
           'ITGA1','SLC6A12','SLIT3','VIM','COL3A1',
           'COL1A2','CALD1','CA2','COL1A1','CHSY1',
           'COX7A1','UNC5B','ADCY4','SDC2')

# gamma

genes <- c('GSTA1','GSTA2','ANPEP','IL32','PRSS3',
           'CTRB2','CPA2','AKR1C3','CTRB1','KLK1',
           'PRSS1','BCAT1','PLA2G1B','GATM','CPA1',
           'CPB1','CEL','SDC4','PNLIP','MGST1',
           'LDHB','CD24','GDF15','CD44','KRT18',
           'CLDN4','CLDN7','RAB11FIP1')

# macrophage
genes <- c('PTGS1','BTK','HPGDS','KIT','CPA3',
           'NCF4','LCP2','RUNX3','BATF','CSF2RB',
           'LAPTM5','TPSB2','GRAP2','PRKCB','TPSAB1',
           'ALOX5AP','RAB27B','HCLS1','ARHGDIB','PTPRC',
           'SLC45A3','CD37','IL4R','RAC2','ABCC4',
           'VIM','SLC18A2','GMPR')

# mast
genes <- c('ITGB2','C1QC','RUNX3','CD86','CCR1',
           'HCK','C1QB','LCP2','PIK3R5','BCL2A1',
           'C1QA','CSF1R','FPR3','DOCK2','LAPTM5',
           'VSIG4','TYROBP','HCLS1','FCGR2A','PTPRC','ACP5','SDS','SIRPA','ALOX5AP','HLA-DRA','TNFRSF1B','WIPF1','VIM','CD74','LITAF')
# quiescent_stellate
genes <- c('WNT16','PTPRZ1','EDNRB','UCN2','GNG11',
           'ITGB3','IFI16','SOX2','HEPH','ZEB2',
           'NGFR','PLA2G4A','SLIT2','COL9A3','A2M',
           'MEF2C','MMP14','NLGN4X','ABCB4','RUNX2',
           'NRXN3','CRYAB','NRXN1','COL4A2','ASPA',
           'SEMA3C','TIMP3','ITGA6','TUBB6','WIPF1',
           'LPCAT2','BVES','PLAUR','AP1S2','WWTR1',
           'VIM','PPP2R2B','CALD1','ITGB8','CD9',
           'ETS1','ITPK1','DAG1','MPDZ')

# schwann
genes <- c('ISL1','CXXC4','GAD2','KCNMB2','PPY',
           'NEUROD1','ABCC8','PAX6','ETV1','CNTN1',
           'PGR','PTPRN','CAMK2B','TSPAN7','GRIA3',
           'UCHL1','TTR','GNG4','CARD11','FEV',
           'CPE','AMOTL1','ID4','NPNT','ENTPD2',
           'PCLO','S100A10','SEMA3E','SLC6A4','CHRM3',
           'INPP5F','PRKACB','ID2','PTPRN2','ABCB1',
           'TUBB2A','SLC40A1','GCNT3','AQP3','ALDH1A1',
           'TUSC3','FGFR1','BHLHE41','CLU','RAP1GAP',
           'FOXA2','MARCKS')

getHyperPvalue(genes, genesets)

# compareplot

# cellenrich
# fisher
# intersect

# pancreas 0.3 mean
int <- data.frame(
  val = c(
    46, 10, 0,
    76, 1, 0,
    25, 4, 0,
    17, 3, 1,
    21, 5, 0,
    70, 6, 1,
    59, 3, 1,
    30, 11, 0,
    52, 10, 0,
    69, 1, 1,
    26, 10, 3,
    78, 8, 2,
    52, 7, 0
  ),
  method = rep(c('CellEnrich','Fisher','Overlap')),
  group = rep(c('acinar', 'activated_stellate', 'alpha', 'beta', 'delta', 'ductal', 'endothelial', 'epsilon', 'gamma', 'macrophage', 'mast', 'quiescent_stellate','schwann'),each = 3)
)

ggplot(int, aes(x = group, y= val, fill = method)) +
  geom_bar(stat = 'identity', position = 'dodge') +
  xlab('Group') +
  ylab('Number of Pathway') +
  ggtitle(label = 'pancreas OddRatio 0.3 mean')


# pancreas 0.3 median

int <- data.frame(
  val = c(
   61, 10, 2,
   92, 1, 0,
   30, 4, 1,
   27, 3, 1,
   4, 5,0,
   66, 6, 0,
   79, 3, 2,
   35, 11, 0,
   40, 10, 0,
   89, 1, 1,
   29, 10, 2,
   90, 8, 2,
   67, 7, 0
  ),
  method = rep(c('CellEnrich','Fisher','Overlap')),
  group = rep(c('acinar', 'activated_stellate', 'alpha', 'beta', 'delta', 'ductal', 'endothelial', 'epsilon', 'gamma', 'macrophage', 'mast', 'quiescent_stellate','schwann'),each = 3)
)

ggplot(int, aes(x = group, y= val, fill = method)) +
  geom_bar(stat = 'identity', position = 'dodge') +
  xlab('Group') +
  ylab('Number of Pathway') +
  ggtitle(label = 'pancreas OddRatio 0.3 median')


hsa <- keggList('hsa')

hsa <- keggLink('pathway', 'hsa')
hsa <- unname(hsa)
hsa <- sort(unique(hsa))
hsa <- sapply(hsa, function(i){gsub('path:','',i)})
hsa <- unname(hsa)

genesets <- list()
for(i in 1:length(hsa)){
  res <- keggGet(hsa[i])[[1]]
  genes <- unname(sapply(res$GENE[grep(';', res$GENE)], function(i){strsplit(i, '; ')[[1]][1]}))
  genesets[[i]] <- genes
  names(genesets)[i] <- gsub(' - Homo sapiens (human)', '', res$NAME)
}

names(genesets) <- unname(sapply(names(genesets), function(i){strsplit(i, ' - ')[[1]][1]}))
grep('Amoe',names(genesets))



