# CellEnrich 

Pathway Enrichment Analysis and Visualization for Single-cell Data

<img src="images/ce_options_run_2024.PNG" width=400 height=270> <img src="images/ce_scatter_bar_2024.PNG" width=400 height=270> 
<img src="images/ce_dyntab_all_2024.PNG" width=400 height=250> <img src="images/pbmc_biplot_or.png" width=400 height=250>  
<img src="images/ce_ga_map_2024.PNG" width=400 height=270> <img src="images/ce_heatmap_pdf.PNG" width=400 height=270>  

## Installation

CellEnrich [manual](https://github.com/noobCoding/CellEnrich/blob/master/CellEnrich_manual.pdf) is available.

NOTE: on a fresh installation, users may need to install some required interpreter compilers for the system to install other R packages:
* C++ compiler 
* gfortran compiler (FYI: tips for [MAC](https://github.com/fxcoudert/gfortran-for-macOS/releases) users or [other OS](https://fortran-lang.org/learn/os_setup/install_gfortran/))
* Seurat >= 5.0.0 is REQUIRED
<br />  

Install Dependent Packages:
*	Some packages must be installed from sources that need compilation and a proper version of [RTools](https://cran.r-project.org/bin/windows/Rtools/).
*	It is recommended to install devtools, and BiocManager packages first before installing the following packages in Github/Bioconductor (not in CRAN).

```R
# install required packages
install.packages('Seurat') # RStudio may need a RESTART for Seurat v5.0.1 to be activated
install.packages('remotes')
install.packages('waiter')
install.packages('farver')
remotes::install_github('vqv/ggbiplot')

# install CellEnrich
remotes::install_github('noobCoding/CellEnrich')
```
<br /> 


## Load essential libraries
```R
library(CellEnrich)
library(Seurat)
```
<br />

## Example with PBMC_3K data 

```R
# Download data, if not downloaded
download.file('https://github.com/noobcoding/CellEnrich/raw/master/data/pbmcData.RData','pbmcData.RData', mode = 'wb')
download.file('https://github.com/noobcoding/CellEnrich/raw/master/data/pbmcClustInfo.RData','pbmcClustInfo.RData', mode = 'wb')
download.file('https://github.com/noobcoding/CellEnrich/raw/master/data/Human_Reactome.RData', 'Human_Reactome.RData', mode = 'wb')

# Load data
load("pbmcData.RData")
load("pbmcClustInfo.RData")

CountData <- pbmcData
GroupInfo <- pbmcClustInfo

# CellEnrich uses normalized count data as input
CountData <- NormalizeData(CountData)

# This will run CellEnrich
CellEnrich(CountData, GroupInfo)
```
<br />  

## Large datasets for testing
Two datasets for testing CellEnrich are too big for hosting on GitHub so you can directly download them at [Zenodo link](https://zenodo.org/records/13891393) including:

* 'GBM_sub' is the Glioblastoma data from [link](https://www.nature.com/articles/s41586-023-06036-1), which contains human HFC(Highly functionally connected) and LFC(Low functionally connected) Glioblastoma cells with added information of GRIA2 expression.
 
```R
# Glioblastoma data can be directly downloaded using the Zenodo link above!
download.file('https://github.com/noobcoding/CellEnrich/raw/master/data/Human_WikiPathways.RData', 'Human_WikiPathways.RData', mode = 'wb')

gbm_sub <- readRDS(file = "GBM_sub.rds")
gbm_sub$counts<-NormalizeData(gbm_sub$counts)

CellEnrich(gbm_sub$counts, gbm_sub$class)
```
<br />  

* 'PD_dat' is the Parkinson's disease data from [link](https://www.nature.com/articles/s41593-022-01061-1), which contains human dopamine cell cluster information about Parkinson's disease vs. Control.
```R
# Parkinson's disease data can be directly downloaded using the Zenodo link above!
download.file('https://github.com/noobcoding/CellEnrich/raw/master/data/Human_Reactome.RData', 'Human_Reactome.RData', mode = 'wb')

PD_dat <- readRDS("PD_dat.rds")
PD_dat$count <- NormalizeData(PD_dat$count)

CellEnrich(PD_dat$count,PD_dat$type)
```
<br />  

## Example with Alzheimer's data 

```R
# Download Alzheimer's data to the working directory
download.file('https://github.com/noobcoding/CellEnrich/raw/master/data/Alzheimer_Counts_sampled.RDS','Alzheimer_Counts_sampled.RDS', mode = 'wb')
download.file('https://github.com/noobcoding/CellEnrich/raw/master/data/Alzheimer_CellType_sampled.RDS','Alzheimer_CellType_sampled.RDS', mode = 'wb')
download.file('https://github.com/noobcoding/CellEnrich/raw/master/data/Human_Reactome.RData', 'Human_Reactome.RData', mode = 'wb')

GroupInfo <- readRDS("Alzheimer_CellType_sampled.RDS")
CountData <- readRDS("Alzheimer_Counts_sampled.RDS")

CountData <- NormalizeData(CountData)

CellEnrich(CountData, GroupInfo)
```
<br /> 

## Dependency

* [R](https://cran.r-project.org/src/base/R-4/) - >= 4.2.0
* [dplyr](https://github.com/tidyverse/dplyr) - 0.8.5
* [DT](https://github.com/rstudio/DT) - 0.13
* [farver](https://cran.r-project.org/web/packages/farver/) - 2.0.3
* [ggplot2](https://github.com/tidyverse/ggplot2) - 3.3.0
* [ggrepel](https://github.com/slowkow/ggrepel) - 0.8.2
* [highcharter](https://github.com/jbkunst/highcharter) - 0.7.0.9001
* [htmltools](https://github.com/rstudio/htmltools) - 0.4.0
* [magrittr](https://github.com/tidyverse/magrittr) - 1.5
* [Rtsne](https://github.com/jkrijthe/Rtsne) - 0.15
* [scales](https://github.com/r-lib/scales) - 1.1.0
* [scMerge](https://github.com/SydneyBioX/scMerge) - 1.5.0
* [scran](https://git.bioconductor.org/packages/scran) - 1.14.6
* [Seurat](https://github.com/satijalab/seurat) - 5.0.1
* [shiny](https://github.com/rstudio/shiny) - 1.4.0.2
* [shinyFeedback](https://github.com/merlinoa/shinyFeedback) - 0.2.0
* [shinyjs](https://github.com/daattali/shinyjs) - 1.1
* [shinymaterial](https://github.com/ericrayanderson/shinymaterial) - 1.0.1
* [SingleCellExperiment](https://git.bioconductor.org/packages/SingleCellExperiment) - 1.8.0
* [sortable](https://github.com/rstudio/sortable) - 0.4.2
* [uwot](https://github.com/jlmelville/uwot) - 0.1.8
* [waiter](https://github.com/JohnCoene/waiter) - 0.1.1.9000* 

## Authors
* Hai Nguyen *hainct@unist.ac.kr* -- [@noobCoding](http://github.com/noobCoding)
* Jinhwan Kim 
* Prof. Dougu Nam *dougnam@unist.ac.kr* 

## License
This project is [MIT](https://opensource.org/licenses/MIT) licensed

