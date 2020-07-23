# CellEnrich 

Pathway enrichment analysis / visualize for Single Cell Data

Online manual is available in [here](https://jhk0530.github.io/CellEnrich/)

## :wrench: Install

```R
if(!require(remotes)){
  install.packages('remotes') # install devtools if not installed.
}
remotes::install_github('vqv/ggbiplot')
remotes::install_github('unistbig/cellenrich')
library(CellEnrich)
```

## :ship: Example run with PBMC data 

```R
# download minimal data to working directory

download.file('https://github.com/jhk0530/CellEnrich/blob/master/pbmcData.RData?raw=true','pbmcData.RData', mode = 'wb')
download.file('https://github.com/jhk0530/CellEnrich/blob/master/pbmcClustinfo.RData?raw=true','pbmcClustinfo.RData', mode = 'wb')
download.file('https://github.com/jhk0530/CellEnrich/blob/master/humanKEGG.RData?raw=true', 'humanKEGG.RData', mode = 'wb')

# load library and data

library(CellEnrich)
load("pbmcClustInfo.RData")
load("pbmcData.RData")


CountData <- pbmcData
GroupInfo <- pbmcClustInfo
rm(pbmcData, pbmcClustInfo) # remove .

# Run cellenrich
CellEnrich(CountData, GroupInfo)

```

## :paperclip: Dependency

* [dplyr](https://github.com/tidyverse/dplyr) - 0.8.5
* [DT](https://github.com/rstudio/DT) - 0.13
* [ggplot2](https://github.com/tidyverse/ggplot2) - 3.3.0
* [ggrepel](https://github.com/slowkow/ggrepel) - 0.8.2
* [highcharter](https://github.com/jbkunst/highcharter) - 0.7.0.9001
* [htmltools](https://github.com/rstudio/htmltools) - 0.4.0
* [magrittr](https://github.com/tidyverse/magrittr) - 1.5
* [Rtsne](https://github.com/jkrijthe/Rtsne) - 0.15
* [scales](https://github.com/r-lib/scales) - 1.1.0
* [scran](https://git.bioconductor.org/packages/scran) - 1.14.6
* [Seurat](https://github.com/satijalab/seurat) - 3.2.0
* [shiny](https://github.com/rstudio/shiny) - 1.4.0.2
* [shinyFeedback](https://github.com/merlinoa/shinyFeedback) - 0.2.0
* [shinyjs](https://github.com/daattali/shinyjs) - 1.1
* [shinymaterial](https://github.com/ericrayanderson/shinymaterial) - 1.0.1
* [SingleCellExperiment](https://git.bioconductor.org/packages/SingleCellExperiment) - 1.8.0
* [sortable](https://github.com/rstudio/sortable) - 0.4.2
* [uwot](https://github.com/jlmelville/uwot) - 0.1.8
* [waiter](https://github.com/JohnCoene/waiter) - 0.1.1.9000

## :blush: Authors
* :octocat: Jinhwan Kim [@jhk0530](http://github.com/jhk0530)
* prof. Dougu Nam *dougnam@unist.ac.kr*

## :memo: License
This project is [MIT](https://opensource.org/licenses/MIT) licensed

