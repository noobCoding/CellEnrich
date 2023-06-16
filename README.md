# CellEnrich 

Pathway enrichment analysis / visualize for Single Cell Data

Online manual is available [here](https://github.com/noobCoding/CellEnrich/wiki)

<img src="data/figures/scatter-area.png"> 
<img src="data/figures/severe-freq.png"> 

## :wrench: Install

```R
if(!require(remotes)){
  install.packages('remotes') 
}
remotes::install_github('vqv/ggbiplot')
install.packages('farver') # install 'farver' if not installed.
remotes::install_github('noobCoding/CellEnrich')
library(CellEnrich)
```

## Example run with Alzheimer data 

```R
# download minimal data to working directory
download.file('https://github.com/noobcoding/CellEnrich/raw/master/data/Alzheimer_Counts_sampled.RDS','Alzheimer_Counts_sampled.RDS', mode = 'wb')
download.file('https://github.com/noobcoding/CellEnrich/raw/master/data/Alzheimer_CellType_sampled.RDS','Alzheimer_CellType_sampled.RDS', mode = 'wb')
download.file('https://github.com/noobcoding/CellEnrich/raw/master/data/WikiPathways_2021_Human.RData', 'WikiPathways_2021_Human.RData', mode = 'wb')

# load library and data
library(CellEnrich)
readRDS("Alzheimer_CellType_sampled.RDS")
readRDS("Alzheimer_Counts_sampled.RDS")

CountData <- Alzheimer_Counts_sampled
GroupInfo <- Alzheimer_CellType_sampled
rm(Alzheimer_Counts_sampled, Alzheimer_CellType_sampled) # remove 

# Run cellenrich
CellEnrich(CountData, GroupInfo)

```

## Dependency

* [dplyr](https://github.com/tidyverse/dplyr) - 0.8.5
* [DT](https://github.com/rstudio/DT) - 0.13
* [faver](https://cran.r-project.org/web/packages/farver/) - 2.0.3
* [ggplot2](https://github.com/tidyverse/ggplot2) - 3.3.0
* [ggrepel](https://github.com/slowkow/ggrepel) - 0.8.2
* [highcharter](https://github.com/jbkunst/highcharter) - 0.7.0.9001
* [htmltools](https://github.com/rstudio/htmltools) - 0.4.0
* [magrittr](https://github.com/tidyverse/magrittr) - 1.5
* [Rtsne](https://github.com/jkrijthe/Rtsne) - 0.15
* [scales](https://github.com/r-lib/scales) - 1.1.0
* [scMerge](https://github.com/SydneyBioX/scMerge) - 1.5.0
* [scran](https://git.bioconductor.org/packages/scran) - 1.14.6
* [Seurat](https://github.com/satijalab/seurat) - 3.2.0
* [shiny](https://github.com/rstudio/shiny) - 1.4.0.2
* [shinyFeedback](https://github.com/merlinoa/shinyFeedback) - 0.2.0
* [shinyjs](https://github.com/daattali/shinyjs) - 1.1
* [shinymaterial](https://github.com/ericrayanderson/shinymaterial) - 1.0.1
* [SingleCellExperiment](https://git.bioconductor.org/packages/SingleCellExperiment) - 1.8.0
* [Slingshot](https://github.com/kstreet13/slingshot) - 1.9.1
* [sortable](https://github.com/rstudio/sortable) - 0.4.2
* [uwot](https://github.com/jlmelville/uwot) - 0.1.8
* [waiter](https://github.com/JohnCoene/waiter) - 0.1.1.9000* 

## Authors
* Dr. Hai Nguyen *hainct@unist.ac.kr* -- [@noobCoding](http://github.com/noobCoding)
* Jinhwan Kim [@jhk0530](http://github.com/jhk0530)
* Prof. Dougu Nam *dougnam@unist.ac.kr* 

## License
This project is [MIT](https://opensource.org/licenses/MIT) licensed

