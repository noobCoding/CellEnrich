# CellEnrich 

Pathway enrichment analysis/visualization for Single Cell Data

<img src="images/LPS_exp_freq.png"> 
<img src="images/LPS_exp_biplot.png" width="1024"> 

## Install

```R
if(!require(remotes)){
  install.packages('remotes') 
}
remotes::install_github('vqv/ggbiplot')
install.packages('waiter')
install.packages('farver')
remotes::install_github('noobCoding/CellEnrich')
```

## Example with primary mouse dendritic cells (DCs) stimulated with lipopolysaccharide (LPS)

```R
# download minimal data to the working directory
download.file('https://github.com/noobcoding/CellEnrich/raw/master/data/LPS_exp.rds','LPS_exp.rds', mode = 'wb')
download.file('https://github.com/noobcoding/CellEnrich/raw/master/data/WikiPathways_2019_Mouse.RData', 'WikiPathways_2019_Mouse.RData', mode = 'wb')

# Load library and data
library(CellEnrich)

LPS_exp <- readRDS("LPS_exp.rds")
CountData <- LPS_exp$counts
GroupInfo <- LPS_exp$meta
 
# Run cellenrich
CellEnrich(CountData, GroupInfo)

```
## Example with PBMC_3K data 

The online manual and example with PBMC 3K data are available [here](https://github.com/noobCoding/CellEnrich/wiki)

## Example with Alzheimer's data 

The online manual and example with Alzheimer's data are available [here](https://noobcoding.github.io/CellEnrich/)


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

