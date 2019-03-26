# FCSGenerator2
Shiny App meant to generate highly customizable FCS files from scratch or use cytometry data to generate "ground truth" files which can then be used as a reference when benchmarking alrogithms.
 
	
## Requirements
  * software: R(Version 3.5.1), Rstudio(optional)
  * R packages: shiny, shinydashboard, shinyjs, flowCore, ggplot2, d3r, reshape2, ggridges, truncnorm, sunburstR, heatmaply, shinyHeatmaply
  
## Quick installation guide

  1. Run the following command in R/RStudio:
```
install.packages(c("shiny", "shinyjs", "shinydashboard", "ggplot2", "d3r", "reshape2", "ggridges", "truncnorm", "sunburstR", "heatmaply", "shinyHeatmaply"))
source("https://bioconductor.org/biocLite.R")
biocLite("Biobase")
biocLite("flowCore")
```
  >You may be asked to reload your environment, if so, accept.
  
  2. Run the next commands:
```
library("devtools")
install_github("isambens/fcsgenerator2")
```

  
## Launching the shiny application

  Run the following commands in R/RStudio:
```
library("FCSGenerator2")
FCSGenerator2.run()
```