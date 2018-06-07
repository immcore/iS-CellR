# iS-CellR
iS-CellR (Interactive platform for Single-cell RNAseq) is a web-based Shiny application designed to provide a comprehensive analysis of single-cell RNA sequencing data. iS-CellR provides a fast method for filtering and normalization of raw data, dimensionality reductions (linear and non-linear) to identify cell types clusters, differential gene expression analysis to locate markers, and inter-/intra-sample heterogeneity analysis. iS-CellR integrates the Seurat package with Shiny's reactive programming framework and interactive visualization using plotly library. iS-CellR runs on any modern web browser and provides access to powerful R libraries through a graphical user interface. Each session of iS-CellR allows the user to share, reproduce and archive results without requiring programming skills.

## Getting started

iS-CellR pipeline overview is illustrated in the figure. iS-CellR is organized into a seven-step process for complete scRNA-seq analysis.

<img src=iS-CellR_workflow.png height="800">


## How to get strated and load data:

## Prerequisite
```{r}
1). Install R (v >= 3.2)
Download and install R from http://cran.us.r-project.org/  

2). Recommended : Install R Studio
Download and install RStudio Desktop from http://rstudio.org/download/desktop

3). Install the "devtools" package from Hadley Wickham
install.packages(“devtools”) # inside R console

4). Install Shiny package
install.packages(“shiny”) # inside R console
```

## Using R
```{r}
install.packages("devtools")
devtools::install_github("rstudio/shiny", dependencies=FALSE)

runUrl('https://github.com/immcore/iS-CellR/archive/master.zip')
# or
shiny::runGitHub("iS-CellR", "immcore")
# or
runApp(/fullpath/iS-CellR) # provide full path to iS-CellR folder
```


## Using Docker
To run container use the command below:

```sh
docker run --rm -p 3838:3838 immcore/is-cellr 
```

After that check with your browser at addresses plus the port 3838 : http://0.0.0.0:3838/

## Troubleshooting
1). Seurat installaton
```{r}
Warning: Installed Rcpp (0.12.12) different from Rcpp used to build dplyr (0.12.11).
Please reinstall dplyr to avoid random crashes or undefined behavior. 

install.packages("dplyr", type = "source")
library(dplyr)
```

If you are using R >= 3.4 on a Mac, the latest versions of R use Clang 4.0.0 and GNU Fortran 6.1 which you will need to successfully compile and install Seurat.
```{r}
curl -O http://r.research.att.com/libs/gfortran-4.8.2-darwin13.tar.bz2
sudo tar fvxz gfortran-4.8.2-darwin13.tar.bz2 -C /
```
##
__*It is important that the input file should follow the same format as descibed below.*__

To get started, please load in a CSV/TSV separated values file. The file should contain the count data where:

1. Columns are Cells
2. Rows are Genes

## Demo data files:

There is one matrix file for Malignant dataset from Tirosh et al., 2016 (Main file for the analysis) and genes file (For STEP: Discriminating marker genes). 

### Count matrix: Maligant50.csv
### Genes file: Genes_MITF_AXL.txt
