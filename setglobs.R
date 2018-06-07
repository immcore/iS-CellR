rm(list = ls())
#.rs.restartR() # Restart R session
# check if pkgs are installed already, if not, install automatically:
source("installPkgsR.R")

# The max limit is 10GB when the app is run locally and 5MB when run from the server.
if(Sys.getenv('SHINY_PORT') == "") options(shiny.maxRequestSize=10000*1024^2)

Logged = FALSE
Uname <- "Anonymous"
values <- reactiveValues(authentication = Uname)

countFile <- reactiveValues(val=NULL)
GeneNames <- reactiveValues(val=NULL)

Featuresobject <- reactiveValues(val=NULL)
GeneListglob <- reactiveValues(val=NULL)
DownloadPlot <- reactiveValues(val=NULL)
CoExprValue <- reactiveValues(val=NULL)
GenesAbsent <- reactiveValues(val=NULL)
QCfilter <-  reactiveValues(val=NULL)
SummaryInfo <- reactiveValues(val=NULL)
PCAClustGlob <- reactiveValues(val=NULL)
PCAClusters <- reactiveValues(val=NULL)
mode <- reactiveValues(n=1,m=1,l=1)
seuratObject <- reactiveValues(val=NULL)
countFileQC <- reactiveValues(val=NULL)
tSNEmatrix <- reactiveValues(val=NULL)
tSNEClusters <- reactiveValues(val=NULL)
tSNE3D <- reactiveValues(val=NULL)
tSNEObj <- reactiveValues(val=NULL)
ClusterLabInfo <- reactiveValues(val=NULL)
dfcluster.ids <- reactiveValues(val=NULL)


options(shiny.maxRequestSize = 6000*1024^2)

source("helpers.R") # Load all the code needed to show feedback on a button click

## options for knitting/rendering rmarkdown chunks
knitr::opts_chunk$set(echo = FALSE, comment = NA, cache = FALSE,
                      message = FALSE, warning = FALSE)

## function to render .md files to html
inclMD <- function(path)
  markdown::markdownToHTML(path, fragment.only = TRUE, options = "", stylesheet = "")

inclRmd <- function(path, r_env = parent.frame()) {
  paste(readLines(path, warn = FALSE), collapse = '\n') %>%
    knitr::knit2html(text = ., fragment.only = TRUE, envir = r_env,  options = "",
                     stylesheet = "") %>%
    gsub("&lt;!--/html_preserve--&gt;","",.) %>%  ## knitr adds this
    gsub("&lt;!--html_preserve--&gt;","",.) %>%   ## knitr adds this
    HTML
  # %>% tagList(., getdeps())
}

## make html table
make_table <- function(dat, width = "50%")
  knitr::kable(dat, align = "c", format = "html",
               table.attr = paste0("class='table table-condensed table-hover' style='width:", width, ";'"))

fixUploadedFilesNames <- function(x) {
  if (is.null(x)) {
    return()
  }
  
  oldNames = x$datapath
  newNames = file.path(dirname(x$datapath),
                       x$name)
  file.rename(from = oldNames, to = newNames)
  x$datapath <- newNames
  x
}

# Seurat steps
selectSteps <- list("Quality control and cell filtering", "Gene variability across single cells", "Linear dimensional reduction (PCA)", "Non-linear dimensional reduction (tSNE)", "Differentially expressed genes", "Discriminating marker genes")

# Packages 
pkgs <- c("threejs_0.3.1","apputils_0.5.1","Seurat_2.2.1","rmarkdown_1.9",
 "stringr_1.3.0","data.table_1.10.4-3","ggrepel_0.7.0","cowplot_0.9.2",
 "gridExtra_2.3","formattable_0.2.0.1","ggjoy_0.4.0","ggridges_0.4.1",
 "psych_1.7.8","evaluate_0.10.1","ggthemes_3.4.0","igraph_1.1.2","pryr_0.1.4",
 "plotly_4.7.1","ggplot2_2.2.1","Matrix_1.2-12","nycflights13_0.2.2","knitcitations_1.0.8",
 "kableExtra_0.7.0","knitr_1.20","dplyr_0.7.4","usethis_1.3.0","devtools_1.13.5.9000",
 "DT_0.4","shinyjs_1.0","shinycssloaders_0.2.0","shinydashboard_0.6.1.9000","shinyBS_0.61",
 "shiny_1.0.5.9000")
