#!/usr/bin/Rscript

# Check and install if the package is not installed and then load them into the R session.

instpkg <- function(pkg,repo){
    new.pkg <- pkg[!(pkg %in% installed.packages()[, "Package"])]
    if (length(new.pkg) && repo == 'CRAN') {
        install.packages(new.pkg, dependencies = TRUE, repos='https://cloud.r-project.org/')
    }
    sapply(pkg, require, character.only = TRUE)
}


# CRAN R packages
CRANpkgs <- c("shinyBS", "shinydashboard", "shinycssloaders", "shinyjs", "DT", "devtools", "dplyr", "knitr", "kableExtra", "knitcitations", "nycflights13", 
	"Matrix", "plotly", "pryr", "igraph", "ggthemes", "evaluate", "psych", "ggjoy", "formattable", "gridExtra", "cowplot", "ggrepel", "data.table", "stringr", "rmarkdown")
instpkg(CRANpkgs, "CRAN")

# check if Dev Shiny installed
pkg <- "shiny"
new.pkg <- pkg[!(pkg %in% installed.packages()[, "Package"])]
if (length(new.pkg)) {
    devtools::install_github("rstudio/shiny", dependencies=FALSE)
    }
sapply(pkg, require, character.only = TRUE)

# check if ggplot2 installed
pkg <- "ggplot2"
new.pkg <- pkg[!(pkg %in% installed.packages()[, "Package"])]
if (length(new.pkg)) {
    devtools::install_github("hadley/ggplot2", dependencies=FALSE)
    }
sapply(pkg, require, character.only = TRUE)

# check if Seurat installed
pkg <- "Seurat"
new.pkg <- pkg[!(pkg %in% installed.packages()[, "Package"])]
if (length(new.pkg)) {
    devtools::install_github("satijalab/seurat", ref = "develop", dependencies=TRUE)
    }
sapply(pkg, require, character.only = TRUE)

# check if leonawicz/apputils installed
pkg <- "apputils"
new.pkg <- pkg[!(pkg %in% installed.packages()[, "Package"])]
if (length(new.pkg)) {
    devtools::install_github("leonawicz/apputils", dependencies=TRUE)
    }
sapply(pkg, require, character.only = TRUE)

# check if Seurat installed
pkg <- "threejs"
new.pkg <- pkg[!(pkg %in% installed.packages()[, "Package"])]
if (length(new.pkg)) {
    devtools::install_github("bwlewis/rthreejs", dependencies=FALSE)
    }
sapply(pkg, require, character.only = TRUE)


