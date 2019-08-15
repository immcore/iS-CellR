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
CRANpkgs <- c("shinyBS", "shinydashboard", "shinydashboardPlus", "shinyFiles", "shinyWidgets", "shinyalert", "htmltools", "shinycssloaders", "shinyjs", "DT", "devtools", "dplyr", "knitr", "kableExtra", "knitcitations", "nycflights13", 
	"Matrix", "plotly", "reticulate", "pryr", "tools", "igraph", "heatmaply", "data.table", "ggthemes", "evaluate", "psych", "ggjoy", "formattable", "gridExtra", "cowplot", "ggrepel", "data.table", "stringr", "rmarkdown")
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

# check if dashboardthemes installed
pkg <- "dashboardthemes"
new.pkg <- pkg[!(pkg %in% installed.packages()[, "Package"])]
if (length(new.pkg)) {
  devtools::install_github("nik01010/dashboardthemes", dependencies=FALSE)
}
sapply(pkg, require, character.only = TRUE)

# check if MAST installed
pkg <- "MAST"
new.pkg <- pkg[!(pkg %in% installed.packages()[, "Package"])]
if (length(new.pkg)) {
    source("https://bioconductor.org/biocLite.R")
    biocLite("MAST")
}
sapply(pkg, require, character.only = TRUE)

# check if UMAP installed
#pkg <- "umap"
#new.pkg <- pkg[!(pkg %in% installed.packages()[, "Package"])]
#if (length(new.pkg)) {
#    reticulate::py_install(c("numpy", "scipy", "scikit-learn", "numba", "umap-learn"))
#}
#sapply(pkg, require, character.only = TRUE)

# check if FI-tSNE installed
FItSNEpath <- reactiveValues(val=NULL)
#FItSNEpath$val <- Sys.which("fast_tsne")[[1]]
FItSNEbin <- "FIt-SNE-master/bin"
if(!file.exists(paste0(FItSNEbin,"/fast_tsne"))){
    #system("git clone https://github.com/FFTW/fftw3.git")
    #system("sh fftw3/bootstrap.sh")
    #system("fftw3/configure --prefix=fftw3")
    #system("make")
    #system("make install")
    #system("wget https://github.com/KlugerLab/FIt-SNE/archive/ec25f1b36598a2d21869d10a258ac366a12f0b05.zip")
    system("unzip FIt-SNE-master.zip")
    #system("rm -r v1.0.0.zip")
    #system("git clone https://github.com/KlugerLab/FIt-SNE.git")
    system("g++ -std=c++11 -O3 FIt-SNE-1.0.0/src/sptree.cpp FIt-SNE-1.0.0/src/tsne.cpp FIt-SNE-1.0.0/src/nbodyfft.cpp -o FIt-SNE-1.0.0/bin/fast_tsne -pthread -lfftw3 -lm")
    FItSNEpath$val <- as.character("FIt-SNE-master/bin")
} else {
    FItSNEpath$val <- as.character("FIt-SNE-master/bin")
}


