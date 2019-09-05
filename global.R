rm(list = ls())

#.rs.restartR() # Restart R session
# check if pkgs are installed already, if not, install automatically:
source("installPkgsR.R")
source("SwitchButton.R")
#options(shiny.autoreload = 100)
source("helpers.R") # Load all the code needed to show feedback on a button click

selectSteps <- list("Quality control and cell filtering", "Gene variability across single cells", "Linear dimensional reduction (PCA)", "Non-linear dimensional reduction (UMAP/tSNE)", "Differentially expressed genes", "Discriminating marker genes")

# The max limit is 10GB when the app is run locally and 5MB when run from the server.
options(shiny.maxRequestSize=10000*1024^2)
