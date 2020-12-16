load.lib<-c("devtools","vegan", "dada2","phyloseq","gridExtra","RColorBrewer",
            "viridis","ggpubr","broom","decontam","indicspecies","lemon","ggalluvial",
            "pheatmap","metagenomeSeq","qiime2R","metagMisc","limma", "Biostrings","DESeq2",
            "edgeR", "caret", "UpSetR", 
            "remotes","readxl","googlesheets","haven", "readr", "rio", "Hmisc", "sqldf", "jsonlite", 
            "XML", "httr", "quantmod", "tidyquant", "rvest", "dplyr", "purrr", "tidyr", 
            "magrittr", "validate", "testthat", "data.table", "stringr", "lubridate", "zoo", "editR", 
            "knitr", "officer", "listviewer", "DT", "ggplot2", "ggiraph", "dygraphs", "googleVis", 
            "metricsgraphics", "sf", "leaflet", "ggmap", "tmap", "tmaptools", "mapsapi", 
            "tidycensus", "glue", "rga", "RSiteCatalyst", "roxygen2", "shiny", "flexdashboard", "openxlsx", 
            "gmodels", "janitor", "car", "rcdimple", "foreach", "scales", "plotly", "highcharter", "profvis", 
            "tidytext", "diffobj", "Prophet", "feather", "fst", "googleAuthR", "cloudyR")


install.lib<-load.lib[!load.lib %in% installed.packages()]
for(lib in install.lib) install.packages(lib,dependencies=TRUE)
sapply(load.lib,require,character=TRUE)

install.packages('vegan')
install_github("microbiome/microbiome") 
devtools::install_github("vmikk/metagMisc")
devtools::install_github("jbisanz/qiime2R")
devtools::install_github("trestletech/shinyAce")
devtools::install_github("swarm-lab/editR")
if (!requireNamespace("BiocManager", quietly=TRUE))
  install.packages("BiocManager")
BiocManager::install("metagenomeSeq")
BiocManager::install("edgeR")
BiocManager::install("DESeq2")
BiocManager::install("decontam")
BiocManager::install("dada2")
BiocManager::install("limma")
BiocManager::install("Biostrings")
install.packages("gt")
