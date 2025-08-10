install_packages <- function() {
  install.packages("tidyverse")
  install.packages("Seurat")
  install.packages("hdf5r") # h5 file handling
  install.packages("patchwork")
  install.packages("BiocManager")
  
  BiocManager::install("speckle")
  BiocManager::install("SingleCellExperiment")
  BiocManager::install("CellBench")
  BiocManager::install("limma")

  BiocManager::install("scater")
  install.packages("devtools")

  BiocManager::install("edgeR")
  install.packages("statmod")
  
  install.packages("ggpubr")

}

load_packages <- function() {
  library(Seurat)
  library(hdf5r)
  library(patchwork)
  
  library(speckle)
  library(SingleCellExperiment)
  library(CellBench)
  library(limma)
  # library(ggplot2)
  library(scater)
  library(patchwork)
  library(edgeR)
  library(statmod)
  
  # plotting
  library(ggpubr)
  
  # Load last to make sure any conflicts are overriden by the tidyverse
  library(tidyverse)
}
