install_packages <- function() {
  install.packages("Seurat")
  install.packages("hdf5r")
  install.packages("patchwork")
  
  install.packages("tidyverse")

}

load_packages <- function() {
  library(Seurat)
  library(hdf5r)
  library(patchwork)
  
  library(tidyverse)
}
