install_packages <- function() {
  install.packages("tidyverse")
  install.packages("Seurat")
  install.packages("hdf5r") # h5 file handling
  install.packages("patchwork")
  
  install.packages("ggpubr")

}

load_packages <- function() {
  library(Seurat)
  library(hdf5r)
  library(patchwork)
  
  # plotting
  library(ggpubr)
  
  # Load last to make sure any conflicts are overriden by the tidyverse
  library(tidyverse)
}
