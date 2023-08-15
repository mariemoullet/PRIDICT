if (!require("BiocManager", quietly = TRUE))
   install.packages("BiocManager")

BiocManager::install("ensembldb")
BiocManager::install("EnsDb.Hsapiens.v86")
BiocManager::install("org.Hs.eg.db")

install.packages("tidyverse")
install.packages("spgs")
install.packages("janitor")
install.packages("digest")
