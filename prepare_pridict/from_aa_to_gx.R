#!/usr/bin/env Rscript

# Uncomment for Rstudio
# setwd(dirname(rstudioapi::getActiveDocumentContext()$path))


# Download packages -------------------------------------------------------


# if (!require("BiocManager", quietly = TRUE))
#   install.packages("BiocManager")

# BiocManager::install("ensembldb")
# BiocManager::install("EnsDb.Hsapiens.v86")
# BiocManager::install("org.Hs.eg.db")

# Load packages -----------------------------------------------------------

library(tidyverse)
library(ensembldb)
library(AnnotationDbi)
library(EnsDb.Hsapiens.v86)
library(org.Hs.eg.db)

# Load data ---------------------------------------------------------------

input_df <- read_csv("input_from_aa.csv")

# Define functions --------------------------------------------------------


get_gnm_from_prt <- function(gene, aa, chrom){
  
  aa <- as.numeric(aa)
  
  edbx <- filter(EnsDb.Hsapiens.v86, filter = ~ seq_name == chrom)
  
  prt <- IRanges(start = aa, end = aa, names = gene)
  
  gnm <- proteinToGenome(prt, edbx)
  
  # -1 to convert to UCSC indexing
  return(gnm[[gene]]@ranges@start - 1)
  
}

# Get genomic positions ---------------------------------------------------

input_df$codon_start <- apply(input_df, 1, function(x) get_gnm_from_prt(x["id"], # use ensembl id  
                                                                        x["aa"], 
                                                                        x["chrom"]))
input_df$codon_end <- input_df$codon_start + 2

input_df$seq_start <- input_df$codon_start - 100
input_df$seq_end <- input_df$codon_end + 100

write_csv(input_df, "intermediate/input_gx_location.csv")






