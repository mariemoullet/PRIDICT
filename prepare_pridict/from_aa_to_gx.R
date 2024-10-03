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
                                                                        x["chrom"] 
                                                                        ))
# Identify rows where codon_start contains NA values
na_rows <- sapply(input_df$codon_start, function(x) any(is.na(x)))

# If there are any rows with missing values, print the warning
if (any(na_rows)) {
  # Print a warning for the dropped rows (id and aa)
  cat(sprintf("Dropping %d rows:\n", sum(na_rows)))
  
  # Print the specific rows that will be dropped
  dropped_rows <- input_df[na_rows, c("id", "aa")]
  apply(dropped_rows, 1, function(row) {
    cat(sprintf("id: %s, aa: %s\n", row["id"], row["aa"]))
  })
}

# Filter the input_df to remove the rows with missing codon_start values
input_df <- input_df[!na_rows, ]
# Check for multiple values in codon_start and print them
multiple_values_rows <- sapply(input_df$codon_start, function(x) length(x) > 1)

if (any(multiple_values_rows)) {
  cat("The following rows have multiple values (corresponding to different exons) in codon_start:\n")
  apply(input_df[multiple_values_rows, ], 1, function(row) {
    codon_values <- paste(as.character(row$codon_start), collapse = ", ")
    cat(sprintf("id: %s, aa: %s, codon_start: %s\n", row["id"], row["aa"], codon_values))
  })
  cat('Selecting the first exon \n')
}

# If codon_start has multiple values, take the first value only
input_df$codon_start <- sapply(input_df$codon_start, function(x) as.numeric(x[1]))
            
input_df$codon_end <- input_df$codon_start + 2

input_df$seq_start <- input_df$codon_start - 100
input_df$seq_end <- input_df$codon_end + 100

# 
# # Function to adjust codon start and end based on strand
# swap_start_end_on_negative_strand <- function(strand, codon_start, codon_end, seq_start, seq_end) {
#   if (strand == "-") {
#     # Flip start and end for the negative strand
#     return(data.frame(new_values_codon_start = codon_end, 
#                       new_values_codon_end = codon_start,
#                       new_values_seq_start = seq_end, 
#                       new_values_seq_end = seq_start))
#   } else {
#     return(data.frame(new_values_codon_start = codon_start, 
#                       new_values_codon_end = codon_end,
#                       new_values_seq_start = seq_start,
#                       new_values_seq_end = seq_end))
#   }
# }
# 
# 
# # Apply function to each row of the dataframe
# output_df <- input_df %>%
#   rowwise() %>%
#   mutate(new_values = list(swap_start_end_on_negative_strand(strand, codon_start, codon_end, seq_start, seq_end))) %>%
#   unnest(cols = c(new_values)) %>%
#   ungroup() %>%
#   as.data.frame() %>%
#   dplyr::select(gene, id, chrom, aa, strand, 
#                 new_values_codon_start, new_values_codon_end, new_values_seq_start, new_values_seq_end) %>%
#   dplyr::rename(codon_start = new_values_codon_start, 
#          codon_end = new_values_codon_end,
#          seq_start = new_values_seq_start, 
#          seq_end = new_values_seq_end)


write_csv(input_df, "intermediate/input_gx_location.csv")






