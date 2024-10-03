# Load packages -----------------------------------------------------------
library(tidyverse)
library(ensembldb)
library(AnnotationDbi)
library(EnsDb.Hsapiens.v86)
library(org.Hs.eg.db)

# Get input gene name from command line argument ---------------------------
args <- commandArgs(trailingOnly = TRUE)  # Retrieve command line arguments
input_name <- args[1]  # Get the first argument as the input gene name

# Define database ---------------------------------------------------------
# Extract the mapping between gene symbols, Ensembl gene IDs, and protein IDs
gene_protein_map <- select(EnsDb.Hsapiens.v86,
                           columns = c("SYMBOL", "GENEID", "PROTEINID"),
                           keys = keys(EnsDb.Hsapiens.v86, keytype = "GENEID"),
                           keytype = "GENEID")

# Check input format (symbol or Ensembl ID) -------------------------------
# If input_name is a valid Ensembl gene ID or gene symbol, proceed with filtering

is_ensembl <- grepl("^ENSG", input_name)  # Check if input is Ensembl ID

if (is_ensembl) {
  gene_col <- "GENEID"
} else {
  gene_col <- "SYMBOL"
}

# Return output -----------------------------------------------------------
# Filter the data by gene symbol or Ensembl gene ID and print the result
result <- gene_protein_map %>%
  dplyr::filter(!!sym(gene_col) == input_name)

# Print result
if (nrow(result) > 0) {
  print(result)
} else {
  cat("No matching gene or protein ID found for", input_name, "\n")
}
