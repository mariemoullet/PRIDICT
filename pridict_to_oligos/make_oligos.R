#!/usr/bin/env Rscript

library(tidyverse)
library(spgs)
library(janitor)
library(digest)

# Uncomment below if running in R studio: 
# setwd(dirname(rstudioapi::getActiveDocumentContext()$path))

# filter argument is TRUE (default) when we want to limit pegRNAs to 200 bp
args = commandArgs(trailingOnly=TRUE)

if (args[1] == TRUE){
    guides_to_keep = 50 # we keep the top 50 and within that keep the top 3 which are <= 200 bp
} else {
    guides_to_keep = 3
}

output_file <- args[2]

# Define functions --------------------------------------------------------


rc <- function(x) {toupper(spgs::reverseComplement(x))}

wrangle_pridict <- function(df, target_mut, keep_top = guides_to_keep){
  # df is pridict output csv file
  # get top 3 ranking pegRNAs 
  df <- df %>%
    arrange(rank) %>%
    head(n = keep_top) %>%
    mutate(Name = paste(target_mut, rank, sep = "_"))
  
  # get oligos 
  df <- df %>% select(Name, protospacer_sequence, extension_oligo_fw)
  
  # remove constant regions for extension oligo
  df <- df %>%
    mutate(extension_oligo_fw = str_remove(extension_oligo_fw, "GTGC"))
  
  # rename columns to match Jonas syntax 
  df <- df %>%
    rename(Spacer = protospacer_sequence, 
           Extension = extension_oligo_fw)
  
  return(df)
}

make_oligos <- function(x, design = "epegRNA") {
  # There are three design option available, classic (for transfection, normal spacer), lentiviral (using the enhanced spacer), and epegRNA (both for transfection and lentivirus, enhanced spacer)
  if(design == "classic") {
    output_df <- mutate(x,
                        ExF = paste0("gtgc", Extension),
                        ExR = paste0("aaaa", rc(Extension)),
                        SpF = paste0("cacc", Spacer, "gtttt"),
                        SpR = paste0("ctctaaaac", rc(Spacer))) %>%
      select(Name, ExF, ExR, SpF, SpR) %>%
      pivot_longer(cols = c(ExF, ExR, SpF, SpR), names_to = "Primer", values_to = "Oligo_sequence") %>%
      unite("Primer_name", c(Name, Primer), sep = "_")
  }
  else if(design == "epegRNA") {
    output_df <- mutate(x,
                        ExF = paste0(toupper("gtgc"), Extension),
                        ExR = paste0(toupper("cgcg"), rc(Extension)),
                        SpF = paste0(toupper("cacc"), Spacer, toupper("gttta")),
                        SpR = paste0(toupper("ctcttaaac"), rc(Spacer))
    ) %>%
      select(Name, Extension, Spacer, ExF, ExR, SpF, SpR) %>%
      pivot_longer(cols = c(Extension, Spacer, ExF, ExR, SpF, SpR), names_to = "Primer", values_to = "Oligo_sequence") %>%
      unite("Primer_name", c(Name, Primer), sep = "_")
  }
  else if(design == "lentiviral") {
    output_df <- mutate(x,
                        ExF = paste0("gtgc", Extension, "tttttttaagc"),
                        ExR = paste0("gtaggcttaaaaaaa", rc(Extension)),
                        SpF = paste0("cacc", Spacer, "gttta"),
                        SpR = paste0("ctcttaaac", rc(Spacer))
    ) %>%
      select(Name, Extension, Spacer, ExF, ExR, SpF, SpR) %>%
      pivot_longer(cols = c(Extension, Spacer, ExF, ExR, SpF, SpR), names_to = "Primer", values_to = "Oligo_sequence") %>%
      unite("Primer_name", c(Name, Primer), sep = "_")
    
  } else {print("Design should be classic, lentiviral or epegRNA")}
}

get_pegrna_from_predict <- function(target_mut) {
  # Primer_name is based on pridict rank
  # oligo name is based on unique primers 
  df <- read_csv(paste0("PRIDICT/predictions/", target_mut, "_pegRNA_Pridict_full.csv"), 
                 show_col_types = FALSE) %>%
    clean_names() %>%
    wrangle_pridict(target_mut) %>%
    make_oligos() %>%
    mutate(length = nchar(Oligo_sequence)) 
  
  return(df)
}

make_pegrna_table <- function(df) {
  df2 <- df %>%
    select(-length) %>%
    separate(col = Primer_name, 
             into = c("gene", "mutation", "pridict_rank", "primer_cat"), 
             sep = "_") %>%
    # only keep forward primers 
    filter(str_ends(primer_cat, "F")) %>%
    # create two columns : ExF and SpF
    pivot_wider(names_from = primer_cat, values_from = Oligo_sequence) %>%
    # using improved scaffold FW from Jonas protocol
    mutate(scaffold = "agagCTATGCTGGAAACAGCATAGCAAGTTTAAATAAGGCTAGTCCGTTATCAACTTGAAAAAGTGGCACCGAGTCG") %>%
    # assemble pegRNA 
    unite("pegRNA", c(SpF, scaffold, ExF), sep = "") %>%
    mutate(`for Gibson assembly` = paste("TATCTTGTGGAAAGGACGAAA", 
                                         pegRNA, 
                                         "CGCGGTTCTATCTAGTTACGCGT", 
                                         sep = "")) %>%
    mutate(length = nchar(pegRNA)) %>%
    mutate(`length with GA overhang` = nchar(`for Gibson assembly`))
  
  return(df2)
}

add_genomic_context <- function(df) {
    # Read the intermediate file
    intermediate_file <- read.csv('intermediate/pridict_input.csv', stringsAsFactors = FALSE)

    # Function to wrangle mutations
    rm_mut <- function(x) {
    x <- gsub("\\(./", "", x)
    x <- gsub("\\)", "", x)
    return(x)
    }

    # Apply the function to editseq column
    intermediate_file$editseq <- sapply(intermediate_file$editseq, rm_mut)

    # Rename editseq column to genomic_context
    colnames(intermediate_file)[colnames(intermediate_file) == "editseq"] <- "genomic_context"  

    # Create sequence_name column in df
    df$sequence_name <- paste(df$gene, df$mutation, sep = "_")

    # Join df with intermediate_file based on sequence_name
    df <- df %>%
    left_join(intermediate_file, by = "sequence_name")

    # Drop sequence_name column
    df <- df %>%
    select(-sequence_name)

    # Rename editseq column

    return(df)
}

# Run ---------------------------------------------------------------------

# Get all available mutations:
muts <- list.files(path = "PRIDICT/predictions/", 
                  pattern = "csv") 

muts <- str_remove(muts, "_pegRNA_Pridict_full.csv")

pridict_list <- list()

for (m in muts) {
  pridict_list[[m]] <- get_pegrna_from_predict(m)
} 

oligos <- bind_rows(pridict_list) %>%
  filter(! grepl("Extension", Primer_name)) %>%
  filter(! grepl("Spacer", Primer_name)) 

# table with full pegRNA
pegrna_table <- make_pegrna_table(oligos) 

# optional : %>%
    # remove primers > 200 bp
if (args[1] == TRUE){
  
  pegrna_table <- pegrna_table %>%
    filter(`length with GA overhang` <= 200) %>%
    group_by(gene, mutation) %>%
    slice_min(as.numeric(pridict_rank), n = 3) 

}
    
# Add genomic context
pegrna_table <- add_genomic_context(pegrna_table)
  
# add unique identifier
pegrna_table$id <- apply(pegrna_table, 1, function(x) digest(x, algo = "crc32"))

write_csv(pegrna_table, output_file)








