from utils import * 
import pandas as pd
from ucsc.api import Sequence 

# read file prepared in from_aa_to_gx.R 
df = pd.read_csv("intermediate/input_gx_location.csv") 

# get codon 
df["codon"] = df.apply(lambda row: get_sequence(row['chrom'],
                                                row['codon_start'], 
                                                row['codon_end']+1,
                                                row['strand']
                                                ), axis=1)

df["sequence"] = df.apply(lambda row: get_sequence(row['chrom'],
                                                   row['seq_start'], 
                                                   row['seq_end']+1,
                                                   row['strand']
                                                   ), axis=1)

# columns need to be 
# gene codon position_label sequence
df = df.rename(columns = {"aa" : "position_label"})
df = df[["gene", "codon", "position_label", "sequence"]]

# add brackets around codon
df["sequence"] = df["sequence"].apply(add_brackets)

# Convert to pridict input format 
# columns are editseq and sequence_name
# where each sequence corresponds to new aa generated from single aa change
output_df = pd.concat([make_pridict_input(row) for _, row in df.iterrows()], ignore_index=True)

output_df.to_csv("intermediate/pridict_input.csv", index = False)


