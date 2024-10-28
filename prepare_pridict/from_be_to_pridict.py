import pandas as pd
from ucsc.api import Sequence
import sys
import numpy as np
from tqdm import tqdm
from utils import *

tqdm.pandas()

# Load data
be_in = pd.read_csv("input_from_be.csv")

# Wrangle into long format with one row per gene, chromosome, position

df = be_in[["Gene", "chr", "Edit_Location"]].dropna()

df = df.rename(columns={"chr": "chromosome", "Gene": "gene"})

df["chromosome"] = df["chromosome"].astype(int)

df["position"] = df["Edit_Location"].str.split(";")
df = df.drop(columns="Edit_Location").explode("position")

# focus on single nt changes for now
single_nt = df[~df.position.str.contains("-")].copy()

# drop chromosome information to keep only position
single_nt['position'] = single_nt['position'].str.split(':').str[-1].astype(int)
single_nt["start_seq"] = single_nt.position - 100
single_nt["end_seq"] = single_nt.position + 100


# get DNA sequence
single_nt["sequence"] = single_nt.apply(
    lambda row: get_sequence(row["chromosome"], row["start_seq"], row["end_seq"], '+'),
    axis=1,
)

# Make edited sequences
# (columns called A and B for pandas explode and then renamed)
single_nt["A"] = single_nt["sequence"].apply(create_editseq)
single_nt["B"] = single_nt.apply(
    lambda row: create_hvgs(row["sequence"], row["chromosome"], row["position"]), axis=1
)

single_nt = single_nt.explode(list("AB"))
single_nt = single_nt.rename(columns={"A": "editseq", "B": "hvgs"})
single_nt = single_nt.reset_index()

single_nt["vep_out"] = single_nt["hvgs"].progress_apply(get_vep)
vep_df = single_nt["vep_out"].str.split("_", expand=True)

vep_df.columns = ["aa_change", "codon_change", "aa_position"]
vep_df["sequence_name"] = vep_df.apply(
    lambda row: get_sequence_label(row["aa_change"], row["aa_position"]), axis=1
)

single_nt.drop(columns="vep_out", inplace=True)
single_nt = single_nt.join(vep_df)

# Output for pridict
#### to get gene name for R wrangling, edit sequence name to include gene name
intermediate = single_nt[["editseq", "sequence_name", "gene"]].copy()

intermediate["sequence_name"] = (
    intermediate["gene"] + "_" + intermediate["sequence_name"]
)

# Make sure sequence name is unique:
# add suffix to duplicate values to make them unique
# use groupby() and cumcount() to add an index that represents the counts of each non-unique element

data = intermediate["sequence_name"]

# use groupby() and cumcount() to add an index that represents the counts of each non-unique element
counts = data.groupby(data).cumcount() + 1
mask = data.duplicated(keep=False)
tmp = pd.DataFrame({"seq": data, "counter": counts.astype(str)})
tmp["seq_out"] = np.where(mask, tmp["seq"] + "mut" + tmp["counter"], tmp["seq"])

intermediate["sequence_name"] = tmp["seq_out"]

intermediate.to_csv("intermediate/pridict_input.csv", index=False)