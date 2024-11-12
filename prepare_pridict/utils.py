import pandas as pd
from ucsc.api import Sequence
import sys, requests, json
import re
import numpy as np


#################################
#### Wrangle input sequence #####
#################################

BASE_URL = 'http://api.genome.ucsc.edu'

def raiseExceptionOfRequest(response):
    if response.get('statusCode') == 400:
        raise NotFoundException('Something went wrong, ' + response.get('error'))

    if response.get('statusCode') == 403:
        raise NotAllowedException('The Requested Resource is not allowed to be accessed' + response.get('error'))

    if response.get('error') is not None:
        raise NotFoundException('An error happened, ' + response.get('error'))

    pass

class Sequence:
    '''
    adapted from 
    https://github.com/Eyadhamza/UCSC-Genomic-REST-Api-Wrapper/blob/main/ucsc/api.py#L327
    '''
    requestUrl = ''
    requestParams = {}
    responseKey = ''

    def __init__(self, genome, chrom, dna=None, hub=None, track=None, start=None, end=None,**kwargs):
        self.end = end
        self.start = start
        self.track = track
        self.hub = hub
        self.chrom = chrom
        self.dna = dna
        self.genome = genome

    @classmethod
    def get(cls, genome, chrom, hubUrl=None, start=None, end=None, revComp=0):
        cls.requestUrl = BASE_URL + '/getData/sequence'
        cls.requestParams = {'hubUrl': hubUrl, 'genome': genome, 'chrom': chrom, 'start': start, 'end': end, 'revComp' : revComp}
        response = requests.get(cls.requestUrl, cls.requestParams).json()
        raiseExceptionOfRequest(response)
        return Sequence(**response)

def get_sequence(chromosome, start_seq, end_seq, strand):
    """
    get genomic DNA sequence from chromosome number, start and end position
    """
    if (strand == '+'):
        rev_comp = 0
    else:
        rev_comp = 1
    
    sequence = Sequence.get(
        genome="hg38", chrom=chromosome, start=start_seq, end=end_seq, revComp=rev_comp
    )
    
    return sequence.dna

def add_brackets(seq):
    """
    for amino acid version of script, add brackets around amino acid of interest
    """
    return seq[:100] + "(" + seq[100:103] + ")" + seq[103:]


def modify_string(row):
    """
    modify the string based on character position
    """
    edit = ""
    for c in range(len(row["Original Codon"])):
        if c == row["Nucleotide Index"] - 1:  # -1 for zero indexing
            edit = edit + "(" + row["Change"].replace(">", "/") + ")"
        else:
            edit += row["Original Codon"][c]

    return edit

def replace_substring(text, replacement, pattern=r"\((.*?)\)"):
    return re.sub(pattern, replacement, text)

def replace_bycolumn(row):
    return replace_substring(row["sequence"], row["pridinput"])

def get_sequence_label(aa_change, aa_position):
    old_aa = aa_change.split("/")[0]

    if "/" in aa_change:
        new_aa = aa_change.split("/")[1]
    else:
        new_aa = ""

    return old_aa + str(aa_position) + new_aa


def create_editseq(seq):
    """
    Replace position 100 with (N/M) where M is all other possible single nt substitutions

    """

    N = seq[99]

    new_sequences = []

    for m in ["A", "T", "C", "G"]:
        if m != N:
            edit = seq[:99] + "(" + N + "/" + m + ")" + seq[100:]
            new_sequences.append(edit)

    return new_sequences

####################################################
###### Functions from AA_Codon_Investigator.py #####
####################################################

class Aminoacid(object):
    def __init__(self, aa, letter, codon_list):
        self.aa = aa
        self.letter = letter
        self.codon_list = codon_list


class Codon(object):
    def __init__(self, codon, aa):
        self.aa = aa
        self.codon = codon


aa_dict = {
    "PHE": "F",
    "LEU": "L",
    "ILE": "I",
    "MET": "M",
    "VAL": "V",
    "SER": "S",
    "PRO": "P",
    "THR": "T",
    "ALA": "A",
    "TYR": "Y",
    "STOP": "STOP",
    "HIS": "H",
    "GLN": "Q",
    "ASN": "N",
    "LYS": "K",
    "ASP": "D",
    "GLU": "E",
    "CYS": "C",
    "TRP": "W",
    "ARG": "R",
    "GLY": "G",
}

aa_df = pd.DataFrame.from_dict(aa_dict, orient="index")
aa_df.columns = ["single_letter"]

codon_dict = {
    "PHE": ["TTT", "TTC"],
    "LEU": ["TTA", "TTG", "CTT", "CTC", "CTA", "CTG"],
    "ILE": ["ATT", "ATC", "ATA"],
    "MET": ["ATG"],
    "VAL": ["GTT", "GTC", "GTA", "GTG"],
    "SER": ["TCT", "TCC", "TCA", "TCG", "AGT", "AGC"],
    "PRO": ["CCT", "CCC", "CCA", "CCG"],
    "THR": ["ACT", "ACC", "ACA", "ACG"],
    "ALA": ["GCT", "GCC", "GCA", "GCG"],
    "TYR": ["TAT", "TAC"],
    "STOP": ["TAA", "TAG", "TGA"],
    "HIS": ["CAT", "CAC"],
    "GLN": ["CAA", "CAG"],
    "ASN": ["AAT", "AAC"],
    "LYS": ["AAA", "AAG"],
    "ASP": ["GAT", "GAC"],
    "GLU": ["GAA", "GAG"],
    "CYS": ["TGT", "TGC"],
    "TRP": ["TGG"],
    "ARG": ["CGT", "CGC", "CGA", "CGG", "AGA", "AGG"],
    "GLY": ["GGT", "GGC", "GGA", "GGG"],
}


aa_info = {}

for aa, letter in aa_dict.items():
    obj = Aminoacid(aa=aa, letter=letter, codon_list=codon_dict[aa])
    aa_info[aa] = obj

codon_info = {}

for aa, codon_list in codon_dict.items():
    for codon in codon_list:
        obj = Codon(codon=codon, aa=aa)
        codon_info[codon] = obj


def single_change_from_seq(seq):
    """
    :param seq: codon of interest : this is what will be modified
    :return: A DataFrame and all possible amino acids that can be created with single nucleotide change ***from the starting codon***
    """

    # Get amino acid from codon
    seq = seq.upper()
    aa = [key for key, value in codon_dict.items() if seq in value][0]

    # Take Amino acid object
    obj = aa_info[aa.upper()]

    # Specify the columns for the DataFrame (will be created for result)
    df_columns = ["Original Codon", "Nucleotide Index", "Change", "New AA", "New Codon"]

    lines = []
    for codon in obj.codon_list:
        # For each codon of the input amino acid
        for k in range(3):
            cd = list(codon)

            # For each nucleotide change for each position in the codon
            for n in ["T", "C", "A", "G"]:
                cd[k] = n
                jcd = "".join(cd)

                # Look if there will be a codon representing any other amino acid
                # if codon_info[jcd].aa != aa: this line removes synonymous information
                if jcd != seq:  # just filter for redundant changes ie dont include A>A
                    lines.append(
                        {
                            "Original Codon": codon,
                            "Nucleotide Index": k + 1,
                            "Change": list(codon)[k] + ">" + n,
                            "New AA": codon_info[jcd].aa,
                            "New Codon": jcd,
                        }
                    )

    new_aa_df = pd.DataFrame(columns=df_columns)
    for l in lines:
        new_aa_df.loc[len(new_aa_df)] = l
    new_aa_df = new_aa_df.reindex(
        columns=["Original Codon", "Nucleotide Index", "Change", "New Codon", "New AA"]
    )

    # Only keep original codon which matches input
    new_aa_df = new_aa_df[new_aa_df["Original Codon"] == seq]

    # Add back amino acid
    new_aa_df["Original AA"] = aa

    # Prepare pridict intput (sequence before edit / sequence after edit)
    # new_aa_df["pridinput"] = "(" + new_aa_df["Original Codon"] + "/" + new_aa_df["New Codon"] + ")"

    # update so single nucleotide change:
    new_aa_df["pridinput"] = new_aa_df.apply(modify_string, axis=1)

    # reorder columns
    new_aa_df = pd.concat(
        [new_aa_df["Original AA"], new_aa_df.drop("Original AA", axis=1)], axis=1
    )
    new_aa_df.reset_index(inplace=True, drop=True)

    return new_aa_df

###########################################
####### Preparing pridict input ###########
###########################################

def make_pridict_input(x):
    codon = x["codon"]
    sequence = x["sequence"]
    pos_label = str(x["position_label"])
    gene = x["gene"]
    
    # create a table of all possible single nt changes giving new aa
    subs_df = single_change_from_seq(codon)
    
    # make a table of input sequences
    input_df = subs_df[["pridinput", "Change", "Original AA", "New AA"]]

    # get single letter change:
    input_df = input_df.join(
        aa_df.rename(columns={"single_letter": "aa_1"}), on="Original AA"
    )
    input_df = input_df.join(
        aa_df.rename(columns={"single_letter": "aa_2"}), on="New AA"
    )

    # add sequence info
    input_df["sequence"] = sequence
    input_df["input_seq"] = input_df.apply(replace_bycolumn, axis=1)

    # wrangle to match pridict batch mode syntax
    # need to write a csv file that has two columns [editseq, sequence_name]
    # pridict_df = input_df[["input_seq", "Change"]]
    # pridict_df.columns = ["editseq", "sequence_name"]
    # pridict_df["sequence_name"] = pridict_df["sequence_name"].str.replace(">", pos_label)

    pridict_df = input_df[["input_seq", "aa_1", "aa_2"]].copy()
    pridict_df["sequence_name"] = pridict_df["aa_1"] + pos_label + pridict_df["aa_2"]
    pridict_df = pridict_df[["input_seq", "sequence_name"]]
    pridict_df.columns = ["editseq", "sequence_name"]

    # make sure sequence names are unique:
    # check for duplicate values
    duplicates = pridict_df["sequence_name"].duplicated(keep=False)

    # add suffix to duplicate values to make them unique
    suffix = duplicates.groupby(duplicates).cumsum().astype(str)

    pridict_df["sequence_name"] = np.where(
        duplicates,
        pridict_df["sequence_name"] + "mut" + suffix,
        pridict_df["sequence_name"],
    )
    
    # add gene to label
    pridict_df["sequence_name"] = gene + "_" + pridict_df["sequence_name"]

    return pridict_df


###########################################
### VEP interface for BEstiamte output ####
###########################################

def create_hvgs(seq, chromosome, position):
    """
    create HVGS notation of all possible single nt substitutions at position of interest
    format: chromosome number:.g position input nt > new nt

    """

    hvgs_list = []

    N = seq[99]

    for m in ["A", "C", "T", "G"]:
        if m != N:
            hvgs = str(chromosome) + ":g." + str(position) + N + ">" + m
            hvgs_list.append(hvgs)

    return hvgs_list


def get_vep(hvgs_in):
    """
    Using VEP API, extract information on codon, aa change and aa position
    """

    server = "https://rest.ensembl.org"
    ext = "/vep/human/hgvs/" + hvgs_in + "?"

    r = requests.get(server + ext, headers={"Content-Type": "application/json"})

    if not r.ok:
        r.raise_for_status()
        sys.exit()

    decoded = r.json()

    if "amino_acids" in decoded[0]["transcript_consequences"][0].keys():
        aa_change = decoded[0]["transcript_consequences"][0]["amino_acids"]
        aa_pos = decoded[0]["transcript_consequences"][0]["protein_start"]
        codon_change = decoded[0]["transcript_consequences"][0]["codons"]

        if decoded[0]["most_severe_consequence"] == "synonymous_variant":
            aa_change = aa_change + "/" + aa_change

    elif "3_prime_UTR" in decoded[0]["most_severe_consequence"]:
        aa_change = "3'UTR"
        aa_pos = re.findall(r'\d+', hvgs_in[hvgs_in.find("g.")+2 : len(hvgs_in)])[0] # get genomic position from HVGS in to discriminate between different 3'UTR changes
        codon_change = ""
        
    elif "splice" in decoded[0]["most_severe_consequence"]:
        aa_change = "splicedonor"
        aa_pos = re.findall(r'\d+', hvgs_in[hvgs_in.find("g.")+2 : len(hvgs_in)])[0] # get genomic position from HVGS in to discriminate between different 3'UTR changes
        codon_change = ""
        
    elif "intron" in decoded[0]["most_severe_consequence"]:
        aa_change = "intron"
        aa_pos = re.findall(r'\d+', hvgs_in[hvgs_in.find("g.")+2 : len(hvgs_in)])[0] # get genomic position from HVGS in to discriminate between different 3'UTR changes
        codon_change = ""
        
    else:
        aa_change = decoded[0]["most_severe_consequence"]
        aa_pos = re.findall(r'\d+', hvgs_in[hvgs_in.find("g.")+2 : len(hvgs_in)])[0]
        codon_change = '?'
        

    vep_out = aa_change + "_" + codon_change + "_" + str(aa_pos)

    return vep_out
    
