# Pridict pipeline

This repo provides scripts to run the pridict algorithm [Mathis et al, 2023](https://www.nature.com/articles/s41587-022-01613-7) automatically for mulitple 
inputs. See documentation [here](./PRIDICT/README.md) for explanations on the PRIDICT model itself. 

## Setting up the pipeline

In order to run the pipeline, create a master folder in which you will keep all the scripts / analysis from this repository. 

Within this folder, create a conda environment with pridict based on
[github](https://github.com/uzh-dqbm-cmi/PRIDICT) instructions :

```
# clone PRIDICT repository
git clone https://github.com/mariemoullet/PRIDICT.git

# navigate into the PRIDICT specific folder within the repository
cd PRIDICT/PRIDICT

# create conda environment and install dependencies for PRIDICT (only has to be done before first run/install)

# use pridict_linux for linux machine or pridict_mac for a macbook
conda env create -f pridict_linux.yml # for linux machine, pridict_mac.yml for mac

# activate the created environment
conda activate pridict

	### ONLY FOR M1 Mac you need to additionally run the following conda install command (tensorflow): 
	conda install conda-forge::tensorflow
	###


# Now install packages to run in batch mode (withou manually providing sequences)
pip install -r requirements.txt

```

The batch generation also depends on R packages. To set up the R packages, install the following packages from 
bioconductor:

```
if (!require("BiocManager", quietly = TRUE))
   install.packages("BiocManager")

BiocManager::install("ensembldb")
BiocManager::install("EnsDb.Hsapiens.v86")
BiocManager::install("org.Hs.eg.db")

install.packages(tidyverse)
install.packages(spgs)
install.packages(janitor)
install.packages(digest)

```

or simply run 

```
Rscript packages.R
```


There are 2 ways to run this pipeline :

## 1. All possible amino acid changes resulting from a single nucleotide change in a specific codon in multiple proteins 

If you have multiple codons/proteins of interest, you can automatically generate pegRNAs for every potential single nucleotide change at your codon of interest which would result in an amino acid change. If there are synonymous mutations, these will be labelled as mutx eg C797Smut1, C797Smut2.

To run this version, create an input file called `input_gene_aa.csv` in the pipeline directory where columns are 
* gene : entrez symbol eg EGFR
* id :ENSEMBL MANE protein id eg ENSP00000275493
* aa : amino acid position eg 797


## 2. All possible single nucleotide substitutions at a specific genomic position 

This is part is designed to use prime editing to engineer mutations at sites targeted by base editors. The pipeline expects an input format with at least the following columns:
* Gene : gene symbol eg EGFR
* chr : chromsome number eg 7
* Edit_Location : genomic locations of possible edits, separated by ;



## Executing the pipeline 

To run the pipeline, in command line run 

```
bash run_batch_mode.sh <intput_filename> [--filter=<filter>] [--output=<output_filename.csv>]
```

The input arguments are
* input file : either input_from_aa.csv or input_from_be.csv as described above
* filter : optional boolean argument : should the file pegRNAs (include Gibson Assembly homology arms) be <= 200 bp? (default is true)
* output file name: optional, default is oligos.csv (nb include extension in argument)

The output returns a file (by default called oligos.csv) which includes:
* gene and mutation for the edit encoded by the pegRNA
* pridict_rank : rank from the pridict algorithm 
* pegRNA : the pegRNA sequence
* for Gibson assembly : pegRNA with added Gibson assembly arms (format is TATCTTGTGGAAAGGACGAAA + pegRNA + GCGCGGTTCTATCTAGTTACGCGT)
* length / length with GA overhang : pegRNA oligo length with or without GA overhang

