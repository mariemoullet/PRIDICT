# PRIDICT: PRIme editing guide RNA preDICTion 

![PRIDICT logo](pridict_logo.png)

For accessing Supplementary Files, click [here](https://github.com/uzh-dqbm-cmi/PRIDICT/tree/supplementary_files).

Repository containing `python` package for running trained `PRIDICT` (PRIme editing guide RNA preDICTion) models. `prieml` package includes modules to setup and run `PRIDICT` models for predicting `prime editing efficiency and product purity`.

To run `PRIDICT` online, see our [webapp](https://pridict.it/).

--------------------------

### Installation using Anaconda (Linux and Mac OS) 🐍
!!! `PRIDICT` can only be installed on `Linux` and `Mac OS` since `ViennaRNA` package is not available for `Windows` !!!

The easiest way to install and manage Python packages on various OS platforms is through [Anaconda](https://docs.anaconda.com/anaconda/install/). Once installed, any package (even if not available on Anaconda channel) could be installed using pip. 

* Install [Anaconda](https://docs.anaconda.com/anaconda/install/).
* Start a terminal and run:
    ```shell
    # clone PRIDICT repository
    git clone https://github.com/uzh-dqbm-cmi/PRIDICT.git
    # navigate into repository
    cd PRIDICT
    # create conda environment and install dependencies for PRIDICT (only has to be done before first run/install)
    conda env create -f pridict.yml
    # activate the created environment
    conda activate pridict
    # run desired PRIDICT command (manual or batch mode, described below)
    # for running batch mode, first create batch file and copy it into 'input' folder
    # batch results will be stored in 'predictions' folder
    python pridict_pegRNA_design.py batch --input-fname batch_example_file.csv --output-fname batchseqs
    ```

* `PRIDICT` environment only has to be installed once. When already installed, follow the following commands to use `PRIDICT` again:
    ```shell
    # open Terminal/Command Line
    # navigate into repository
    cd PRIDICT
    # activate the created environment
    conda activate pridict
    # run desired PRIDICT command (manual or batch mode, described below)
    # for running batch mode, first create batch file and copy it into 'input' folder
    # batch results will be stored in 'predictions' folder
    python pridict_pegRNA_design.py batch --input-fname batch_example_file.csv --output-fname batchseqs
    ```

--------------------------

### Running PRIDICT in 'manual' mode:
  ####  Required:
  -  `--sequence-name`: name of the sequene (i.e. unique id for the sequence)
  -  `--sequence`: target sequence to edit in quotes (format: `"xxxxxxxxx(a/g)xxxxxxxxxx"`; minimum of 100 bases up and downstream of brackets are needed)
  ####  Optional:
  -  `--output-dir`: output directory where results are dumped on disk (default: "./predictions")
  -  `--use-5folds`: Use all 5-folds trained models. Default is to use fold-1 model
  -  `--cores`: Number of cores to use for multiprocessing. Default value 0 uses all available cores.
  -  `--nicking`: Additionally, design nicking guides for edit (PE3) with DeepSpCas9 prediction.
  -  `--ngsprimer`: Additionally, design NGS primers for edit based on Primer3 design.
```shell

python pridict_pegRNA_design.py manual --sequence-name seq1 --sequence 'GCCTGGAGGTGTCTGGGTCCCTCCCCCACCCGACTACTTCACTCTCTGTCCTCTCTGCCCAGGAGCCCAGGATGTGCGAGTTCAAGTGGCTACGGCCGA(G/C)GTGCGAGGCCAGCTCGGGGGCACCGTGGAGCTGCCGTGCCACCTGCTGCCACCTGTTCCTGGACTGTACATCTCCCTGGTGACCTGGCAGCGCCCAGATGCACCTGCGAACCACCAGAATGTGGCCGC'
``` 
--------------------------

### Running in batch mode:
  ####  Required:
  -  `--input-fname`: input file name - name of csv file that has two columns [`editseq`, `sequence_name`]. See `batch_template.csv` in the `./input` folder
  ####  Optional:
  -  `--input-dir` : directory where the input csv file is found on disk
  -  `--output-dir`: directory on disk where to dump results (default: "./predictions")
  -  `--output-fname`: output filename used for the saved results
  -  `--combine-results`: Compile all results in one dataframe
  -  `--use-5folds`: Use all 5-folds trained models. Default is to use fold-1 model
  -  `--cores`: Number of cores to use for multiprocessing. Default value 0 uses all available cores.
  -  `--nicking`: Additionally, design nicking guides for edit (PE3) with DeepSpCas9 prediction.
  -  `--ngsprimer`: Additionally, design NGS primers for edit based on Primer3 design.
```shell

 python pridict_pegRNA_design.py batch --input-fname batch_example_file.csv --output-fname batchseqs

``` 
--------------------------

### Citation

If you find our work is useful in your research, please cite the following paper:

> @article {Mathis et al.,  
	author = {Mathis, Nicolas and Allam, Ahmed and Kissling, Lucas and  Marquart, Kim Fabiano and Schmidheini, Lukas and Solari, Cristina and Balázs, Zsolt and Krauthammer, Michael and Schwank, Gerald},  
	title = {Predicting prime editing efficiency and product purity by deep learning},  
	year = {2023},  
	doi = {},  
	URL = {},  
	eprint = {},  
	journal = {}  
}
