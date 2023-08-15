#!/bin/bash

# clean up existing pridict output files 
if ls PRIDICT/predictions/*.csv >/dev/null 2>&1; then
    rm PRIDICT/predictions/*.csv
fi

# clean up existing intermediate files
if ls intermediate/*.csv >/dev/null 2>&1; then
    rm intermediate/*.csv
fi

# activate conda env for pridict
source activate pridict

# Set default values for optional arguments
FILTER="TRUE"
OUTPUT="oligos.csv"

# Check if the correct number of arguments were provided
if [ $# -lt 1 ]; then
  echo "Missing filename"
  echo "Usage: $0 <filename> [--filter=<bool>] [--output=<output_filename.csv>]"
  exit 1
fi

# Parse optional arguments
for arg in "$@"; do
  case "$arg" in
    --filter=*)
      FILTER="${arg#*=}"
      shift
      ;;
    --output=*)
      OUTPUT="${arg#*=}"
      shift
      ;;
    *)
      ;;
  esac
done


# Check if the file exists
# if [ ! -f "$1" ]; then
#   echo "Error: File does not exist"
#   exit 1
# fi

# Prepare pridict input:
# Determine which Python script to run based on the file name
case "$1" in
  input_from_aa.csv)
    echo "Retrieving sequences..."
    Rscript prepare_pridict/from_aa_to_gx.R 
    python prepare_pridict/from_aa_to_pridict.py $1
    ;;
  input_from_be.csv)
    echo "Retrieving sequences..."
    python prepare_pridict/from_be_to_pridict.py $1
    ;;
  *)
    echo "Error: Invalid file name"
    exit 1
    ;;
esac

# Run pridict  
echo "Running pridict..."
python PRIDICT/pridict_pegRNA_design.py batch --input-fname ../../intermediate/pridict_input.csv --output-fname batchseqs


# Convert pridict output to oligos for order
# run the R script with the filter argument
Rscript pridict_to_oligos/make_oligos.R $FILTER $OUTPUT


# Exit successfully
echo "Done!"
exit 0
