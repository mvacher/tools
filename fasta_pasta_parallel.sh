#!/bin/bash

FASTA_DIRECTORY="/dd_groupdata/bernard/fasta/"
OUTPUT_DIRECTORY="/dd_groupdata/bernard/outputs/"


process() { 
  echo "Processing file=$1 "
  output_file=`basename $1`	
  python fasta_pasta.py -f "$1" -i output.hits -o "$OUTPUT_DIRECTORY""$output_file"
}

export -f process
export OUTPUT_DIRECTORY

# Get the list of fasta file
list=$(find $FASTA_DIRECTORY*.fa)
./parallel process ::: "${list[@]}"

#Post processing: 
