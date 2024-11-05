#!/bin/bash

# Process start/stop coordinates to only use the "5'-region" of the gene
awk -F $'\t' '
!/^#/ { 
    # Extract ID from the first field in $9 of the gff
    split($9, info, ";");
    sub(/^ID=/, "", info[1]);  # Remove "ID=" prefix
    ID = info[1];              # Set ID value

    #Calculate coordinates based on strand orientation
    if ($7 == "-") {
        start = ($5 - 10 >= 0) ? $5 - 10 : 0;
        end = $5 + 27;
    } else {
        start = ($4 - 27 >= 0) ? $4 - 27 : 0;
        end = $4 + 10;
    }
    
    # output in bed format, keeping sequence sample name, gene name and strand
    print $1 "\t" start "\t" end "\t" ID "\t" $7
}' "$1" 1> "$2"
