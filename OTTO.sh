#!/bin/bash

# Initialize variables for options
ASOfile=""
cores=1
mismatches=2
output="OTTOput"
kmer_db="dedup.db"
idx_tree="out"

################################################################## specify help message! #################################################################
# Function to display help message
usage() {
    echo "Usage: $0 -a <ASOfile> -c <cores> -m <mismatches> -o <output> -k <kmer_db> -i <idx_tree>"
    exit 1
}

# Parse command line arguments
while getopts ":a:c:m:o:k:i:" opt; do
  case $opt in
    a)
      ASOfile=$OPTARG
      ;;
    c)
      cores=$OPTARG
      ;;
    m)
      mismatches=$OPTARG
      ;;
    o)
      output=$OPTARG
      ;;
    k)
      kmer_db=$OPTARG
      ;;
    i)
      idx_tree=$OPTARG
      ;;
    \?)
      echo "Invalid option: -$OPTARG" >&2
      usage
      ;;
    :)
      echo "Option -$OPTARG requires an argument." >&2
      usage
      ;;
  esac
done

# Check if ASOfile is provided
if [ -z "$ASOfile" ]; then
  echo "ASOfile is required."
  usage
fi

# check if ASO-file exists
if [ ! -f "$ASOfile" ]; then
  echo "$ASOfile does not exist!"
  usage
fi

# check if kmer-DB exists
if [ ! -d "$kmer_db" ]; then
  echo "$kmer_db does not exist!"
  usage
fi

# check if the centroid index hash-tree directory is there
if [ ! -d "$idx_tree" ]; then
  echo "$idx_tree does not exist!"
  usage
fi

# check if output directory is already present
if [ ! -d "$output" ]; then
  mkdir $output
else
  echo "$output already exists! Delete it and retry."
  exit 1
fi

# permute ASO sequences to generate mismatching sequences
# and if necessary extend shorter ASO sequences to 12 nucleotides
python scripts/ASO_permutation.py -a $ASOfile \
                                  -m $mismatches \
                                  -o $output

# query the kmer-DB for matching sequences and return
# hit centroids
for kmerlength in $output/*permutation_fasta.txt; do
  nopath=${kmerlength#*/}
  k=${nopath%permutation_fasta.txt};
  if [ ! -f "${kmer_db}/${k}.db" ]; then
    echo "there is no database for k-mer length ${k} in directory ${kmer_db}!"
    echo "use the OTTO snakemake-pipeline to create the neccessary database."
    exit 1
  fi;
done

for kmerlength in $output/*permutation_fasta.txt; do
  nopath=${kmerlength#*/}
  k=${nopath%permutation_fasta.txt};
  scripts/kmer-db_stdout_noInfo/kmer-db new2all -multisample-fasta \
                                              -sparse \
                                              -t 1 "${kmer_db}/${k}.db" ${kmerlength} \
                                              /dev/null ;
done | \

# match the hit centroids to aligned 5'-regions and return
# sample, gene and COG
python scripts/expand_matches.py $idx_tree \
                                    $cores \
                                    $output > ${output}/results_expanded.txt
