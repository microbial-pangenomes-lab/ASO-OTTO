# ASO-OTTO
ASO-OTTO is a k-mer based tool to that finds targets and off-targets of antisense oligonucleotides (ASOs) in large genome collections. It comes with a Snakemake pipeline to create searchable databases from your own fasta files and relies on the fast k-mer lookup of Kmer-db to handle the search requests.


# Installation

Currently, ASO-OTTO needs to be cloned from github and the submodule (a forked and modified version of) Kmer-db needs to be initialized.

clone repository:
```
git clone --recursive https://github.com/microbial-pangenomes-lab/ASO-OTTO.git
```

Follow the steps listed in the forked [Kmer-db](https://github.com/haneubau/kmer-db_stdout_noInfo) repository to compile.


To use the ASO-OTTO pipeline, you will need to install Snakemake as well. It can be [installed](https://snakemake.readthedocs.io/en/stable/getting_started/installation.html) by following the instructions on their website.

Basic instructions for installing Snakemake using conda:
```
conda create -c conda-forge -c bioconda -n snakemake snakemake 
```  
This installs just the basic snakemake environment. If you want to run the pipeline with this environment you will always have to specify `--use-conda` in the snakemake command!

If you also want to use the visualization script of ASO-OTTO, we recommend to create another conda environment with:
```
conda create -c conda-forge -c bioconda -n OTTO numpy=1.26.2 pandas=2.1.4 seaborn=0.13.2 matplotlib=3.8.4
```

After that you should be good to go!

# Quick-start

## running the pipeline
When creating your own searchable dataset you should first have a look in the config file of the pipeline with a text editor of your choice:

```
cd ASO-OTTO/OTTO_snakemake
nano config/config.yaml
```
Here you will have to change some variables, depending on what you want to do and the specifics of you computer.
Most importantly you will have to tell the pipeline whether you want to download the Diamond database for eggNOG and the Kraken2 DB.
That can be done by changing the variables `eggNOG_DB: True`/`Kraken2_DB: True` to `eggNOG_DB: False`/`Kraken2_DB: False`.
If you already have those databases on your system you can skip this step!
The download paths can be specified with the variables`pteggDB` and `ptkrakDB`. If you already have the databases, you still need to set the path to the databases!

If you work on a system that can handle ~50GB of RAM being used, you can keep the variable on `fastmode: True` and eggNOG will load the entire Diamond DB into memory when annotating protein sequences. Otherwise you need to switch it to `fastmode: False`. 


After entering the path to your fasta files with `input_file: "data/fasta/"` and specifying the file suffix of your fasta files with `file_suffix: ".fasta"`, you can start the pipeline with:

```
conda activate snakemake
snakemake --cores X --use-conda
```

This will create the databases in `ASO-OTTO/OTTO_snakemake/output` as `kX.db`, depending on how many sizes of k you specified in the config file.
The hash tree-like folder is created as `ASO-OTTO/OTTO_snakemake/output/IDtree` and is needed as input for ASO-OTTO.

With that you are set to run ASO-OTTO on your own dataset!

## running ASO-OTTO
after preparing everything and converting your genomes into a searchable database you can run ASO-OTTO with:
```
cd ..
bash OTTO.sh -a MyASO.fasta -c 1 -m 2 -o MyOutputFolder -k OTTO_snakemake/output -i OTTO_snakemake/output/IDtree
```

This will create the `results_expanded.txt`, which contains the actual output, `idx_ASO_target.txt` containing your ASOs and their ID and `kXpermutation.fasta` + `kXpermutation_fasta.txt` depending on what lengths your queried ASOs are.

If your `MyASO.fasta` only contained one ASO, you can directly use the convenience script to visualize it, otherwise you'll need to split the `results_expanded.txt` by the `orig_ID` column.
```
conda activate OTTO
python scripts/OTTO_viz_new.py -i MyOutputFolder/results_expanded.txt -coi MyCOGOfInterest -t OTTO_snakemake/output/full_taxonomy.tsv 
```

That should give you a broad overview of what species/genes were targeted by your ASO!

# Prerequisites
ASO-OTTO uses:
```
Kmer-db
```

ASO-OTTO pipeline uses (handled by snakemake):
```
Biopython
numpy
pandas
eggNOG
prodigal
Kraken2
Taxonomy-ranks
bedtools
```

OTTO_viz_new.py also requires (have to be installed manually):
```
pandas
numpy
matplotlib 
seaborn
```
