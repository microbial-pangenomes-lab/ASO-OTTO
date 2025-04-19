import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns
import numpy as np
import argparse

def get_options():
    description = "Visualize OTTO outputs"
    parser = argparse.ArgumentParser(description=description)

    parser.add_argument("-i", "--infile",
                        required=True,
                        help = "file containing the output from OTTO"
                        )

    parser.add_argument("-coi", "--COG_of_interest",
                        type=str,
                        default=None,
                        help = "COG to be displayed at all times"
                        )

    parser.add_argument("-t", "--taxfile",
                        required=True,
                        help = "file containing the taxonomic data for the dataset"
                        )

    parser.add_argument("-tr", "--tax_rank",
                        default="species",
                        help = """taxonomic rank to be displayed.
                        choose from superkingdom, phylum, class, order, family, genus, species."""
                        )

    parser.add_argument("-s", "--min_samples",
                        type=int,
                        default=1000,
                        help = """minimum number of samples needed
                        for a species to be displayed. (default: %(default)d)"""
                        )

    parser.add_argument("-c", "--top_COGs",
                        type=int,
                        default=10,
                        help = """number of most targeted COGs to display.
                        The COG of interest will be added regardless of the
                        number of its hits! (default: %(default)d)"""
                        )

    parser.add_argument("-on", "--out_name", default="sample", help = "name to be used as out-file prefix.")

    parser.add_argument("-fmt", "--out_format",
                                                default="svg",
                                                help = """format of out-file.
                                                (formats offered by plt.savefig() --> svg, pdf, png...)"""
                                                )
    # parser.add_argument("-v", action='count',
    #                     default=0,
    #                     help='Increase verbosity level')

    # parser.add_argument('--version', action='version',
    #                   version='%(prog)s '+__version__)

    return parser.parse_args()

def main():
    args = get_options()

    infile=args.infile
    taxfile=args.taxfile
    tax_rank=args.tax_rank
    min_samp=args.min_samples
    nr_COG=args.top_COGs
    coi=args.COG_of_interest
    ofname=args.out_name
    fmt=args.out_format
    #load OTTO-results and the taxonomy file
    print(coi)
    data=pd.read_csv(infile,sep=',')
    taxa=pd.read_csv(taxfile,sep='\t')
    #create DF for desired taxonomic rank
    tax=taxa[["sample_ID",tax_rank]]

    taxdat=pd.merge(data,tax,on="sample_ID",how="left")
    taxdat["sage"]=taxdat["sample_ID"] + ":" + taxdat["gene_ID"]

    taxcount=pd.DataFrame(tax[tax_rank].value_counts())
    taxcount=taxcount[taxcount["count"] >= min_samp]

    dedupdict={}
    COGdict={}
    total_hits_df=pd.DataFrame(index=list(taxcount.index))
    mima=data["mismatch"].nunique()
    #deduplicate the OTTO results. there will be duplicate hits if the provided ASOs are shorter than 12nt!
    for x in range(0,mima):
        dedupdict[str(x)]=taxdat[taxdat["mismatch"] == x].drop(["gene_ID","ASO"],axis=1).drop_duplicates()

        COGdict[str(x)]=set(pd.DataFrame(dedupdict[str(x)]["COG"].value_counts().head(nr_COG)).index)
        COGdict[str(x)].add(coi)
        # keep specified taxonomic rank with more than "min_samp" samples
        df_filtered = dedupdict[str(x)][dedupdict[str(x)][tax_rank].isin(list(taxcount.index))][[tax_rank,"sample_ID"]].drop_duplicates()
        # Count total hits per taxonomic rank
        hits_per_tax = df_filtered.groupby(tax_rank).size()
        # Calculate percentage of total hits per species for this mismatch
        total_hits_percentage = (hits_per_tax / taxcount.sort_index()["count"]) * 100
        # Store results in the total_hits_df
        total_hits_df[f"mismatch_{x}_total_hits"] = total_hits_percentage

    cog_dict={}
    for i, df in dedupdict.items():
        # Initialize cog_hits_df with species as index
        cog_dict[i] = pd.DataFrame(index=list(taxcount.index))
        # Initialize a set to keep track of processed COGs for each iteration
        processed_COGs = COGdict[str(i)]
        # Initialize a DataFrame to hold the sum of COG percentages per species
        perc_per_tax = pd.Series(0, index=cog_dict[i].index)
        # Loop through each specified COG and calculate percentages for each species

        for COG in processed_COGs:
            # Filter for the current COG and species
            df_filtered = df[(df[tax_rank].isin(list(taxcount.index))) & (df['COG'] == f'{COG}')]
            # Count the number of samples hit in this COG for each species
            cog_hits_per_tax = df_filtered.groupby(tax_rank)['sample_ID'].count()
            # Calculate the percentage of samples within a species hit in this COG
            cog_hits_percentage = (cog_hits_per_tax / taxcount.sort_index()["count"]) * 100
            # Store the percentage in the cog_hits_df
            cog_dict[i][f'mismatch_{i}_{COG}'] = cog_hits_percentage
            # Add to the total percentage sum for the current species
            perc_per_tax += cog_hits_percentage.fillna(0)

    heatmap_data = pd.concat([total_hits_df["mismatch_0_total_hits"], cog_dict[str(0)]], axis=1)
    if mima > 1:
        for x in range(1,mima):
            heatmap_data = pd.concat([heatmap_data,total_hits_df[f"mismatch_{x}_total_hits"], cog_dict[str(x)]],axis=1)
    #light grey=COG was not hit in this taxonomic rank, grey=less than 1 percent of samples were hit in this taxonomic rank
    cmap1 = plt.get_cmap('viridis').copy()
    cmap1.set_bad('xkcd:light grey')
    cmap1.set_under('xkcd:grey')

    plt.figure(figsize=(15, 10))
    ax = sns.heatmap(heatmap_data, cmap=cmap1, square=True, annot=False, fmt="g", linewidths=.5, yticklabels=True, xticklabels=True, vmin=0.01, vmax=100)
    plt.ylabel(tax_rank)
    ax.tick_params(top=True, bottom=False, labeltop=True, labelbottom=False)
    plt.xticks(rotation=90, fontsize=14) #, ha='left')
    plt.yticks(fontsize=14, fontstyle='italic')
    plt.tight_layout()
    # plt.colorbar(label="percentage")
    plt.savefig(f"{ofname}_hits.{fmt}", format=f"{fmt}", transparent=True)

if __name__ == "__main__":
    main()
