#this script look for imrinted genes from reference's paper and extract data from scRNA
#BIOBIX UGENT

import pandas as pd
import glob

# Read the list of imprinted genes (Tine bulk dataset) 
with open('bc_imprinted_genes.txt', 'r') as file:
    imprinted_genes = [line.strip() for line in file]

# Create an empty DataFrame to store concatenated data
concatenated_data = pd.DataFrame()

# Find all files matching the pattern "*_LAF_cnv_cell.txt"
input_files = glob.glob("*_LAF_cnv_cell.txt")

# Loop through each input file and filter the data
for file_name in input_files:
    # Read the genotype data from the file
    genotype_data = pd.read_csv(file_name, sep='\t')
    
    # Filter the genotype data to include only the rows where the gene is in the list of imprinted genes
    imprinted_genotype_data = genotype_data[genotype_data['Gene'].isin(imprinted_genes)]
    
    # Concatenate the filtered data to the existing DataFrame
    concatenated_data = pd.concat([concatenated_data, imprinted_genotype_data])

# Write the concatenated imprinted genotype data to a single output file
concatenated_data.to_csv('imprinted_genes.tsv', sep='\t', index=False)
