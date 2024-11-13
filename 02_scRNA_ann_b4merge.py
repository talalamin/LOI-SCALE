#This  scripts is for second step  for post-processing and cratin of appropraite file system.
# It is used to merge annovar annotation  with pileup parsed data.
#It create two new columns at the end of file for Gene and known rsid inforamtion from vcf file.
#BiOBIX UGENT


# Read gene annotation file and create a dictionary of chr_loc to gene mapping with values
def create_gene_dict(gene_file):
    gene_dict = {}
    with open(gene_file, 'r') as gene_fh:
        for line in gene_fh:
            chr, loc, gene = line.strip().split('\t')
            gene_dict[(chr, loc)] = gene
    return gene_dict

# Rsid file and create a dictionary of chr_loc to rsid mapping with values
def create_rsid_dict(rsid_file):
    rsid_dict = {}
    with open(rsid_file, 'r') as rsid_fh:
        for line in rsid_fh:
            chr, loc, _, _, rsid = line.strip().split('\t')
            rsid_dict[(chr, loc)] = rsid
    return rsid_dict

# Update genotype file with gene names and rsids
def update_genotype_file(genotype_file, gene_dict, rsid_dict, output_file):
    with open(genotype_file, 'r') as genotype_fh, open(output_file, 'w') as output_fh:
        for idx, line in enumerate(genotype_fh):
            if idx == 0: 
                output_fh.write(line.strip() + '\tGene\tRsid\n') 
            else:
                data = line.strip().split('\t')
                chr, loc = data[0], data[1]
                gene = gene_dict.get((chr, loc), '_')  # Get gene name, place a dasf '_' if gene name not found
                rsid = rsid_dict.get((chr, loc), '.')  # Get rsid, place a dot '.' if rsid is not found
                output_fh.write('\t'.join(data + [gene, rsid]) + '\n')  # Write updated line to output file

# Read prefixes from Sample_list.txt
prefixes_file = 'sample_list.txt'
with open(prefixes_file, 'r') as prefixes_fh:
    prefixes = prefixes_fh.read().splitlines()

# Process each prefix
for prefix in prefixes:
    genotype_file = f'{prefix}.parse25.txt'
    gene_file = f'{prefix}.gene.txt'
    rsid_file = f'{prefix}.rsid.txt'
    output_file = f'{prefix}_genotype_ann.txt'

    # Create dictionaries for gene names and rsids with values
    gene_dict = create_gene_dict(gene_file)
    rsid_dict = create_rsid_dict(rsid_file)

    # Update genotype file with gene names and rsids
    update_genotype_file(genotype_file, gene_dict, rsid_dict, output_file)
