#This script add cellular annotation to genotype file on basis of RGID
#BIOBIX UGENT
import pandas as pd
import glob

def write_header(output_file, columns):
    output_file.write('\t'.join(columns) + '\n')

chunk_size = 15000

# Find all files matching the pattern "*_geno_LAF_cnv.txt"
input_files = glob.glob("*_geno_LAF_cnv.txt")
cell_types_bc = 'cell_types_bc.tsv'

for input_file in input_files:
    # Extract the base file name without the extension and split it by '_'
    base_name = input_file.split('_')[0]

    # Create the ouptut file name by adding '_LAF_cnv_cell.txt' to the base file name
    output_file_name = f"{base_name}_LAF_cnv_cell.txt"
    output_header_written = False

    file1_iter = pd.read_csv(input_file, sep='\t', iterator=True, chunksize=chunk_size)
    file2_iter = pd.read_csv(cell_types_bc, sep='\t', iterator=True, chunksize=chunk_size)

    with open(output_file_name, 'w') as output_file:
        for file1_chunk in file1_iter:
            for file2_chunk in file2_iter:
                merged_chunk = pd.merge(file1_chunk, file2_chunk, on='RGID', how='left')

                if not output_header_written:
                    write_header(output_file, merged_chunk.columns)
                    output_header_written = True

                merged_chunk.to_csv(output_file, header=False, index=False, sep='\t', mode='a')
                break

            file2_iter = pd.read_csv(cell_types_bc, sep='\t', iterator=True, chunksize=chunk_size)

print("Processing completed!")
