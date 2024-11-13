#This script finds CNV from reference paper, break them into intervals and match it with current files
#BIOBIX UGENT
import glob

def match_with_cnv(input_file_name):
    output_file_name = input_file_name.replace("_geno_LAF.txt", "_geno_LAF_cnv.txt")

    with open(input_file_name, "r") as f1, open("Reference_cnv.txt", "r") as f2, open(output_file_name, "w") as fout:
        # Read and write the headers
        header1 = f1.readline().strip().split('\t')

        # Write the header from the genotype file and add the 'CNV_Status' column
        fout.write('\t'.join(header1) + '\tCNV_Status\n')

        
        next(f2)

        # Store CNV intervals with statuses in a dictionary
        cnv_intervals = {}
        for line in f2:
            cols = line.strip().split('\t')
            contig, interval, cnv, status = cols[0], cols[1], cols[2], cols[4] 
            if interval != 'Position':  
                start, stop = map(int, interval.split('-'))
                positions = range(start, stop + 1)
                cnv_intervals.setdefault(contig, {}).update({pos: status for pos in positions})

        # Loop through the genotype data, adding CNV status while retaining the original structure
        for line in f1:
            cols = line.strip().split('\t')
            contig, position = cols[0], int(cols[1])

            # Look for CNV status in the dictionary
            if contig in cnv_intervals and position in cnv_intervals[contig]:
                status = cnv_intervals[contig][position]
                # Check if status is 'gain' or 'loss', and add it to the line
                if status in ['gain', 'loss']:
                    fout.write('\t'.join(cols) + f'\t{status}\n')
                else:
                    # Add a dash for CNV status if status is not 'gain' or 'loss'
                    fout.write('\t'.join(cols) + '\t-\n')
            else:
                # Add a dash for CNV status if position not found
                fout.write('\t'.join(cols) + '\t-\n')

# Find all files matching the pattern "*_geno_LAF.txt"
input_files = glob.glob("*_geno_LAF.txt")

for input_file_name in input_files:
    match_with_cnv(input_file_name)
