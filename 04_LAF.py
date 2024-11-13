#This script calcultes Lesser allele fraction
#BIOBIX UGENT

import glob

def calculate_laf(ref_count, alt_count):
    # Avoid division by zero by checking total count before calculating LAF
    total_count = ref_count + alt_count
    if total_count == 0:
        return None  # Return None if total count is zero

    # Calculate the LAF as the ratio of the smaller count to the total count
    return min(ref_count, alt_count) / total_count

# Find all files matching the pattern "*_output.txt"
input_files = glob.glob("*_output.txt")

for input_file_name in input_files:
    output_file_name = input_file_name.replace("_output.txt", "_geno_LAF.txt")

    with open(input_file_name, "r") as input_file, open(output_file_name, "w") as output_file:
        # Read the header line
        header = input_file.readline().strip()
        output_file.write(header + "\tLAF\n")  # Write the header line for the output file, including the LAF column

        # Process each subsequent line
        for line in input_file:
            # Split the line by tabs
            fields = line.strip().split("\t")

            # Skip lines that might contain header information or are empty
            if len(fields) < 12:
                continue

            # Extract the necessary fields
            chr_val = fields[0]
            loc_val = fields[1]
            ref_val = fields[2]
            alt_val = fields[3]

            # Handling the conversion to int with error checking
            try:
                ref_count = int(fields[4])
                alt_count = int(fields[5])
            except ValueError:
                # Handle the case where conversion fails
                continue  # Skip this line if conversion fails

            genotype = fields[6]
            RGID = fields[7]
            second_alt = fields[8]
            second_alt_count = fields[9]
            gene = fields[10]
            rsid = fields[11]

            # Calculate LAF
            laf = calculate_laf(ref_count, alt_count)

            if laf is None:
                continue 

            # Write the output line with the LAF value appended
            output_file.write(f"{chr_val}\t{loc_val}\t{ref_val}\t{alt_val}\t{ref_count}\t{alt_count}\t{genotype}\t{RGID}\t{second_alt}\t{second_alt_count}\t{gene}\t{rsid}\t{laf}\n")
