#This script is for 3rd step of the pipeline. where mpilup base files are properly parsed.
#BIOBIX UGENT
import csv
import logging
from collections import defaultdict

def determine_genotype(ref_base, alt_base):
    if not alt_base:
        return f"Homozygous (0/{ref_base})"
    else:
        return f"Heterozygous ({ref_base}/{alt_base})"

def process_input(rgid_prefix):
    input_filename = f"{rgid_prefix}_genotype_ann.txt"
    output_filename = f"{rgid_prefix}_output.txt"

    with open(input_filename, "r") as input_file, open(output_filename, "w") as output_file:
        reader = csv.DictReader(input_file, delimiter="\t")
        header = reader.fieldnames

        output_file.write("chr\tloc\tref\talt\trefCount\taltCount\tgenotype\tRGID\tsecondAlt\tsecondAltCount\tGene\tRsid\n")

        for line in reader:
            counts = defaultdict(int)
            for base in ['A', 'T', 'C', 'G']:
                counts[base] = int(line[base])

            total_count = sum(counts.values())

            if total_count == 0:
                logging.info(f"Skipping position {line['chr']}:{line['loc']} as it has no counts")
                continue

            ref_base = line['ref']
            alt_base = ""
            alt_count = 0

            for base in counts:
                if counts[base] > alt_count and base != ref_base:
                    alt_base = base
                    alt_count = counts[base]

            if ref_base not in counts or counts[ref_base] == 0:
                ref_base = "0"
                ref_count = 0
                genotype = determine_genotype(alt_base, "")
            else:
                ref_count = counts[ref_base]
                genotype = determine_genotype(ref_base, alt_base)

            if not alt_base:
                alt_count = 0
                alt_base = "0"

            second_alt_base = ""
            second_alt_count = 0

            for base in counts:
                if counts[base] > second_alt_count and base != ref_base and base != alt_base:
                    second_alt_base = base
                    second_alt_count = counts[base]

            if not second_alt_base:
                second_alt_count = 0
                second_alt_base = "0"

            output_line = [line['chr'], line['loc'], ref_base, alt_base, ref_count, alt_count, genotype, rgid_prefix, second_alt_base, second_alt_count, line['Gene'], line['Rsid']]
            output_file.write("\t".join(map(str, output_line)) + "\n")
            logging.info(f"Writing output line: {output_line}")

            print(f"Processed: {input_filename}")

def main():
    sample_list_filename = 'sample_list.txt'

    logging.basicConfig(level=logging.INFO, format="%(message)s")

    with open(sample_list_filename, 'r') as sample_file:
        sample_prefixes = [line.strip() for line in sample_file if line.strip()]

    for rgid_prefix in sample_prefixes:
        process_input(rgid_prefix)

if __name__ == "__main__":
    main()
