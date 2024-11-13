#!/bin/bash

#Variant calling for full-length scRNA data (QC passed fastq files)
#TALAL AMIN: BIOBIX (UGENT)

##############################################################################################
##############################################################################################
########################### STEP 1: DEFINE VARIABLES and PARAMETERS ##########################
##############################################################################################
##############################################################################################

# 1: Define variables AND PARAMETERS
genomeDir=~/tools/ASE_ScRNA_Breast/GRCh38/
reference=~/tools/ASE_ScRNA_Breast/GRCh38/GRCh38.fa
knownSites=~/tools/ASE_ScRNA_Breast/data_breastcancer_fastq/dbsnp/Homo_sapiens_assembly38.dbsnp138.vcf
knownSitesCompressed=~/tools/ASE_ScRNA_Breast/data_breastcancer_fastq/dbsnp/dpsnp138.vcf.gz
pile2base=~/tools/ASE_ScRNA_Breast/data_breastcancer_fastq/scripts/pileup2base_no_strand.pl
samples=( $(cat sample_list.txt) )
gatk_options="--java-options \"-XX:GCTimeLimit=50 -XX:GCHeapFreeLimit=10 -XX:+PrintFlagsFinal -XX:ConcGCThreads=1 -XX:ParallelGCThreads=1 -XX:CICompilerCount=2  -XX:+PrintGCDetails -Xloggc:gc_log.log -Xmx15G -Djava.io.tmpdir=/data/tmp\""

##############################################################################################
##############################################################################################
########################### STEP 2: PREPARE REFERENCE GENOME #################################
##############################################################################################
##############################################################################################

# 2: Generate the genome using STAR
STAR --runMode genomeGenerate --genomeDir $genomeDir --genomeFastaFiles $reference --sjdbGTFfile ~/tools/ASE_ScRNA_Breast/GRCh38/gencode.v44.primary_assembly.annotation.gtf --sjdbOverhang 99 --runThreadN 12

##############################################################################################
##############################################################################################
########################### STEP 3: PRE-PROCESSING & CONVERSION ##############################
##############################################################################################
##############################################################################################

# 3.1: Align reads to the genome and output genecounts in parallel
parallel -j 4 "STAR --genomeDir \"$genomeDir\" --readFilesIn {}_1.fastq {}_2.fastq --outFilterIntronMotifs RemoveNoncanonicalUnannotated --outSAMmapqUnique 60 --outFileNamePrefix {} --runThreadN 8" ::: "${samples[@]}"

# 3.2: Convert .sam to .bam in parallel
parallel -j 5 "samtools view -bS {}Aligned.out.sam > {}Aligned.out.bam" ::: "${samples[@]}"

# 3.3: Filter .bam files in parallel
parallel -j 5 "perl ~/tools/ASE_ScRNA_Breast/data_breastcancer_fastq/scripts/filter_sam_v2.pl {}Aligned.out.bam {}Aligned.out.filtered.sam && samtools view -bS {}Aligned.out.filtered.sam > {}Aligned.out.filtered.bam" ::: "${samples[@]}"

# 3.4:  Sort .bam files
parallel -j 5 "picard SortSam INPUT={}Aligned.out.filtered.bam OUTPUT={}Aligned.out.filtered.sorted.bam SORT_ORDER=coordinate" ::: "${samples[@]}"

# 3.5: Split reads with N in cigar
parallel -j 5 "gatk SplitNCigarReads -R $reference -I {}Aligned.out.filtered.sorted.bam -O {}Aligned.out.filtered.sorted.split.bam" ::: "${samples[@]}"

# 3.6: Add read group, and index
for samp in "${samples[@]}"; do
    # Get the current directory name and store it in a variable
    dir_name=$(basename "$PWD")
    # Use the variable in the picard command
    picard AddOrReplaceReadGroups INPUT="${samp}Aligned.out.filtered.sorted.split.bam" OUTPUT="${samp}Aligned.out.filtered.sorted.rg.bam" RGID="${samp}" RGLB=transcriptomic RGPL=ILLUMINA RGPU=machine RGSM="$dir_name"
    samtools index "${samp}Aligned.out.filtered.sorted.rg.bam"

    # 3.7: append to rg.bam.list
    echo "${samp}Aligned.out.filtered.sorted.rg.bam" >> rg.bam.list
done

# Step 3.8: Mark duplicates
parallel -j 5 "picard MarkDuplicates INPUT={}Aligned.out.filtered.sorted.rg.bam  OUTPUT={}Aligned.out.filtered.sorted.dedup.rg.bam  METRICS_FILE={}Aligned.dedup.metrics.txt" ::: "${samples[@]}"

# Step 3.9: Base Recalibrator
parallel -j 5 "gatk BaseRecalibrator -I {}Aligned.out.filtered.sorted.dedup.rg.bam  -R $reference --known-sites $knownSites -O {}recal_data.table" ::: "${samples[@]}"

# Step 3.10: Base quality score recalibration
parallel -j 5 "gatk ApplyBQSR -R $reference -I {}Aligned.out.filtered.sorted.dedup.rg.bam --bqsr-recal-file {}recal_data.table -O {}Aligned.out.filtered.sorted.dedup.bqsr.rg.bam" ::: "${samples[@]}"

##############################################################################################
##############################################################################################
########################### STEP 4: CALL VARIANTS ON SINGLE FILES ############################
##############################################################################################
##############################################################################################

# Step 4.1: Create VCF file from BAM file
parallel -j 5 "bcftools mpileup -Ob -o {}.bcf -f $reference {}Aligned.out.filtered.sorted.dedup.bqsr.rg.bam" ::: "${samples[@]}"
parallel -j 5 "bcftools call -vmO z -o {}.vcf.gz {}.bcf" ::: "${samples[@]}"
parallel -j 5 "bcftools index -t {}.vcf.gz" ::: "${samples[@]}"


##############################################################################################
##############################################################################################
########################### STEP 5: ANNOTATION AND MPILEUP ###################################
##############################################################################################
##############################################################################################
# Step 5.1: Annotate VCF file with rsid
parallel -j 5 "bcftools annotate -c ID -a $knownSitesCompressed {}.vcf.gz > {}.rs.vcf.gz" ::: "${samples[@]}"

#Step 5.2: Extract rsid from rs file
parallel -j 5 "bcftools query -f '%CHROM\t%POS\t%REF\t%ALT\t%ID\n' {}.rs.vcf.gz > {}.rsid.txt" ::: "${samples[@]}"

#Step 5.3: annotate with gene names using annovar
parallel -j 5 "perl ~/tools/annovar/table_annovar.pl {}.vcf.gz ~/tools/annovar/humandb/ -buildver hg38 -out {}.gene -remove -protocol refGene -operation g -nastring . -vcfinput" ::: "${samples[@]}"

#Step 5.4: extract gene names
parallel -j 5 "bcftools query -f '%CHROM\t%POS\t%INFO/Gene.refGene\n' {}.gene.hg38_multianno.vcf > {}.gene.txt" ::: "${samples[@]}"
EOF
# Step 5.5: samtools mpileup
echo "Samtools mpileup in action..."
parallel -j 5 "samtools mpileup -q 25 {}Aligned.out.filtered.sorted.dedup.bqsr.rg.bam  --reference $reference -o {}.q25.mpileup" ::: "${samples[@]}"
echo "Samtools mpileup calling completed."

# Step 5.6: mpileup to base
parallel -j 5 "perl $pile2base {}.q25.mpileup 25 {}.parse25.txt" ::: "${samples[@]}"

##############################################################################################
##############################################################################################
########################### STEP 6:Post processing files ###########################################
##############################################################################################
##############################################################################################

#Step 6.1: Add gene and rsid information to base file
python ~/tools/ASE_ScRNA_Breast/data_breastcancer_fastq/scripts/02_scRNA_ann_b4merge.py

#Step 6.2: Parse base file to generate genotyping information
python ~/tools/ASE_ScRNA_Breast/data_breastcancer_fastq/scripts/03_scRNA_parser.py

#Step 6.3: Calculate Lesser Allele fraction of cells
echo "Calculating LAF of individual cells.."
python ~/tools/ASE_ScRNA_Breast/data_breastcancer_fastq/scripts/04_LAF.py
echo "LAF calculation of individual cell is finished :-) "

#Step 6.4: Find CNV information from reference paper and add in file
echo "Matching CNV from reference paper.."
python ~/tools/ASE_ScRNA_Breast/data_breastcancer_fastq/scripts/05_CNV.py
echo "CNV information from reference paper is merged...:-)"

#Step 6.5: Add cells annotation by checking RGIDs
#Note to self: there is an alternative version for code so please double check
echo "Adding cell annotations by matching RGIDs..."
python ~/tools/ASE_ScRNA_Breast/data_breastcancer_fastq/scripts/06_RGID_cell.py
echo "Cell annotations are complete :-)"

#Step 6.6: Find list of imprinted genes and store results
echo "Finding loss of imprinting..."
python ~/tools/ASE_ScRNA_Breast/data_breastcancer_fastq/scripts/07_imprinted_genes.py
echo "All imprinted genes are in the new file. Good Luck myself :-) "

echo "All scripts executed."


##############################################################################################
##############################################################################################
########################### CALL VARIANTS ON MERGED FILE #####################################
##############################################################################################
##############################################################################################

merged_bam="merged.bam"
#samtools merge "$merged_bam" "${samples[@]/%/Aligned.out.filtered.sorted.dedup.rg.bam}"

# Create an array of merged samples
merged_samples=("$merged_bam")

# Continue with remaining steps on the merged BAM
for merged_sample in "${merged_samples[@]}"; do
    #Step 15: Create VCF file from merged BAM
    bcftools mpileup -Ob -o "${merged_sample}.bcf" -f "$reference" "$merged_sample"
    bcftools call -vmO z -o "${merged_sample}.vcf.gz" "${merged_sample}.bcf"
    bcftools index -t "${merged_sample}.vcf.gz"

    # Step 17: samtools mpileup
    samtools mpileup -q 25 "$merged_sample" --reference "$reference" -o "${merged_sample}.q25.mpileup"
    
    # Step 19: mpileup to base
    perl "$pile2base" merged.bam.q25.mpileup 25 merged.25.out
   
    #Step 20: merge parser
    python ../merge_parser.py
  
done
