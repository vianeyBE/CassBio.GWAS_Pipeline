#!/bin/bash

# Initial requeriments
input_vcf=$1
output=$2

# BGI filtering code
# vcftools --gzvcf ${input_vcf} --max-missing 0.90 --maf 0.05 --recode --recode-INFO-all --out ${output}

# GATK filtering code
vcftools --gzvcf ${input_vcf} --max-missing 0.90 --maf 0.05 --minDP 6 --maxDP 40 --minQ 30 --recode --recode-INFO-all --out ${output}

# Another filtering code
# vcftools --gzvcf ${input_vcf} --max-missing 0.90 --maf 0.05 --minDP 10 --maxDP 40 --min-meanDP 10 --max-meanDP 40 --minQ 30 --recode --recode-INFO-all --out ${output}