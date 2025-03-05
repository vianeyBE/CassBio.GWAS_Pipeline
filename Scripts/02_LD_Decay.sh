#!/bin/bash

# Check if all arguments are provided
if [[ $# -ne 6 ]]; then
    echo "Usage: $0 <dirIn> <vcfFile> <dirOut> <prefix> <dist> <mode>"
    echo "Example: $0 /path/to/vcfs mydata.vcf /path/to/output myprefix 5000 by_chr"
    exit 1
fi

# Paths to tools (adjust these paths to your installation)
PopLDdecay="/datas3/Cassava/Cassava_basics/Tools/Programs/PopLDdecay/bin/PopLDdecay"
PlotOnePop="/datas3/Cassava/Cassava_basics/Tools/Programs/PopLDdecay/bin/Plot_OnePop.pl"

# Input arguments
dirIn=$1
vcfFile=$2
dirOut=$3
prefix=$4
dist=$5
mode=$6

# Create necessary directories
mkdir -p ${dirOut}
mkdir -p ${dirIn}/split_chr

# Start log (optional but recommended)
exec > >(tee ${dirOut}/ld_decay_analysis_${prefix}.log) 2>&1

# Print summary of parameters
echo "Starting LD Decay Analysis"
echo "Input directory: ${dirIn}"
echo "VCF file: ${vcfFile}"
echo "Output directory: ${dirOut}"
echo "Prefix: ${prefix}"
echo "Max distance (bp): ${dist}"
echo "Mode: ${mode}"

# Function for whole-genome analysis
run_whole_genome() {
    echo "Running LD Decay for the whole genome"

    ${PopLDdecay} -InVCF ${dirIn}/${vcfFile} \
                  -OutStat ${dirOut}/${prefix}.PopLDdecay.MaxDist_${dist}.stat \
                  -MaxDist ${dist}

    gzip -f ${dirOut}/${prefix}.PopLDdecay.MaxDist_${dist}.stat

    perl ${PlotOnePop} -inFile ${dirOut}/${prefix}.PopLDdecay.MaxDist_${dist}.stat.gz \
                       -output ${dirOut}/${prefix}.PopLDdecay.MaxDist_${dist}.plot

    echo "Whole genome analysis completed"
}

# Function for per-chromosome analysis
run_by_chr() {
    echo "Running LD Decay per chromosome"

    for chr in {01..18}; do
        chr_id="chr${chr}"
        echo "Processing $chr_id"

        vcftools --vcf ${dirIn}/${vcfFile} \
                 --chr ${chr_id} --recode --out ${dirIn}/split_chr/${prefix}.${chr_id}

        ${PopLDdecay} -InVCF ${dirIn}/split_chr/${prefix}.${chr_id}.recode.vcf \
                      -OutStat ${dirOut}/${prefix}.${chr_id}.PopLDdecay.MaxDist_${dist}.stat \
                      -MaxDist ${dist}

        gzip -f ${dirOut}/${prefix}.${chr_id}.PopLDdecay.MaxDist_${dist}.stat

        perl ${PlotOnePop} -inFile ${dirOut}/${prefix}.${chr_id}.PopLDdecay.MaxDist_${dist}.stat.gz \
                           -output ${dirOut}/${prefix}.${chr_id}.MaxDist_${dist}.plot

        echo "Chromosome $chr_id completed"
    done

    echo "Per-chromosome analysis completed"
}

# Run the appropriate analysis based on the mode
if [[ $mode == "whole_genome" ]]; then
    run_whole_genome
elif [[ $mode == "by_chr" ]]; then
    run_by_chr
else
    echo "ERROR: Mode '$mode' is not recognized. Use 'by_chr' or 'whole_genome'."
    exit 1
fi

echo "LD Decay Analysis completed. Results are saved in: ${dirOut}"
