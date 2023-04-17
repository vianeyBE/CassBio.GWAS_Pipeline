# Cassava Bioinformatics Platform: GWAS Pipeline
Pipeline for Genome-wide association studies in Cassava

**Authors**: Vianey Barrera-Enriquez and Camilo E. Sanchez

The pipeline has six main steps:

1. Quality control (*in progress*)
2. Imputation (*in progress*)
3. LD Decay
4. Population structure and covariable selection (*in progress*)
5. GWAS analysis using GAPIT3
6. Annotation of results (*in progress*)
7. Boxplot of significant markers: genotypes vs. phenotype
8. Customizable Manhattan plots (*in progress*)

## 1. Quality control

*In progress*

```sh
inputFile = "/path/to/inputFile.vcf.gz"
prefix = "gwas_results"

# Get lisf of samples from VCF 
bcftools query -l ${inputFile} > ${prefix}.list.samples.txt

# Get stats
bcftools stats ${inputFile} > ${prefix}.stats.txt

# Calculate allele frequency 
vcftools --gzvcf ${inputFile} --freq2 --out ${prefix}.frequency

# Stats mean depth per individual
vcftools --gzvcf ${inputFile} --depth --out ${prefix}.depth

# Mean depth per site
vcftools --gzvcf ${inputFile} --site-mean-depth --out ${prefix}.depth_site

# Site quality
vcftools --gzvcf ${inputFile} --site-quality --out ${prefix}.quality

# Missing data per individual
vcftools --gzvcf ${inputFile} --missing-indv --out ${prefix}.missing

# Missing data per site
vcftools --gzvcf ${inputFile} --missing-site --out ${prefix}.missing_site

# Heterozygosity and inbreeding coefficients per individual 
vcftools --gzvcf ${inputFile} --het --out ${prefix}.heterozygosity

# Hardy-Weinberg p-value 
vcftools --gzvcf ${inputFile} --hardy --out ${prefix}.hwe
```

## 2. Imputation

*In progress*

## 3. LD Decay

*In progress*

## 4. Population structure and covariable selection

*In progress*

## 5. GWAS: GAPIT

This R script runs a Genome-Wide Association Study (GWAS) analysis using the GAPIT3 R package. It saves the results of each trait in an individual folder.

### Usage

```R
GAPIT3(phenofile, genofile, wdir, trait_list = NULL)
```

### Arguments

- `phenofile`: measurements or BLUPs of the traits. Rows are individuals, and columns are traits/variables.
- `genofile`: genotype data in hapmap format.
- `wdir`: working directory. Here, it will create the folders for each trait.
- `trait_list`: a vector with the trait's name to do GWAS. The default is `NULL`, which will use the full phenotype dataset.

### Details

This function loads the required packages and data, then performs a GWAS analysis using the `GAPIT` function from the GAPIT3 package. It loops through each trait in `trait_list`, creating a new folder for each trait in the working directory and saving the results of the GWAS analysis for that trait in that folder.

### Examples

```R
GAPIT3(phenofile = "my_phenotypes.csv", genofile = "my_genotypes.hmp", wdir = "my_results_folder", trait_list = c("Trait1", "Trait2", "Trait3"))
```

This will perform a GWAS analysis using the phenotype data in ; my_phenotypes.csv and the genotype data in my_genotypes.hmp. It will create a new folder for each trait in the trait_list vector in the my_results_folder directory.

### Output
The function will create a folder for each trait in the trait_list vector in the working directory. In each folder, it will save the results of the GWAS analysis for that trait. For more information about GAPIT outputs, please check the official GAPIT documentation. 

### Dependencies
- `devtools`
- `GAPIT3`

## 6. Annotation of results

This code annotates the gen containing the significat SNPs from the GAPIT results. Additional, it retrieves the closest genes downstream and upstream.

### Arguments:
- `pat`: Location path of the files.
- `wdyw`: Gene feature to annotate.
- `Mdir`: Name of the directory that contains the GAPIT results.
- `mod`: Names of the models to filter.
- `gff3`: gff3 file from the genome version used for alignment
- `annotationFile`: annotation details of the genes. txt file from the genome version used for alignment. 

### Output
A single CSV file containing relevant gene information plus SNPs' P-value, trait, model and effect. 

### Dependencies
- `tidyverse`

## 7. Boxplot: Genotype vs Phenotype

This R function generates a boxplot for a given SNP. The function takes as inputs a CSV file with the phenotype values and a list of SNPs. The function can also receive an optional CSV file to add extra information about the samples to the plot (categories, family, ect.) The output of the function is a PDF file with the plot.

### Usage

```R
GWAS_Boxplot(outputname, dir, phenofile, genofile, snp_list_file, labelfile = NULL)
```

### Arguments
- `outputname`: (required) character string with the base name for the output file.
- `dir`: (required) character string with the directory where the output file will be saved.
- `phenofile`: (required) character string with the name of the phenotype file in tabular format. The first column should contain the sample names, and the rest of the columns should contain the phenotypes.
- `genofile`: (required) character string with the name of the genotype file in hapmap format.
- `snp_list_file`: (required) character string with the name of the CSV file with three columns:
    - Column 01 - name: SNPS, list of SNPS to plot, name should be the same as in the geno data.
    - Column 02 - name: trait, name of the trait as in the pheno data.
    - Column 03 - name: xlabel, name of the trait to be included as a label.
- `labelfile`: (optional) character string with the name of the CSV file with two columns:
    - Column 01 - name: Taxa, sample names.
    - Column 02 - name: label, label or category to add to the plot.


### Details

The function loads all the necessary packages and data files. It then checks for the presence of an optional file with sample labels and prepares the data for the boxplot. The function generates a boxplot for each SNP specified in the input CSV file. It uses ggplot2 to generate the plot, with the genotypes on the x-axis and the phenotype on the y-axis. The function adds a label to the x-axis to indicate the trait being plotted. If a labelfile is provided, it adds extra information about the samples, coloring the data by the levels in the provided file.

### Example

```R

# Generate boxplot without extra labels
GWAS_Boxplot("outputname", ".path/to/save/plots/", "phenotype.csv", "genotype.hmp", "snp_list.csv")

# Generate boxplot with extra labels
GWAS_Boxplot("outputname", ".path/to/save/plots/", "phenotype.csv", "genotype.hmp", "snp_list.csv", "labelfile.csv")

```

### Output
A single PDF file containing the boxplot of the SNPs. 

### Dependencies
- `tidyverse`
- `tibble`
- `dplyr`
- `janitor`
- `ggplot2`
- `Biostrings`
- `hrbrthemes`
- `forcats`
- `ggsignif`
- `RColorBrewer`

## 8. Manhattan plots

*In progress*

---
## Citation

- Petr Danecek, James K Bonfield, Jennifer Liddle, John Marshall, Valeriu Ohan, Martin O Pollard, Andrew Whitwham, Thomas Keane, Shane A McCarthy, Robert M Davies, Heng Li, Twelve years of SAMtools and BCFtools, GigaScience, Volume 10, Issue 2, February 2021, giab008, https://doi.org/10.1093/gigascience/giab008

- Petr Danecek, Adam Auton, Goncalo Abecasis, Cornelis A. Albers, Eric Banks, Mark A. DePristo, Robert E. Handsaker, Gerton Lunter, Gabor T. Marth, Stephen T. Sherry, Gilean McVean, Richard Durbin, 1000 Genomes Project Analysis Group, The variant call format and VCFtools, Bioinformatics, Volume 27, Issue 15, August 2011, Pages 2156–2158, https://doi.org/10.1093/bioinformatics/btr330

- B L Browning, X Tian, Y Zhou, and S R Browning (2021) Fast two-stage phasing of large-scale sequence data. Am J Hum Genet 108(10):1880-1890. doi:10.1016/j.ajhg.2021.08.005

- B L Browning, Y Zhou, and S R Browning (2018). A one-penny imputed genome from next generation reference panels. Am J Hum Genet 103(3):338-348. doi:10.1016/j.ajhg.2018.07.015

- Chi Zhang, Shan-Shan Dong, Jun-Yang Xu, Wei-Ming He, Tie-Lin Yang, PopLDdecay: a fast and effective tool for linkage disequilibrium decay analysis based on variant call format files, Bioinformatics, Volume 35, Issue 10, May 2019, Pages 1786–1788, https://doi.org/10.1093/bioinformatics/bty875

- Thibaut Jombart, adegenet: a R package for the multivariate analysis of genetic markers, Bioinformatics, Volume 24, Issue 11, June 2008, Pages 1403–1405, https://doi.org/10.1093/bioinformatics/btn129

- Xiuwen Zheng, David Levine, Jess Shen, Stephanie M. Gogarten, Cathy Laurie, Bruce S. Weir, A high-performance computing toolset for relatedness and principal component analysis of SNP data, Bioinformatics, Volume 28, Issue 24, December 2012, Pages 3326–3328, https://doi.org/10.1093/bioinformatics/bts606

- Wang J., Zhang Z., GAPIT Version 3: Boosting Power and Accuracy for Genomic Association and Prediction, Genomics, Proteomics & Bioinformatics (2021), doi: https://doi.org/10.1016/j.gpb.2021.08.005.

- Huang M, Liu X, Zhou Y, Summers RM, Zhang Z. BLINK: A package for the next level of genome-wide association studies with both individuals and markers in the millions. Gigascience. https://doi.org/10.1093/gigascience/giy154.

- Liu X., Huang M., Fan B., Buckler E. S., Zhang Z., 2016 Iterative Usage of Fixed and Random Effect Models for Powerful and Efficient Genome-Wide Association Studies. PLoS Genet. 12: e1005767. https://doi.org/10.1371/journal.pgen.1005767.

- David M. Goodstein, Shengqiang Shu, Russell Howson, Rochak Neupane, Richard D. Hayes, Joni Fazo, Therese Mitros, William Dirks, Uffe Hellsten, Nicholas Putnam, and Daniel S. Rokhsar, Phytozome: a comparative platform for green plant genomics, Nucleic Acids Res. 2012 40 (D1): D1178-D1186

- R Core Team (2021). R: A language and environment for statistical computing. R Foundation for Statistical Computing, Vienna, Austria. URL https://www.R-project.org/.

## Contact
For questions or feedback about this pipeline, please contact Vianey Barrera-Enriquez at v.barrera@cgiar.org.
