# Cassava Bioinformatics Platform: GWAS Pipeline

Pipeline for Genome-wide Association Studies in Cassava

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

### Usage

*In progress*

### Arguments

*In progress*

### Details

*In progress*

### Examples

```sh

# Initial requeriments
inputFile <- "/path/to/inputFile.vcf.gz"
prefix <- "gwas_results"

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

### Output

*In progress*

### Dependencies

*In progress*

## 2. Imputation

*In progress*

## 3. LD Decay

*In progress*

## 4. Population structure and covariable selection

This script contains the customized functions to perform four different multivariate methods: (1) Non-Metric Multidimensional Scaling (NDMS), (2) Multidimensional Scaling (MDS), (3) Principal Component Analysis (PCA), and (4) Discriminant Analysis of Principal Components (DAPC).

### Usage

```R

# NDMS, MDS, and PCA
NDMS(dir, phenofile, dist = "bray", groups = F)
MDS(dir, phenofile, dist = "gower", groups = F)
PCA(dir, type, groups = F, phenofile = NULL, labels = NULL, genofile = NULL, vcf = NULL, gds = NULL, PC.retain = F)

# DAPC
library(adegenet)
grp <- find.clusters(my_genind)
dapc <- dapc(my_genind, grp$grp)

```

### Arguments

1. For NDMS:
- `dir`: Directory where data is located.
- `phenofile`: A database of genotypes/individuals (rows) and their traits (columns). First column must be genotypes/individuals names.
- `dist`: Dissimilarity index/measure to use (default = "bray").
- `groups`: Boolean value indicating whether the data includes different treatments/groups (default = F). If `TRUE` is selected then the second column of the `phenofile` must be the groups/treatments.

2. For MDS:
- `dir`: Directory where data is located.
- `phenofile`: A database of genotypes/individuals (rows) and their traits (columns). First column must be genotypes/individuals names.
- `dist`: Dissimilarity index/measure to use (default = "gower").
- `groups`: Boolean value indicating whether the data includes different treatments/groups (default = F). If `TRUE` is selected then the second column of the `phenofile` must be the groups/treatments.

3. For PCA:
- `dir`: Directory where data is located.
- `type`: A character string indicating wheter data is genotypic ("geno") of phenotypic ("pheno").
- `groups`: Boolean value indicating whether the data includes different treatments/groups (default = F). If `TRUE` is selected then the second column of the `phenofile` must be the groups/treatments.
- `phenofile`: A database of genotypes/individuals (rows) and their traits (columns). First column must be genotypes/individuals names.
- `genofile`: An object of class SNPGDSFileClass (GDS file read with the `snpgdsOpen` function from `SNPRelate` package).
- `labels`: When provide a `genofile` and `groups` argument is `TRUE`, please provide a dataframe with genotypes/individuals in the first column and groups/treatment data in the second column.
- `gds`: A Genomic Data Structures (GDS) file (a reformatted VCF file with the `snpgdsVCF2GDS` function from `vcfR` package).
- `vcf`: A Variant Call Format (VCF) file containing DNA polymorphism data such as SNPs, insertions, deletions and structural variants, together with rich annotations.
- `PC.retain`: Boolean value indicating whether to analyze how many PCs retain (default = F).

4. For DAPC:
- `my_genind`: An object of class genind (vcf file read with the `vcfR2genind` function from `vcfR` package).

### Details

This series of functions retrieves the required pacakges and data to perform four multivariate methods, potentially used to study population structure and select covariable for next analyses. Especifically, this script includes the functions for Non-Metric Multidimensional Scaling (`NDMS`), Multidimensional Scaling (`MDS`), Principal Component Analysis (`PCA`), and Discriminant Analysis of Principal Components (`DAPC`). The first two functions are built only for phenotypic data while the PCA can handle both data type (phenotypic or genotypic), and the DAPC only uses genotypic data. For all functions interactive plots (made using `plotly` package) in html format. The DAPC part is not built as a function give the nature of the functions used in which prompts are necessary to continue. Please run this methods line by line.

### Examples

```R

# NDMS example(s)
# Set arguments
# Phenotypic with groups
dir <- "D:/OneDrive - CGIAR/Cassava_Bioinformatics_Team/02_CTS_Drought_Family/01_Phenotype_Preliminar_Analysis/"
phenofile <- read.csv(paste0(dir, "Prueba.csv"), header = T)
groups <- T

# Phenotypic without groups
dir <- "D:/OneDrive - CGIAR/Cassava_Bioinformatics_Team/02_CTS_Drought_Family/01_Phenotype_Preliminar_Analysis/"
phenofile <- read.csv(paste0(dir, "Prueba.csv"), header = T)
phenofile <- phenofile[-2]

# Run function
NDMS(dir, phenofile, dist, groups)



# MDS example(s)
# Set arguments
# Phenotypic with groups
dir <- "D:/OneDrive - CGIAR/Cassava_Bioinformatics_Team/02_CTS_Drought_Family/01_Phenotype_Preliminar_Analysis/"
phenofile <- read.csv(paste0(dir, "Prueba.csv"), header = T)
groups <- T

# Phenotypic without groups
dir <- "D:/OneDrive - CGIAR/Cassava_Bioinformatics_Team/02_CTS_Drought_Family/01_Phenotype_Preliminar_Analysis/"
phenofile <- read.csv(paste0(dir, "Prueba.csv"), header = T)
phenofile <- phenofile[-2]

# Run function
MDS(dir, phenofile, dist, groups)



# PCA example(s)
# Set arguments
# Phenotypic with groups
dir <- "D:/OneDrive - CGIAR/Cassava_Bioinformatics_Team/02_CTS_Drought_Family/01_Phenotype_Preliminar_Analysis/"
phenofile <- read.csv(paste0(dir, "Prueba.csv"), header = T)
groups <- T

# Phenotypic without groups
dir <- "D:/OneDrive - CGIAR/Cassava_Bioinformatics_Team/02_CTS_Drought_Family/01_Phenotype_Preliminar_Analysis/"
phenofile <- read.csv(paste0(dir, "Prueba.csv"), header = T)
phenofile <- phenofile[-2]

# Genotypic with groups
dir <- "D:/OneDrive - CGIAR/Cassava_Bioinformatics_Team/03_GWAS_PPD_Populations/02_PCA/"
labels <- read.csv(paste0(dir, "GWAS_PPD.labels.csv"))
vcf <- "GWAS_PPD.snps.filter_info.missing_0.10.imputation.vcf.gz"
type <- "Geno"
groups <- T

dir <- "D:/OneDrive - CGIAR/Cassava_Bioinformatics_Team/01_ACWP_F2_Phenotype/01_Population_Structure/"
labels <- read.csv(paste0(dir, "AM1588_labels.csv"))
vcf <- "AM1588_MAP.miss0.05.recode.vcf"
type <- "Geno"
groups <- T

# Run function
PCA(dir, phenofile, genofile, gds, vcf, type = "Geno", groups = T, PC.retain = F)



# DAPC example(s)
# Set arguments
dir <- "D:/OneDrive - CGIAR/Cassava_Bioinformatics_Team/01_ACWP_F2_Phenotype/01_Population_Structure/"
vcf <- "AM1588_MAP.miss0.05.recode.vcf"

dir <- "D:/OneDrive - CGIAR/Cassava_Bioinformatics_Team/03_GWAS_PPD_Populations/02_PCA/"
vcf <- "GWAS_PPD.snps.filter_info.missing_0.10.imputation.vcf.gz"

# Load libraries
library(adegenet)
library(grDevices)
library(vcfR)
library(tidyverse)
library(plotly)
library(htmlwidgets)

# Load VCF files and convert them into a genind object
setwd(dir)
vcf <- read.vcfR(vcf, verbose = T)
my_genind <- vcfR2genind(vcf)
my_genind

# Identify clusters
# Shows a graph with the accumulated variance explained by the eigenvalues of the PCA
grp <- find.clusters(my_genind)

# Performs the DAPC
dapc <- dapc(my_genind, grp$grp)

```

### Output

Interactive plots (made using `plotly` package) in html format.

### Dependencies

- `vegan`
- `tidyverse`
- `plotly`
- `htmlwidgets`
- `SNPRelate`
- `gdsfmt`
- `psych`
- `adegenet`
- `grDevices`
- `vcfR`

## 5. GWAS: GAPIT

This R script runs a Genome-Wide Association Study (GWAS) analysis using the GAPIT3 R package. It saves the results of each trait in an individual folder.

### Usage

```R

GAPIT3(phenofile, genofile, wdir, trait_list = NULL)

```

### Arguments

- `phenofile`: A database of genotypes/individuals (rows) and the trait's measurements or BLUPs (columns).
- `genofile`: Genotype data in hapmap format.
- `wdir`: A working directory where the function will create the folders for each trait.
- `trait_list`: A vector with the trait's name. The default is `NULL`, which will use the full phenotype dataset.

### Details

This function loads the required packages and data, then performs a GWAS analysis using the `GAPIT` function from the `GAPIT3` package. It loops through each trait in `trait_list`, creating a new folder for each trait in the working directory and saving the results of the GWAS analysis for that trait in that folder.

### Examples

```R

GAPIT3(phenofile = "my_phenotypes.csv", genofile = "my_genotypes.hmp", wdir = "my_results_folder", trait_list = c("Trait1", "Trait2", "Trait3"))

```

This will perform a GWAS analysis using the phenotype data in `my_phenotypes.csv` and the genotype data in `my_genotypes.hmp`. It will create a new folder for each trait in the `trait_list` vector in the `my_results_folder` directory.

### Output

The function will create a folder for each trait in the `trait_list` vector in the working directory. In each folder, it will save the results of the GWAS analysis for that trait. For more information about `GAPIT` function outputs, please check the official `GAPIT3` package documentation.

### Dependencies

- `devtools`
- `GAPIT3`

## 6. Annotation of results

This code annotates the gen containing the significat SNPs from the `GAPIT3` results. Additional, it retrieves the closest genes downstream and upstream.

### Usage

```R

GWAS_Annotation(Mdir, pat, mod, wdyw)

```

### Arguments

- `Wdir`: Name of the directory that contains the GAPIT results. For example: home/user/folder.
- `Ddir`: Directory where is located the annotation files (annot, GFF files).
- `pat`: Enter the path of file names to look for. For example: GAPIT.Association.GWAS_Results. The path must finish with a point (.).
- `mod`: Enter the model(s) of interest (Options: BLINK, GLM, MLM, FarmCPU).
- `wdyw`: Enter what are you looking for to annotate (Options: CDS, five_prime_UTR, gene, mRNA, three_prime_UTR).
- `annot`: Annotation details of the genes. txt file from the genome version used for alignment.
- `GFF`: gff3 file from the genome version used for alignment.
- `version`: (Optional) You can choose between the genome of reference version 6.1 or 8.1 (Options: 6.1 or 8.1. Default = 6.1).

### Details

*In progress*

### Examples

```R

# Set arguments
Wdir <- "D:/OneDrive - CGIAR/Cassava_Bioinformatics_Team/03_GWAS_PPD_Populations/04_GWAS/GAPIT_Results/"
Ddir <- "D:/OneDrive - CGIAR/Cassava_Bioinformatics_Team/00_Data/"
pat <- "GAPIT.Association.GWAS_Results."
mod <- c("BLINK", "FarmCPU", "MLM")
wdyw <- "gene"

# Run function
GWAS_Annotation(Wdir, pat, mod, wdyw)

```

### Output

A single CSV file containing relevant gene information plus SNPs' P-values, traits, models, and effects.

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
    - Column 01 - Name: SNPS. List of SNPS to plot. The name should be the same as in the `genofile` data.
    - Column 02 - Name: trait. Name of the trait as in the `phenofile` data.
    - Column 03 - Name: xlabel. Name of the trait to be included as a label.
- `labelfile`: (optional) character string with the name of the CSV file with two columns:
    - Column 01 - Name: Taxa. Sample names.
    - Column 02 - Name: label. Label or category to add to the plot.

### Details

The function loads all the necessary packages and data files. It then checks for the presence of an optional file with sample labels and prepares the data for the boxplot. The function generates a boxplot for each SNP specified in the input CSV file. It uses `ggplot2` package to generate the plot, with the genotypes on the x-axis and the phenotype on the y-axis. The function adds a label to the x-axis to indicate the trait being plotted. If `labelfile` is provided, it adds extra information about the samples, coloring the data by the levels in the provided file.

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

This R function generates a Manhattan plot(s) for a given set of GWAS results. The function scans whithin directories and takes as inputs the CSV file that contain the results of previous GWAS analyses. The output of the function is a plot in html format made using `plotly` package.

### Usage

```R

Manhattan(Mdir, pat, mod, wtd)

```

### Arguments

- `Mdir`: Name of the directory that contains the GAPIT results. For example: home/user/folder.
- `pat`: Enter the path of file names to look for whithin directories. For example: QTL_LOD_Intervals. The path must finish with a point (.).
- `mod`: Enter the model(s) of interest. Options: BLINK, GLM, MLM, FarmCPU.
- `wtd`: How many traits do you want to plot. Options: One, Several, All.
- `colors`: (Optional) Colors of the chromosomes in Manhattan plots. If you want to change the colors, provide 2 or more. It can be color names o hex codes (Default: grey and skyblue).

### Details

*In progress*

### Examples

```R

# Set arguments
Mdir <- "D:/OneDrive - CGIAR/Cassava_Bioinformatics_Team/01_ACWP_F1_Metabolomics/07_GWAS"
pat <- "GAPIT.Association.GWAS_Results."
mod <- c("BLINK", "FarmCPU", "MLM")
wtd <- "One"

# Run function
Manhattan(Mdir, pat, mod, wtd)

```

### Output

Manhattan plots for the traits(s) selected.

### Dependencies

- `tidyverse`
- `ggtext`
- `plotly`

---
## References

- Danecek, P., Bonfield, J. K., Liddle, J., Marshall, J., Ohan, V., Pollard, M. O., Whitwham, A., Keane, T., McCarthy, S. A., Davies, R. M., and Li, H. (2021). Twelve years of SAMtools and BCFtools. GigaScience, 10(2): giab008. doi: 10.1093/gigascience/giab008.

- Danecek, P., Auton, A., Abecasis, G., Albers, C. A., Banks, E., DePristo, M. A., Handsaker, R. E., Lunter, G., Marth, G. T., Stephen T. Sherry, McVean, G., Durbin, R., and 1000 Genomes Project Analysis Group. (2011). The variant call format and VCFtools. Bioinformatics, 27(15): 2156 - 2158. doi: 10.1093/bioinformatics/btr330.

- Browning, B. L., Tian, X., Zhou, Y., and Browning, S. R. (2021) Fast two-stage phasing of large-scale sequence data. The American Journal of Human Genetics, 108(10): 1880 - 1890. doi: 10.1016/j.ajhg.2021.08.005.

- Browning, B. L., Zhou, Y., and Browning, S. R. (2018). A one-penny imputed genome from next generation reference panels. The American Journal of Human Genetics, 103(3):338-348. doi: 10.1016/j.ajhg.2018.07.015.

- Zhang, C., Dong, S. S., Xu, J. Y., He, W. M., and Yang, T. L (2019). PopLDdecay: a fast and effective tool for linkage disequilibrium decay analysis based on variant call format files. Bioinformatics, 35(10): 1786 - 1788, doi: 10.1093/bioinformatics/bty875.

- Jombart, T. (2008). adegenet: a R package for the multivariate analysis of genetic markers, Bioinformatics, 24(11): 1403 - 1405. doi: 10.1093/bioinformatics/btn129.

- Zheng, X., Levine, D., Shen, J., Gogarten, S. M., Laurie, C., and Weir, B. C. (2012). A high-performance computing toolset for relatedness and principal component analysis of SNP data. Bioinformatics, 28(24): 3326 - 3328. doi: 10.1093/bioinformatics/bts606.

- Wang, J., and Zhang Z. (2021). GAPIT Version 3: Boosting power and accuracy for genomic association and prediction. Genomics Proteomics Bioinformatics, 19(4): 629 - 640. doi: 10.1016/j.gpb.2021.08.005.

- Huang, M., Liu, X., Zhou, Y., Summers, R. M., and Zhang, Z. BLINK: A package for the next level of genome-wide association studies with both individuals and markers in the millions. GigaScience, 8: 2018, 1 - 12.doi: 10.1093/gigascience/giy154.

- Liu X., Huang M., Fan B., Buckler E. S., and Zhang Z. (2016). Iterative usage of fixed and random effect models for powerful and efficient Genome-Wide Association Studies. PLoS Genetics, 12(3): e1005767. doi: 10.1371/journal.pgen.1005767.

- Goodstein, D. M., Shu, S., Howson, R., Neupane, R., Hayes, R. D., Fazo, J., Mitros, T., Dirks, W., Hellsten, U., Putnam, N., and Rokhsar, D. S. (2012). Phytozome: a comparative platform for green plant genomics. Nucleic Acids Research, 40(D1): D1178 - D1186. doi: 10.1093/nar/gkr944.

- R Core Team. (2021). R: A language and environment for statistical computing. R Foundation for Statistical Computing, Vienna, Austria. URL: https://www.R-project.org/.

## Contact

For questions or feedback about this pipeline, please contact Vianey Barrera-Enriquez at v.barrera@cgiar.org.
