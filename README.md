# Cassava Bioinformatics Platform: GWAS Pipeline

Pipeline for Genome-Wide Association Studies (GWAS) in Cassava (Manihot esculenta)

## Description

This repository contains a modular pipeline for performing genome-wide association studies (GWAS) in cassava (Manihot esculenta). It covers the workflow from genotype quality control, linkage disequilibrium analysis, population structure analysis, and GWAS analysis, to functional annotation of significant markers.

## Authors

- **Vianey Barrera-Enriquez** - [CGIAR]
- **Camilo E. Sanchez** - [CGIAR]

## Requirements

### Technologies and tools

- R (>= 4.1.0)
- GAPIT3
- VCFtools
- BCFtools
- PopLDdecay
- adegenet, vegan, SNPRelate, plotly, tidyverse, etc.

### Installation

```bash

git clone https://github.com/your-username/cassava-gwas-pipeline.git

```

## ⚙️ Workflow Overview
| 1 | Quality control | ✅ Completed |
| 2 | LD decay analysis | 🔧 In progress |
| 3 | LD pruning | 🔧 In progress |
| 4 | Population structure analysis | ✅ Completed |
| 5 | GWAS - GAPIT3 | ✅ Completed |
| 6 | GWAS - EMMAX | 🔧 In progress |
| 7 | Functional annotation | 🔧 In progress |
| 8 | Marker validation boxplots | ✅ Completed |




## 1. Quality control

### Description

A simple Bash script for filtering Variant Call Format (VCF) files using VCFtools. The script applies customizable quality control filters to genetic variant data. The script currently includes three different filtering strategies (commented and active), allowing users to select the most appropriate filter set for their project.

### Arguments

- `input_vcf`: Path to VCF file to filter.
- `output`: Prefix to the output VCF filtered file.

### Usage

```sh

bash 01_QC.sh input_vcf output

```

### Example

```sh

bash 01_QC.sh raw_variants.vcf.gz filtered_variants

```

### Dependencies

- `VCFtools`: (https://vcftools.github.io/)








## 2. LD decay

### Description

*In progress*

### Usage

```sh

bash 02_LD_Decay.sh

```

### Arguments

- `XXX`: 
- `XXX`: 

### Dependencies

- `XXX`
- `XXX`







## 3. LD pruning

*In progress*








## 4. Population structure and covariable selection

This script contains the customized functions to perform four different multivariate methods: (1) Non-Metric Multidimensional Scaling (NDMS), (2) Multidimensional Scaling (MDS), (3) Principal Component Analysis (PCA), and (4) Discriminant Analysis of Principal Components (DAPC). This series of functions retrieves the required pacakges and data to perform four multivariate methods, potentially used to study population structure and select covariable for next analyses. Especifically, this script includes the functions for Non-Metric Multidimensional Scaling (`NDMS`), Multidimensional Scaling (`MDS`), Principal Component Analysis (`PCA`), and Discriminant Analysis of Principal Components (`DAPC`). The first two functions are built only for phenotypic data while the PCA can handle both data type (phenotypic or genotypic), and the DAPC only uses genotypic data. For all functions interactive plots (made using `plotly` package) in html format. The DAPC part is not built as a function give the nature of the functions used in which prompts are necessary to continue. Please run this methods line by line.

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

GWAS_Annotation(Mdir, pat, mod, wdyw, annot, GFF)

```

### Arguments

- `pat`: Location path of the files.
- `wdyw`: Gene feature to annotate.
- `Mdir`: Name of the directory that contains the `GAPIT3` results.
- `mod`: Names of the models to filter (options: )
- `gff3`: gff3 file from the genome version used for alignment.
- `annotationFile`: Annotation details of the genes. txt file from the genome version used for alignment.
- `version`: (Options: 6.1 or 8.1. Default = 6.1).

### Details

*In progress*

### Examples

*In progress*

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

### Dependencies

- `tidyverse`
- `ggtext`
- `plotly`

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
