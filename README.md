# Cassava Bioinformatics Platform: GWAS Pipeline 🧬

A modular and reproducible pipeline for performing genome-wide association studies (GWAS) in Manihot esculenta (cassava), covering the complete analysis workflow from genotype quality control to the functional annotation of significant markers.

This pipeline is designed to streamline complex GWAS analyses while ensuring flexibility and scalability for high-throughput plant breeding projects. Each module is independent and well-documented, enabling easy customization or integration into broader genomic workflows.





## Authors 🙋

For questions or feedback about this pipeline, please contact:

- **Vianey Barrera-Enriquez** - [CGIAR]: v.barrera@cgiar.org
- **Camilo E. Sanchez** - [CGIAR]: c.e.sanchez@cgiar.org





## Workflow overview 🚀
| 1 | Quality control - Filters SNPs and samples based on quality metrics
| 2 | Linkgae desequilibirum (LD) decay analysis - Estimates and visualizes LD decay across the genome
| 3 | Linkgae desequilibirum (LD) pruning - Removes SNPs in high LD to reduce redundancy
| 4 | Population structure analysis - Performs PCA to assess population stratification
| 5 | GWAS - GAPIT3 - Runs GWAS using models from the GAPIT3 package
| 6 | GWAS - EMMAX - Performs GWAS using the EMMAX mixed model
| 7 | Functional annotation - Annotates significant SNPs with gene-level information
| 8 | Marker validation boxplots - Generates boxplots to visualize genotype–phenotype effects





## Module descriptions 🧩
Each module contains a README.md describing:
- Description: A brief explanation of what the module does.
- Arguments: Parameters or flags used when running the module.
- Usage: Command-line or script usage.
- Example: Example invocation with test data.
- Dependencies: Required packages, tools, or environments.





## Installation and dependencies 🔧

```bash

git clone https://github.com/vianeyBE/cassava-gwas-pipeline.git

```





## Technologies and tools 🔧

- R (>= 4.1.0): GAPIT, adegenet, vegan, SNPRelate, plotly, tidyverse, etc.
- Bash
- Perl
- VCFtools
- BCFtools
- PopLDdecay
- Snakemake
- EMMAX





## 1. Quality control

### Description

A simple Bash script for filtering Variant Call Format (VCF) files using VCFtools. The script applies customizable quality control filters to genetic variant data. The script currently includes three different filtering strategies (commented and active), allowing users to select the most appropriate filter set for their project.

### Arguments

- `input_vcf`: Path to VCF file to filter.
- `output`: Prefix to the output VCF filtered file.

### Usage

```sh

bash 01_QC.sh <input_vcf> <output>

```

### Example

```sh

bash 01_QC.sh raw_variants.vcf.gz filtered_variants

```

### Dependencies

- `VCFtools`: https://vcftools.github.io/





## 2. LD decay

### Description

This script performs Linkage Disequilibrium (LD) Decay analysis using PopLDdecay. It supports two modes:
- Whole genome: Calculates LD decay across the entire genome.
- Per chromosome: Splits the VCF file by chromosome and calculates LD decay for each chromosome separately.

### Arguments

- `dirIn`: Directory containing the input VCF file.
- `vcfFile`: Name of the input VCF file.
- `dirOut`: Directory where output files will be saved.
- `prefix`: Prefix for output files.
- `dist`: Maximum distance (in base pairs) for pairwise LD calculation.
- `mode`: Analysis mode — either ''whole_genome'' or ''by_chr''.

### Usage

```sh

bash 02_LD_Decay.sh <dirIn> <vcfFile> <dirOut> <prefix> <dist> <mode>

```

### Example

```sh

bash 02_LD_Decay.sh /path/to/vcf gs.vcf /path/to/output gs_2023 10000 whole_genome

```

### Dependencies

- `PopLDdecay`: https://github.com/BGI-shenzhen/PopLDdecay
- `vcftools`: https://vcftools.github.io/
- `perl`: For the PopLDdecay plotting script







## 3. LD pruning

### Description

XXX

### Arguments

- `XXX`: 

### Usage

```sh

bash 

```

### Example

```sh

bash 

```

### Dependencies

- `XXX`: 








## 4. Population structure and covariable selection

XXXXX

### Arguments

- `dir`: Directory where data is located.
- `data`: If phenofile: A csv file of genotypes (rows) and their traits (columns). First column must be genotypes. If genofile: A Variant Call Format (VCF) file.
- `labels`: Dataframe with genotypes in the first column and groups data in the second column.
- `PC.retain`: Boolean value indicating whether to analyze how many PCs retain (default = F).

### Usage

```R

PCA(dir = <dir>, data = <data>, labels = <labels>, PC.retain = <PC.retain>)

```

### Example

```R

PCA(dir = "/path/to/vcf/", data = "snps_filter.vcf", labels = "snps_filter_labels.csv", PC.retain = FALSE)

```

### Dependencies

- `SNPRelate`
- `tidyverse`
- `psych`
- `tools`









## 5. GWAS: GAPIT

This R script runs a Genome-Wide Association Study (GWAS) analysis using the GAPIT3 R package. It saves the results of each trait in an individual folder. This function loads the required packages and data, then performs a GWAS analysis using the `GAPIT` function from the `GAPIT3` package. It loops through each trait in `trait_list`, creating a new folder for each trait in the working directory and saving the results of the GWAS analysis for that trait in that folder.

### Arguments

- `phenofile`: A database of genotypes/individuals (rows) and the trait's measurements or BLUPs (columns).
- `genofile`: Genotype data in hapmap format.
- `wdir`: A working directory where the function will create the folders for each trait.
- `trait_list`: A vector with the trait's name. The default is `NULL`, which will use the full phenotype dataset.

### Usage

```R

GAPIT3(phenofile, genofile, wdir, trait_list = NULL)

```

### Examples

```R

GAPIT3(phenofile = "my_phenotypes.csv", genofile = "my_genotypes.hmp", wdir = "my_results_folder", trait_list = c("Trait1", "Trait2", "Trait3"))

```

### Dependencies

- `devtools`
- `GAPIT3`






## 6. GWAS: EMMAX

XXX

### Arguments

- `XXX`: 

### Usage

```bash

conda 

```

### Examples

```bash

conda

```

### Dependencies

- `XXX`







## 6. Annotation of results

This code annotates the gen containing the significat SNPs from the `GAPIT3` results. Additional, it retrieves the closest genes downstream and upstream.

### Arguments

- `pat`: Location path of the files.
- `wdyw`: Gene feature to annotate.
- `Mdir`: Name of the directory that contains the `GAPIT3` results.
- `mod`: Names of the models to filter (options: )
- `gff3`: gff3 file from the genome version used for alignment.
- `annotationFile`: Annotation details of the genes. txt file from the genome version used for alignment.
- `version`: (Options: 6.1 or 8.1. Default = 6.1).

### Usage

```R

GWAS_Annotation(Mdir, pat, mod, wdyw, annot, GFF)

```

### Examples

```R

GWAS_Annotation(Mdir, pat, mod, wdyw, annot, GFF)

```

### Dependencies

- `XXX`









## 7. Boxplot: Genotype vs Phenotype

This R function generates a boxplot for a given SNP. The function takes as inputs a CSV file with the phenotype values and a list of SNPs. The function can also receive an optional CSV file to add extra information about the samples to the plot (categories, family, ect.) The output of the function is a PDF file with the plot. The function loads all the necessary packages and data files. It then checks for the presence of an optional file with sample labels and prepares the data for the boxplot. The function generates a boxplot for each SNP specified in the input CSV file. It uses `ggplot2` package to generate the plot, with the genotypes on the x-axis and the phenotype on the y-axis. The function adds a label to the x-axis to indicate the trait being plotted. If `labelfile` is provided, it adds extra information about the samples, coloring the data by the levels in the provided file.

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

### Usage

```R

GWAS_Boxplot(outputname, dir, phenofile, genofile, snp_list_file, labelfile = NULL)

``` 

### Example

```R

# Generate boxplot without extra labels
GWAS_Boxplot("outputname", ".path/to/save/plots/", "phenotype.csv", "genotype.hmp", "snp_list.csv")

# Generate boxplot with extra labels
GWAS_Boxplot("outputname", ".path/to/save/plots/", "phenotype.csv", "genotype.hmp", "snp_list.csv", "labelfile.csv")

```

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






## References 📚

- Browning, B. L., Tian, X., Zhou, Y., and Browning, S. R. 2021. Fast two-stage phasing of large-scale sequence data. American Journal of Human Genetics, 108(10):1880–1890. https://doi.org/10.1016/j.ajhg.2021.08.005
- Browning, B. L., Zhou, Y., and Browning, S. R. 2018. A one-penny imputed genome from next generation reference panels. American Journal of Human Genetics, 103(3):338–348. https://doi.org/10.1016/j.ajhg.2018.07.015
- Danecek, P., Auton, A., Abecasis, G., Albers, C. A., Banks, E., DePristo, M. A., Handsaker, R. E., Lunter, G., Marth, G. T., Sherry, S. T., McVean, G., Durbin, R., and 1000 Genomes Project Analysis Group. 2011. The variant call format and VCFtools. Bioinformatics, 27(15):2156–2158. https://doi.org/10.1093/bioinformatics/btr330
- Danecek, P., Bonfield, J. K., Liddle, J., Marshall, J., Ohan, V., Pollard, M. O., Whitwham, A., Keane, T., McCarthy, S. A., Davies, R. M., and Li, H. 2021. Twelve years of SAMtools and BCFtools. GigaScience, 10(2):giab008. https://doi.org/10.1093/gigascience/giab008
- Goodstein, D. M., Shu, S., Howson, R., Neupane, R., Hayes, R. D., Fazo, J., Mitros, T., Dirks, W., Hellsten, U., Putnam, N., and Rokhsar, D. S. 2012. Phytozome: a comparative platform for green plant genomics. Nucleic Acids Research, 40(D1):D1178–D1186.
- Huang, M., Liu, X., Zhou, Y., Summers, R. M., and Zhang, Z. 2019. BLINK: A package for the next level of genome-wide association studies with both individuals and markers in the millions. GigaScience. https://doi.org/10.1093/gigascience/giy154
- Jombart, T. 2008. adegenet: a R package for the multivariate analysis of genetic markers. Bioinformatics, 24(11):1403–1405. https://doi.org/10.1093/bioinformatics/btn129
- Liu, X., Huang, M., Fan, B., Buckler, E. S., and Zhang, Z. 2016. Iterative usage of fixed and random effect models for powerful and efficient genome-wide association studies. PLoS Genetics, 12:e1005767. https://doi.org/10.1371/journal.pgen.1005767
- R Core Team. 2021. R: A language and environment for statistical computing. R Foundation for Statistical Computing, Vienna, Austria. https://www.R-project.org/
- Wang, J., and Zhang, Z. 2021. GAPIT Version 3: Boosting power and accuracy for genomic association and prediction. Genomics, Proteomics & Bioinformatics. https://doi.org/10.1016/j.gpb.2021.08.005
- Zhang, C., Dong, S. S., Xu, J. Y., He, W. M., and Yang, T. L. 2019. PopLDdecay: a fast and effective tool for linkage disequilibrium decay analysis based on variant call format files. Bioinformatics, 35(10):1786–1788. https://doi.org/10.1093/bioinformatics/bty875
- Zheng, X., Levine, D., Shen, J., Gogarten, S. M., Laurie, C., and Weir, B. S. 2012. A high-performance computing toolset for relatedness and principal component analysis of SNP data. Bioinformatics, 28(24):3326–3328. https://doi.org/10.1093/bioinformatics/bts606