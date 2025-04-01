# Cassava Bioinformatics Platform: GWAS Pipeline 🧬

A modular and reproducible pipeline for performing genome-wide association studies (GWAS) in Manihot esculenta (cassava), covering the complete analysis workflow from genotype quality control to the functional annotation of significant markers.

This pipeline is designed to streamline complex GWAS analyses while ensuring flexibility and scalability for high-throughput plant breeding projects. Each module is independent and well-documented, enabling easy customization or integration into broader genomic workflows.






## Authors 🙋

For questions or feedback about this pipeline, please contact:

- **Vianey Barrera-Enriquez** - [CGIAR]: v.barrera@cgiar.org
- **Camilo E. Sanchez** - [CGIAR]: c.e.sanchez@cgiar.org





## Workflow overview 🚀

**| 1 | Quality control** 🧹 Filters SNPs and samples based on quality metrics

**| 2 | Linkage desequilibirum (LD) decay analysis** 📉 Estimates and visualizes LD decay across the genome

**| 3 | Linkage desequilibirum (LD) pruning** ✂️ Removes SNPs in high LD to reduce redundancy

**| 4 | Population structure analysis** 🧭 Performs PCA to assess population stratification

**| 5 | GWAS - GAPIT3** 📊 Runs GWAS using models from the GAPIT3 package

**| 6 | GWAS - EMMAX** 🧮 Performs GWAS using the EMMAX mixed model

**| 7 | Functional annotation** 🧾 Annotates significant SNPs with gene-level information

**| 8 | Marker validation boxplots** 📈 Generates boxplots to visualize genotype–phenotype effects





## Module descriptions 🧩

Each module contains a README.md describing:

- **Description**: A brief explanation of what the module does
- **Input arguments**: Parameters or flags used when running the module
- **Example usage**: Command-line or script usage example
- **Example output structure**: Example of the outputs produced by the module
- **Dependencies**: Required packages, tools, or environments





## Installation and dependencies 🛠️

```bash

git clone https://github.com/vianeyBE/cassava-gwas-pipeline.git

```

Tools and softwares to install before use this pipeline:

- **`R (>= 4.1.0)`**: **`GAPIT`**, **`adegenet`**, **`vegan`**, **`SNPRelate`**, **`tidyverse`**, etc
- **`Bash`**
- **`Perl`**
- **`VCFtools`**
- **`BCFtools`**
- **`PopLDdecay`**
- **`Snakemake`**
- **`EMMAX`**





## 1. Quality control 🧹

### Description

This module applies quality control filters to genotype data in VCF format using `VCFtools`. It provides a simple and customizable Bash script that supports different filtering strategies tailored for different pipelines. You can activate the desired strategy by commenting/uncommenting the relevant lines in the script.

This step ensures that only high-quality and informative SNPs are retained for downstream analysis.

### Input arguments

- **`input_vcf`**: Path to the input VCF file (can be compressed with **`.gz`**)
- **`output`**: Prefix for the output file(s) generated after filtering. The output will be a **`.recode.vcf`** file

### Example usage

``` sh

bash 01_GWAS_QC.sh <input_vcf> <output>

```

### Example output structure

``` markdown

01_quality_control/
├── raw_variants.vcf.gz          # Input file
└── filtered_variants.recode.vcf # Output after filtering

```

### Dependencies

This module requires:

- **`VCFtools`**: A robust toolkit for manipulating and filtering VCF files.

Ensure `vcftools` is accessible in your environment path. You can install it using `conda`:

``` sh

conda install -c bioconda vcftools

```





## 2. LD decay 📉

### Description

This module performs Linkage Disequilibrium (LD) decay analysis using `PopLDdecay`, a fast and effective tool designed to calculate pairwise LD between SNPs and visualize the decay of LD with increasing genomic distance.

The script supports two flexible modes:
- **Whole-genome analysis**: Computes LD decay across the entire genome from a single VCF file.
- **Per-chromosome analysis**: Splits the input VCF by chromosome and performs LD decay analysis on each separately, enabling finer resolution and scalability.

All steps are logged and output files are compressed and plotted automatically.

### Input arguments

- **`dirIn`**: Directory containing the input VCF file.
- **`vcfFile`**: Name of the input VCF file.
- **`dirOut`**: Directory where output files will be saved.
- **`prefix`**: Prefix for output files.
- **`dist`**: Maximum distance (in base pairs) for pairwise LD calculation.
- **`mode`**: Analysis mode — either ''whole_genome'' or ''by_chr''.

### Example sage

```sh

bash 02_LD_Decay.sh <dirIn> <vcfFile> <dirOut> <prefix> <dist> <mode>

```

### Example output structure

``` markdown

02_ld_decay/
├── split_chr/                             # Temporary VCFs per chromosome
│   ├── cassava_ld.chr01.recode.vcf
│   └── ...
├── cassava_ld.chr01.PopLDdecay.MaxDist_10000.stat.gz
├── cassava_ld.chr01.MaxDist_10000.plot
├── cassava_ld.chr02.PopLDdecay.MaxDist_10000.stat.gz
├── cassava_ld.chr02.MaxDist_10000.plot
├── ...
└── ld_decay_analysis_cassava_ld.log       # Complete run log

```

### Dependencies

This module requires the following tools:

- **`PopLDdecay`**: For LD calculation and plotting
- **`VCFtools`**: Used in per-chromosome mode to subset the VCF by chromosome
- **`Perl`**: Required for running the bundled Plot_OnePop.pl plotting script from PopLDdecay

To install `PopLDdecay`:

```sh

git clone https://github.com/BGI-shenzhen/PopLDdecay.git
cd PopLDdecay
make

```





## 3. LD pruning ✂️

### Description

This module performs linkage disequilibrium (LD) pruning to reduce the number of highly correlated (non-informative) SNPs from the input genotype file (VCF). LD pruning enhances GWAS accuracy by reducing multicollinearity and improves time-efficiency in downstream modeling.

The workflow is implemented using Snakemake and is designed to operate chromosome-wise. It includes format conversion, pruning using PLINK’s `--indep-pairwise` method, and output conversion to `HapMap` format for `GWAS` compatibility.

### Input arguments

The pruning parameters are set via hard-coded values in the Snakefile:

- **`path`**: Path to the directory where is located the input VCF file 
- **`file`**: Name of the input VCF file
- **`window`**: The size (in kilobases) of the sliding window used by PLINK to calculate LD between SNPs
- **`step`**: The number of SNPs to shift the window forward at each iteration
- **`r²`**: The LD threshold above which one of two highly correlated SNPs is removed

You can configure until which step of the workflow you want to reach. For example:

``` yaml

rule all:
    input:
        outdir = directory(f"{path}/LD/{window}_{step}_{r2}_{get_basename(file)}/05_hapmap")

```

Tool-specific paths and thread count are controlled via **`config/config.yaml`**:

``` yaml

PLINK:
  path: /path/to/plink
  filtering:
    threads: <threads>

tassel:
  path: /path/to/tassel-5/run_pipeline.pl
  threads: <threads>

```

### Example usage

Ensure the paths and filenames are correctly configured in the **`Snakefile`**. If everything is correctly configurated, this command will tell you: 

``` sh

snakemake --np 

```

Then execute the workflow with:

``` sh

snakemake --cores <cores>

```

### Example output structure

``` markdown

LD/
└── 10000_1_0.25_GATK2.snps.vcf.gz/ # This name is given based on the Snakefile inputs
    ├── 01_plink/
    ├── 02_filtering/
    ├── 03_prune/
    ├── 04_vcf/
    └── 05_hapmap/

```

### Dependencies

This module depends on the following tools:

- **`PLINK`**: For LD pruning and format conversion
- **`TASSEL5`**: For VCF to `HapMap` conversion
- **`Snakemake`**: For workflow orchestration

Ensure the following are available and properly configured:

- `plink` binary in your environment
- `run_pipeline.pl` from `TASSEL`
- `config/config.yaml`** and `Snakemake` version ≥ 6.0





## 4. Population structure 🎯

This module performs **Principal Component Analysis (PCA)** to assess **population structure** in both phenotypic and genotypic datasets. It supports:
- **Phenotypic PCA** using standard `.csv` input, allowing exploratory analysis and variance decomposition
- **Genotypic PCA** using `.vcf` files, with conversion to `GDS` format via `SNPRelate`, suitable for large-scale SNP datasets
- Optional analysis for **retaining the most informative principal components**, using methods like `parallel analysis`, `Velicer’s MAP`, and `scree plots`
- Generation of **interactive or static plots** with or without group labels to visualize stratification.

This module helps detect population substructure and adjust for confounding effects in the downstream `GWAS`

### Input arguments

- **`dir`**: Directory where the input data file is located. Used as the working directory during execution.
- **`data`**: For **phenotypic data**: A `.csv` file where rows are genotypes and columns are traits. First column must be genotype IDs. For **genotypic data**: A `.vcf` file with SNP data.
- **`labels`** (optional): A `.csv` file with genotype IDs in the first column and group/class information in the second column (used for coloring plots).
- **`PC.retain`** (optional; default = `FALSE`): Boolean value indicating whether to perform an in-depth analysis to determine the optimal number of PCs to retain using multiple statistical methods.

### Example usage

This command will run **PCA** on a `.vcf` file, include sample labels for grouping in the plot, and skip the PC-retention diagnostics.

``` R

PCA(
  dir = "/path/to/vcf/",
  data = "data.vcf",
  labels = "data_labels.csv",
  PC.retain = FALSE
)

```

### Example output structure

``` markdown

04_population_structure/
├── snps_filter.gds                 # Converted GDS file from VCF
├── PCA_snps_filter.jpg             # PCA plot (static image)
├── GWAS_PCA.html                   # PCA plot (interactive, if phenotypic)
├── snps_filter.PCA_SNPRelate.csv   # Eigenvalues / explained variance
├── PCA_summary.txt                 # PC retention statistics (if enabled)
└── lib/                            # Support files for HTML visualization

```

### Dependencies

This module uses several R packages for data handling, analysis, and visualization:

- **`SNPRelate`**: Performs efficient PCA on large genotype datasets
- **`tidyverse`**: Data wrangling and plotting
- **`psych`**: Parallel analysis and factor retention
- **`tools`**: File and extension handling
- **`plotly`** (optional): Interactive plots for phenotypic PCA
- **`ggplot2`**: High-quality static plots

To install the dependencies:

``` R

install.packages(c("tidyverse", "psych", "tools", "plotly"))
if (!requireNamespace("BiocManager", quietly = TRUE)) install.packages("BiocManager")
BiocManager::install("SNPRelate")

```





## 5. GWAS: GAPIT 📊

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






## 6. GWAS: EMMAX 🧮

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







## 6. Functional annotation 🧾

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









## 7. Boxplot: Genotype vs Phenotype 📈

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