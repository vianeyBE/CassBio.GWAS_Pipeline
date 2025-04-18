# Cassava Bioinformatics Platform: GWAS Pipeline 🧬

A modular, scalable, and reproducible pipeline for Genome-Wide Association Studies (GWAS) in Manihot esculenta (cassava). This end-to-end workflow spans from genotype quality control and population structure analysis to GWAS modeling and functional annotation.

This pipeline is designed to streamline complex GWAS analyses while ensuring flexibility and scalability for high-throughput plant breeding projects. Each module is independent and well-documented, enabling easy customization or integration into broader genomic workflows.

**Highlights** ✨:
- Modular design: Run individual components independently
- Built for plant breeding: Tailored filters and annotation
- Scalable: Snakemake-driven parallel execution
- Compatible: GAPIT, EMMAX, TASSEL, PLINK, and more






## Authors 🙋

For questions or feedback about this pipeline, please contact:

- **Camilo E. Sánchez-Sarria** - [CGIAR]: c.e.sanchez@cgiar.org
- **Vianey Barrera-Enriquez** - [CGIAR]: v.barrera@cgiar.org





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

To clone the repository:

```bash

# Clone the repository
git clone https://github.com/vianeyBE/cassava-gwas-pipeline.git
cd cassava-gwas-pipeline/

# OPTIONAL: Set up a conda environment
conda create -n cassava_gwas_env snakemake plink vcftools bcftools perl -c bioconda -c conda-forge
conda activate cassava_gwas_env

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

This module applies quality control filters to genotype data in `.vcf` format using `VCFtools`. It provides a simple and customizable `Bash` script that supports different filtering strategies tailored for different pipelines. You can activate the desired strategy by commenting/uncommenting the relevant lines in the script

This step ensures that only high-quality and informative SNPs are retained for downstream analysis

### Input arguments

- **`input_vcf`**: Path to the input `.vcf` file (can be compressed with `.gz`)
- **`output`**: Prefix for the output file(s) generated after filtering. The output will be a `.recode.vcf` file

### Example usage

Structure to execute the command:

``` sh

bash 01_GWAS_QC.sh <input_vcf> <output>

```

This command will apply the currently active filter set (GATK-style by default) to the `raw_variants.vcf.gz` file and produce a filtered VCF named `filtered_variants.recode.vcf`

``` sh

bash 01_QC.sh raw_variants.vcf.gz filtered_variants

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

### Example usage

Structure of to execture the command

``` sh

bash 02_LD_Decay.sh <dirIn> <vcfFile> <dirOut> <prefix> <dist> <mode>

```

This will split `mydata.vcf` into 18 chromosomes, compute LD decay up to `10 kb` for each chromosome, and generate corresponding `.stat.gz` files and LD decay plots in the specified output directory.

```sh

bash 02_LD_Decay.sh /path/to/vcf mydata.vcf /path/to/output cassava_ld 10000 by_chr

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

The pruning parameters are set via hard-coded values in the `Snakefile`:

- **`path`**: Path to the directory where is located the input `.vcf` file 
- **`file`**: Name of the input VCF file
- **`window`**: The size (in kilobases) of the sliding window used by PLINK to calculate LD between SNPs
- **`step`**: The number of SNPs to shift the window forward at each iteration
- **`r²`**: The LD threshold above which one of two highly correlated SNPs is removed
- **`r²`**: A mandatory flag by `Snakemake` in which you define until which step of the pipeline you want to reach

Tool-specific paths and thread count are controlled via `config/config.yaml`:

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

First, please configure initial parameters and the final output in the `Snakefile`. For example as:

``` yaml

# Define objects
path = "/path/to/vcf"
file = "mydata.vcf.gz"

window = 10000
step = 10
r2 = 0.25

# Include modular rules for different stages of the workflow
include: "rules/common.smk"
include: "rules/plink.smk"

# Define the final output of the workflow
rule all:
    input:
        outdir = directory(f"{path}/LD/{window}_{step}_{r2}_{get_basename(file)}/04_vcf")

```

Ensure the paths and filenames are correctly configured in the `Snakefile`. If everything is correctly configured, this command will display the planned workflow without executing it:

``` sh

snakemake --np

```

Then execute the workflow with:

``` sh

snakemake --cores 50 

```

After executing the previous commands. You will filter `mydata.vcf.gz` using `PLINK` using a window size of `'10000`, step of `10` and r² trehsold of 0.25 and you will obtain a `.vcf` as an output

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

### Description

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

### Description

This module performs **Genome-Wide Association Studies (GWAS)** using the `GAPIT3` `R package`. It supports multiple GWAS models such as `MLMM`, `BLINK`, or `FarmCPU`, and allows batch analysis of multiple traits.

The script reads phenotypic and genotypic data, prepares the input format for `GAPIT`, and executes **GWAS** independently for each trait listed. Each trait's results are organized in dedicated output folders within the working directory. The function also provides detailed logging and quality control summaries (e.g. number of markers and samples).

This implementation is ideal for high-throughput association studies across diverse traits and environments.

### Input arguments

- **`phenofile`**: A `.csv` file containing phenotype data. Rows represent individuals, columns represent traits. The first column must be individual/genotype IDs
- **`genofile`**: Genotypic data in HapMap format (`.hmp.txt`)
- **`dir`**: Working directory where output folders for each trait will be created
- **`models`**: A character vector specifying one or more GWAS models to use (e.g., `MLMM`, `Blink`, `FarmCPU`).
- **`trait_list`** (optional): A vector of trait names to analyze. If `NULL`, all traits in the phenotype file will be included.

### Example usage

This command will perform GWAS for the three traits using the specified models, and output results into separate folders inside `my_gwas_output/`

```R
GAPIT3(
  phenofile = "my_phenotypes.csv",
  genofile = "my_genotypes.hmp.txt",
  dir = "my_gwas_output/",
  models = c("MLMM", "Blink", "FarmCPU"),
  trait_list = c("PlantHeight", "Yield", "DiseaseScore")
)

```

### Example output structure

Each folder contains model-specific outputs including **Manhattan plots**, **QQ plots**, and result tables for significant associations

``` markdown

05_gwas_gapit/
├── PlantHeight/
│   ├── GAPIT.MLMM.Plot.*
│   ├── GAPIT.MLMM.Result.csv
│   └── ...
├── Yield/
│   ├── GAPIT.FarmCPU.Result.csv
│   └── ...
├── DiseaseScore/
│   ├── GAPIT.Blink.Plot.*
│   └── ...
└── GAPIT_run.log

```

### Dependencies

This module requires the following `R` packages:

- **`GAPIT3`**: Core **GWAS** models and visualization.
- **`devtools`**: To install `GAPIT3` from GitHub.

To install them:

``` R

install.packages("devtools")
devtools::install_github("jiabowang/GAPIT3", force = TRUE)

```

Make sure you have `GAPIT3` installed and loaded before running this module





## 6. GWAS: EMMAX 🧮

### Description

This module performs **Genome-Wide Association Studies (GWAS)** using the **EMMAX** software, a high-performance tool for mixed-model association analysis that accounts for population structure and relatedness. Implemented as a **Snakemake workflow**, this module automates all steps from:
- Preparing input genotype and phenotype data
- Generating a kinship matrix
- Running EMMAX for each trait
- Calculating **PVE (percent variance explained)** for significant SNPs
- Producing **Manhattan and QQ plots**
- Saving trait-specific results in structured folders

It supports **multiple traits** and applies statistical corrections (Bonferroni or FDR) to identify significant SNPs

### Input arguments

- **`geno`**: A genotype file (`.vcf`)
- **`pheno`**: Phenotype table with columns: `Taxa` and `trait(s)`
- **`prefix`**: Base name for outputs
- **`method`**: Significance threshold correction method (`Bonferroni` or `FDR`)
- **`alpha`**: Threshold significance level (p.e., alpha = 0.05)

### Example usage

This will launch the workflow, process all traits in the phenotype file, and save EMMAX results and plots in the specified working directory

```bash

snakemake --cores 16 --resources mem_mb=5150

```

### Example output structure

Each trait is processed separately, and visual outputs are saved under `results/`

``` markdown

08_gwas_emmax/
├── kinship_matrix.kinf                         # Kinship matrix
├── trait1/
│   ├── results/
│   │   ├── EMMAX_results_trait1.csv           # Final filtered SNP results with PVE
│   │   ├── EMMAX_Manhattan_trait1.pdf         # Manhattan plot
│   │   └── EMMAX_QQ_trait1.pdf                # QQ plot
│   ├── trait1.ps                               # Association stats from EMMAX
│   ├── trait1.reml                             # REML variance estimates
│   └── freq.frq                                # Minor allele frequencies
├── trait2/
│   └── ...
└── logs/

```

### Dependencies

This module depends on the following tools:

- **`PLINK`**: For genotype preprocessing and conversion from VCF to binary format (`.bed` - `.bim` - `.fam`) compatible with EMMAX
- **`EMMAX`**: Efficient Mixed-Model Association tool for GWAS
- **`emmax-kin`**: Companion tool to compute kinship matrices required by EMMAX
- **`R`**: For post-GWAS statistical analysis and high-quality visualization of Manhattan and QQ plots
- **`Snakemake`**: For workflow orchestration and scalable execution across traits

Ensure the following are available and properly configured:

- The `plink` binary (v1.9 or above) is accessible in your system path or specified in `config.yaml`
- EMMAX executables: `emmax` and `emmax-kin` are correctly pointed to in the config
- `R (≥ 4.0)` with the following packages installed:

``` R 

install.packages("tidyverse")

```

- `plot_gwas_results.R` script is present in the working directory or path specified in the workflow
- A valid `config.yaml` file with paths, filenames, models, and thresholds
- `Snakemake` version ≥ 6.0 installed and available in your environment

To install Snakemake (if needed):

``` bash 

conda install -c bioconda snakemake

```





## 6. Functional annotation 🧾

### Description

This module performs **functional annotation of significant SNPs** identified through GWAS using the GAPIT3 module. It annotates the genes that directly overlap with each SNP, as well as the closest **upstream and downstream genes** within a defined genomic window (±10 kb)

The function supports **recursive search across trait-specific folders**, filters SNPs based on significance threshold (**Bonferroni correction**), and annotates results using genome-specific GFF3 and functional annotation tables. It is compatible with both **cassava genome v6.1 and v8.1**

The final output includes SNP location, nearby gene information, gene ontology (GO), annotation summary, and effect sizes from the **GWAS**

### Input arguments

- **`version`**: Reference genome version used for alignment and annotation. Options: `6.1` or `8.1` (default = `6.1`)
- **`recursive`**: Logical. If `TRUE`, the function searches subdirectories recursively to locate **GWAS** result files (e.g., from multiple traits)
- **`Wdir`**: Path to the working directory containing **GAPIT3** results (e.g., trait-specific folders or a summary `.csv` file)
- **`name`**: If `recursive = TRUE`: Name pattern of the summary result files to look for. If `recursive = FALSE`: File name of a merged summary result file (e.g., `Summary_GWAS_Results.csv`)
- **`mod`**: A character vector specifying which GWAS models to include (e.g., `MLMM`, `Blink`, `FarmCPU`)
- **`wdyw`**: Feature type to annotate from the GFF3 file. Options: `gene`, `CDS`, `mRNA`, `five_prime_UTR`, `three_prime_UTR`

### Example usage

This command will annotate significant SNPs (after **Bonferroni correction**), fetch nearby genes, and generate a clean annotation report saved in the working directory.

``` R

GWAS_Annotation(
  version = "8.1",
  recursive = FALSE,
  Wdir = "/path/to/GAPIT_results/",
  name = "Summary_GWAS_Results.csv",
  mod = c("BLINK", "FarmCPU", "MLMM"),
  wdyw = "gene"
)

```

### Example output structure

``` markdown

06_annotation/
├── GWAS_Annotation.csv         # Final annotated SNP-gene associations
├── GAPIT_trait1/               # GAPIT3 output folders (input)
├── GAPIT_trait2/               # ...
├── Mesculenta_v6.1.gene.gff3   # GFF3 file (used internally)
├── Mesculenta_v6.1.annotation_info.txt
└── (log and temp files)

```

### Dependencies

This module uses `base` and `tidyverse` `R` packages:

- **`tidyverse`**: Data wrangling and filtering
- **`Base`**: For string and file management

To install the dependencies:

``` R

install.packages("tidyverse")

```

Ensure your annotation files (**GFF3** and **annotation table**) are formatted for the specified genome version and are available in the working directory or updated within the function paths





## 7. Boxplot: Genotype vs Phenotype 📈

### Description

This module generates **publication-ready boxplots** to visualize the effect of individual SNPs on phenotypic traits. Each boxplot compares phenotypic values (y-axis) across **genotype groups (x-axis)** for SNPs of interest. The plots help validate and interpret **marker-trait associations** identified through GWAS

The function accepts:
- A **phenotype file** (`.csv` file)
- A **genotype file** (`HapMap` format)
- A **list of SNP–trait pairs** to plot
- Optionally, a **label file** to color samples by categories (e.g. family, treatment, population group)

It supports **recursive search** for multiple traits/models, automatically handles IUPAC allele encoding, and outputs:

- A **multi-page PDF** with all boxplots
- An **Excel summary** of the plotted data

### Input arguments
- **`prefix`**: Output prefix. Used for naming the PDF and Excel files
- **`dir`**: Directory containing all input files and where outputs will be saved
- **`phenofile`**: 
  `.csv` file with phenotype data:
    1. Sample IDs (Taxa)
    2. Other columns: Traits.
- **`genofile`**: Genotype data in HapMap format (`.hmp.txt`)
- **`snpList`**:
  Either:
  - A `.csv` file with 3 columns:
    1. `SNPs`: SNP names as in genotype file
    2. `trait`: Trait name in phenotype file
    3. `model`: GWAS model name
  OR a **search pattern** if `recursive = TRUE`
- **`recursive`**: Logical. If `TRUE`, the function searches subdirectories for SNP-trait result files
- `**labelfile`** (optional):
  A `.csv` file with 2 columns:
    1. `Taxa`: Sample names
    2. `label`: Grouping variable (e.g., population, treatment)
- **`order`** (optional): Character vector specifying the desired plot order of levels in `labelfile`

### Example usage

```R

# Without group labels
GWAS_Boxplot(
  prefix = "ACWP_F2_Results",
  dir = "/path/to/data/",
  phenofile = "phenotypes.csv",
  genofile = "genotypes.hmp.txt",
  snpList = "significant_snps.csv",
  recursive = FALSE
)

# With group labels and custom order
GWAS_Boxplot(
  prefix = "ACWP_F2_Results",
  dir = "/path/to/data/",
  phenofile = "phenotypes.csv",
  genofile = "genotypes.hmp.txt",
  snpList = "significant_snps.csv",
  recursive = FALSE,
  labelfile = "sample_labels.csv",
  order = c("S", "IS", "I", "IR", "R")
)


``` 

### Example output structure

``` markdown

07_boxplots/
├── ACWP_F2_Results.SNPs_boxplot.pdf     # PDF file with one boxplot per SNP
├── ACWP_F2_Results.trait_x_snp.xlsx     # Excel summary of genotype–phenotype values
├── phenotypes.csv                       # Input phenotype data
├── genotypes.hmp.txt                    # Input genotype data (hapmap format)
├── significant_snps.csv                # List of SNPs to plot
└── sample_labels.csv                   # Optional: sample groupings

```

### Dependencies

This module requires several R packages:
- **`tidyverse`**: Data wrangling and visualization
- **`janitor`**: Data cleaning
- **`Biostrings`**: For SNP encoding
- **`hrbrthemes`**: Improved plot aesthetics
- **`ggsignif`**: Significance bars on plots
- **`RColorBrewer`**: Custom color palettes
- **`openxlsx`**: Excel export

To install the dependencies:

``` R

install.packages(c("tidyverse", "janitor", "hrbrthemes", "ggsignif", "RColorBrewer", "openxlsx"))
if (!requireNamespace("BiocManager", quietly = TRUE)) install.packages("BiocManager")
BiocManager::install("Biostrings")

```





## Citing this pipeline 📌

If you use this pipeline in your research or publication, please cite it as:

> Authors: Sánchez-Sarria, C. E., and Barrera-Enríquez, V  
> Title: *GWAS modular pipeline for Cassava: A reproducible and scalable workflow for pre-processing steps, GWAS and annotation in Cassava (Manihot esculenta)*
> GitHub repository: https://github.com/vianeyBE/CassBio.GWAS_Pipeline  
> Version: v1.0  
> Year: 2025

Alternatively, cite this repository using the following BibTeX entry:

``` bibtex

@misc{cassava_gwas_pipeline,
  author       = {Sánchez-Sarria, C. E., and Barrera-Enríquez, V},
  title        = {GWAS modular pipeline for Cassava},
  year         = 2025,
  version      = {v1.0},
  url          = {https://github.com/vianeyBE/CassBio.GWAS_Pipeline},
  note         = {GWAS modular pipeline for Cassava}
}

```





## License 📄

This pipeline is released under the MIT License (see LICENSE file)





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