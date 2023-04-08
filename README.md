# Cassava Bioinformatics Platform: GWAS_Pipeline
Pipeline for Genome-wide association studies in Cassava

**Autor**: Vianey Barrera-Enriquez and Camilo E. Sanchez

The pipeline has four main steps:

1. Quality control (*in progress*)
2. GWAS Analysis using GAPIT
3. Annotation of results (*in progress*)
4. Boxplot of significant markers genotypes vs. phenotype  

## 2. GWAS: GAPIT

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
- devtools
- GAPIT3

## 3. Annotation

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

## 3. Boxplot: Genotype vs Phenotype

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

# generate boxplot without extra labels
GWAS_Boxplot("outputname", ".path/to/save/plots/", "phenotype.csv", "genotype.hmp", "snp_list.csv")

# generate boxplot with extra labels
GWAS_Boxplot("outputname", ".path/to/save/plots/", "phenotype.csv", "genotype.hmp", "snp_list.csv", "labelfile.csv")

```

### Output
A single PDF file containing the boxplot of the SNPs. 

### Dependencies
- tidyverse
- tibble
- dplyr
- janitor
- ggplot2
- Biostrings
- hrbrthemes
- forcats
- ggsignif
- RColorBrewer

---

## Contact
For questions or feedback about this script, please contact Vianey Barrera-Enriquez at vbarrera@cgiar.org.
