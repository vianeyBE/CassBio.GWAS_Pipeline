# Cassava_Bioinformatics_Platform.GWAS_Pipeline
Pipeline for Genome-wide association studies in Cassava

Autor: Vianey Barrera-Enriquez and Camilo E. Sanchez

The pipeline has four main steps:

1. Quality control (*in progress*)
2. GWAS Analysis using GAPIT
3. Annotation of results
4. Boxplot of significant markers genotypes vs. phenotype  

## 2. GWAS: GAPIT


### GAPIT3 R script documentation

This R script runs a Genome-Wide Association Study (GWAS) analysis using the GAPIT3 R package. It saves the results of each trait in an individual folder.

### Usage

```R
GAPIT3(phenofile, genofile, wdir, trait_list = NULL)`
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

This will perform a GWAS analysis using the phenotype data in my_phenotypes.csv and the genotype data in my_genotypes.hmp. It will create a new folder for each trait in the trait_list vector in the my_results_folder directory.

### Output
The function will create a folder for each trait in the trait_list vector in the working directory. In each folder, it will save the results of the GWAS analysis for that trait. For more information about GAPIT outputs, please check the official GAPIT documentation. 

### Dependencies
- devtools
- GAPIT3

## License
This script is licensed under the MIT License. See the LICENSE file for details.

## Contact
For questions or feedback about this script, please contact Vianey Barrera-Enriquez at vbarrera@cgiar.org.





