# Short name: GWAS analysis.
# Description: Run GWAS analysis using GAPIT3 R package.
# Output: It saves the results of each trait in a individual folder.
#
# Author: Vianey Barrera-Enriquez (vpbarrerae@gmail.com)
#
# Arguments:
# phenofile: Measurements or BLUPs of the traits. Row are individuals and columns are traits/variables
# genofile: Genotype data on hapmap format
# models: Vector with the model or models to evaluate
# dir: Working directory. Here it will create the folders for each trait
# trait_list: Vector with the trait's name to do GWAS. Default NULL



##### To do #####
# Everything is good!



# 0: Function init -------------------------------------------------------------
GAPIT3 <- function(phenofile, genofile, dir, models, trait_list = NULL){
  
  
  
  # 1: Load packages and data --------------------------------------------------
  # Install packages if needed
  if (!require(devtools)) install.packages(devtools)
  devtools::install_github("jiabowang/GAPIT3", force = T)
  
  # Load packages
  library(devtools)
  library(GAPIT)
  
  # Load data
  myY2 <- read.csv(paste0(dir, phenofile))
  myG <- read.delim(paste0(dir, genofile), head = F)
  taxacol <- names(myY2)[1] <- "Taxa"
  
  # Determines if use complete dataset or not
  if (is.null(trait_list)){
    
    # Informative message
    message("Using full phenotype dataset")
    
    # 
    myY <- myY2
    trait_list <- names(myY)[-1]
    
  } else {
    
    # Informative message
    message("Subsetting according provided list")
    
    # 
    myY <- myY2[, c(taxacol, trait_list)]
    
  }
  
  # Informative messages
  # About phenofile
  message("Number of traits in phenofile: ", dim(myY)[2] - 1)
  message("Number of samples in phenofile: ", dim(myY)[1])
  
  # About genofile
  message("Number of samples in genofile: ", dim(myG)[2] - 11)
  message("Number of SNPs in genofile: ", dim(myG)[1] - 1)
  
  
  
  # 2: GAPIT -------------------------------------------------------------------
  # GWAS
  for (trait in trait_list){
    
    # Creation and setting working directories
    setwd(dir)
    dir.create(trait)
    setwd(paste(dir, trait, sep = ""))
    
    # Informative message
    message("Running GWAS on trait: ", trait)
    
    # Actually run the GAPIT function to GWAS
    myGAPIT <- GAPIT(
      
      Y = myY[, c(taxacol, trait)], # Phenotypic data
      G = myG, # Genotypic data
      PCA.total = 10, # Principal components to retain
      model = models, # Model(s) to use
      Multiple_analysis = T, # Allow multiple analyses
      Geno.View.output = T, # Visualization of geno data
      Phenotype.View = T, # Visualization of pheno data
      Model.selection = T, # Select a model
      Random.model = T, # Do not include random effects
      kinship.cluster = "average", # Clustering method
      kinship.group = "Mean", # Type of grouping of indvs
      Inter.Plot = F # Interactive plots
      
    )
    
    # Informative message
    message("Finished GWAS on trait: ", trait)
    setwd(dir)
    
  }
  
  # 3: Function ends -----------------------------------------------------------
  message("Done!")
  
}



###### Example(s) ######
# Set arguments
# phenofile <- "Traits_noparents.csv"
# genofile <- "GATK_noparents.hmp.txt"
# models <- c("MLMM", "Blink", "FarmCPU")
# dir <- "/datas3/Cassava/Cassava_Analysis/Cassava_BGI_VCF/2024-07-17_group7/GWAS/GAPIT/"
# trait_list <- NULL

# Run function
# GAPIT3(phenofile, genofile, wdir, trait_list = NULL)
