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
# wdir: Working directory. Here it will create the folders for each trait
# trait_list: Vector with the trait's name to do GWAS. Default NULL



##### To do #####
# Everything is good!



# 0: Function init -------------------------------------------------------------
GAPIT3 <- function(phenofile, genofile, wdir, trait_list = NULL){
  
  
  
  # 1: Load packages and data --------------------------------------------------
  if (!require(devtools)) install.packages(devtools)
  devtools::install_github("jiabowang/GAPIT3", force = T)
  
  library(devtools)
  library(GAPIT)

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
    myY <- myY2[, c(taxacol, selection)]
    
  }
  
  # Informative messages
  message("Number of traits in phenofile: ", dim(myY)[2] - 1)
  message("Number of samples in phenofile: ", dim(myY)[1])
  message("Number of samples in genofile: ", dim(myG)[2] - 11)
  
  
  
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
      
      Y = myY[, c(taxacol, trait)],
      G = myG,
      PCA.total = 10,
      model = models,
      Multiple_analysis = T,
      # Geno.View.output = F,
      # Phenotype.View = F,
      Model.selection = T,
      Random.model = F,
      # kinship.cluster = "average", 
      # kinship.group = "Mean",
      # Inter.Plot = F
      
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
# phenofile <- "cbsv_pheno_filter.csv"
# genofile <- "10000_50_0.8_03_group6_filter.hmp.txt"
# models <- c("Blink", "FarmCPU")
# dir <- "D:/OneDrive - CGIAR/00_CassavaBioinformaticsPlatform/04_CBSD_Group6/08_LD_Test/"
# trait_list <- NULL

# Run function
# GAPIT3(phenofile, genofile, wdir, trait_list = NULL)
