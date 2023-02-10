#
# Author: Vianey Barrera-Enriquez (vpbarrerae@gmail.com)
# Run GWAS analysis using GAPIT3 R package. It saves the results of each trait
# in a individual folder. 
# phenofile: measurements or BLUPs of the traits. Row are individuals and Columns are traits/variables
# genofile: genotype data on hapmap format
# wdir: working directory. Here it will create the folders for each trait
# trait_list: vector with the trait's name to do GWAS. Default NULL

GAPIT3 <- function(phenofile, genofile, wdir, trait_list=NULL){

  # 1: Load Packages and Data --------------------------------------------------
  library(GAPIT3)
  
  myY2  <- read.csv(phenofile) 
  myG <- read.delim(genofile, head = FALSE)
  
  taxacol <- names(myY2)[1] <- "Taxa"
  
  
  if (is.null(trait_list)){
    message("Using full phenotype dataset")
    myY <- myY2
    trait_list <- names(myY)[-1]
  }
  else {
    message("Subsetting according provided list")
    myY <- myY2[,c(taxacol, selection)]
  }
  
  message("Number of Trait: ", dim(myY)[2]-1)
  message("Number of Samples: ", dim(myY)[1]-1)
  
  # 2: GAPIT -------------------------------------------------------------------

  for ( trait in trait_list ){
    dir.create(wdir)
    setwd(wdir)
    dir.create(trait)
    setwd(paste(wdir, trait, sep = "" ))
    message("Running GWAS on trait: ", trait )
    
    myGAPIT <- GAPIT(
      Y = myY[, c(taxacol, trait)],
      G = myG,
      PCA.total=10,
      model = c("GLM","MLM","FarmCPU","Blink"),
      Multiple_analysis = TRUE,
      #Geno.View.output = FALSE,
      #Phenotype.View = FALSE,
      Model.selection = TRUE,
      Random.model = FALSE,
      #kinship.cluster = "average", 
      #kinship.group = "Mean",
      #Inter.Plot = FALSE
    )
    
    message("Finished GWAS on trait: ", trait )
    setwd(wdir)
  }
  
  message("Done!")
}






