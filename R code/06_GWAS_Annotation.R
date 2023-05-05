# Short name: Annotation for GWAS results
# Description (1): Get the annotation of the gen containing the SNP using the GWAS GAPIT results
# Description (2): Retrieves the closest genes down and up stream
# Output: Basename of the output file
# 
# Authors: Camilo E. Sanchez (c.e.sanchez@cgiar.org) and Vianey Barrera-Enriquez (vpbarrera@gmail.com)
#
# Arguments:
# Wdir: Name of the directory that contains the GAPIT results. For example: home/user/folder.
# Ddir: Directory where is located the annotation files (annot, GFF files).
# pat: Enter the path of file names to look for. For example: GAPIT.Association.GWAS_Results.
# mod: Enter the model(s) of interest (Options: BLINK, GLM, MLM, FarmCPU).
# wdyw: Enter what are you looking for to annotate (Options: CDS, five_prime_UTR, gene, mRNA, three_prime_UTR).
# annot: Annotation details of the genes. txt file from the genome version used for alignment.
# GFF: gff3 file from the genome version used for alignment.
# version: (Optional) You can choose between the genome of reference version 6.1 or 8.1 (Options: 6.1 or 8.1. Default = 6.1).



###### To do ######
#



# 0: Function init -------------------------------------------------------------

GWAS_Annotation <- function(Wdir, pat, mod, wdyw, annot, GFF, version){
  
  
  
  # 1: Load all the directories, info, and data --------------------------------
  
  # Load packages
  if (!require(tidyverse)) install.packages(tidyverse)
  library(tidyverse)
  
  # Load annotation data and re-organize it
  
  message("Loading required files to do the annotation")
  
  if (version == 6.1) {
    
    # Set working directory where files are located
    setwd(Ddir)
    
    annot <- read.delim(annot, header = F) %>%
      rename(ID = 1, Locus = 2, Trans = 3, Peptide = 4, GO = 10, AT.name = 12, AT.define = 13) %>%
      select(ID, Locus, Trans, Peptide, GO, AT.name, AT.define)
    
    GFF <- read.delim(GFF, header = F, comment.char = "#") %>%
      rename(Chr = V1, What = V3, Start = V4, End = V5) %>%
      tidyr::separate(col = V9, into = c("ID", "na"), sep = ";") %>%
      tidyr::separate(col = ID, into = c("na2", "na3"), sep = "=") %>%
      tidyr::separate(col = na3, into = c("Name", "na4"), sep = ".v") %>%
      select(Chr, What, Start, End, Name) %>%
      mutate(Chr = recode(Chr, Chromosome01 = 1, Chromosome02 = 2, Chromosome03 = 3, Chromosome04 = 4,
                          Chromosome05 = 5, Chromosome06 =  6, Chromosome07 = 7, Chromosome08 = 8, 
                          Chromosome09 = 9, Chromosome10 = 10, Chromosome11 = 11, Chromosome12 = 12,
                          Chromosome13 = 13, Chromosome14 = 14, Chromosome15 = 15, Chromosome16 = 16, 
                          Chromosome17 = 17, Chromosome18 = 18)) %>%
      na.omit() %>% 
      dplyr::filter(What %in% wdyw)
    
  } else {
    
    # Set working directory where files are located
    setwd(Ddir)
    
    annot <- read.delim(annot, header = T) %>%
      rename(ID = 1, Locus = 2, Trans = 3, Peptide = 4, GO = 10, AT.name = 11, AT.define = 12) %>%
      select(ID, Locus, Trans, Peptide, GO, AT.name, AT.define)
    
    GFF <- read.delim(GFF, header = F, comment.char = "#") %>%
      rename(Chr = V1, What = V3, Start = V4, End = V5) %>%
      tidyr::separate(col = V9, into = c("ID", "na"), sep = ";") %>%
      tidyr::separate(col = ID, into = c("na2", "na3"), sep = "=") %>%
      tidyr::separate(col = na3, into = c("Name", "na4"), sep = ".v") %>%
      select(Chr, What, Start, End, Name) %>%
      mutate(Chr = recode(Chr, Chromosome01 = 1, Chromosome02 = 2, Chromosome03 = 3, Chromosome04 = 4,
                          Chromosome05 = 5, Chromosome06 =  6, Chromosome07 = 7, Chromosome08 = 8, 
                          Chromosome09 = 9, Chromosome10 = 10, Chromosome11 = 11, Chromosome12 = 12,
                          Chromosome13 = 13, Chromosome14 = 14, Chromosome15 = 15, Chromosome16 = 16, 
                          Chromosome17 = 17, Chromosome18 = 18)) %>%
      na.omit() %>% 
      dplyr::filter(What %in% wdyw)
    
  }
  
  
  
  # 2: Find all the CSVs with the results and filter them ----------------------
  
  message("Getting list of CSV files...")
  
  # Get the names of the files
  setwd(Wdir)
  names <- list.files(path = Wdir, pattern = paste0(pat, "."), all.files = F, full.names = F, recursive = T)
  
  message("Reading GWAS files...")
  message(paste("GWAS files found:", length(names)))
  
  # Create an empty list and the variable to iterate
  # Read, rename, modify and save all the csv files in the list
  csv.l <- list()
  p <- 0
  
  # Loop to find GWAS files
  for (i in 1:length(names)){
    
    # Iterator
    p <- p + 1
    
    # Database handling
    csv.l[[i]] <- read.csv(paste0(Wdir, "/", names[p])) %>%
      mutate(traits = paste0("/", names[p]), Start = Pos - 10000, End = Pos + 10000) %>%
      tidyr::separate(col = traits, into = c("na", "trait"), sep = paste0(pat, ".")) %>%
      tidyr::separate(col = trait, into = c("Model", "Trait", "na2"), sep = "\\.") %>%
      dplyr::rename(PValue = P.value, Nobs = nobs) %>%
      select(SNP, Chr, Pos, Start, End, PValue, MAF, Nobs, Effect, Model, Trait) %>%
      filter(PValue <= (0.05/length(SNP)))
    
    # Progress bar
    cat('\r', i, ' files processed |', rep('=', i / 2), ifelse(i == length(names), '|\n', '>'), sep = '')
    
  }
  
  # Merge the data frames of the list in a single data frame and filter by model type
  GWAS <- bind_rows(csv.l)
  GWAS <- filter(GWAS, Model %in% mod)
  
  message(paste("There are", dim(GWAS)[1], "SNPs after filtering"))
  
  
  
  # 3: Match the results with the gene annotation database ---------------------
  
  # Creates a conditional in which if there is at least one result of GWAS models, the annotation continues
  if (dim(GWAS)[1] > 0){
    
    message("Retriving annotation...")
    
    # Creates the objects to perform the loop
    p <- 0
    gensL <- list()
    
    # Loop to annotate each results
    for (i in 1:nrow(GWAS)){
      
      # Iterator
      p <- p + 1
      
      # Creates the databases to filter the data inside, upstream, and downstream
      data <- dplyr::filter(GFF, GWAS$Chr[p] == Chr)
      data.der <- dplyr::filter(data, GWAS$Pos[p] <= Start)
      data.izq <- dplyr::filter(data, GWAS$Pos[p] >= End)
      
      # Filter and merge the created databases
      gensL[[i]] <- rbind(
        
        # Genes in a 10.000 window 
        dplyr::filter(data, GWAS$Pos[p] >= Start & GWAS$Pos[p] <= End) %>%
          mutate(SNP = GWAS$SNP[p],
                 SNP.Pos = GWAS$Pos[p],
                 SNP.Start = GWAS$Start[p],
                 SNP.End = GWAS$End[p],
                 SNP.Location = "Inside",
                 Distance = "---",
                 PValue = GWAS$PValue[p],
                 MAF = GWAS$MAF[p],
                 Effect = GWAS$Effect[p],
                 Model = GWAS$Model[p],
                 Trait = GWAS$Trait[p]),
        
        # Looking genes downstream 
        data.der[which.min(abs(GWAS$Pos[p] - data.der$Start)),] %>%
          mutate(SNP = GWAS$SNP[p],
                 SNP.Pos = GWAS$Pos[p],
                 SNP.Start = GWAS$Start[p],
                 SNP.End = GWAS$End[p],
                 SNP.Location = "Downstream",
                 Distance = min(abs(GWAS$Pos[p] - data.der$Start)),
                 PValue = GWAS$PValue[p],
                 MAF = GWAS$MAF[p],
                 Effect = GWAS$Effect[p],
                 Model = GWAS$Model[p],
                 Trait = GWAS$Trait[p]),
        
        # Looking genes upstream
        data.izq[which.min(abs(GWAS$Pos[p] - data.izq$End)),] %>%
          mutate(SNP = GWAS$SNP[p],
                 SNP.Pos = GWAS$Pos[p],
                 SNP.Start = GWAS$Start[p],
                 SNP.End = GWAS$End[p],
                 SNP.Location = "Upstream",
                 Distance = min(abs(GWAS$Pos[p] - data.izq$End)),
                 PValue = GWAS$PValue[p],
                 MAF = GWAS$MAF[p],
                 Effect = GWAS$Effect[p],
                 Model = GWAS$Model[p],
                 Trait = GWAS$Trait[p])
        
      )
      
      # Convert the values of the column into characters to avoid issues
      gensL[[i]]$Distance <- as.character(gensL[[i]]$Distance)
      
      # Progress bar
      cat('\r', i, ' files processed |', rep('=', i / 2), ifelse(i == nrow(GWAS), '|\n', '>'), sep = '')
      
    }
    
    # Merge the data frames of the list in a single data frame
    gensLD <- bind_rows(gensL)
    
    message(paste("There are", dim(gensLD)[1], "genes to annotate"))
    
    
    
    # 4: Formatting dataframe --------------------------------------------------
    
    message("Obtaining gene information...")
    
    # Creates the objects to perform the loop
    p <- 0
    gensLCc <- list()
    
    # Loop
    for (i in 1:nrow(gensLD)){
      
      # Iterator
      p <- p + 1
      
      # 
      gensLCc[[i]] <- filter(annot, 
                              annot$Locus == gensLD$Name[p] |
                               annot$Trans == gensLD$Name[p] |
                                annot$Peptide == gensLD$Name[p]) %>%
        mutate(Chr = gensLD$Chr[p],
               Effect = gensLD$Effect[p],
               Model = gensLD$Model[p],
               Trait = gensLD$Trait[p],
               PValue = gensLD$PValue[p],
               MAF = gensLD$MAF[p],
               Gen.Start = gensLD$Start[p],
               Gen.End = gensLD$End[p],
               Distance = gensLD$Distance[p],
               SNP = gensLD$SNP[p],
               SNP.Pos = gensLD$SNP.Pos[p],
               SNP.Location = gensLD$SNP.Location[p],
               SNP.Start = gensLD$SNP.Start[p],
               SNP.End = gensLD$SNP.End[p])
      
      # Progress bar
      cat('\r', i, ' files processed |', rep('=', i / 10), ifelse(i == nrow(gensLD), '|\n', '>'), sep = '')
      
    }
    
    # Merge the data frames of the list in a single data frame and modify it
    SNP_annotation <- bind_rows(gensLCc) %>%
      select(SNP, Chr, Model, PValue, MAF, Trait, Locus, Gen.Start, Gen.End, GO, AT.name, AT.define,
             Effect, SNP.Pos, SNP.Start, SNP.End, SNP.Location, Distance)
    
    
    
    # 5: Save output -----------------------------------------------------------
    
    message("Saving output file: 'GWAS_Annotation.csv' in working directory")
    
    write.csv(SNP_annotation, file = "GWAS_Annotation.csv", quote = F, row.names = F)
    
    message("Done! ")
    
  } else {
    
    message("No SNPs to annotate")
    
  }
  
}



###### Example(s) ######
# Set arguments
# Wdir <- "D:/OneDrive - CGIAR/Cassava_Bioinformatics_Team/03_GWAS_PPD_Populations/04_GWAS/GAPIT_Results/"
# Ddir <- "D:/OneDrive - CGIAR/Cassava_Bioinformatics_Team/00_Data/"
# pat <- "GAPIT.Association.GWAS_Results"
# mod <- c("BLINK", "FarmCPU", "MLM")
# wdyw <- "gene"
# annot <- "Mesculenta_305_v6.1.annotation_info.txt"
# GFF <- "Mesculenta_305_v6.1/Mesculenta_305_v6.1.gene.gff3"
# version <- "6.1"

# Run function
# GWAS_Annotation(Wdir, pat, mod, wdyw)
