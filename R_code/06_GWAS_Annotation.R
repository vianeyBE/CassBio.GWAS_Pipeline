# Short name: Annotation for GWAS results
# Description (1): Get the annotation of the gen containing the SNP using the GWAS GAPIT results
# Description (2): Retrieves the closest genes down and up stream
# Output: Basename of the output file
# 
# Authors: Camilo E. Sanchez (c.e.sanchez@cgiar.org) and Vianey Barrera-Enriquez (vpbarrera@gmail.com)
#
# Arguments:
# Ddir: Directory where are located annotation files (annot and gff files) (For example: home/user/folder).
# annot: Annotation details of the genes. txt file from the genome version used for alignment.
# gff: gff3 file from the genome version used for alignment.
# version: You can choose among the different reference genome versions (Options: 6.1 or 8.1).
# Wdir: Directory that contains the GAPIT results (Example: home/user/folder).
# name: The path of file names to look for. (Example: QTL_LOD_Intervals).
# mod: The model(s) of interest (Options: BLINK, GLM, MLM, FarmCPU).
# wdyw: What are you looking for to annotate (Options: CDS, five_prime_UTR, gene, mRNA, three_prime_UTR).



###### To do ######
# Everything is good!



# 0: Function init -------------------------------------------------------------

GWAS_Annotation <- function(Ddir, annot, gff, version, Wdir, name, mod, wdyw){
  
  
  
  # 1: Load all the directories, info, and data --------------------------------
  
  # Load packages
  if (!require(tidyverse)) install.packages(tidyverse)
  
  library(tidyverse)
  
  # Load annotation data and re-organize it
  
  if (version == 6.1){
    
    message("Loading required files to do the annotation")
    
    # Set working directory where files are located
    setwd(Ddir)
    
    annot <- read.delim(annot, header = F) %>%
      rename(ID = 1, Locus = 2, Trans = 3, Peptide = 4, GO = 10, AT.name = 12, AT.define = 13) %>%
      select(ID, Locus, Trans, Peptide, GO, AT.name, AT.define)
    
    gff <- read.delim(gff, header = F, comment.char = "#") %>%
      rename(Chr = V1, What = V3, Start = V4, End = V5) %>%
      tidyr::separate(col = V9, into = c("ID", "na"), sep = ";", extra = "drop") %>%
      tidyr::separate(col = ID, into = c("na2", "na3"), sep = "=", extra = "drop") %>%
      tidyr::separate(col = na3, into = c("Name", "na4"), sep = ".v", extra = "drop") %>%
      select(Chr, What, Start, End, Name) %>%
      dplyr::filter(grepl("Chromosome", Chr)) %>%
      mutate(Chr = recode(Chr, Chromosome01 = 1, Chromosome02 = 2, Chromosome03 = 3, Chromosome04 = 4,
                          Chromosome05 = 5, Chromosome06 =  6, Chromosome07 = 7, Chromosome08 = 8, 
                          Chromosome09 = 9, Chromosome10 = 10, Chromosome11 = 11, Chromosome12 = 12,
                          Chromosome13 = 13, Chromosome14 = 14, Chromosome15 = 15, Chromosome16 = 16, 
                          Chromosome17 = 17, Chromosome18 = 18)) %>%
      dplyr::filter(What %in% wdyw)
    
  } else {
    
    message("Loading required files to do the annotation")
    
    # Set working directory where files are located
    setwd(Ddir)
    
    annot <- read.delim(annot, header = T) %>%
      rename(ID = 1, Locus = 2, Trans = 3, Peptide = 4, GO = 10, AT.name = 11, AT.define = 12) %>%
      select(ID, Locus, Trans, Peptide, GO, AT.name, AT.define)
    
    gff <- read.delim(gff, header = F, comment.char = "#") %>%
      rename(Chr = V1, What = V3, Start = V4, End = V5) %>%
      tidyr::separate(col = V9, into = c("ID", "na"), sep = ";", extra = "drop") %>%
      tidyr::separate(col = ID, into = c("na2", "na3"), sep = "=", extra = "drop") %>%
      tidyr::separate(col = na3, into = c("Name", "na4"), sep = ".v", extra = "drop") %>%
      select(Chr, What, Start, End, Name) %>%
      dplyr::filter(grepl("Chromosome", Chr)) %>%
      mutate(Chr = recode(Chr, Chromosome01 = 1, Chromosome02 = 2, Chromosome03 = 3, Chromosome04 = 4,
                          Chromosome05 = 5, Chromosome06 =  6, Chromosome07 = 7, Chromosome08 = 8, 
                          Chromosome09 = 9, Chromosome10 = 10, Chromosome11 = 11, Chromosome12 = 12,
                          Chromosome13 = 13, Chromosome14 = 14, Chromosome15 = 15, Chromosome16 = 16, 
                          Chromosome17 = 17, Chromosome18 = 18)) %>%
      dplyr::filter(What %in% wdyw)
    
  }
  
  
  
  # 2: Find all the CSVs with the results and filter them ----------------------
  
  # Set working directory
  setwd(Wdir)
  
  # Get the names of the files
  message("Getting list of CSV files...")
  names <- list.files(path = Wdir, pattern = paste0(name, "."), all.files = F, full.names = F, recursive = T)
  
  # Informative messages
  message(paste("GWAS files found:", length(names)))
  message("Reading GWAS files...")
  
  # Create an empty list and the variable to iterate
  # Read, rename, modify and save all the csv files in the list
  csv_L <- list()
  
  # Loop to find GWAS files
  for (i in 1:length(names)){
    
    # Get the name of the trait and the model
    tr <- str_split_1(paste0(names[i]), "/")
    int <- str_split_1(tr[2], paste0(name, "."))
    mdl <- str_split_1(int[2], paste0(".", tr[1]))
    
    # Database handling
    csv_L[[i]] <- read.csv(paste0(Wdir, names[i])) %>%
      mutate(Trait = tr[1], Model = mdl[1], Start = Pos - 10000, End = Pos + 10000) %>%
      select(SNP, Chr, Pos, Start, End, P.value, MAF, Effect, Model, Trait) %>%
      filter(P.value <= (0.05/length(SNP)))
    
    # Progress bar
    cat('\r', i, ' files processed |', rep('=', i / 2), ifelse(i == length(names), '|\n', '>'), sep = '')
    
  }
  
  # Merge the data frames of the list in a single data frame and filter by model type
  GWAS <- bind_rows(csv_L)
  GWAS <- filter(GWAS, Model %in% mod)
  
  message(paste("There are", dim(GWAS)[1], "SNPs after filtering"))
  
  
  
  # Conditional: If there is at least one result of the GWAS models, the annotation continues
  if (dim(GWAS)[1] > 0){
    
    
    
    # 3: Match the results with the gene annotation database -------------------
    
    message("Retriving annotation...")
    
    # Creates an empty list to store the data
    gff_L <- list()
    
    # Filters and merges the gff3 data
    for (i in 1:nrow(GWAS)){
      
      # Filter the gff3 file per chromosome
      data <- dplyr::filter(gff, GWAS$Chr[i] == Chr)
      data.der <- dplyr::filter(data, GWAS$Pos[i] <= Start)
      data.izq <- dplyr::filter(data, GWAS$Pos[i] >= End)
      
      # Combine the different filtered data
      gff_L[[i]] <- rbind(
        
        # Genes in a 10.000 window 
        dplyr::filter(data, GWAS$Pos[i] >= Start & GWAS$Pos[i] <= End) %>%
          mutate(SNP = GWAS$SNP[i], SNP.Pos = GWAS$Pos[i], 
                 SNP.Start = GWAS$Start[i], SNP.End = GWAS$End[i], 
                 SNP.Location = "Inside", Distance = "---",
                 Pvalue = GWAS$P.value[i], MAF = GWAS$MAF[i], Effect = GWAS$Effect[i],
                 Model = GWAS$Model[i], Trait = GWAS$Trait[i]),
        
        # Looking genes down stream 
        data.der[which.min(abs(GWAS$Pos[i] - data.der$Start)),] %>%
          mutate(SNP = GWAS$SNP[i], SNP.Pos = GWAS$Pos[i], 
                 SNP.Start = GWAS$Start[i], SNP.End = GWAS$End[i], 
                 SNP.Location = "Right", Distance = min(abs(GWAS$Pos[i] - data.der$Start)),
                 Pvalue = GWAS$P.value[i], MAF = GWAS$MAF[i], Effect = GWAS$Effect[i],
                 Model = GWAS$Model[i], Trait = GWAS$Trait[i]),
        
        # Looking genes up stream
        data.izq[which.min(abs(GWAS$Pos[i] - data.izq$End)),] %>%
          mutate(SNP = GWAS$SNP[i], SNP.Pos = GWAS$Pos[i],
                 SNP.Start = GWAS$Start[i], SNP.End = GWAS$End[i],
                 SNP.Location = "Left", Distance = min(abs(GWAS$Pos[i] - data.izq$End)),
                 Pvalue = GWAS$P.value[i], MAF = GWAS$MAF[i], Effect = GWAS$Effect[i],
                 Model = GWAS$Model[i], Trait = GWAS$Trait[i])
        
      )
      
      # Convert the values of the column into characters to avoid issues
      gff_L[[i]]$Distance <- as.character(gff_L[[i]]$Distance)
      
      # Progress bar
      cat('\r', i, ' files processed |', rep('=', i / 3), ifelse(i == nrow(GWAS), '|\n', '>'), sep = '')
      
    }
    
    # Merge the data frames of the list in a single data frame
    gff_LM <- bind_rows(gff_L)
    
    message(paste("There are", dim(gff_LM)[1], "genes to annotate"))
    
    
    
    # 4: Formatting dataframe --------------------------------------------------
    
    message("Obtaining gene information...")
    
    # Creates an empty list to store the data
    annot_L <- list()
    
    # Filters and merges the annot data
    for (i in 1:nrow(gff_LM)){
      
      # Filter the annotation data and creates new columns based on previous data
      annot_L[[i]] <- filter(annot, 
                             Locus == gff_LM$Name[i] | 
                               Trans == gff_LM$Name[i] | 
                               Peptide == gff_LM$Name[i]) %>%
        mutate(Chr = gff_LM$Chr[i], Effect = gff_LM$Effect[i], Pvalue = gff_LM$Pvalue[i],
               MAF = gff_LM$MAF[i], Model = gff_LM$Model[i], Trait = gff_LM$Trait[i],
               Gen.Start = gff_LM$Start[i], Gen.End = gff_LM$End[i], Distance = gff_LM$Distance[i],
               SNP = gff_LM$SNP[i], SNP.Pos = gff_LM$SNP.Pos[i], SNP.Location = gff_LM$SNP.Location[i],
               SNP.Start = gff_LM$SNP.Start[i], SNP.End = gff_LM$SNP.End[i])
      
      # Progress bar
      cat('\r', i, ' files processed |', rep('=', i / 4), ifelse(i == nrow(gff_LM), '|\n', '>'), sep = '')
      
    }
    
    # Merge the data frames of the list in a single data frame and modify it
    SNP_annotation <- bind_rows(annot_L) %>%
      select(SNP, Chr, Model, Trait, Pvalue, Effect, MAF, Locus, Gen.Start, Gen.End, AT.define,
             SNP.Pos, SNP.Location, Distance, SNP.Start, SNP.End)
    
    
    
    # 5: Save output -----------------------------------------------------------
    
    message("Saving output file: 'GWAS_Annotation.csv' in working directory")
    
    write.csv(SNP_annotation, file = "GWAS_Annotation.csv", quote = F, row.names = F)
    
    message("Done!")
    
  } else {
    
    message("No SNPs to annotate")
    
  }
  
}



# Example(s) -------------------------------------------------------------------
# Set arguments
# Ddir <- "D:/OneDrive - CGIAR/Cassava_Bioinformatics_Team/00_Data/"
# annot <- "Mesculenta_305_v6.1/Mesculenta_305_v6.1.annotation_info.txt"
# gff <- "Mesculenta_305_v6.1/Mesculenta_305_v6.1.gene.gff3"
# version <- "6.1"
# Wdir <- "D:/OneDrive - CGIAR/Cassava_Bioinformatics_Team/03_GWAS_PPD_Populations/04_GWAS/GAPIT_Results/"
# name <- "GAPIT.Association.GWAS_Results"
# mod <- c("BLINK", "FarmCPU", "MLM")
# wdyw <- "gene"



# Run function -----------------------------------------------------------------
# GWAS_Annotation(Ddir, annot, gff, version, Wdir, name, mod, wdyw)
