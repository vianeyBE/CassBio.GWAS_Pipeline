# Short name: Annotation for GWAS results
# Description (1): Get the annotation of the gen containing the SNP using the GWAS GAPIT results
# Description (2): Retrieves the closest genes down and up stream
# Output: Basename of the output file
# 
# Authors: Camilo E. Sanchez (c.e.sanchez@cgiar.org) and Vianey Barrera-Enriquez (vpbarrera@gmail.com)
#
# Arguments:
# Mdir: Name of the directory that contains the GAPIT results. For example: home/user/folder.
# Ddir: 
# version: 
# pat: Enter the path of file names to look for. For example: QTL_LOD_Intervals
# mod: Enter the model(s) of interest. Options: BLINK, GLM, MLM, FarmCPU
# wdyw: Enter what are you looking for to annotate. Options: CDS, five_prime_UTR, gene, mRNA, three_prime_UTR.



###### To do ######
# 1: Convert to function



# 0: Function init -------------------------------------------------------------

GWAS_Annotation <- function(Mdir, Ddir, version, pat, mod, wdyw){
  
  
  
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
  
  
}








# 3: Match the results with the gene annotation database -----------------------

# Creates a conditional
# If there is at least one result of the GWAS models, the annotation is done
if (dim(GWAS)[1] > 0){
  
  message("Retriving annotation...")
  
  # Creates an empty list to store the data
  gensL <- list()
  
  # Loop
  for (i in 1:nrow(GWAS)){
    
    # 
    data <- dplyr::filter(GFF, GWAS$Chr[i] == chr)
    data.der <- dplyr::filter(data, GWAS$Pos[i] <= start)
    data.izq <- dplyr::filter(data, GWAS$Pos[i] >= end)
    
    # 
    gensL[[i]] <- rbind(
      
    # Genes in a 10.000 window 
    dplyr::filter(data, GWAS$Pos[i] >= start & GWAS$Pos[i] <= end) %>%
      mutate(SNP = GWAS$SNP[i],
             SNP.pos = GWAS$Pos[i],
             SNP.start = GWAS$start[i],
             SNP.end = GWAS$end[i],
             SNP.Location = "Inside",
             Distance = "---",
             p_value = GWAS$P.value[i],
             MAF = GWAS$MAF[i],
             effect = GWAS$Effect[i],
             model = GWAS$model[i],
             trait = GWAS$trait[i]),
    
    # Looking genes down stream 
    data.der[which.min(abs(GWAS$Pos[i] - data.der$start)),] %>%
      mutate(SNP = GWAS$SNP[i],
             SNP.pos = GWAS$Pos[i],
             SNP.start = GWAS$start[i],
             SNP.end = GWAS$end[i],
             SNP.Location = "Right",
             Distance = min(abs(GWAS$Pos[i] - data.der$start)),
             p_value = GWAS$P.value[i],
             MAF = GWAS$MAF[i],
             effect = GWAS$Effect[i],
             model = GWAS$model[i],
             trait = GWAS$trait[i]),
             
    # Looking genes up stream
    data.izq[which.min(abs(GWAS$Pos[i] - data.izq$end)),] %>%
      mutate(SNP = GWAS$SNP[i],
             SNP.pos = GWAS$Pos[i],
             SNP.start = GWAS$start[i],
             SNP.end = GWAS$end[i],
             SNP.Location = "Left",
             Distance = min(abs(GWAS$Pos[i] - data.izq$end)),
             p_value = GWAS$P.value[i],
             MAF = GWAS$MAF[i],
             effect = GWAS$Effect[i],
             model = GWAS$model[i],
             trait = GWAS$trait[i])
    )
    
    # Convert the values of the column into characters to avoid issues
    gensL[[i]]$Distance <- as.character(gensL[[i]]$Distance)
    
    # Progress bar
    cat('\r', i, ' files processed |', rep('=', i / 3), ifelse(i == nrow(GWAS), '|\n', '>'), sep = '')
    
    }
  
  # Merge the data frames of the list in a single data frame
  gensLD <- bind_rows(gensL)
  
  message(paste("There are", dim(gensLD)[1], "genes to annotate"))
  
  
  
  # 4: Formatting dataframe ----------------------------------------------------
  
  message("Obtaining gene information...")
    
  # Creates an empty list to store the data
  gensLCc <- list()
  
  # Loop
  for (i in 1:nrow(gensLD)){
    
    # 
    gensLCc[[i]] <- filter(annot,
                            annot$Gen1 == gensLD$name[i] | 
                             annot$Gen2 == gensLD$name[i] | 
                              annot$Gen3 == gensLD$name[i]) %>%
      mutate(chr = gensLD$chr[i],
             effect = gensLD$effect[i],
             model = gensLD$model[i],
             trait = gensLD$trait[i],
             p_value = gensLD$p_value[i],
             MAF = gensLD$MAF[i],
             gen_start = gensLD$start[i],
             gen_end = gensLD$end[i],
             distance = gensLD$Distance[i],
             SNP = gensLD$SNP[i],
             SNP_pos = gensLD$SNP.pos[i],
             SNP_location = gensLD$SNP.Location[i],
             SNP_start = gensLD$start[i],
             SNP_end = gensLD$SNP.end[i])
    
    # Progress bar
    cat('\r', i, ' files processed |', rep('=', i / 4), ifelse(i == nrow(gensLD), '|\n', '>'), sep = '')
    
    }
    
    # Merge the data frames of the list in a single data frame and modify it
    SNP_annotation <- bind_rows(gensLCc) %>%
      select(SNP, chr, model, p_value, MAF, trait, Gen1, name, gen_start, gen_end, GO, NPI, effect,
             SNP_pos, SNP_start, SNP_end, SNP_location, distance) %>%
      rename(gen_name = Gen1, gen_name_extend = name)
    
    
    
    # 5: Save output -----------------------------------------------------------
    
    message("Saving output file: 'GWAS_Annotation.csv' in working directory")
    
    write.csv(SNP_annotation, file = "GWAS_Annotation.csv", quote = F, row.names = F)
    
    message("Done!")
    
  } else {
    
    message("No SNPs to annotate")
    
    }


###### Example(s) ######
# Set arguments
 Ddir <- "D:/OneDrive - CGIAR/Cassava_Bioinformatics_Team/00_Data/"
 Wdir <- "D:/OneDrive - CGIAR/Cassava_Bioinformatics_Team/03_GWAS_PPD_Populations/04_GWAS/GAPIT_Results/"
 version <- "6.1"
 name <- "GAPIT.Association.GWAS_Results"
 mod <- c("BLINK", "FarmCPU", "MLM")
 wdyw <- "gene"

# Run function
# GWAS_Annotation(Mdir, Ddir, version, pat, mod, wdyw)
