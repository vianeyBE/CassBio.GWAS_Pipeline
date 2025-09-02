# Short name: Annotation for GWAS results
# Description (1): Get the annotation of the SNP using the GWAS GAPIT results
# Description (2): Retrieves the closest bp down and up stream
# Output: Basename of the output file
# 
# Authors: Camilo E. Sanchez (c.e.sanchez@cgiar.org)
# Authors: Vianey Barrera-Enriquez (vpbarrera@gmail.com)
#
# Arguments:
# version: Choose among the reference genome versions (Options: 6.1 or 8.1)
# recursive: If the functions look file by file or one file already merged (Options: TRUE/FALSE)
# pb: Distance in pair bases
# dir: Directory that contains the GAPIT results (Example: home/user/folder)
# file_name: The path of file names to look for. (Example: QTL_LOD_Intervals)
# mod: The model(s) of interest (Options: BLINK, GLM, MLM, FarmCPU, etc)
# wdyw: What to annotate (Options: CDS, five_prime_UTR, gene, mRNA, three_prime_UTR)



###### To do ######
# Everything is ok!



# 0: Function init -------------------------------------------------------------

GWAS_Annotation <- function(version, recursive, pb, dir, file_name, mod, wdyw){
  
  
  
  # 1: Load all the directories, info, and data --------------------------------
  
  # Load packages
  if (!require(tidyverse)) install.packages(tidyverse)
  
  library(tidyverse)
  
  # Load annotation data and re-organize it
  if (version == 6.1){
    
    # Informative message
    message("Loading required files to do the annotation")
    
    # Set objects
    Ddir <- 
      "D:/OneDrive - CGIAR/00_BioInf_Platform/00_Basics/03_Server_Bioinfo-CO/00_ReferenceGenome/Mesculenta_v6.1/"
    annot <- "Mesculenta_305_v6.1.annotation_info.txt"
    gff <- "Mesculenta_305_v6.1.gene.gff3"
    
    # Read and modify the files
    annot <- read.delim(paste0(Ddir, annot), header = F) %>%
      rename(ID = 1, Name = 2, Trans = 3, Peptide = 4, GO = 10, AT.name = 12, Annotation = 13) %>%
      select(ID, Name, GO, Annotation)
    
    gff <- read.delim(paste0(Ddir, gff), header = F, comment.char = "#") %>%
      rename(Chr = V1, What = V3, Start = V4, End = V5, Strand = V7) %>%
      tidyr::separate(col = V9, into = c("ID", "na"), sep = ";", extra = "drop") %>%
      tidyr::separate(col = ID, into = c("na2", "na3"), sep = "=", extra = "drop") %>%
      tidyr::separate(col = na3, into = c("Name", "na4"), sep = ".v", extra = "drop") %>%
      select(Chr, What, Start, End, Name, Strand) %>%
      dplyr::filter(grepl("Chromosome", Chr)) %>%
      mutate(Chr = recode(Chr, Chromosome01 = 1, Chromosome02 = 2, Chromosome03 = 3, Chromosome04 = 4,
                          Chromosome05 = 5, Chromosome06 =  6, Chromosome07 = 7, Chromosome08 = 8, 
                          Chromosome09 = 9, Chromosome10 = 10, Chromosome11 = 11, Chromosome12 = 12,
                          Chromosome13 = 13, Chromosome14 = 14, Chromosome15 = 15, Chromosome16 = 16, 
                          Chromosome17 = 17, Chromosome18 = 18)) %>%
      dplyr::filter(What %in% wdyw)
    
    # Informative message
    message("Annotation files loaded succesfully")
    
  } else {
    
    # Informative message
    message("Loading required files to do the annotation")
    
    # Set objects
    Ddir <- 
      "D:/OneDrive - CGIAR/00_BioInf_Platform/00_Basics/03_Server_Bioinfo-CO/00_ReferenceGenome/Mesculenta_v8.1/"
    annot <- "Mesculenta_671_v8.1.annotation_info.txt"
    gff <- "Mesculenta_671_v8.1.gene.gff3"
    
    # Read and modify the files
    annot <- read.delim(paste0(Ddir, annot), header = T) %>%
      rename(ID = 1, Name = 2, GO = 10, Annotation = 12) %>%
      select(Name, GO, Annotation)
    
    gff <- read.delim(paste0(Ddir, gff), header = F, comment.char = "#") %>%
      rename(Chr = V1, What = V3, Start = V4, End = V5, Strand = V7) %>%
      tidyr::separate(col = V9, into = c("ID", "na"), sep = ";", extra = "drop") %>%
      tidyr::separate(col = ID, into = c("na2", "na3"), sep = "=", extra = "drop") %>%
      tidyr::separate(col = na3, into = c("Name", "na4"), sep = ".v", extra = "drop") %>%
      select(Chr, What, Start, End, Name, Strand) %>%
      dplyr::filter(grepl("Chromosome", Chr)) %>%
      mutate(Chr = recode(Chr, Chromosome01 = 1, Chromosome02 = 2, Chromosome03 = 3, Chromosome04 = 4,
                          Chromosome05 = 5, Chromosome06 =  6, Chromosome07 = 7, Chromosome08 = 8, 
                          Chromosome09 = 9, Chromosome10 = 10, Chromosome11 = 11, Chromosome12 = 12,
                          Chromosome13 = 13, Chromosome14 = 14, Chromosome15 = 15, Chromosome16 = 16, 
                          Chromosome17 = 17, Chromosome18 = 18)) %>%
      dplyr::filter(What %in% wdyw)
    
    # Informative message
    message("Annotation files loaded succesfully")
    
  }
  
  
  
  # 2: Read GWAS results and filter them ---------------------------------------
  
  # Set working directory
  setwd(dir)
  
  # Read GWAS results files
  if (recursive == T){
    
    # Get the names of the files
    message("Getting list of CSV files...")
    names <- list.files(path = dir, pattern = paste0(file_name, "."), 
                        all.files = F, full.names = F, recursive = T)
    
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
      int <- str_split_1(tr[2], paste0(file_name, "."))
      mdl <- str_split_1(int[2], paste0(".", tr[1]))
      
      # Database handling
      csv_L[[i]] <- read.csv(paste0(dir, names[i])) %>%
        mutate(Trait = tr[1], Model = mdl[1], Start = Pos - pb, End = Pos + pb) %>%
        select(SNP, Chr, Pos, Start, End, P.value, MAF, Effect, PVE, 
               Model, Trait) %>%
        filter(P.value <= (0.05/length(SNP)))
      
      # Progress bar
      cat('\r', i, ' files processed |', rep('=', i / 2), 
          ifelse(i == length(names), '|\n', '>'), sep = '')
      
    }
    
    # Merge dataframes of the list in a single data frame and filter by model type
    GWAS <- bind_rows(csv_L)
    
    message(paste("There are", dim(GWAS)[1], "SNPs after filtering"))
    
  } else {
    
    GWAS <- read.csv(paste0(dir, file_name)) %>%
      select(Trait, Model, SNP, Chr, Pos, PValue, MAF, Effect, PVE)
    
  }
  
  
  
  # Conditional: If there is at least one SNP, the annotation continues
  if (dim(GWAS)[1] > 0){
    
    
    
    # 3: Match the results with the gene annotation database -------------------
    # Informative message
    message("Retriving annotation...")
    
    # Creates an empty list to store the data
    gff_L <- list()
    
    # Filters and merges the gff3 data
    for (i in 1:nrow(GWAS)){
      
      # SNP info
      snp_chr <- GWAS$Chr[i]
      snp_pos <- GWAS$Pos[i]
      
      # Filter all the data by chr and position less than pb
      gff_L[[i]] <- dplyr::filter(gff, Chr == snp_chr) %>%
        mutate(Dist_SNP_to_Start = abs(Start - snp_pos),
               Dist_SNP_to_End = abs(End - snp_pos)) %>%
        filter(Dist_SNP_to_Start <= pb | Dist_SNP_to_End <= pb) %>%
        filter(!duplicated(Name)) %>%
        mutate(Trait = GWAS$Trait[i], Model = GWAS$Model[i], SNP = GWAS$SNP[i],
               SNP_Chr = GWAS$Chr[i], SNP_Pos = GWAS$Pos[i],  
               Pvalue = GWAS$PValue[i], MAF = GWAS$MAF[i],
               Effect = GWAS$Effect[i], PVE = GWAS$PVE[i]) %>%
        select(Trait, Model, SNP, SNP_Chr, SNP_Pos, Pvalue, MAF, Effect, PVE, 
               Name, Start, End, Strand, Dist_SNP_to_Start, Dist_SNP_to_End) %>%
        mutate(Location = case_when(
          SNP_Pos > Start & SNP_Pos < End ~ "Inside",
          SNP_Pos < Start & SNP_Pos < End ~ "Upstream",
          SNP_Pos > Start & SNP_Pos > End ~ "Downstream",
          TRUE ~ NA_character_))
      
      # Progress bar
      cat('\r', i, ' SNPs processed |', rep('=', i / 3), 
          ifelse(i == nrow(GWAS), '|\n', '>'), sep = '')
      
    }
    
    # Merge the data frames of the list in a single data frame
    gff_LM <- bind_rows(gff_L)
    
    message(paste("There are", dim(gff_LM)[1], "genes to annotate"))
    
    
    
    # 4: Formatting dataframe --------------------------------------------------
    # Informative message
    message("Obtaining gene information...")
    
    # Left join between genetic data and annotation 
    SNP_annotation <- left_join(gff_LM, annot, by = "Name",
                                relationship = "many-to-many") %>%
      mutate(ID = paste0(Trait, Model, SNP, PVE, Name, Annotation)) %>%
      #filter(Annotation != "") %>%
      filter(!duplicated(ID)) %>%
      select(-ID)
    
    # Informative message
    message(paste("There are", dim(SNP_annotation)[1], "genes annotated"))
    
    
    
    # 5: Save output -----------------------------------------------------------
    # Informative message
    message("Saving output file: ", 
            paste0(gsub(".csv", "", file_name)),
            "_GWAS_Annotation.csv in working directory")
    
    # Actually saves it!
    write.csv(SNP_annotation, 
              file = paste0(gsub(".csv", "", file_name), "_GWAS_Annotation.csv"), 
              row.names = F)
    
    # Informative
    message("Done!")
    
  } else {
    
    message("No SNPs to annotate")
    
  }
  
}



# Example(s) -------------------------------------------------------------------
# Set arguments
 version <- "8.1"
 recursive <- "F"
 pb <- 10000
 dir <- "D:/OneDrive - CGIAR/00_BioInf_Platform/01_ACWP/03_F2/02_F2_Pheno/12_MaleSterility/Annotation/"
 file_name <- "GWAS_results_MS.csv"
 wdyw <- "gene"



# Run function -----------------------------------------------------------------
# GWAS_Annotation(version, recursive, pb, dir, file_name, mod, wdyw)
