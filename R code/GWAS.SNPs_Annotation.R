# Author: Camilo E. Sanchez (c.e.sanchez@cgiar.org) Vianey Barrera-Enriquez (vpbarrera@gmail.com)
# Get the annotation of the gen containing the SNP using the GWAS GAPIT results
# Additionally, it retrieve the closest  gene down and up stream.
# It requires: 
# Mdir: the directory that contains the GAPIT results separated by traits
# output: basename of the output file
# gff3: gff3 annotation format. It contains gene ID, transcript ID, GO ID, star, end
# annotationFile: text anotation file containing additional description of the genes

### ADD p value % MAF from GWAS RESULTS #######
# Configure the requirements
pat <- as.character()
wdyw <- as.character()
Mdir <- as.character()
mod <- as.character()

# Configure them manually on code
pat <- "GAPIT.Association.GWAS_Results."
wdyw <- "gene"
Mdir <- "D:/OneDrive - CGIAR/Cassava_Bioinformatics_Team/03_GWAS_PPD_Populations/04_GWAS/GAPIT_Results"
mod <- c("BLINK", "FarmCPU", "MLM")

# Configure them using a prompt
if (rlang::is_empty(pat)) {
  pat <- readline(prompt =
  "Enter the path of file names to looking for (p.e., QTL_LOD_Intervals. [Must finish with a point (.)]): ")
}

if (rlang::is_empty(wdyw)) {
  wdyw <- readline(prompt =
  "Enter what are you looking for to anotate (options: CDS, five_prime_UTR, gene, mRNA, three_prime_UTR): ")
}

if (rlang::is_empty(Mdir)) {
  Mdir <- readline(prompt = "Enter the working directory (p.e., home/user/folder): ")
}

if (rlang::is_empty(mod)) {
  Mdir <- readline(prompt = "Enter the model(s) of interest (options: BLINK, GLM, MLM, FarmCPU): ")
}



# 1: Load all the directories, info, and data ----------------------------------

# Load packages
if (!require(tidyverse)) install.packages(tidyverse)
library(tidyverse)

# Load data and re-organize it
setwd("D:/OneDrive - CGIAR/Cassava_Bioinformatics_Team/00_Data/Mesculenta_305_v6.1")

message("Loading required files to do the annotation")

annot <- read.delim("Mesculenta_305_v6.1.annotation_info.txt", header = F) %>%
  rename(ID = V1, Gen1 = V2, Gen2 = V3, Gen3 = V4, PF = V5, PTHR = V6, KOG = V7, na = V8, K0 = V9,
         GO = V10, AT = V11, NPI = V12, name = V13) %>%
  select(ID, Gen1, Gen2, Gen3, GO, NPI, name)

GFF <- read.delim("Mesculenta_305_v6.1.gene.gff3", header = F, comment.char = "#") %>%
  rename(chr = V1, phyto = V2, what = V3, start = V4, end = V5, na = V6, sign = V7, NPI = V8, sep = V9) %>%
  tidyr::separate(col = sep, into = c("ID", "na"), sep = ";") %>%
  separate(col = ID, into = c("na2", "na3"), sep = "=") %>%
  separate(col = na3, into = c("name", "na4"), sep = ".v" ) %>%
  select(chr, what, start, end, name) %>%
  mutate(chr = recode(chr, Chromosome01 = 1, Chromosome02 = 2, Chromosome03 = 3, Chromosome04 = 4,
                      Chromosome05 = 5, Chromosome06 =  6, Chromosome07 = 7, Chromosome08 = 8, 
                      Chromosome09 = 9, Chromosome10 = 10, Chromosome11 = 11, Chromosome12 = 12,
                      Chromosome13 = 13, Chromosome14 = 14, Chromosome15 = 15, Chromosome16 = 16, 
                      Chromosome17 = 17, Chromosome18 = 18)) %>%
  na.omit() %>% 
  dplyr::filter(what %in% wdyw)



# 2: Find all the CSVs with the results and filter them ------------------------

message("Getting list of CSV files...")

# Get the names of the files
setwd(Mdir)
names <- list.files(path = Mdir, pattern = pat, all.files = F, full.names = F, recursive = T)

message(paste("GWAS files found:", length(names)))

# Create an empty list and the variable to iterate
# Read, rename, modify and save all the csv files in the list
message("Reading GWAS files...")

# Creates the objects to perform the loop
csv.l <- list()
p <- 0

# Loop
for (i in 1:length(names)){
  p <- p + 1
  
  name.F <- gsub(" / ", "/", paste(Mdir, "/", names[p]))
  name.T <- gsub(paste(Mdir), "", paste(name.F))
  
  csv.l[[i]] <- read.csv(name.F) %>%
    mutate(traits = name.T, start = Pos - 10000, end = Pos + 10000) %>%
    tidyr::separate(col = traits, into = c("na", "trait"), sep = paste(pat)) %>%
    tidyr::separate(col = trait, into = c("model", "trait", "na2"), sep = "\\.") %>%
    select(SNP, Chr, Pos, start, end, P.value, MAF, nobs, Effect, model, trait) %>%
    filter(P.value <= (0.05/length(SNP)))
  
  # Progress bar
  cat('\r', i, ' files processed |', rep('=', i / 2), ifelse(i == length(names), '|\n', '>'), sep = '')
}

# Merge the data frames of the list in a single data frame and filter by model type
GWAS <- bind_rows(csv.l)
GWAS <- filter(GWAS, model %in% mod)

message(paste("There are", dim(GWAS)[1], "SNPs after filtering" ))

# Delete unnecessary objects
rm(csv.l, i, p, name.F, name.T)



# 3: Match the results with the gene annotation database -----------------------

# Creates a conditional
# If there is at least one result of the GWAS models, the annotation is done
if (dim(GWAS)[1] > 0){
  
  message("Retriving Annotation...")
  
  # Creates the objects to perform the loop
  p <- 0
  gensL <- list()
  
  # Loop
  for (i in 1:nrow(GWAS)){
    p <- p + 1
    
    data <- dplyr::filter(GFF, GWAS$Chr[p] == chr)
    data.der <- dplyr::filter(data, GWAS$Pos[p] <= start)
    data.izq <- dplyr::filter(data, GWAS$Pos[p] >= end)
    
    gensL[[i]] <- rbind(
    # Genes in a 10.000 window 
    dplyr::filter(data, GWAS$Pos[p] >= start & GWAS$Pos[p] <= end) %>%
      mutate(SNP = GWAS$SNP[p],
             SNP.pos = GWAS$Pos[p],
             SNP.start = GWAS$start[p],
             SNP.end = GWAS$end[p],
             SNP.Location = "Inside",
             Distance = "---",
             p_value = GWAS$P.value[p],
             MAF = GWAS$MAF[p],
             effect = GWAS$Effect[p],
             model = GWAS$model[p],
             trait = GWAS$trait[p]),
    
    # Looking genes down stream 
    data.der[which.min(abs(GWAS$Pos[p] - data.der$start)),] %>%
      mutate(SNP = GWAS$SNP[p],
             SNP.pos = GWAS$Pos[p],
             SNP.start = GWAS$start[p],
             SNP.end = GWAS$end[p],
             SNP.Location = "Right",
             Distance = min(abs(GWAS$Pos[p] - data.der$start)),
             p_value = GWAS$P.value[p],
             MAF = GWAS$MAF[p],
             effect = GWAS$Effect[p],
             model = GWAS$model[p],
             trait = GWAS$trait[p]),
             
    # Looking genes up stream
    data.izq[which.min(abs(GWAS$Pos[p] - data.izq$end)),] %>%
      mutate(SNP = GWAS$SNP[p],
             SNP.pos = GWAS$Pos[p],
             SNP.start = GWAS$start[p],
             SNP.end = GWAS$end[p],
             SNP.Location = "Left",
             Distance = min(abs(GWAS$Pos[p] - data.izq$end)),
             p_value = GWAS$P.value[p],
             MAF = GWAS$MAF[p],
             effect = GWAS$Effect[p],
             model = GWAS$model[p],
             trait = GWAS$trait[p])
    )
    
    # Convert the values of the column into characters to avoid issues
    gensL[[i]]$Distance <- as.character(gensL[[i]]$Distance)
    
    # Progress bar
    cat('\r', i, ' files processed |', rep('=', i / 3), ifelse(i == nrow(GWAS), '|\n', '>'), sep = '')
    }
  
  # Merge the data frames of the list in a single data frame
  gensLD <- bind_rows(gensL)
  
  message(paste("There are", dim(gensLD)[1], "genes to annotate"))
  
  # Delete unnecessary objects
  rm(gensL, i, p, data, data.der, data.izq)
  
  
  
  # 4: Formatting dataframe ----------------------------------------------------
  
  message("Obtaining gene information...")
    
  # Creates the objects to perform the loop
  p <- 0
  gensLCc <- list()
  
  # Loop
  for (i in 1:nrow(gensLD)){
    p <- p + 1
    gensLCc[[i]] <- filter(annot, annot$Gen1 == gensLD$name[p] | annot$Gen2 == gensLD$name[p] |
                             annot$Gen3 == gensLD$name[p]) %>%
      mutate(chr = gensLD$chr[p],
             effect = gensLD$effect[p],
             model = gensLD$model[p],
             trait = gensLD$trait[p],
             p_value = gensLD$p_value[p],
             MAF = gensLD$MAF[p],
             gen_start = gensLD$start[p],
             gen_end = gensLD$end[p],
             distance = gensLD$Distance[p],
             SNP = gensLD$SNP[p],
             SNP_pos = gensLD$SNP.pos[p],
             SNP_location = gensLD$SNP.Location[p],
             SNP_start = gensLD$start[p],
             SNP_end = gensLD$SNP.end[p])
    
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
    
    message("Done! ")
    
  } else {
    message("No SNPs to annotate")
    }
