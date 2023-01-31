#
# Author: Camilo E. Sanches (c.e.sanchez@cgiar.org) Vianey Barrera-Enriquez (vpbarrera@gmail.com)
# Get the annotation of the gen containing the SNP using the GWAS GAPIT results
# Additionally, it retrieve the closest  gene down and up stream.
# It requires: 
# Mdir: the directory that contains the GAPIT results separated by traits
# output: basename of the output file
# gff3: gff3 annotation format. It contains gene ID, transcript ID, GO ID, star, end
# annotationFile: text anotation file containing additional description of the genes

annotation <- function(Mdir, output, gff3, annotationFile,
                       pat = "GAPIT.Association.GWAS_Results.", 
                       mod = c("BLINK", "FarmCPU", "MLM"),
                       gen = c("gene")){
  
  # 1: Load all the directories, info, and data ----------------------------------
  
  # Load packages
  library(dplyr)
  library(tidyverse)
  
  
  # Load data and re-organize it
  annot <- read.delim(annotationFile, header = F) %>%
    rename(ID = V1, Gen1 = V2, Gen2 = V3, Gen3 = V4, PF = V5, PTHR = V6, KOG = V7, na = V8, K0 = V9,
           GO = V10, AT = V11, NPI = V12, name = V13) %>%
    select(ID, Gen1, Gen2, Gen3, GO, NPI, name)
  
  GFF <- read.delim(gff3, header = F, comment.char = "#") %>%
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
    dplyr::filter(what %in% gen)
  
  # 2: Find all the CSVs with the results and filter them ------------------------
  message("Getting list of CSV files...")
  
  # Get the names of the files
  names <- list.files(path = Mdir, pattern = pat, all.files = F, full.names = F, recursive = T)
  message(paste("GWAS files found:", length(names)))
  
  # Create an empty list and the variable to iterate
  # Read, rename, modify and save all the csv files in the list 
  csv.l <- list()
  p <- 0
  
  message
  
  for (i in 1:length(names)){
    
    p <- p + 1
    name.F <- gsub("/ ", "/", paste(Mdir, names[p]))
    name.T <- gsub(paste(Mdir), "", paste(name.F))
    db <- read.csv(name.F) %>% mutate(traits = name.T, start = Pos - 10000, end = Pos + 10000) %>%
      tidyr::separate(col = traits, into = c("na", "trait"), sep = paste(pat)) %>%
      tidyr::separate(col = trait, into = c("model", "trait", "na2"), sep = "\\.") %>%
      select(SNP, Chr, Pos, start, end, P.value, MAF, nobs, H.B.P.Value, Effect, model, trait) %>%
      filter(P.value <= (0.05/length(SNP))) # p adj value = 0.001/no.markers
    csv.l[[i]] <- db
    
    # Progress bar
    cat('\r',
        i,
        ' files processed |',
        rep('=', i / 4),
        ifelse(i == length(names), '|\n',  '>'), sep = '')
  }
  
  # Merge the data frames of the list in a single data frame and delete the GLM results
  GWAS <- Reduce(function(...) merge(..., all = T), csv.l)
  GWAS <- filter(GWAS, model %in% mod)
  
  message(paste("There are", dim(GWAS)[1], "SNPs after filtering" ))
  
  # 3: Match the results with the gene annotation database -----------------------
  
  message("Retriving Annotation...")
  
  p <- 0
  gensL <- list()
  
  for (i in 1:nrow(GWAS)){
    p <- p + 1
    
    data <- dplyr::filter(GFF, GWAS$Chr[p] == chr)
    data.der <- dplyr::filter(data, GWAS$Pos[p] <= start)
    data.izq <- dplyr::filter(data, GWAS$Pos[p] >= end)
    
    gensL[[i]] <- rbind(
      # Genes in a 10.000 window 
      dplyr::filter(data, GWAS$Pos[p] >= start & GWAS$Pos[p] <= end) %>%
        mutate(SNP = GWAS$SNP[p], effect = GWAS$Effect[p], SNP.pos = GWAS$Pos[p],
               SNP.start = GWAS$start[p], SNP.end = GWAS$end[p], 
               SNP.Location = "Inside", Distance = "---",
               model = GWAS$model[p], trait = GWAS$trait[p]),
      # Looking genes down stream 
      data.der[which.min(abs(GWAS$Pos[p] - data.der$start)),] %>%
        mutate(SNP = GWAS$SNP[p], effect = GWAS$Effect[p], SNP.pos = GWAS$Pos[p],
               SNP.start = GWAS$start[p], SNP.end = GWAS$end[p],
               SNP.Location = "Right", Distance = min(abs(GWAS$Pos[p] - data.der$start)),
               model = GWAS$model[p], trait = GWAS$trait[p]),
      # Looking genes up stream
      data.izq[which.min(abs(GWAS$Pos[p] - data.izq$end)),] %>%
        mutate(SNP = GWAS$SNP[p], effect = GWAS$Effect[p], SNP.pos = GWAS$Pos[p],
               SNP.start = GWAS$start[p], SNP.end = GWAS$end[p],
               SNP.Location = "Left", Distance = min(abs(GWAS$Pos[p] - data.izq$end)),
               model = GWAS$model[p], trait = GWAS$trait[p])
    )
  }
  
  # 4: Formatting dataframe ----------------------------------------------------
  
  # Merge the data frames of the list in a single data frame
  gensLD <- Reduce(function(...) merge(..., all = T), gensL)
  
  # 
  p <- 0
  gensLCc <- list()
  
  for (i in 1:nrow(gensLD)){
    p <- p + 1
    gensLC <- filter(annot, annot$Gen1 == gensLD$name[p] | annot$Gen2 == gensLD$name[p] | 
                       annot$Gen3 == gensLD$name[p]) %>%
      mutate(chr = gensLD$chr[p], 
             gen.start = gensLD$start[p], gen.end = gensLD$end[p], 
             SNP = gensLD$SNP[p], effect = gensLD$effect[p], 
             SNP.pos = gensLD$SNP.pos[p], 
             SNP.start = gensLD$start[p], SNP.end = gensLD$SNP.end[p], 
             SNP.Location = gensLD$SNP.Location[p], distance = gensLD$Distance[p], 
             model = gensLD$model[p], trait = gensLD$trait[p])
    gensLCc[[i]] <- gensLC
  }
  
  # Merge the data frames of the list in a single data frame and modify it
  gensF <- Reduce(function(...) merge(..., all = T), gensLCc) %>% 
    select(SNP, model, trait, chr, Gen1, name, gen.start, gen.end, GO, NPI, effect, SNP.pos, SNP.start,
           SNP.end, SNP.Location, distance) %>%
    rename(gen.name = Gen1, gen.name.extend = name)
  
  # 5: Save output ---------------------------------------------------------------
  
  message(paste("Saving output file: ", output, ".GWAS_Annotation.txt", sep = ""))
  
  write.table(gensF, file = paste(output, ".GWAS_Annotation.txt", sep = ""), 
              sep = "\t", quote = FALSE, row.names = FALSE)
  
  # Remove objects that will not be used
  rm(csv.l, data, data.der, data.izq, db, i, gensL, gensLC, gensLCc, gensLD, name.F, name.T, names, p)
  
  message("Done! ")
  
}

