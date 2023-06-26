# Short name: GWAS boxplots
# Description: Plot a boxplot for a SNP: x is genotype, y is phenotype. It is possible to add extra information about the samples using the labelfile
# Output: A single PDF file containing the boxplot of the SNPs
#
# Author: Vianey Barrera-Enriquez (vpbarrerae@gmail.com)
#
# Arguments:
# outputname: Name that will have the output
# dir: Directory to save output
# phenofile: Phenotype data in tabular format. Column with samples name should be names 'Taxa'
# genofile: Genotype data in hapmap format
# snp_list_file: CSV file with three columns:
#                Column 01: SNPS. List of SNPS to plot, name should be the same as in the geno data
#                Column 02: trait. Name of the trait as in the pheno data
#                Column 03: xlabel. Name of the trait to be included as label
# labelfile: CSV file with two columns (optional)
#            Column 01: taxa. Sample names
#            Column 02: label. Label or category to add at the plot



##### To do #####
# 1: Change order and color of labels
# 2: Add significance for pairs



# 0: Function init -------------------------------------------------------------

GWAS_Boxplot <- function(outputname, dir, phenofile, genofile, snp_list_file, labelfile = NULL){
  
  
  
  # 1: Load all the directories, info, and data --------------------------------
  
  # Load packages
  if (!require(tibble)) install.packages(tibble)
  if (!require(dplyr)) install.packages(dplyr)
  if (!require(janitor)) install.packages(janitor)
  if (!require(ggplot2)) install.packages(ggplot2)
  if (!require(Biostrings)) install.packages(Biostrings)
  if (!require(hrbrthemes)) install.packages(hrbrthemes)
  if (!require(forcats)) install.packages(forcats)
  if (!require(ggsignif)) install.packages(ggsignif)
  if (!require(RColorBrewer)) install.packages(RColorBrewer)
  
  library(tibble)
  library(dplyr)
  library(janitor)
  library(ggplot2)
  library(Biostrings)
  library(hrbrthemes)
  library(forcats)
  library(ggsignif)
  library(RColorBrewer)
  
  # Load files
  pheno <- read.csv(phenofile, head = T)
  myG <- read.delim(genofile, head = F)
  geno <- myG %>% column_to_rownames(var = 'V1') %>% janitor::row_to_names(1)
  snp_list <- read.csv(snp_list_file, header = T)
  
  message("Traits: ", paste(unique(snp_list$trait), collapse = ', '))
  message("SNPs to Plot: ", length(unique(snp_list$SNPS)))
  
  dir.create(dir, showWarnings = F)
  
  # 
  iupac <- Biostrings::IUPAC_CODE_MAP
  iupac["A"]<- "AA"
  iupac["C"]<- "CC"
  iupac["G"]<- "GG"
  iupac["T"]<- "TT"
  
  
  
  # 2:     ---------------------------------------------------------------------
  
  if (!is.null(labelfile)){
    
    
    
    # 2.1: Labels --------------------------------------------------------------
    
    message("Labels were provided, they will be included in the boxplot")
    
    # 
    label <- read.csv(labelfile, header = T,
                      col.names = c('Taxa', 'Levels'),
                      colClasses = c('character', 'factor'))
    
    message("Labels to colour data: ", paste(levels(label$Levels), collapse = ', '))
    
    # 
    palette <- brewer.pal(n = length(levels(label$Levels)), name = "Dark2")
    # palette <-  c('#219ebc', '#023047', '#fb8500')
    
    
    
    # 2.2: Prepare data for boxplot with labels --------------------------------
    
    matrix_GT <- label
    
    # 
    pdf(paste(dir, outputname, ".SNPs_boxplot.pdf", sep = ''), onefile = T)
    
    # 
    for (i in 1:dim(snp_list)[1]){
      
      # 
      snpname <- snp_list$SNPS[i]
      trait <- snp_list$trait[i]
      x_label <- snp_list$xlabel[i]
      
      # Join label and genotype data
      matrix_GT <- merge(label, t(geno[snpname, -c(1:10) ]), by.x = 'Taxa', by.y = 'row.names')
      
      # Join label-genotype and pheno data
      data <- merge(matrix_GT, pheno[,c("Taxa", trait)], by = 'Taxa')
      
      # Replace IUPAC to nucleotides
      data$snp <-  as.factor(iupac[data[,snpname]])
      
      # List of comparisons
      comp <- as.list(data.frame((matrix(combn(unique(data$snp), 2), ncol = 3))))
      
      
      # 2.3 Boxplot ------------------------------------------------------------
      
      plot <- data %>% ggplot(aes(x = snp, y = get(trait))) +
        geom_boxplot(fill = c("#FC9AA2", "#FFDAC1", "#B5EAD7")[1:length(levels(data$snp))],
                     alpha = 0.45) + # Boxplot
        geom_jitter(aes(color = Levels, shape = Levels), size = 2, height = 0, width = 0.2,
                   alpha = 1) + # Add Phenotype Data
        geom_signif(comparisons = comp, map_signif_level = T) + # Significance
        labs(x = snpname, y = x_label) + # Labs names
        scale_color_manual(values = palette) +
        theme_minimal() + # Aesthetic 
        theme(legend.position = "bottom", legend.title = element_blank()) 
        
      print(plot)
      
      # Progress bar
      cat('\r', i, ' SNPs processed |', rep('=', i/2 ),
          ifelse(i == dim(snp_list)[1], '|\n',  '>'), sep = '')
      
    }
    
    dev.off()
    
  } else {
    
    message("No labels provided")
    
    # 3.1 Prepare data without labels ------------------------------------------
    
    # 
    pdf(paste(dir, outputname, ".SNPs_boxplot.pdf", sep = ''), onefile = T)
    
    # 
    for (i in 1:dim(snp_list)[1]){
      
      # 
      snpname <- snp_list$SNPS[i]
      trait <- snp_list$trait[i]
      x_label <- snp_list$xlabel[i]
      
      # Join pheno and geno data
      data <- merge(pheno, t(geno[snpname, -c(1:10)]), by.x = 'Taxa', by.y = 'row.names')
      data$snp <- as.factor(Biostrings::IUPAC_CODE_MAP[data[, snpname]])
      
      # comp = as.list(data.frame((combn(unique(data$snp),2))))
      comp <- combn(unique(as.character(data$snp)), 2, simplify = F)
      
      
      
      # 3.2 Boxplot ------------------------------------------------------------
      
      # 
      n <- length(levels(data$snp))
      
      # 
      plot <- data %>% ggplot(aes(x = snp, y = get(trait))) +
        geom_boxplot(fill = c("#FC9AA2", "#FFDAC1", "#B5EAD7")[1:length(levels(data$snp))],
                     alpha = 0.45) + # Boxplot
        geom_jitter(size = 2, height = 0, width = 0.2, alpha = 1) + # Points
        geom_signif(comparisons = comp, map_signif_level = T) + # Significance
        labs(x = snpname, y = x_label) + # Lab names
        theme_minimal() + theme(legend.position = "none") # Aesthetic  
      
      print(plot)
      
      # Progress bar
      cat('\r', i, ' SNPs processed |', rep('=', i/2),
          ifelse(i == dim(snp_list)[1], '|\n',  '>'), sep = '')
      
    }
    
    dev.off()
    
  }
  
  message("Done!")
  
}



###### Example(s) ######
# Set arguments
# outputname <- 
# dir <- 
# phenofile <- 
# genofile <- 
# snp_list_file <-  
# labelfile <- 

# Run function
# GWAS_Boxplot(outputname, dir, phenofile, genofile, snp_list_file, labelfile)
