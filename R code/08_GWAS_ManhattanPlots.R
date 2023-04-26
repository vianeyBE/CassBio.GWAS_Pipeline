# Short name: Manhattan plot for GWAS
# Description: Make a customizable Manhattan plot the GWAS GAPIT results
# Output: Manhattan plots for one, several, or all traits
#
# Authors: Camilo E. Sanchez (c.e.sanchez@cgiar.org) and Vianey Barrera-Enriquez (vpbarrera@gmail.com)
#
# Arguments: 
# pat: The location path of the files
# mod: Names of the models to filter
# Mdir: Name of the directory that contains the GAPIT results



# 1: # Configure the initial requirements --------------------------------------
# Manually using the console
message("Enter the path of file names to looking for\n\n",
        "For example: QTL_LOD_Intervals. The path must finish with a point (.)\n\n",
        "Finish with two tabs")
pat <- scan(what = character(), n = 1)

message("Enter the model(s) of interest\n\n",
        "Options: BLINK, GLM, MLM, FarmCPU\n\n",
        "Finish with two tabs")
mod <- scan(what = character(), n = 4)

message("Enter the working directory\n\n",
        "For example: home/user/folder\n\n",
        "Finish with two tabs")
Mdir <- scan(what = character(), n = 1)

# Set as default
if (rlang::is_empty(pat)) {
  pat <- "GAPIT.Association.GWAS_Results."
  } 
if (rlang::is_empty(mod)) {
  mod <- c("BLINK", "FarmCPU", "MLM")
  }
if (rlang::is_empty(Mdir)) {
  Mdir <- "D:/OneDrive - CGIAR/Cassava_Bioinformatics_Team/01_ACWP_F1_Metabolomics/07_GWAS"
  }



# 2: Load all the directories, info, and data ----------------------------------

# Load packages
if (!require(tidyverse)) install.packages(tidyverse)
if (!require(ggtext)) install.packages(ggtext)

library(tidyverse)
library(ggtext)



# 3: Find all the CSVs with the results and organize them ----------------------

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
  
  dframe <- read.csv(name.F) %>%
    mutate(rename = name.T) %>%
    tidyr::separate(col = rename, into = c("batch", "data"), sep = paste(pat))
  
  if (gdata::startsWith(dframe[1,10], "Far")) {
    dframe <- dframe %>% tidyr::separate(col = data, into = c("na", "na2", "trait"), sep = ".m") %>%
      tidyr::separate(col = trait, into = c("trait", "na3"), sep = ".csv") %>%
      mutate(model = "FarmCPU") %>%
      select(SNP, Chr, Pos, P.value, MAF, nobs, Effect, model, trait)
  } else {
    dframe <- dframe %>% tidyr::separate(col = data, into = c("model", "trait"), sep = ".m") %>%
      tidyr::separate(col = trait, into = c("trait", "na"), sep = ".csv") %>%
      select(SNP, Chr, Pos, P.value, MAF, nobs, Effect, model, trait)
  }
  
  csv.l[[i]] <- dframe

  # Progress bar
  cat('\r', i, ' files processed |', rep('=', i / 4), ifelse(i == length(names), '|\n', '>'), sep = '')
}

# Merge the data frames of the list in a single data frame
message("All GWAS files were processed")

GWAS <- bind_rows(csv.l)
GWAS <- GWAS %>% na.omit() %>% dplyr::filter(model %in% mod)
GWAS$Pos<- as.numeric(GWAS$Pos)

message(paste("There are:\n",
              dim(GWAS)[1], "SNPs to plot\n",
              length(table(GWAS$trait)), "traits\n",
              length(table(GWAS$model)), "different models (specified in the 'mod' object)"))

# Delete unnecessary objects
rm(i, p, name.F, name.T, dframe)



# 4: Preparing the data and plot -----------------------------------------------

wtd <- readline(prompt = message("How many traits do you want to plot (One/Several/All)?:"))

message("The color of the chromosomes in the Manhattan plots are default to 'grey' and 'skyblue'\n\n",
        "If you want to change the colors, provide 2 or more. (It can be color names o hex codes)\n\n",
        "Finish with two tabs")

colors <- scan(what = character())

if (rlang::is_empty(colors)) {
  colors <- c("grey", "skyblue")
  } 

if (wtd == "One") {
  
  tr <- readline(prompt = message(paste(unique(GWAS$trait), collapse = "\n"), "\n\n",
                                  "Choose a trait from the list above and paste its name (exact match): "))
  
  data <- GWAS %>% dplyr::filter(trait == tr)
  
  message(paste("You chose the following trait:", tr, "\n\n",
                "For that trait, there are:\n",
                dim(data)[1], "SNPs to plot\n"))
  
  data_cum <- data %>% group_by(Chr) %>%
    summarise(max_bp = max(Pos)) %>% 
    mutate(bp_add = lag(cumsum(max_bp), default = 0)) %>% 
    select(Chr, bp_add)
  
  data <- data %>% inner_join(data_cum, by = "Chr") %>% 
    mutate(bp_cum = Pos + bp_add)
  
  br <- max(-log10(0.001/dim(data)[1]), -log10(min(data$P.value))) / 5
  
  axis_set <- data %>% group_by(Chr) %>%
    summarise(min = min(Pos), max = max(Pos)) %>%
    inner_join(data_cum, by = "Chr") %>%
    mutate(center = ((max - min) /2) + bp_add)
  
  message("Making the plot. It can take a few seconds. Please be patient.")
  
  # Manhattan plot
  ggplot(data, aes(bp_cum, -log10(P.value), color = as_factor(Chr),
                   size = -log10(P.value), shape = model)) +
    geom_hline(yintercept = -log10(0.001/dim(data)[1]), color = "black", linetype = "dashed") +
    geom_point(alpha = 0.75) +
    scale_color_manual(values = rep(c("grey", "skyblue"), 18)) +
    scale_y_continuous(expand = c(0, 0), limits = c(0, br * 5.1),
                       breaks = c(round(br, 2), round(br * 2, 2), round(br * 3, 2), round(br * 4, 2),
                                  round(br * 5, 2))) +
    scale_x_continuous(expand = c(0, 0), label = axis_set$Chr, breaks = axis_set$center) +
    guides(color = "none", size = "none") +
    labs(y = "-log<sub>10</sub>(p)", x = "Chromosome", title = paste("Trait: ", tr), shape = "Model") +
    theme_classic() +
    theme(legend.position = "bottom", axis.title.y = element_markdown(),
          panel.grid.major.x = element_blank(), panel.grid.minor.x = element_blank(),
          legend.title = element_text(colour = "black", size = 12),
          legend.text = element_text(colour = "black", size = 12),
          plot.title = element_text(size = 16, color = "black"),
          axis.text = element_text(size = 12, color = "black"),
          axis.title = element_text(size = 14, color = "black"),
          axis.title.x = element_text(margin = margin(t = 10)))
  
  } else {
    if (wtd == "Several") {
    
    message(paste(unique(GWAS$trait), collapse = "\n"), "\n\n",
            "Choose the traits from the list above and paste their names in the next function\n\n",
            "Finish with two tabs")
    
    tr <- scan(what = character())
    
    data <- GWAS %>% dplyr::filter(trait %in% tr)
    
    message(paste("You chose the following traits:", paste(tr, collapse = "  ")))
    
    data_cum <- data %>% group_by(Chr) %>%
      summarise(max_bp = max(Pos)) %>% 
      mutate(bp_add = lag(cumsum(max_bp), default = 0)) %>% 
      select(Chr, bp_add)
    
    data <- data %>% inner_join(data_cum, by = "Chr") %>% 
      mutate(bp_cum = Pos + bp_add)
    
    br <- max(-log10(0.001/dim(data)[1]), -log10(min(data$P.value))) / 5
    
    axis_set <- data %>% group_by(Chr) %>%
      summarise(min = min(Pos), max = max(Pos)) %>%
      inner_join(data_cum, by = "Chr") %>%
      mutate(center = ((max - min) /2) + bp_add)
    
    message("Making the plots. It can take a few seconds. Please be patient.")
    
    # Manhattan plot
    ggplot(data, aes(bp_cum, -log10(P.value), color = as_factor(Chr),
                     size = -log10(P.value))) +
      geom_hline(yintercept = -log10(0.001/dim(data)[1]), color = "black", linetype = "dashed") +
      geom_point(alpha = 0.75) +
      scale_color_manual(values = rep(colors, 18)) +
      scale_y_continuous(expand = c(0, 0), limits = c(0, br * 5.1),
                         breaks = c(round(br, 2), round(br * 2, 2), round(br * 3, 2), round(br * 4, 2),
                                    round(br * 5, 2))) +
      scale_x_continuous(expand = c(0, 0), label = axis_set$Chr, breaks = axis_set$center) +
      guides(color = "none", size = "none") +
      labs(y = "-log<sub>10</sub>(p)", x = "Chromosome") +
      theme_classic() +
      facet_wrap(~ trait, ncol = 1) +
      theme(legend.position = "bottom", axis.title.y = element_markdown(),
            panel.grid.major.x = element_blank(), panel.grid.minor.x = element_blank(),
            legend.title = element_text(colour = "black", size = 12),
            legend.text = element_text(colour = "black", size = 12),
            plot.title = element_text(size = 16, color = "black"),
            axis.text = element_text(size = 12, color = "black"),
            axis.title = element_text(size = 14, color = "black"),
            axis.title.x = element_text(margin = margin(t = 10)),
            strip.text = element_text(size = 12, color = "black", face = "bold"),
            strip.background = element_rect(color = "black", fill = "white"))
    
    } else {
      if (wtd == "All") {
      
      data <- GWAS
      
      message(paste(dim(data)[1], "SNPs and", length(unique(data$trait)), "traits will be plotted"))
      
      data_cum <- data %>% group_by(Chr) %>%
        summarise(max_bp = max(Pos)) %>% 
        mutate(bp_add = lag(cumsum(max_bp), default = 0)) %>% 
        select(Chr, bp_add)
      
      data <- data %>% inner_join(data_cum, by = "Chr") %>% 
        mutate(bp_cum = Pos + bp_add)
      
      br <- max(-log10(0.001/dim(data)[1]), -log10(min(data$P.value))) / 5
      
      axis_set <- data %>% group_by(Chr) %>%
        summarise(min = min(Pos), max = max(Pos)) %>%
        inner_join(data_cum, by = "Chr") %>%
        mutate(center = ((max - min) /2) + bp_add)
      
      message("Making the plot. It can take a few seconds. Please be patient.")
      
      # Manhattan plot
      ggplot(data, aes(bp_cum, -log10(P.value), color = as_factor(Chr),
                       size = -log10(P.value), shape = model)) +
        geom_hline(yintercept = -log10(0.001/dim(data)[1]), color = "black", linetype = "dashed") +
        geom_point(alpha = 0.75) +
        scale_color_manual(values = rep(c("grey", "skyblue"), 18)) +
        scale_y_continuous(expand = c(0, 0), limits = c(0, br * 5.1),
                           breaks = c(round(br, 2), round(br * 2, 2), round(br * 3, 2), round(br * 4, 2),
                                      round(br * 5, 2))) +
        scale_x_continuous(expand = c(0, 0), label = axis_set$Chr, breaks = axis_set$center) +
        guides(color = "none", size = "none") +
        labs(y = "-log<sub>10</sub>(p)", x = "Chromosome", title = "All traits", shape = "Model") +
        theme_classic() +
        theme(legend.position = "bottom", axis.title.y = element_markdown(),
              panel.grid.major.x = element_blank(), panel.grid.minor.x = element_blank(),
              legend.title = element_text(colour = "black", size = 12),
              legend.text = element_text(colour = "black", size = 12),
              plot.title = element_text(size = 16, color = "black"),
              axis.text = element_text(size = 12, color = "black"),
              axis.title = element_text(size = 14, color = "black"),
              axis.title.x = element_text(margin = margin(t = 10)))
      
      } else {
      message("You can only choose the following options: 'One', 'Several', or 'All' (exact match)")
    }
  } 
}
