# Short name: Multivariate methods for GWAS covariable selection
# Description: It can perform several multivariate methods (NDMS, or MDS) to identify covariates for GWAS 
# Output(s): Interactive html plots
#
# Authors: Camilo E. Sanchez (c.e.sanchez@cgiar.org) and Vianey Barrera-Enriquez (vpbarrera@gmail.com)
#
# Arguments: 
# 1. For NDMS:
# dir: Directory where data is located.
# phenofile: A database of genotypes/individuals (rows) and their traits (columns). First column must be genotypes/individuals names.
# dist: Dissimilarity index/measure to use (default = "bray").
# groups: Boolean value indicating whether the data includes different treatments/groups (default = F). If `TRUE` is selected then the second column of the `phenofile` must be the groups/treatments.
#
# 2. For MDS:
# dir: Directory where data is located.
# phenofile: A database of genotypes/individuals (rows) and their traits (columns). First column must be genotypes/individuals names.
# dist: Dissimilarity index/measure to use (default = "gower").
# groups: Boolean value indicating whether the data includes different treatments/groups (default = F). If `TRUE` is selected then the second column of the `phenofile` must be the groups/treatments.



##### To do #####
# Everything is good!



# 1: Non-metric multidimensional scaling (NDMS) --------------------------------

NDMS <- function(dir, phenophile, dist = "bray", groups = F){
  
  # Set working directory
  setwd(dir)
  
  # Load libraries
  library(vegan)
  library(tidyverse)
  library(plotly)
  library(htmlwidgets)
  
  # Conditional for color case
  if (groups == T){
    
    # Database handling
    names <- phenofile[1]
    tr <- phenofile[2]
    phenofile <- phenofile[3:dim(phenofile)[2]]
    
    # Peform the NDMS
    NDMS <- metaMDS(phenofile, distance = dist, binary = F, autotransform = T)
    NDMS.pd <- as.data.frame(scores(MDS, display = c("sites")))
    
    # Merge MDS data with their names
    NDMS.f <- cbind(names, tr, NDMS.pd) 
    NDMS.f <- NDMS.f %>% rename(Label = 1, Treat = 2)
    
    # Plot MDS
    fig <- plot_ly(data = NDMS.f, x = ~ NMDS1, y = ~ NMDS2, color = ~ as.factor(Treat), type = 'scatter',
                   mode = 'markers', symbol = ~ as.factor(Treat), text = ~ Label) %>%
      layout(legend = list(title = list(text = '<b> Groups </b>'), orientation = 'h'))
    
    # Save the plot
    saveWidget(fig, "GWAS_NDMS.html", selfcontained = F, libdir = "lib")
    saveWidget(as_widget(fig), "GWAS_NDMS.html")
    
    } else {
      
      # Database handling
      names <- phenofile[1]
      phenofile <- phenofile[-1]
      
      # Perform NDMS
      NDMS <- metaMDS(phenofile, distance = dist, binary = F, autotransform = T)
      NDMS.pd <- as.data.frame(scores(NDMS, display = c("sites")))
      
      # Merge MDS data with their names
      NDMS.f <- cbind(names, NDMS.pd)
      NDMS.f <- NDMS.f %>% rename(Label = 1)
      
      # Plot MDS
      fig <- plot_ly(data = NDMS.f, x = ~ NMDS1, y = ~ NMDS2, type = 'scatter',
                     mode = 'markers', text = ~ Label)
      
      # Save the plot
      saveWidget(fig, "GWAS_NDMS.html", selfcontained = F, libdir = "lib")
      saveWidget(as_widget(fig), "GWAS_NDMS.html")
      
    }
  
}



# 1.1: NDMS example(s) ---------------------------------------------------------
# Phenotypic with groups
# dir <- "D:/OneDrive - CGIAR/Cassava_Bioinformatics_Team/02_CTS_Drought_Family/01_Phenotype_Preliminar_Analysis/"
# phenofile <- read.csv(paste0(dir, "Prueba.csv"), header = T)
# groups <- T

# Phenotypic without groups
# dir <- "D:/OneDrive - CGIAR/Cassava_Bioinformatics_Team/02_CTS_Drought_Family/01_Phenotype_Preliminar_Analysis/"
# phenofile <- read.csv(paste0(dir, "Prueba.csv"), header = T)
# phenofile <- phenofile[-2]



# Run function -----------------------------------------------------------------
# NDMS(dir, phenofile, dist, groups)



# 2: Multidimensional scaling (MDS) --------------------------------------------

MDS <- function(dir, phenofile, dist = "gower", groups = F){
  
  # Set working directory
  setwd(dir)
  
  # Load libraries
  library(vegan)
  library(tidyverse)
  library(plotly)
  library(htmlwidgets)
  
  # Conditional for color case
  if (groups == T){
    
    # Database handling
    names <- phenofile[1]
    tr <- phenofile[2]
    phenofile <- phenofile[3:dim(phenofile)[2]]
    
    # Performs MDS, calculates cumulative variance
    diss <- vegdist(phenofile, method = dist) # dist function can also be used
    MDS <- cmdscale(diss, eig = T, x.ret = T)
    MDS.var <- round(MDS$eig / (sum(MDS$eig) * 100), 3)
    
    # Merge NDMS data with their names
    MDS.f <- cbind(names, tr, MDS$points)
    MDS.f <- MDS.f %>% rename(Label = 1, Treat = 2, MDS1 = 3, MDS2 = 4)
    
    # Plot
    fig <- plot_ly(data = MDS.f, x = ~ MDS1, y = ~ MDS2, color = ~ as.factor(Treat), type = 'scatter',
                   mode = 'markers', symbol = ~ as.factor(Treat), text = ~ Label) %>%
      layout(legend = list(title = list(text = '<b> Groups </b>'), orientation = 'h'),
             xaxis = list(title = paste("MDS1 - ", MDS.var[1], "%", sep = "")), 
             yaxis = list(title = paste("MDS2 - ", MDS.var[2], "%", sep = "")))
    
    # Save the plot
    saveWidget(fig, "GWAS_MDS.html", selfcontained = F, libdir = "lib")
    saveWidget(as_widget(fig), "GWAS_MDS.html")
    
  } else {
    
    # Database handling
    names <- phenofile[1]
    phenofile <- phenofile[2:dim(phenofile)[2]]
    
    # Performs MDS, calculates cumulative variance
    diss <- vegdist(phenofile, method = dist) # dist function can also be used
    MDS <- cmdscale(diss, eig = T, x.ret = T)
    MDS.var <- round(MDS$eig / (sum(MDS$eig) * 100), 3)
    
    # Merge NDMS data with their names
    MDS.f <- cbind(names, MDS$points)
    MDS.f <- MDS.f %>% rename(Label = 1, MDS1 = 2, MDS2 = 3)
    
    # Plot
    fig <- plot_ly(data = MDS.f, x = ~ MDS1, y = ~ MDS2, type = 'scatter',
                   mode = 'markers', text = ~ Label) %>%
      layout(xaxis = list(title = paste("MDS1 - ", MDS.var[1], "%", sep = "")), 
             yaxis = list(title = paste("MDS2 - ", MDS.var[2], "%", sep = "")))
    
    # Save the plot
    saveWidget(fig, "GWAS_MDS.html", selfcontained = F, libdir = "lib")
    saveWidget(as_widget(fig), "GWAS_MDS.html")
    
  }
  
}



# 2.1: MDS example(s) ----------------------------------------------------------
# Phenotypic with groups
# dir <- "D:/OneDrive - CGIAR/Cassava_Bioinformatics_Team/02_CTS_Drought_Family/01_Phenotype_Preliminar_Analysis/"
# phenofile <- read.csv(paste0(dir, "Prueba.csv"), header = T)
# groups <- T

# Phenotypic without groups
# dir <- "D:/OneDrive - CGIAR/Cassava_Bioinformatics_Team/02_CTS_Drought_Family/01_Phenotype_Preliminar_Analysis/"
# phenofile <- read.csv(paste0(dir, "Prueba.csv"), header = T)
# phenofile <- phenofile[-2]



# 2.2: Run MDS function --------------------------------------------------------
# MDS(dir, phenofile, dist, groups)
