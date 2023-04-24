# Short name: Multivariate methods for GWAS covariable selection
# Description: It can perform several multivariate methods (MDS, PCA, DAPC) to identify covariates for GWAS 
# Output(s): It depends on the chosen method. For example:
#            - NDMS: An interactive html plot
#            - MDS: An interactive html plot
#            - PCA:
#            - DAPC:
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
#
# 3. For PCA:
# dir: Directory where data is located.
# type: A character string indicating wheter data is genotypic ("geno") of phenotypic ("pheno").
# groups: Boolean value indicating whether the data includes different treatments/groups (default = F). If `TRUE` is selected then the second column of the `phenofile` must be the groups/treatments.
# phenofile: A database of genotypes/individuals (rows) and their traits (columns). First column must be genotypes/individuals names.
# genofile: An object of class SNPGDSFileClass (GDS file read with the `snpgdsOpen` function from `SNPRelate` package).
# labels: When provide a `genofile` and `groups` argument is `TRUE`, please provide a dataframe with genotypes/individuals in the first column and groups/treatment data in the second column.
# gds: A Genomic Data Structures (GDS) file (a reformatted VCF file with the `snpgdsVCF2GDS` function from `vcfR` package).
# vcf: A Variant Call Format (VCF) file containing DNA polymorphism data such as SNPs, insertions, deletions and structural variants, together with rich annotations.
# PC.retain: Boolean value indicating whether to analyze how many PCs retain (default = F).
#
# 4. For DAPC:
# my_genind: An object of class genind (vcf file read with the `vcfR2genind` function from `vcfR` package).



####### To do ####### 
# 1. Change gds/vcf part



# Test arguments
dir <- "D:/OneDrive - CGIAR/Cassava_Bioinformatics_Team/02_CTS_Drought_Family/01_Phenotype_Preliminar_Analysis/"
dir <- "D:/OneDrive - CGIAR/Cassava_Bioinformatics_Team/03_GWAS_PPD_Populations/02_PCA/"
dist <- c("gower")
type <- "Geno"

# Phenotypic test data
setwd(dir)
phenofile <- read.csv("Prueba.csv", header = T) # Add colors
phenofile <- phenofile[-2] # No color

# Genotypic test data
setwd(dir)
labels <- read.csv('GWAS_PPD.labels.csv')
vcf <- 'GWAS_PPD.snps.filter_info.missing_0.10.imputation.vcf.gz'
gds <- 'GWAS_PPD.gds'

dir <- "D:/OneDrive - CGIAR/Cassava_Bioinformatics_Team/01_ACWP_F2_Phenotype/01_Population_Structure"
setwd(dir)
vcf <- "AM1588_MAP.miss0.05.recode.vcf"




# 1: Non-metric multidimensional scaling (NDMS) --------------------------------

NDMS <- function(dir, phenophile, dist, groups = F){
  
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

NDMS(dir, phenofile, dist, groups)



# 2: Multidimensional scaling (MDS) --------------------------------------------

MDS <- function(dir, phenofile, dist, groups = F){
  
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

MDS(dir, phenofile, dist, groups)



# 3: Principal component analysis (PCA) ----------------------------------------

PCA <- function(dir, type, groups = T, phenofile = NULL, genofile = NULL, labels = NULL, gds = NULL,
                vcf = NULL, PC.retain = F){
  
  # Set working directory
  setwd(dir)
  
  # Load libraries
  library(SNPRelate)
  library(gdsfmt)
  library(tidyverse)
  library(psych)
  library(plotly)
  library(htmlwidgets)
  
  # Conditional to determine if the data is genotypic of phenotypic
  if (type == "Pheno"){
    
    message("Selected data option: Phenotypic")
    
    # Conditional for groups / no groups database handling
    if (groups == T){
      
      # Database handling for color/groups
      names <- phenofile[1]
      tr <- phenofile[2]
      phenofile <- phenofile[3:dim(phenofile)[2]]
      
    } else {
      
      # Database handling for non color/groups
      names <- phenofile[1]
      phenofile <- phenofile[3:dim(phenofile)[2]]
      
    }
    
    # Function continues to the PCA calculation itself
    
    # Performs a 'normal' PCA
    PCA <- prcomp(phenofile)
    
    # Creates the table to draw scree plots and analyze the number of PCs to retain
    PC <- data.frame(PC = 1:dim(phenofile)[2], eigen = PCA$sdev^2,
                     var = ((sqrt(PCA$sdev^2))^2 / (sum(var(phenofile)))) * 100) %>%
      mutate(var.cum = cumsum(var))
    
    # Conditional for PCs retaining analysis 
    if (PC.retain == T){
      
      message("Principal components retaining analyses in process")
      
      # Establish a cutoff and obtain the genetic correlation matrix
      cutoff <- 1 / length(PC$PC) * 100
      
      # Makes the parallel, VSS, Velicer's MAP, and BIC analyses
      FA.p <- fa.parallel(phenofile, fm = "pa", fa = "pc", n.iter = 30, plot = F)
      FAC <- nfactors(phenofile)
      dev.off()
      
      # Obtain the summary statistics
      FAC.t <- FAC[["vss.stats"]] %>% mutate(map = FAC$map)
      
      # What would you get with different ways to retain PCs
      message(paste("There are a total of", length(PC$PC), "PCs \n\nIf you: \n",
                    "(1) Use the 'K1' method, you would retain",
                    dim(dplyr::filter(PC, eigen > 1))[1], "PCs \n",
                    "(2) Accept all the PCs that explain more than one variable’s worth of data, you would retain",
                    dim(dplyr::filter(PC, var > cutoff))[1], "PCs \n",
                    "(3) Retain the PCs that explain at least a 70% of cumulative variance, you would retain",
                    dim(dplyr::filter(PC, var.cum < 70))[1], "PCs \n",
                    "(4) Establish an 80% threshold, you would retain",
                    dim(dplyr::filter(PC, var.cum < 80))[1], "PCs \n",
                    "(5) Perform a Parallal analysis, you would retain",
                    FA.p$ncomp, "PCs \n",
                    "(6) Apply the 'Very Simple Structure criteria', you would retain",
                    which.max(FAC.t$cfit.2), "PCs \n",
                    "(7) Try with the 'Velicer's MAP criteria', you would retain",
                    which.min(na.omit(FAC.t$map)), "PCs \n\n",
                    "Also, you could retain a fixed number of PCs or retain PCs based on the scree plot"))
      
    } else {
      
      message("NOT doing the PCs retaining analyses")
      
    }
    
    # Function continues to the plots
    # Scree plots
    # Individual variance
    plot_ly(data = PC, x = ~ PC, y = ~ var, type = "scatter", mode = "lines+markers", text = ~ round(var, 2),
            line = list(color = "grey")) %>%
      layout(xaxis = list(title = "PC"), yaxis = list(title = "Explained variance (%)"))
    
    # Cumulative variance
    plot_ly(data = PC, x = ~ PC, y = ~ var.cum, type = "scatter", mode = "lines+markers",
            text = ~ round(var.cum, 2), line = list(color = "grey")) %>%
      layout(xaxis = list(title = "PC"), yaxis = list(title = "Cumulative variance (%)"))
    
    # Conditional for groups / no groups PCA plot
    if (groups == T){
      
      # Table to plot the PCA
      dt <- cbind(names, tr, data.frame(PCA[["x"]])) %>%
        dplyr::rename(sample.id = 1, label = 2) %>%
        select(sample.id, label, PC1, PC2)
      
      # Plot the PCA ###################
      fig <- plot_ly(data = dt, x = ~ PC1, y = ~ PC2, color = ~ as.factor(label), type = "scatter",
                     mode = "markers", symbol = ~ label, symbols = c("circle", "x"), text = ~ sample.id,
                     marker = list(size = 6)) %>%
        layout(legend = list(title = list(text = "<b> Groups </b>"), orientation = "h"),
               xaxis = list(title = paste("Dimension 1 - ", round(PC$var)[1], "%")),
               yaxis = list(title = paste("Dimension 2 - ", round(PC$var)[2], "%")))
      
    } else {
      
      # Table to plot the PCA
      dt <- cbind(names, data.frame(PCA[["x"]])) %>%
        dplyr::rename(sample.id = 1) %>%
        select(sample.id, PC1, PC2)
      
      # Plot the PCA
      fig <- plot_ly(data = dt, x = ~ PC1, y = ~ PC2, type = "scatter", mode = "markers",
                     text = ~ sample.id, marker = list(size = 6)) %>%
        layout(xaxis = list(title = paste("Dimension 1 - ", round(PC$var)[1], "%")),
               yaxis = list(title = paste("Dimension 2 - ", round(PC$var)[2], "%")))
      
    }
    
    # Function continues to save the PCA plot
    
    # Save the plot
    saveWidget(fig, "GWAS_PCA.html", selfcontained = F, libdir = "lib")
    saveWidget(as_widget(fig), "GWAS_PCA.html")
    
  } else {
    
    message("Selected data option: Genotypic")
    
    # Conditional to determine how to proceed with the files
    if (is.null(genofile) & is.null(gds)) {
      
      message("VCF file provided...\n\n", "Reading VCF file...\n\n",
              "Reformatting it to a GDS file and then transforming it to a SNPGDSFileClass file")
      
      snpgdsVCF2GDS(vcf, paste0(substring(vcf, 1, nchar(vcf)-7), ".gds"), ignore.chr.prefix = "chromosome")
      genofile <- snpgdsOpen(gds)
      sample_id <- read.gdsn(index.gdsn(genofile, "sample.id"))
      
      if (is.null(genofile) & is.null(vcf)) {
        
        message("GDS file provided...\n\n", "Reading CDS file...\n\n",
                "Transforming it to a SNPGDSFileClass file")
        
        genofile <- snpgdsOpen(gds)
        sample_id <- read.gdsn(index.gdsn(genofile, "sample.id"))
        
      } else {
          
          message("Genofile (SNPGDSFileClass file) is provided...\n\n",
                  "PCA function continues regularly")
        
      }
    }
    
    # Function continues to the PCA calculation itself
    # Performs the SNPRelate's PCA
    PCA <- snpgdsPCA(genofile, autosome.only = T, remove.monosnp = T, need.genmat = T,
                     algorithm = "exact", eigen.method = "DSPEVX")
    
    # Creates the table to draw scree plots and analyze the number of PCs to retain 
    PC <- data.frame(PC = 1:length(PCA$varprop), id = PCA[["sample.id"]],
                     eigen = PCA[["eigenval"]], var = PCA$varprop * 100) %>%
      mutate(var.cum = cumsum(var))
    
    # Conditional for PCs retaining analysis 
    if (PC.retain == T) {
      
      message("Principal components retaining analyses in process")
      
      # Establish a cutoff and obtain the genetic correlation matrix
      cutoff <- 1 / length(PC$PC) * 100
      cor.m <- PCA[['genmat']]
      cor.m <- as.matrix(replace(cor.m, cor.m > 1, 1))
      
      # Makes the parallel, VSS, Velicer's MAP, and BIC analyses
      FA.p <- fa.parallel(cor.m, n.obs = dim(cor.m)[1], fm = "pa", fa = "pc", n.iter = 30, plot = F)
      FAC <- nfactors(cor.m, diagonal = T, fm = "pa", n.obs = dim(cor.m)[1], plot = F)
      
      # Obtain the summary statistics
      FAC.t <- FAC[["vss.stats"]] %>% mutate(map = FAC$map)
      
      # What would you get with different ways to retain PCs
      message(paste("There are a total of", length(PC$PC), "PCs \n\nIf you: \n",
                    "(1) Use the 'K1' method, you would retain",
                    dim(dplyr::filter(PC, eigen > 1))[1], "PCs \n",
                    "(2) Accept all the PCs that explain more than one variable’s worth of data, you would retain",
                    dim(dplyr::filter(PC, var > cutoff))[1], "PCs \n",
                    "(3) Retain the PCs that explain at least a 70% of cumulative variance, you would retain",
                    dim(dplyr::filter(PC, var.cum < 70))[1], "PCs \n",
                    "(4) Establish an 80% threshold, you would retain",
                    dim(dplyr::filter(PC, var.cum < 80))[1], "PCs \n",
                    "(5) Perform a Parallal analysis, you would retain",
                    FA.p$ncomp, "PCs \n",
                    "(6) Apply the 'Very Simple Structure criteria', you would retain",
                    which.max(FAC.t$cfit.2), "PCs \n",
                    "(7) Try with the 'Velicer's MAP criteria', you would retain",
                    which.min(na.omit(FAC.t$map)), "PCs \n\n",
                    "Also, you could retain a fixed number of PCs or retain PCs based on the scree plot"))
      
    } else {
      
      message("NOT doing the PCs retaining analyses")
      
    }
    
    # Function continues to the plots
    # Scree plots
    # Individual variance
    plot_ly(data = PC, x = ~ PC, y = ~ var, type = "scatter", mode = "lines+markers", text = ~ round(var, 2),
            line = list(color = "grey")) %>%
      layout(xaxis = list(title = "PC"), yaxis = list(title = "Explained variance (%)"))
    
    # Cumulative variance
    plot_ly(data = PC, x = ~ PC, y = ~ var.cum, type = "scatter", mode = "lines+markers",
            text = ~ round(var.cum, 2), line = list(color = "grey")) %>%
      layout(xaxis = list(title = "PC"), yaxis = list(title = "Cumulative variance (%)"))
    
    # Table to plot the PCA
    tab <- data.frame(sample.id = PCA$sample.id, stringsAsFactors = F,
                      PC1 = PCA$eigenvect[,1], PC2 = PCA$eigenvect[,2])
    
    if (groups == T){
      
      # Adding groups to the table for plotting
      labels <- labels %>% dplyr::rename(sample.id = 1, groups = 2)
      dt <- tab %>% inner_join(labels, by = "sample.id")
      
      # Plot the PCA
      fig <- plot_ly(data = dt, x = ~ PC1, y = ~ PC2, color = ~ as.factor(groups), type = "scatter",
                     mode = "markers", symbol = ~ as.factor(groups), symbols = c("circle", "x"),
                     text = ~ sample.id, marker = list(size = 6)) %>%
        layout(legend = list(title = list(text = "<b> Groups </b>"), orientation = "h"),
               xaxis = list(title = paste("Dimension 1 - ", round(PC$var)[1], "%")),
               yaxis = list(title = paste("Dimension 2 - ", round(PC$var)[2], "%")))
      
    } else {
      
      # Plot the PCA
      fig <- plot_ly(data = dt, x = ~ PC1, y = ~ PC2, type = "scatter", mode = "markers",
                     text = ~ sample.id, marker = list(size = 6)) %>%
        layout(xaxis = list(title = paste("Dimension 1 - ", round(PC$var)[1], "%")),
               yaxis = list(title = paste("Dimension 2 - ", round(PC$var)[2], "%")))
      
    }
    
    # Save the plot
    saveWidget(fig, "GWAS_PCA.html", selfcontained = F, libdir = "lib")
    saveWidget(as_widget(fig), "GWAS_PCA.html")
    
  }
}


PCA(dir, type, groups = T, phenofile, genofile = NULL, vcf, gds, PC.retain = F)


  
# 4: Discriminant analysis of principal components (DAPC) ----------------------

# Load data
vcf <- read.vcfR('GWAS_PPD.snps.filter_info.missing_0.10.imputation.vcf.gz', verbose = T)
my_genind <- vcfR2genind(vcf)
my_genind

# Load libraries
library(adegenet)
library(grDevices)
library(vcfR)
library(tidyverse)
library(plotly)
library(htmlwidgets)

# Identify clusters
# Shows a graph with the accumulated variance explained by the eigenvalues of the PCA
# All PC were retained (200, actually there is less, about 150)
# The second graph shows the elbow at k = 3 (Number of clusters)
grp <- find.clusters(my_genind)

# Transform the data using PCA and then performs a discriminant Analysis.
# DAPC actually benefit from less PC. Here 80 will be selected
dapc <- dapc(my_genind, grp$grp)

# Plot the DAPC
pdf("GWAS_DAPC_scatter.pdf", width = 10, height = 8) 
scatter(dapc, scree.pca = F, ratio.pca = 0.3, pch = 20, cell = 1, solid = 0.6, cex = 2.5, clab = 0,
        scree.da = F, leg = T, txt.leg = paste("Cluster", 1:3))
dev.off()

set.seed(4)
contrib <- loadingplot(dapc$var.contr, axis = 2, lab.jitter = 1)

groups <- tibble::rownames_to_column(as.data.frame(dapc[["grp"]]), var = "Name")
colnames(groups) <- c('Name','group_dapc')

write.csv(groups, 'GWAS_PPD.dapc_groups.csv', row.names = F, col.names = T)
