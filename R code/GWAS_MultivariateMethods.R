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
#        dir: Name of the directory that contains the data
#        traits: A matrix/database of genotypes/individuals in rows and their traits in in columns
#        dist: Dissimilarity index/measure to use. The default is bray.
#        trata: Boolean value indicating whether the data includes different treatments
# 2. For MDS:
#        dir: Name of the directory that contains the data
#        traits: A matrix/database of genotypes/individuals in rows and their traits in in columns
#        dist: Dissimilarity index/measure to use. The default is bray.
#        trata: Boolean value indicating whether the data includes different treatments
# 3. For PCA:
#        dir: Name of the directory that contains the data
#        type: A character string indicating wheter data is genotypic ("geno") of phenotypic ("pheno")
#        trata: Boolean value indicating whether the data includes different treatments
#        genofile: An object of class SNPGDSFileClass. A SNP - GDS file
#        phenofile: A matrix/database of genotypes/individuals in rows and their traits in in columns
#        vcf: 
#        gds: 
#        PC.retain: 
#        
# 4. For DAPC:
#



####### To do ####### 
# 1: Add function structures for PCA (including plotly - htmal plots)
# 2: Finish head description

# Test arguments
dir <- "D:/OneDrive - CGIAR/Cassava_Bioinformatics_Team/02_CTS_Drought_Family/01_Phenotype_Preliminar_Analysis/"
dist <- c("gower")
trata <- F

# Phenotypic test data
setwd(dir)
traits <- read.csv("Prueba.csv", header = T) # Add colors
traits <- traits[-2] # No color





# 1: Non-metric multidimensional scaling (NDMS) --------------------------------

NDMS <- function(dir, traits, dist, trata = F){
  
  # Set working directory
  setwd(dir)
  
  # Load libraries
  library(vegan)
  library(tidyverse)
  library(plotly)
  library(htmlwidgets)
  
  # Conditional for color case
  if (trata == T){
    
    # Database handling
    names <- traits[1]
    tr <- traits[2]
    traits <- traits[3:dim(traits)[2]]
    
    # Peform the NDMS
    NDMS <- metaMDS(traits, distance = dist, binary = F, autotransform = T)
    NDMS.pd <- as.data.frame(scores(MDS, display = c("sites")))
    
    # Merge MDS data with their names
    NDMS.f <- cbind(names, tr, NDMS.pd) 
    NDMS.f <- NDMS.f %>% rename(Label = 1, Treat = 2)
    
    # Plot MDS
    fig <- plot_ly(data = NDMS.f, x = ~ NMDS1, y = ~ NMDS2, color = ~ as.factor(Treat), type = 'scatter',
                   mode = 'markers', symbol = ~ as.factor(Treat), text = ~ Label) %>%
      layout(legend = list(title = list(text = '<b> Treatment </b>'), orientation = 'h'))
    
    # Save the plot
    saveWidget(fig, "GWAS_NDMS.html", selfcontained = F, libdir = "lib")
    saveWidget(as_widget(fig), "GWAS_NDMS.html")
    
    } else {
      
      # Database handling
      names <- traits[1]
      traits <- traits[-1]
      
      # Perform NDMS
      NDMS <- metaMDS(traits, distance = dist, binary = F, autotransform = T)
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

NDMS(dir, traits, dist, trata)



# 2: Multidimensional scaling (MDS) --------------------------------------------

MDS <- function(dir, traits, dist, trata = F){
  
  # Set working directory
  setwd(dir)
  
  # Load libraries
  library(vegan)
  library(tidyverse)
  library(plotly)
  library(htmlwidgets)
  
  # Conditional for color case
  if (trata == T){
    
    # Database handling
    names <- traits[1]
    tr <- traits[2]
    traits <- traits[3:dim(traits)[2]]
    
    # Performs MDS, calculates cumulative variance
    diss <- vegdist(traits, method = dist) # dist function can also be used
    MDS <- cmdscale(diss, eig = T, x.ret = T)
    MDS.var <- round(MDS$eig / (sum(MDS$eig) * 100), 3)
    
    # Merge NDMS data with their names
    MDS.f <- cbind(names, tr, MDS$points)
    MDS.f <- MDS.f %>% rename(Label = 1, Treat = 2, MDS1 = 3, MDS2 = 4)
    
    # Plot
    fig <- plot_ly(data = MDS.f, x = ~ MDS1, y = ~ MDS2, color = ~ as.factor(Treat), type = 'scatter',
                   mode = 'markers', symbol = ~ as.factor(Treat), text = ~ Label) %>%
      layout(legend = list(title = list(text = '<b> Treatment </b>'), orientation = 'h'),
             xaxis = list(title = paste("MDS1 - ", MDS.var[1], "%", sep = "")), 
             yaxis = list(title = paste("MDS2 - ", MDS.var[2], "%", sep = "")))
    
    # Save the plot
    saveWidget(fig, "GWAS_MDS.html", selfcontained = F, libdir = "lib")
    saveWidget(as_widget(fig), "GWAS_MDS.html")
    
  } else {
    
    # Database handling
    names <- traits[1]
    traits <- traits[2:dim(traits)[2]]
    
    # Performs MDS, calculates cumulative variance
    diss <- vegdist(traits, method = dist) # dist function can also be used
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

MDS(dir, traits, dist, trata)



# 3: Principal component analysis (PCA) ----------------------------------------


# Files names
labels <- read.csv('GWAS_PPD.labels.csv')
vcf <- 'GWAS_PPD.snps.filter_info.missing_0.10.imputation.vcf.gz'
gds <- 'GWAS_PPD.gds'
type <- "Geno"

# Load files
setwd("D:/OneDrive - CGIAR/Cassava_Bioinformatics_Team/03_GWAS_PPD_Populations/02_PCA")




PCA <- function(dir, type, trata, phenofile = NULL, genofile = NULL, vcf = NULL, gds = NULL, PC.retain = T){
  
  # Set working directory
  setwd(dir)
  
  # Load libraries
  library(SNPRelate)
  library(gdsfmt)
  library(tidyverse)
  library(psych)
  library(plotly)
  library(htmlwidgets)
  
}

# Conditional to determine if the data is genotypic of phenotypic
if (type == "Pheno") {
  
  message("Selected data option: Phenotypic")
  
  # Internal conditional to determine if perform or not the PCs retaining
  if (PC.retain == T) {
    
    message("Principal components retaining analyses in process")
    
  } else {
    
    message("NOT doing the PCs retaining analyses")
    
    }
  
} else {
  
  message("Selected data option: Genotypic")
  
  # Conditional within the genotypic data to determine if genofile is empty or necessary to read from dir
  if (is.null(genofile)) {
    
    message("Genofile is NULL\n\nReading vcf and gds files\n\nCreating genofile")
    
    snpgdsVCF2GDS(vcf, gds, ignore.chr.prefix = "chromosome")
    genofile <- snpgdsOpen(gds)
    sample_id <- read.gdsn(index.gdsn(genofile, "sample.id"))
    
  } else { 
    
    message("Genofile is not NULL\n\nPCA function continues regularly")
    
  }
  
  # Continue function
  # Performs the SNPRelate's PCA
  PCA <- snpgdsPCA(genofile, autosome.only = T, remove.monosnp = T, need.genmat = T,
                   algorithm = "exact", eigen.method = "DSPEVX")
  
  # Creates the table to draw scree plots and analyze the number of PCs to retain 
  PC <- data.frame(PC = 1:length(PCA$varprop), id = PCA[["sample.id"]],
                   eigen = PCA[["eigenval"]], var = PCA$varprop * 100) %>%
    mutate(var.cum = cumsum(var))
  
  # PCs retaining part
  if (PC.retain == T) {
    
    message("Principal components retaining analyses in process")
    
    # Establish a cutoff and obtain the genetic correlation matrix
    cutoff <- 1 / length(PC$PC) * 100
    cor.m <- PCA[['genmat']]
    cor.m <- as.matrix(replace(cor.m, cor.m > 1, 1))
    
    # Makes the parallel, VSS, Velicer's MAP, and BIC analyses
    FA.p <- fa.parallel(cor.m, n.obs = dim(cor.m)[1], fm = "pa", fa = "pc", n.iter = 30, plot = F)
    FAC <- nfactors(cor.m, diagonal = T, fm = "pa", n.obs = dim(cor.m)[1])
    
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
  
  }






# Scree plots
ggplot(PC, aes(x = PC, y = eigen)) +
  geom_point(shape = 1) +
  geom_point(data = dplyr::filter(PC, eigen > 1), color = "black") +
  geom_hline(yintercept = 1) +
  labs(y = "Eigen values", x = "PC") +
  theme_classic() +
  theme(axis.title = element_text(size = 14, color = "black"),
        axis.text = element_text(size = 14, color = "black"),
        axis.title.y = element_text(margin = margin(r = 10)),
        axis.title.x = element_text(margin = margin(t = 10))) +
  scale_y_continuous(breaks = c(0, 5, 10, 15, 20)) +
  scale_x_continuous(breaks = c(0, 25, 50, 75, 100, 125, 150))

ggplot(PC, aes(x = PC, y = var)) +
  geom_point(shape = 1) +
  geom_point(data = dplyr::filter(PC, var.cum < 80), aes(x = PC, y = var), color = "black") +
  geom_point(data = dplyr::filter(PC, var.cum < 70), aes(x = PC, y = var), color = "red") +
  geom_hline(yintercept = cutoff) +
  labs(y = "Explained variance (%)", x = "PC") +
  theme_classic() +
  theme(axis.title = element_text(size = 14, color = "black"),
        axis.text = element_text(size = 14, color = "black"),
        axis.title.y = element_text(margin = margin(r = 10)),
        axis.title.x = element_text(margin = margin(t = 10))) +
  scale_y_continuous(breaks = c(0, 2, 4, 6, 8, 10, 12)) +
  scale_x_continuous(breaks = c(0, 25, 50, 75, 100, 125, 150))

ggplot(PC, aes(x = PC, y = var.cum)) +
  geom_point(shape = 1) +
  geom_point(data = dplyr::filter(PC, var.cum < 80), aes(x = PC, y = var.cum), color = "black") +
  geom_point(data = dplyr::filter(PC, var.cum < 70), aes(x = PC, y = var.cum), color = "red") +
  geom_hline(yintercept = max(dplyr::filter(PC, var > cutoff)$var.cum)) +
  labs(y = "Cumulative variance (%)", x = "PC") +
  theme_classic() +
  theme(axis.title = element_text(size = 14, color = "black"),
        axis.text = element_text(size = 14, color = "black"),
        axis.title.y = element_text(margin = margin(r = 10)),
        axis.title.x = element_text(margin = margin(t = 10))) +
  scale_y_continuous(breaks = c(0, 20, 40, 60, 80, 100)) +
  scale_x_continuous(breaks = c(0, 25, 50, 75, 100, 125, 150))



# Table to plot the PCA
tab <- data.frame(sample.id = PCA$sample.id, stringsAsFactors = F,
                  EV1 = PCA$eigenvect[,1], EV2 = PCA$eigenvect[,2]) # The first two eigenvectors
dt <- merge(tab, labels, by.x = 'sample.id', by.y = 'Taxa')

# Plot the PCA
plot_ly(data = dt, x = ~ EV1, y = ~ EV2, color = ~ label, type = 'scatter', mode = 'markers',
        symbol = ~ label, symbols = c('circle', 'x'), text = ~ sample.id, marker = list(size = 6)) %>%
  layout(xaxis = list(title = paste('Dimension 1 - ', round(PC$var)[1], '%')),
         yaxis = list(title = paste('Dimension 2 - ', round(PC$var)[2], '%')))

# Save the plot
saveWidget(fig, "GWAS_PCA.html", selfcontained = F, libdir = "lib")
saveWidget(as_widget(fig), "GWAS_PCA.html")



# 4: Discriminant analysis of principal components (DAPC) ----------------------

vcf <- read.vcfR('GWAS_PPD.snps.filter_info.missing_0.10.imputation.vcf.gz', verbose = T)
my_genind <- vcfR2genind(vcf)
my_genind

DAPC <- function(vcf){}

# Load packages and data
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

pdf("GWAS_DAPC_scatter.pdf", width = 10, height = 8) 
scatter(dapc, scree.pca = F, ratio.pca = 0.3, pch = 20, cell = 1, solid = 0.6, cex = 2.5, clab = 0,
        scree.da = F, leg = T, txt.leg = paste("Cluster", 1:3))
dev.off()

set.seed(4)
contrib <- loadingplot(dapc$var.contr, axis = 2, lab.jitter = 1)

groups <- tibble::rownames_to_column(as.data.frame(dapc[["grp"]]), var = "Name")
colnames(groups) <- c('Name','group_dapc')

write.csv(groups, 'GWAS_PPD.dapc_groups.csv', row.names = F, col.names = T)
