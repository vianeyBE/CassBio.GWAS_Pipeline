# 3. For PCA:
# dir: Directory where data is located.
# phenofile: A database of genotypes/individuals (rows) and their traits (columns). First column must be genotypes/individuals names.
# genofile: An object of class SNPGDSFileClass (GDS file read with the `snpgdsOpen` function from `SNPRelate` package).
# labels: When provide a `genofile` and `groups` argument is `TRUE`, please provide a dataframe with genotypes/individuals in the first column and groups/treatment data in the second column.
# gds: A Genomic Data Structures (GDS) file (a reformatted VCF file with the `snpgdsVCF2GDS` function from `vcfR` package).
# vcf: A Variant Call Format (VCF) file containing DNA polymorphism data such as SNPs, insertions, deletions and structural variants, together with rich annotations.
# type: A character string indicating wheter data is phenotypic ("pheno") or genotypic ("geno") (default = Pheno).
# groups: Boolean value indicating whether the data includes different treatments/groups (default = F). If `TRUE` is selected then the second column of the `phenofile` must be the groups/treatments.
# PC.retain: Boolean value indicating whether to analyze how many PCs retain (default = F).

# 3: Principal component analysis (PCA) ----------------------------------------

PCA <- function(dir, phenofile, genofile, labelfile, gds, vcf, type, PC.retain, prefixVCF, output, num.thread){
  
  # Set working directory
  setwd(dir)
  
  #
  if (!require(SNPRelate)) install.packages(SNPRelate)
  if (!require(gdsfmt)) install.packages(gdsfmt)
  if (!require(tidyverse)) install.packages(tidyverse)
  if (!require(psych)) install.packages(psych)
  if (!require(plotly)) install.packages(plotly)
  if (!require(htmlwidgets)) install.packages(htmlwidgets)
  
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
      
      # Plot the PCA
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
    
    # Informative message
    message("Selected data option: Genotypic")
    
    # Loading files ------------------------------------------------------------
    
    # If vcf is provided
    if (!is.null(vcf)){
      
      message("VCF file provided...\n\n", 
              "Reading VCF file...\n\n",
              "Reformatting it to a GDS file and then transforming it to a SNPGDSFileClass file")
      
      # Paste the path
      gds <- paste0(dir, output, ".gds")
      
      # Read and transform the files
      snpgdsVCF2GDS(vcf, gds, ignore.chr.prefix = prefixVCF, verbose = F)
      genofile <- snpgdsOpen(gds)
      sample_id <- read.gdsn(index.gdsn(genofile, "sample.id"))
      
    } else {
      
      # If gds is provided
      if (!is.null(gds)){
        
        message("GDS file provided...\n\n", 
                "Reading CDS file...\n\n",
                "Transforming it to a SNPGDSFileClass file")
        
        # Read and transform the files
        genofile <- snpgdsOpen(gds)
        sample_id <- read.gdsn(index.gdsn(genofile, "sample.id"))
        
      } else {
        
        # If genofile is provided (in SNPGDSFileClass format)
        if (!is.null(genofile)){
          
          message("Genofile (SNPGDSFileClass file) is provided...\n\n",
                  "PCA function continues regularly")
          
        }
      }
    }
    
    if ( !is.null(labelfile)){
      
      message("Label file provided by the user ")
      groups <- T
      labels <- read.csv(labelfile)
      
    } else {
      
      message("No labels provided ")
      groups <- F
      
    }
    
    message("Files loaded successfully!")
    
    
    # PCA - SNPRelate ----------------------------------------------------------
    PCA <- snpgdsPCA(genofile, autosome.only = T, remove.monosnp = T, need.genmat = T,
                     algorithm = "exact", eigen.method = "DSPEVX", 
                     num.thread = num.thread, verbose = F)
    
    # Saving PCs
    PC <- data.frame(PC = 1:length(PCA$varprop), id = PCA[["sample.id"]],
                     eigen = PCA[["eigenval"]], var = PCA$varprop * 100) %>%
      mutate(var.cum = cumsum(var))
    
    write.csv(PCA[["eigenval"]], paste0(output, '.PC_SNPrelated.csv'), quote = F)
    
    
    # PCs retaining analysis ---------------------------------------------------
    
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
      
      message("Not doing the PCs retaining analysis")
      
    }
    
    
    # Plots --------------------------------------------------------------------
    
    # Scree plot - Individual variance
    ind_var <- plot_ly(data = PC, x = ~ PC, y = ~ var, type = "scatter", 
                       mode = "lines+markers", text = ~ round(var, 2),
                       line = list(color = "grey")) %>% 
      layout(xaxis = list(title = "PC"), yaxis = list(title = "Explained variance (%)"))
    
    htmlwidgets::saveWidget(as_widget(ind_var), paste0(output, ".ind_variance.html"))
    
    # Scree plot -  Cumulative variance
    cum_var <- plot_ly(data = PC, x = ~ PC, y = ~ var.cum, type = "scatter",
                       mode = "lines+markers", text = ~ round(var.cum, 2), 
                       line = list(color = "grey")) %>%
      layout(xaxis = list(title = "PC"), yaxis = list(title = "Cumulative variance (%)"))
    
    htmlwidgets::saveWidget(as_widget(cum_var), paste0(output, ".cum_variance.html"))
    
    # Table to plot the PCA
    tab <- data.frame(sample.id = PCA$sample.id, stringsAsFactors = F,
                      PC1 = PCA$eigenvect[,1], PC2 = PCA$eigenvect[,2])
    
    # PCA - Scatter Plots PC1 vs PC2
    # Adding Groups 
    if (groups == T){
      
      # Labeling
      labels <- labels %>% dplyr::rename(sample.id = 1, groups = 2)
      dt <- tab %>% inner_join(labels, by = "sample.id")
      
      # Plotting
      fig <- plot_ly(data = dt, x = ~ PC1, y = ~ PC2, 
                     color = ~ as.factor(groups), #symbol = ~ as.factor(groups),
                     type = "scatter", mode = "markers", 
                     text = ~ sample.id, marker = list(size = 6)) %>%
        layout(legend = list(title = list(text = "<b> Groups </b>"), orientation = "h"),
               xaxis = list(title = paste("Dimension 1 - ", round(PC$var)[1], "%")),
               yaxis = list(title = paste("Dimension 2 - ", round(PC$var)[2], "%")))
      
    } else { # Without Groups
      
      # Plotting
      fig <- plot_ly(data = tab, x = ~ PC1, y = ~ PC2, 
                     type = "scatter", mode = "markers", 
                     text = ~ sample.id, marker = list(size = 6)) %>%
        layout(xaxis = list(title = paste("Dimension 1 - ", round(PC$var)[1], "%")),
               yaxis = list(title = paste("Dimension 2 - ", round(PC$var)[2], "%")))
      
    }
    
    htmlwidgets::saveWidget(as_widget(fig), paste0(output, ".PCA_SNPRelated.html"))
    
  }
  
}



# 3.1: PCA example(s) ----------------------------------------------------------
# Genotypic without groups
# dir <- "D:/OneDrive - CGIAR/00_CassavaBioinformaticsPlatform/03_PPD/00_Data/00_Geno/"
# labelfile <- NULL
# vcf <- "chrs_PPD_v6_filterGATK.vcf"
# type <- "Geno"
# groups <- F
# output <- "chrs_PPD_v6_filterGATK"
# PC.retain = F
# prefixVCF = "Chromosome"
# num.thread = 4


# 3.2: Run PCA function --------------------------------------------------------
# PCA(dir, phenofile, genofile, gds, vcf, type = "Geno", groups = T, PC.retain = F)