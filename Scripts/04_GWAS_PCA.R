# Short name: PCA
# Description: Principal component analysis for phenotypic and genotypic data
# Output: It saves the results of PCA and its plots in the working folder
#
# Author: Camilo E. Sánchez-Sarria (c.e.sanchez@cgiar.org)
#
# Arguments:
# dir: Directory where data is located
# data: phenofile -> A CSV of genotypes (rows) and traits (columns)
#       genofile -> A Variant Call Format (VCF) file
# labels: CSV with genotypes in first column and grouping data in the second
# ellipses: If TRUE -> Draws ellipses on groups
# PC.retain: If TRUE -> Analyzes how many PCs retain (default = F)




##### To do #####
# 1. Change Phenotypic PCA accordingly




# 0: Function init -------------------------------------------------------------
PCA <- function(dir, data, labels, ellipses, PC.retain){
  
  
  
  
  # 1: Load packages and data --------------------------------------------------
  # Install packages if needed
  if (!require("BiocManager", quietly = T)) install.packages("BiocManager")
  if (!require(SNPRelate)) BiocManager::install("SNPRelate")
  if (!require(tidyverse)) install.packages("tidyverse")
  if (!require(psych)) install.packages("psych")
  if (!require(tools)) install.packages("tools")
  
  # Load libraries
  library(SNPRelate)
  library(tidyverse)
  library(psych)
  library(tools)
  
  
  
  
  # 2: Determines which path continue in the analysis -------------------------- 
  # Set working directory
  setwd(dir)
  
  # Get extension for file name
  ext <- file_ext(paste0(data))
  
  # Conditional to determine if the data is genotypic of phenotypic
  if (ext == "csv"){
    
    
    
    
    # 2.1: Phenotypic PCA ------------------------------------------------------
    # Informative message
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
    
    
    
    
    
    # 2.2: Genotypic PCA -------------------------------------------------------
    # 2.3: Read VCF and transform the data -------------------------------------
    # Informative message
    message(paste0("Selected data option: Genotypic\n\n",
                   "Reading VCF file: ", data, "\n"))
    
    # Reading VCF
    prefixVCF <- readline("Enter the VCF chromosome prefix
                          If no prefix, skip, just enter: ")
    
    # Informative message
    message("Reformatting VCF to a GDS file...\n")
    
    # Extract the prefix
    prefix <- file_path_sans_ext(data)
    
    # Generate a path to store the gds file
    gds <- paste0(dir, prefix, ".gds")
      
    # Transform the VCF to a gds depending on the chr prefix
    if (prefixVCF == "") {
      
      snpgdsVCF2GDS(data, gds, verbose = T)
      
    } else {
      
      snpgdsVCF2GDS(data, gds, ignore.chr.prefix = prefixVCF, verbose = T)
      
    }
    
    # Read the gds file and its index
    genofile <- snpgdsOpen(gds)
    sample_id <- read.gdsn(index.gdsn(genofile, "sample.id"))
    
    
    
    
    # 2.4: Perform the PCA -----------------------------------------------------
    message("Doing PCA...\n")
    PCA <- snpgdsPCA(genofile, autosome.only = F, remove.monosnp = F, 
                     need.genmat = T, algorithm = "exact", 
                     eigen.method = "DSPEVX", num.thread = 4, verbose = F)
    
    # Saving PCs data
    PC <- data.frame(PC = 1:length(PCA$varprop), id = PCA[["sample.id"]],
                     eigen = PCA[["eigenval"]], var = PCA$varprop * 100) %>%
      mutate(var.cum = cumsum(var))
    
    # Table to plot the PCA
    tab <- data.frame(sample.id = PCA$sample.id, stringsAsFactors = F,
                      PC1 = PCA$eigenvect[, 1],PC2 = PCA$eigenvect[, 2])
    
    # Save the output dataframe
    message("Saving PCA output...\n")
    write.csv(PC, file = paste0("PCA_results_", prefix, ".csv"), row.names = F)
    dev.off()
    
    
    
    
    # 2.5: PCs retaining analysis ----------------------------------------------
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
    
    
    
    
    # 2.6: Plots -------------------------------------------------------------
    # PCA plots
    message("Plotting...\n")
    
    # Scree plot - Individual variance
    ind_var <- PC %>% dplyr::filter(!is.nan(var)) %>%
      ggplot(aes(x = as.factor(PC), y = var)) +
      geom_bar(stat = "identity", width = 0.5, fill = "grey", color = "grey") +
      labs(x = "PCs", y = "Explained variance (%)") +
      theme_classic() +
      theme(axis.title = element_text(color = "black", size = 14),
            axis.text  = element_text(color = "black", size = 14))
    
    # Scree plot -  Cumulative variance
    cum_var <- PC %>% dplyr::filter(!is.nan(var)) %>%
      ggplot(aes(x = as.factor(PC), y = var.cum)) +
      geom_bar(stat = "identity", width = 0.5, fill = "grey", color = "grey") +
      labs(x = "PCs", y = "Cumulative variance (%)") +
      theme_classic() +
      theme(axis.title = element_text(color = "black", size = 14),
            axis.text  = element_text(color = "black", size = 14))
    
    
    # Simple PCA: Without groups and without ellipses
    pca_s <- ggplot(tab, aes(x = PC1, y = PC2, label = sample.id)) +
      geom_point(size = 1.5, color = "#91bfdb") +
      labs(x = paste("PC1 (", round(PC$var[1], 2), "%)"),
           y = paste("PC2 (", round(PC$var[2], 2), "%)")) +
      geom_vline(xintercept = 0, linetype = "dashed") +
      geom_hline(yintercept = 0, linetype = "dashed") +
      theme_bw() +
      theme(axis.title = element_text(color = "black", size = 16),
            axis.text = element_text(color = "black", size = 16),
            panel.grid = element_blank())
    
    # If groups
    if (!is.null(labels) == TRUE) {
      
      # Informative message
      message("Label file provided by the user...\n")
      
      # Read csv file
      labels <- read.csv(labels)
      message("Files loaded successfully!\n")
      
      # Labeling
      labels <- labels %>% dplyr::rename(sample.id = 1, groups = 2)
      dt <- tab %>% inner_join(labels, by = "sample.id")
        
      # Define personalized colors and shapes
      n_groups <- length(unique(dt$groups))
      shapes <- rep(0:25, length.out = n_groups)
      colors <- scales::hue_pal()(n_groups)
      pal <- stats::setNames(colors, sort(unique(dt$groups)))
        
    } else {
      
      # Informative message
      message("No labels provided...\n")
      
    }
    
    # PCA plots with groups
    if (ellipses == TRUE) {
    
    # Plotting with groups and ellipses
    #dt <- dt %>% dplyr::filter(groups != "Unassigned")
    pca_e <- ggplot(dt, aes(x = PC1, y = PC2, color = as.factor(groups), 
                            shape = as.factor(groups))) +
      geom_point(size = 2) +
      #stat_ellipse(aes(group = as.factor(groups), fill = as.factor(groups)),
      #             type = "norm", level = 0.95, geom = "polygon", 
      #             alpha = 0.05, linewidth = 0) +
      stat_ellipse(aes(group = as.factor(groups)), type = "norm", 
                   level = 0.95, linewidth = 0.3) +
      scale_shape_manual(values = shapes) +
      scale_color_manual(values = colors) +
      scale_fill_manual(values  = colors) +
      geom_hline(yintercept = 0, linetype = "dashed") +
      geom_vline(xintercept = 0, linetype = "dashed") +
      labs(x = paste0("PC1 (", round(PC$var[1], 2), " %)"),
           y = paste0("PC2 (", round(PC$var[2], 2), " %)")) +
      theme_bw() +
      theme(legend.title = element_blank(), legend.position = "bottom",
            legend.text = element_text(color = "black", size = 14),
            axis.title = element_text(color = "black", size = 14),
            axis.text = element_text(color = "black", size = 14),
            legend.key.size = unit(1.2, "lines"), 
            panel.grid = element_blank(), legend.box = "horizontal") +
      guides(color = guide_legend(nrow = 1, override.aes = list(size = 4)),
             shape = guide_legend(nrow = 1, override.aes = list(size = 4)))
    
    } else {
    
    # Plotting with groups but without ellipses
    #dt <- dt %>% dplyr::filter(groups != "Unassigned")
    pca_g <- ggplot(dt, aes(x = PC1, y = PC2, color = as.factor(groups), 
                            shape = as.factor(groups))) +
      geom_point(size = 2) +
      scale_shape_manual(values = shapes) +
      scale_color_manual(values = colors) +
      geom_hline(yintercept = 0, linetype = "dashed") +
      geom_vline(xintercept = 0, linetype = "dashed") +
      labs(x = paste0("PC1 (", round(PC$var[1], 2), " %)"),
           y = paste0("PC2 (", round(PC$var[2], 2), " %)")) +
      theme_bw() +
      theme(legend.title = element_blank(), legend.position = "bottom",
            legend.text = element_text(color = "black", size = 14),
            axis.title = element_text(color = "black", size = 14),
            axis.text = element_text(color = "black", size = 14),
            legend.key.size = unit(1.2, "lines"), 
            panel.grid = element_blank(), legend.box = "horizontal") +
      guides(color = guide_legend(nrow = 1, override.aes = list(size = 4)),
             shape = guide_legend(nrow = 1, override.aes = list(size = 4)))
        
    }
    
    
    
    
    # 2.2.5: Save the plots ----------------------------------------------------
    # Save the PCA with high resolution (600 dpi)
    pdf(paste0("PCA_", prefix, ".pdf"), onefile = T, width = 10, height = 6)
    
    print(ind_var)
    print(cum_var)
    print(pca_s)
    print(pca_g)
    print(pca_e)
    
    dev.off()
    
  }
  
}



# 3.1: PCA example(s) ----------------------------------------------------------
# Examples
# dir <- "D:/OneDrive - CGIAR/00_BioInf_Platform/09_DiversityPanel/03_Paper/"
# data <- "Diversity.GATK.vcf"
# labels <- "Labels.csv"
# ellipses <- T
# PC.retain <- F



# 3.2: Run PCA function --------------------------------------------------------
# PCA(dir, data, labels, ellipses, PC.retain)
 