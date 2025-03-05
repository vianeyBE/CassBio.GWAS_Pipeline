# Short name: Manhattan plot for GWAS
# Description: Make a customizable Manhattan plot the GWAS GAPIT results
# Output: Manhattan plots for one, several, or all traits
#
# Authors: Camilo E. Sanchez (c.e.sanchez@cgiar.org) and Vianey Barrera-Enriquez (vpbarrera@gmail.com)
#
# Arguments:
# Mdir: Name of the directory that contains the GAPIT results. For example: home/user/folder.
# pat: Enter the path of file names to look for. For example: GAPIT.Association.GWAS_Results.
# mod: Enter the model(s) of interest. Options: BLINK, GLM, MLM, FarmCPU.
# wtd: How many traits do you want to plot (Options: One, Several, All).
# colors: (Optional) Colors of the chromosomes in Manhattan plots.
#         If you want to change the colors, provide 2 or more.
#         It can be color names o hex codes (Default: grey and skyblue).



###### To do ######
# Everything is good!



# 0: Function init -------------------------------------------------------------

Manhattan <- function(Mdir, pat, mod, wtd, colors = c("grey", "skyblue")){
  
  # 1: Load packages -----------------------------------------------------------
  
  if (!require(tidyverse)) install.packages(tidyverse)
  if (!require(ggtext)) install.packages(ggtext)
  if (!require(plotly)) install.packages(plotly)
  if (!require(htmlwidgets)) install.packages(htmlwidgets)
  
  library(tidyverse)
  library(ggtext)
  library(plotly)
  library(htmlwidgets)
  
  
  
  # 2: Find all the CSVs with the results and organize them --------------------
  
  # Informative message
  message("Getting list of CSV files...")
  
  # Get the names of the files
  setwd(Mdir)
  names <- list.files(path = Mdir, pattern = paste0(pat, "."), all.files = F, full.names = F, recursive = T)
  
  # Informative messages
  message(paste0("GWAS files found: ", length(names), "\n\nReading GWAS files..."))
  
  # Create an empty list to store the data the variable to iterate
  csv.l <- list()
  
  # Loop to find, read, rename, modify and save all the csv files in the list
  for (i in 1:length(names)){
    
    # Read and modify each one of the csvs
    dframe <- read.csv(paste0(Mdir, "/", names[i])) %>%
      mutate(rename = paste0("/", names[i])) %>%
      tidyr::separate(col = rename, into = c("batch", "data"), sep = paste0(pat, "."), extra = "drop")
    
    # Conditional given the differences in FarmCPU models
    if (gdata::startsWith(dframe[1,10], "Far")){
      
      # Modified the database to obtain the model and the trait
      dframe <- dframe %>% tidyr::separate(col = data, into = c("na", "na2", "trait"),
                                           sep = ".m", extra = "drop") %>%
        tidyr::separate(col = trait, into = c("trait", "na3"), sep = ".csv", extra = "drop") %>%
        mutate(model = "FarmCPU") %>%
        select(SNP, Chr, Pos, P.value, MAF, nobs, Effect, model, trait)
      
    } else {
      
      # Modified the database to obtain the model and the trait
      dframe <- dframe %>% tidyr::separate(col = data, into = c("model", "trait"), 
                                           sep = ".m", extra = "drop") %>%
        tidyr::separate(col = trait, into = c("trait", "na"), sep = ".csv", extra = "drop") %>%
        select(SNP, Chr, Pos, P.value, MAF, nobs, Effect, model, trait)
      
    }
    
    # Save the dataframe in the list
    csv.l[[i]] <- dframe
    
    # Progress bar
    cat('\r', i, ' files processed |', rep('=', i / 4), ifelse(i == length(names), '|\n', '>'), sep = '')
    
  }
  
  # Informative message
  message("All GWAS files were processed")
  
  # Merge the data frames of the list in a single data frame
  GWAS <- bind_rows(csv.l)
  GWAS <- GWAS %>% na.omit() %>% dplyr::filter(model %in% mod)
  GWAS$Pos<- as.numeric(GWAS$Pos)
  
  # Informative message
  message(paste("There are:\n",
                dim(GWAS)[1], "SNPs to plot\n",
                length(table(GWAS$trait)), "traits\n",
                length(table(GWAS$model)), "different models (specified in the 'mod' object)"))
  
  
  
  # 3: Preparing and plotting the data -----------------------------------------
  
  if (wtd == "One"){
    
    # Prompt to receive the trait
    tr <- readline(prompt = message(paste(unique(GWAS$trait), collapse = "\n"), "\n\n",
                                    "Choose a trait from the list above and paste its name (exact match): "))
    
    # Filter the data by trait
    data <- GWAS %>% dplyr::filter(trait == tr)
    
    # Informative message
    message(paste("You chose the following trait:", tr, "\n\n",
                  "For that trait, there are:\n",
                  dim(data)[1], "SNPs to plot\n"))
    
    # Prepare the data to plot
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
      mutate(center = ((max - min) / 2) + bp_add)
    
    # Informative message
    message("Making the plots. It can take a few seconds. Please be patient.")
    
    # Manhattan plot
    fig_man <- ggplot(data, aes(bp_cum, -log10(P.value), color = as_factor(Chr),
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
    
    # CM plot
    fig_cm <- ggplot(data, aes(bp_cum, Effect, color = as_factor(Chr), size = Effect,
                               shape = model)) +
      geom_point(alpha = 0.75) +
      geom_hline(yintercept = 0, color = "black", linetype = "dashed") +
      scale_color_manual(values = rep(c("grey", "skyblue"), 18)) +
      scale_y_continuous(expand = c(0.1, 0.1), limits = c(min(data$Effect), max(data$Effect))) +
      scale_x_continuous(expand = c(0, 0), label = axis_set$Chr, breaks = axis_set$center) +
      guides(color = "none", size = "none") +
      labs(y = "SNP effect", x = "Chromosome", title = paste("Trait: ", tr), shape = "Model") +
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
    
    if (wtd == "Several"){
      
      # Prompt to receive the traits
      message(paste(unique(GWAS$trait), collapse = "\n"), "\n\n",
              "Choose the traits from the list above and paste their names in the next function\n\n",
              "Finish with two tabs")
      
      tr <- scan(what = character())
      
      # Filter the data by trait
      data <- GWAS %>% dplyr::filter(trait %in% tr)
      
      # Informative message
      message(paste("You chose the following traits:", paste(tr, collapse = "  ")))
      
      # Prepare the data to plot
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
      
      # Informative message
      message("Making the plots. It can take a few seconds. Please be patient.")
      
      # Manhattan plot
      fig_man <- ggplot(data, aes(bp_cum, -log10(P.value), color = as_factor(Chr),
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
      
      # CM plot
      fig_cm <- ggplot(data, aes(bp_cum, Effect, color = as_factor(Chr), size = Effect,
                                 shape = model)) +
        geom_point(alpha = 0.75) +
        geom_hline(yintercept = 0, color = "black", linetype = "dashed") +
        scale_color_manual(values = rep(c("grey", "skyblue"), 18)) +
        scale_y_continuous(expand = c(0.1, 0.1), limits = c(min(data$Effect), max(data$Effect))) +
        scale_x_continuous(expand = c(0, 0), label = axis_set$Chr, breaks = axis_set$center) +
        guides(color = "none", size = "none") +
        labs(y = "SNP effect", x = "Chromosome", title = paste("Trait: ", tr), shape = "Model") +
        theme_classic() +
        facet_wrap(~ trait, ncol = 1) +
        theme(legend.position = "bottom", axis.title.y = element_markdown(),
              panel.grid.major.x = element_blank(), panel.grid.minor.x = element_blank(),
              legend.title = element_text(colour = "black", size = 12),
              legend.text = element_text(colour = "black", size = 12),
              plot.title = element_text(size = 16, color = "black"),
              axis.text = element_text(size = 12, color = "black"),
              axis.title = element_text(size = 14, color = "black"),
              axis.title.x = element_text(margin = margin(t = 10)))
      
    } else {
      
      if (wtd == "All"){
        
        # Copy the data
        data <- GWAS
        
        # Informative message
        message(paste(dim(data)[1], "SNPs and", length(unique(data$trait)), "traits will be plotted"))
        
        # Prepate the data to plot
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
        
        # Informative message
        message("Making the plots. It can take a few seconds. Please be patient.")
        
        # Manhattan plot
        fig_man <- ggplot(data, aes(bp_cum, -log10(P.value), color = as_factor(Chr),
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
        
        # CM plot
        fig_cm <- ggplot(data, aes(bp_cum, Effect, color = as_factor(Chr), size = Effect,
                                   shape = model)) +
          geom_point(alpha = 0.75) +
          geom_hline(yintercept = 0, color = "black", linetype = "dashed") +
          scale_color_manual(values = rep(c("grey", "skyblue"), 18)) +
          scale_y_continuous(expand = c(0.1, 0.1), limits = c(min(data$Effect), max(data$Effect))) +
          scale_x_continuous(expand = c(0, 0), label = axis_set$Chr, breaks = axis_set$center) +
          guides(color = "none", size = "none") +
          labs(y = "SNP effect", x = "Chromosome", title = paste("Trait: ", tr), shape = "Model") +
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
        
        # Inform about an error in the argument
        message("You can only choose the following options: 'One', 'Several', or 'All' (exact match)")
        
      }
    }
  }
  
  
  
  # 4: Saving the plots  -------------------------------------------------------
  
  # Transform the plots in plotly objects
  fig_man <- ggplotly(fig_man)
  fig_cm <- ggplotly(fig_cm)
  
  # Save plots
  saveWidget(fig_man, "GWAS_Manhattan.html", selfcontained = F, libdir = "lib_man")
  saveWidget(as_widget(fig_man), "GWAS_Manhattan.html")
  
  saveWidget(fig_cm, "GWAS_CM.html", selfcontained = F, libdir = "lib_cm")
  saveWidget(as_widget(fig_cm), "GWAS_CM.html")
  
  
  
  # 5: Function ends -----------------------------------------------------------
  
  
  
  message("Done!")
  
}



# Example(s) -------------------------------------------------------------------
# Set arguments
# Mdir <- "D:/OneDrive - CGIAR/Cassava_Bioinformatics_Team/01_ACWP_F1_Metabolomics/07_GWAS"
# pat <- "GAPIT.Association.GWAS_Results"
# mod <- c("BLINK", "FarmCPU", "MLM")
# wtd <- "One"



# Run function -----------------------------------------------------------------
# Manhattan(Mdir, pat, mod, wtd)
