# Short name: Multivariate methods for GWAS covariates
# Description: It can perform several multivariate methods (MDS, PCA, DAPC) to identify covariates for GWAS 
# Output(s): It depends on the chosen method. For example:
#            - MDS:
#            - PCA:
#            - DAPC:
#
# Authors: Camilo E. Sanchez (c.e.sanchez@cgiar.org) and Vianey Barrera-Enriquez (vpbarrera@gmail.com)
#
# Arguments: 
# 



# 1: Load packages and data ----------------------------------------------------
if (!require(SNPRelate)) install.packages(SNPRelate)
if (!require(gdsfmt)) install.packages(gdsfmt)
if (!require(ggplot2)) install.packages(ggplot2ggplot2)
if (!require(ggtree)) install.packages(ggtree)
if (!require(ape)) install.packages(ape)
if (!require(plotly)) install.packages(plotly)
if (!require(tidyverse)) install.packages(tidyverse)

library(SNPRelate)
library(gdsfmt)
library(ggplot2)
library(ggtree)
library(ape)
library(plotly)
library(tidyverse)
library(psych)



####### To do ####### 
# 1: Transform into functions structure for each method
# 2: Interactive graphs (using plotly and saving them as html format)
# 3: For DAPC: Customizable labels
# 4: Finish the script script on the head
# 5: Add the NDMS and MDS parts

annotation <- read.csv('GWAS_PPD.labels.csv')
vcf <- 'GWAS_PPD.snps.filter_info.missing_0.10.imputation.vcf.gz'
gds <- 'GWAS_PPD.gds'
snpgdsVCF2GDS(vcf, gds, ignore.chr.prefix = "chromosome")
genofile <- snpgdsOpen(gds)
sample_id <- read.gdsn(index.gdsn(genofile, "sample.id"))



# 3: PCA -----------------------------------------------------------------------
pca <- snpgdsPCA(genofile, autosome.only = T, remove.monosnp = T, algorithm = "exact",
                 eigen.method = "DSPEVX", need.genmat = T)

# Creates the table to to analyze the number of PCs to retain as well as to draw scree plots
PC <- data.frame(PC = 1:length(pca$varprop), id = pca[["sample.id"]],
                 eigen = pca[["eigenval"]], var = pca$varprop * 100) %>%
  mutate(var.cum = cumsum(var))

# Establish a cutoff
cutoff <- 1 / length(PC$PC) * 100

# Obtain the genetic correlation matrix
cor.m <- pca[['genmat']]
cor.m <- as.matrix(replace(cor.m, cor.m > 1, 1))

# Makes the parallel, VSS, Velicer's MAP, and BIC analyses
fa.p <- fa.parallel(cor.m, n.obs = dim(cor.m)[1], fm = "pa", fa = "pc", n.iter = 30, plot = F)
fac <- nfactors(cor.m, diagonal = T, fm = "pa", n.obs = dim(cor.m)[1])

# Obtain the summary statistics
fac.t <- fac[["vss.stats"]] %>% mutate(map = fac$map)

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
              fa.p$ncomp, "PCs \n",
              "(6) Apply the 'Very Simple Structure criteria', you would retain",
              which.max(fac.t$cfit.2), "PCs \n",
              "(7) Try with the 'Velicer's MAP criteria', you would retain",
              which.min(na.omit(fac.t$map)), "PCs \n\n",
              "Also, you could retain a fixed number of PCs or retain PCs based on the scree plot"))



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
tab <- data.frame(sample.id = pca$sample.id, stringsAsFactors = F,
                  EV1 = pca$eigenvect[,1], EV2 = pca$eigenvect[,2]) # the first two eigenvectors

dt <-  merge(tab, annotation, by.x = 'sample.id', by.y = 'Taxa')

# Plot the PCA
fig <- dt %>% plot_ly(x = ~ EV1, y = ~ EV2, color = ~ label, type = 'scatter', mode = 'markers', 
                      symbol = ~ label, symbols = c('circle', 'x'), text = ~ sample.id, 
                      marker = list(size = 6))
fig <- fig %>% layout(title = 'PCA', 
                      xaxis = list(title = paste('Dimension 1 - ', round(PC$var)[1], '%')), 
                      yaxis = list(title = paste('Dimension 2 - ', round(PC$var)[2], '%')))
fig

# Save the plot
saveWidget(fig, "GWAS_PCA.html", selfcontained = F, libdir = "lib")
htmlwidgets::saveWidget(as_widget(fig), "GWAS_PCA.html")





# 4: DAPC ----------------------------------------------------------------------
# Load packages and data
if (!require(adegenet)) install.packages(adegenet)
if (!require(vcfR)) install.packages(vcfR)
if (!require(tibble)) install.packages(tibble)

library(adegenet)
library(vcfR)
library(tibble)

vcf <- read.vcfR('GWAS_PPD.snps.filter_info.missing_0.10.imputation.vcf.gz', verbose = T)
my_genind <- vcfR2genind(vcf)
my_genind


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
