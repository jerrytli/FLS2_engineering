#-----------------------------------------------------------------------------------------------
# Coaker Lab - Plant Pathology Department UC Davis
# Author: Danielle M. Stevens
# Last Updated: 08/21/2024
# Script Purpose: Run the R scripts for FLS2-flagellin engineering paper
# Inputs: N/A
# Outputs: N/A
#-----------------------------------------------------------------------------------------------

######################################################################
# set path to data
######################################################################

#setwd to where repo was cloned and maintained
library(rstudioapi)
setwd(dirname(rstudioapi::getSourceEditorContext()$path))

######################################################################
# library packages need to load
######################################################################

# packages for plotting phylogenetic tree
library(ggtree)
library(ggtreeExtra)
library(phangorn)
library(ggnewscale)


# packages for loading data
library(xlsx)
library(readr)


# packages for plotting data
library(jsonlite)
library(ggbeeswarm)
library(ggplot2)
library(viridis)
library(patchwork)

# packages for chemical property calculations
library(Peptides)
library(Biostrings)
library(alakazam)


# R version 4.0.4
#install.packages("https://cran.r-project.org/src/contrib/Archive/estimability/estimability_1.4.1.tar.gz", repos = NULL, type = "source")
#install.packages("https://cran.r-project.org/src/contrib/Archive/emmeans/emmeans_1.7.5.tar.gz", repos = NULL, type = "source")
#install.packages("https://cran.r-project.org/src/contrib/Archive/pbkrtest/pbkrtest_0.5-0.1.tar.gz", repos = NULL, type = "source")
#install.packages('car')
#install.packages("factoextra", repos = "https://cran.rstudio.com/")

# packages for PCA anaylysis
library(FactoMineR)
library(factoextra)
library(corrplot)


######################################################################
# scrips to run
######################################################################

# process and characterize fkg22 variant frequency
source("./01_flg22_abundance_and_conservation.R")

# build flagellin phylogenetic tree 
# run the following commands on the command line 
system("mafft --thread 10 --maxiterate 1000 --auto --reorder filC_full_length.fasta > fliC_aligned")
system("iqtree2 s fliC_aligned -st AA -bb 1000 -T AUTO -m TEST")

# plot phylogenetic tree
source("./02_flagellin_tree.R")

# Compare the chemistry of the exposed surface along the LRR for each FLS2 homolog
# determine which chemical property leads to greatest degree of seperation
source("./03_FLS2_rank_properies.R")


# determine base values along the LRR of top two chemistries - bulkiness and charge
source("./04_FLS2_top_properties_heatmap.R")



#########################################################
# funciton - convert DNAstringset attribute to dataframe
#########################################################

# turning AAsequences (fasta) into dataframe -> uses Biostrings functionality/structure
dss2df <- function(dss){
  return(data.frame(width = dss@ranges@width, names = names(dss), seq = as.character(dss), stringsAsFactors = FALSE))
}


# turning AAmultiplesequence alignment into dataframe -> uses Biostrings functionality/structure
aa2df <- function(dss){
  return(data.frame(names = rownames(dss), seq = as.character(dss), stringsAsFactors = FALSE))
}

# turning readAAStringSet alignment into dataframe -> uses Biostrings functionality/structure
aa2df_v2 <- function(dss){
  return(data.frame(names = names(dss), seq = paste(dss), stringsAsFactors = FALSE))
}













