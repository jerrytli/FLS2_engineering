#-----------------------------------------------------------------------------------------------
# Coaker Lab - Plant Pathology Department UC Davis
# Author: Danielle M. Stevens
# Last Updated: 08/21/2024
# Script Purpose: Process and plot flagellin phylogenetic tree
# Inputs Necessary: database Stevens_et_al_MAMP_database.xlsx and tree fliC_aligned.contree
# Outputs: Flagellin Phylogenic tree
#-----------------------------------------------------------------------------------------------


######################################################################
# load flg22 variants from Stevens et al., 2024 and immune outcome data accross multiple studies
######################################################################

import_flg22_hits <- xlsx::read.xlsx("../Stevens_et_al_MAMP_database.xlsx", sheetName = "Sheet1")
import_flg22_hits <- import_flg22_hits[import_flg22_hits$MAMP_Hit %in% c("flg22_consensus"),]
import_flg22_hits <- import_flg22_hits[,c(2:9)]


######################################################################
# Run flagellin hits in bash to build a phylogenetic tree
######################################################################

Flagellin_protein_tree <- read.tree("../Flagellin_tree/clip_version1/filC_aligned.clipkit.contree")
Flagellin_protein_tree <- phangorn::midpoint(Flagellin_protein_tree, node.labels='label')



for (i in 1:length(Flagellin_protein_tree$tip.label)){
  Flagellin_protein_tree$tip.label[i] <- paste(strsplit(Flagellin_protein_tree$tip.label[i], "|", fixed = T)[[1]][1],
                                               strsplit(Flagellin_protein_tree$tip.label[i], "|", fixed = T)[[1]][5],
                                               sep = "|")
}

for (i in 1:length(Flagellin_protein_tree$tip.label)){
  Flagellin_protein_tree$tip.label[i] <- sub("(\\|[^.]+)\\.", "\\1_", Flagellin_protein_tree$tip.label[i])
}


# add info for ROS outcome 
import_flg22_hits["Immunogenicity"] <- "NT"
for (i in 1:nrow(import_flg22_hits)){
  if (import_flg22_hits$MAMP_Sequence[i] %in% flg22_ros_screen_outcome$Sequence == FALSE){
    next()
  }
  if (import_flg22_hits$MAMP_Sequence[i] %in% flg22_ros_screen_outcome$Sequence == TRUE){
    import_flg22_hits$Immunogenicity[i] <- flg22_ros_screen_outcome[flg22_ros_screen_outcome$Sequence %in% import_flg22_hits$MAMP_Sequence[i],2]
  }
} 


# copy over for figure before
all_flagellin_values <- import_flg22_hits


# add in new label to match with protein tree tips
import_flg22_hits["New_label"] <- paste(import_flg22_hits$Protein_Name, import_flg22_hits$File_Name, sep = "|")
import_flg22_hits <- import_flg22_hits[c("New_label", "Percent_Identity", "Genera", "Immunogenicity", "File_Name")]
import_flg22_hits <- import_flg22_hits[import_flg22_hits$New_label %in% Flagellin_protein_tree$tip.label,]


######################################################################
# Run flagellin hits in bash to build a phylogenetic tree
######################################################################

Flagellin_protein_tree <- ggtree(Flagellin_protein_tree, layout = "rectangular", ladderize = T, size = 0.2, linetype = 1, 
                                 aes(color = as.numeric(label) > 70 | isTip), show.legend = FALSE) %<+% import_flg22_hits +
  scale_colour_grey(start = .7, end = .1) +
  theme(legend.position = "bottom") +
  
  new_scale("color") +
  geom_tippoint(aes(color = Genera), size = 0.01, show.legend = FALSE) +
  
  scale_color_manual("Genera", values = Genera_colors) +
  geom_treescale(x = -0.15, y = -1, linesize = 1, family = "Arial", offset = 28) 




# label the nodes - this is for visualization only to determine which nodes to expand and collapse
Flagellin_protein_tree + geom_text2(aes(subset = !isTip, label = node), hjust = -.3, size = 2) 

# save as 5x5.5 inches - pdf
Flagellin_protein_tree %>% 
  #xanthomonas
  ggtree::collapse(5563, 'min') %>% 
  #agro
  ggtree::collapse(7063, 'min') %>%  
  #pseudomonas
  ggtree::collapse(6246, 'max') %>% 
  #ggtree::collapse(8402, 'min') %>%
  ggtree::collapse(8421, 'min') %>% 
  ggtree::collapse(8471, 'min') %>% 
  #ralstonia
  ggtree::collapse(8481, 'max') %>% 
  #pect, erwinina, dickeya
  ggtree::collapse(5565, 'min')



# ---------------------------------------------------------------



all_flagellin_values_flg22 <- subset(all_flagellin_values, all_flagellin_values$MAMP_Hit == "flg22_consensus")


# determine abundance within the dataset and then remove redundant eptiopes
flg22_abudance <- as.data.frame(all_flagellin_values_flg22 %>% group_by(MAMP_Sequence) %>% count(n=n()))
flg22_abudance["%Abundance"] <- 100*(flg22_abudance$nn/flg22_abudance$n)
all_flagellin_values_flg22 <- distinct(all_flagellin_values_flg22, MAMP_Sequence, .keep_all = T)
all_flagellin_values_flg22 <- all_flagellin_values_flg22 %>% inner_join(flg22_abudance)

# determine total mean of flagellin AA siminlarity values
mean_aa_flg22_value <- mean(all_flagellin_values_flg22$Percent_Identity)
mean_by_genera_flg22 <-as.data.frame(all_flagellin_values_flg22 %>% group_by(Genera) %>% summarise(mean = mean(Percent_Identity)))




# distinct(all_flagellin_values, MAMP_Sequence, .keep_all = T)[distinct(all_flagellin_values_flg22, MAMP_Sequence, .keep_all = T)$Genera %in% c("Agrobacterium"),],

# filter out curteobacterium and leifstonia
all_flagellin_values <- all_flagellin_values[!all_flagellin_values$Genera %in% c("Curtobacterium","Leifsonia", "Streptomyces"),]

# determine abundance within the dataset and then remove redundant eptiopes
flg22_abudance <- as.data.frame(all_flagellin_values %>% group_by(Genera, MAMP_Sequence) %>% count(n=n()))
genera_totals <- as.data.frame(flg22_abudance %>% group_by(Genera) %>% summarise(genera_total = sum(nn)))
flg22_abudance <- merge(x = flg22_abudance, y = genera_totals, by = "Genera", all = TRUE)

flg22_abudance["%Abundance"] <- 100*(flg22_abudance$nn/flg22_abudance$genera_total)
all_flagellin_values <- distinct(all_flagellin_values, MAMP_Sequence, .keep_all = T)
all_flagellin_values <- all_flagellin_values %>% inner_join(flg22_abudance)

# determine total mean of flagellin AA siminlarity values
mean_aa_flg22_value <- mean(all_flagellin_values$Percent_Identity)
mean_by_genera_flg22 <-as.data.frame(all_flagellin_values %>% group_by(Genera) %>% summarise(mean = mean(Percent_Identity)))

# color code for Immunogenicity
Immunogenicity_colors <- c("#206ba3","#005b96","#d7573b", "#BFBFBF","#E1AD01")
names(Immunogenicity_colors) <- c("Immunogenic","Weakly Immunogenic","Non-Immunogenic","NT","Unknown")




all_flagellin_values$Genera_f = factor(all_flagellin_values$Genera, levels=c('Agrobacterium','Pseudomonas','Ralstonia','Pectobacterium','Erwinia','Dickeya','Xanthomonas'))
ggplot(all_flagellin_values, 
       aes(x = Percent_Identity, y = `%Abundance`, fill = Immunogenicity)) +
  geom_bar(position="stack", stat = "identity",  width = 4.2) +
  facet_grid(Genera_f ~ .) +
  theme_linedraw() +
  #xlim(20,100) +
  scale_fill_manual("Immunogenicity", values = Immunogenicity_colors) +
  xlab("%AA Similarity to Consensus") +
  theme(legend.position = "none", panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        panel.background = element_blank(),
        axis.text.x = element_text(color = "black", size = 7), axis.text.y = element_text(color = "black", size = 7),
        axis.title.x = element_text(size = 8), axis.title.y = element_text(size =8), strip.text.y = element_text(size =3))

geom_text(aes(label = Weight), vjust = -0.2)





