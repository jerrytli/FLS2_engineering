#-----------------------------------------------------------------------------------------------
# Coaker Lab - Plant Pathology Department UC Davis
# Author: Danielle M. Stevens
# Last Updated: 08/21/2024
# Script Purpose: Process Flg22 data for Jerry's Manuscript
# Inputs Necessary: Dataframe of mined MAMPs from Stevens et al., 2024. PNAS
# Outputs: Plots for frequency and distribution of flg22 variants
#-----------------------------------------------------------------------------------------------


######################################################################
# load which flg22 variants jerry characterized
######################################################################

jerry_flg22 <- data.frame("Short_hand" = c('Atu','Rrh','Rso1',
                                           'Rso2','Dda','Eam',
                                           'Pvi','Xfr','Xor'), 
                          "MAMP_Seq" = c("SRVSSGLRVKSASDNAAYWSIA","QQVTSGFRVGSASDDPTYWSMA","QRLSTGLRVNSAQDDSAAYAAS",
                                         "QRLSTGMRVNSAQDDAAAYASA","ERLSSGLAINSAKDNAAGSGIV","ERLSSGSRITSSKEDAAGQAIS",
                                         "SRLSSGLKVQNARDNVGVLSTI","QQLSSGKRITSFSVDAAGGAIA","QQLSSGKRITSFAVDAAGGAIA"))



######################################################################
# Supplemental Figure 1A - Flg22 variant similarity to consensus in resoect to other studies
######################################################################

# load in table of previously characterized flg22 variants
flg22_ros_screen_outcome <- xlsx:::read.xlsx("./../flg22_variant_past_outcome_data.xlsx", sheetName = "Sheet1")
flg22_ros_screen_outcome <- flg22_ros_screen_outcome[,c(1:3)]
flg22_ros_screen_outcome <- rbind(flg22_ros_screen_outcome, cbind("Sequence" = jerry_flg22$MAMP_Seq, 
                                                                  "Known.Outcome" = rep("Unknown", length(jerry_flg22)), 
                                                                  "Reference" = rep("This Study", length(jerry_flg22))))



Percent_Similarity_to_Consensus <- list()
pull_ref_MAMP <- subset(load_reference_MAMPs_fasta, names == 'flg22_consensus')
for (i in 1:nrow(flg22_ros_screen_outcome)){
  Percent_Similarity_to_Consensus[[i]] <- Biostrings::pid(Biostrings::pairwiseAlignment(pull_ref_MAMP$seq, 
                                                                                        flg22_ros_screen_outcome$Sequence[i], type = "global-local", 
                                                                                        gapOpening = 100, gapExtension = 100, substitutionMatrix = BLOSUM45), type = "PID1")
}

flg22_ros_screen_outcome <- cbind(flg22_ros_screen_outcome, "Percent_Identity" = unlist(Percent_Similarity_to_Consensus))


#  ------- plot AA similarity of jerry's flg22 variants in comparison to the consensus and other flg22 papers -------

ggplot(subset(flg22_ros_screen_outcome, flg22_ros_screen_outcome$Reference != "Parys et al., 2021"), 
       aes(x = 0, y = Percent_Identity, color = Reference)) +
  geom_beeswarm(side = 1, cex = 3, corral.width = 0.05, stat="identity") +
  coord_flip() +
  my_ggplot_theme +
  xlab("Variants Reported") +
  ylab("%AA Similarity to Consensus") +
  theme(panel.grid.major.y = element_line(color = "grey"), 
        panel.border = element_rect(colour = "black", fill = NA, size = 0.7),
        axis.ticks.y = element_blank())



######################################################################
# Supplemental Figure 1B - Flg22 variant frequency using dataset from flg22 
######################################################################


temp_hold_filtered_seq <- filtered_hold_MAMP_seqs[filtered_hold_MAMP_seqs$MAMP_Sequence %in% jerry_flg22$MAMP_Seq,]
occurance_data <- data.frame(temp_hold_filtered_seq %>% group_by(MAMP_Sequence,Genera) %>% summarise(This_Study=n()))
total_flg22_count <- data.frame(filtered_hold_MAMP_seqs %>% group_by(MAMP_Hit,Genera) %>% summarise(Total_Number=n()))
total_flg22_count <- subset(total_flg22_count, total_flg22_count$MAMP_Hit =="flg22_consensus")
occurance_data <- na.omit(total_flg22_count %>% left_join(occurance_data, by = c(Genera = "Genera")))
occurance_data$Genera <- factor(occurance_data$Genera, levels = c("Dickeya","Ralstonia","Xanthomonas",
                                                                  "Erwinia","Agrobacterium","Pseudomonas"))

#  ------- plot frequency of jerry's flg22 variants in comparison to the total dataset in MAMP paper -------

ggplot(reshape2::melt(occurance_data[,2:5]), aes(x=Genera,  y=value, fill=variable)) +
  geom_bar(position="fill", stat="identity") + 
  coord_flip() +
  my_ggplot_theme + 
  geom_text(aes(label = value), position = position_fill(vjust = 0.5, reverse=F)) +
  scale_fill_brewer(palette = "Paired", direction = -1) +
  ylab("Epitope Frequency") +
  theme(axis.text.x = element_blank(), axis.ticks.x = element_blank())



######################################################################
# Supplemental Figure 1C - Flg22 variant freequency using dataset from flg22 
######################################################################

# determine frequency of each species for each mamp variant

# ----------------------------------agro and relatives ----------------------------------
Atu_flg22 <- list()
Atu_variant <- subset(temp_hold_filtered_seq, temp_hold_filtered_seq$MAMP_Sequence == "SRVSSGLRVKSASDNAAYWSIA")
for (i in 1:nrow(Atu_variant)){Atu_flg22[[i]] <- strsplit(Atu_variant$File_Name[i],"_")[[1]][2]}
Atu_flg22 <- table(as.data.frame(unlist(Atu_flg22)))

Rrh_flg22 <- list()
Atu_variant <- subset(temp_hold_filtered_seq, temp_hold_filtered_seq$MAMP_Sequence == "QQVTSGFRVGSASDDPTYWSMA")
for (i in 1:nrow(Atu_variant)){Rrh_flg22[[i]] <- strsplit(Atu_variant$File_Name[i],"_")[[1]][2]}
Rrh_flg22 <- table(as.data.frame(unlist(Rrh_flg22)))

# ---------------------------------- ralstonia ----------------------------------
Rals1_flg22 <- list()
Atu_variant <- subset(temp_hold_filtered_seq, temp_hold_filtered_seq$MAMP_Sequence == "QRLSTGLRVNSAQDDSAAYAAS")
for (i in 1:nrow(Atu_variant)){Rals1_flg22[[i]] <- strsplit(Atu_variant$File_Name[i],"_")[[1]][2]}
Rals1_flg22 <- table(as.data.frame(unlist(Rals1_flg22)))

Rals2_flg22 <- list()
Atu_variant <- subset(temp_hold_filtered_seq, temp_hold_filtered_seq$MAMP_Sequence == "QRLSTGMRVNSAQDDAAAYASA")
for (i in 1:nrow(Atu_variant)){Rals2_flg22[[i]] <- strsplit(Atu_variant$File_Name[i],"_")[[1]][2]}
Rals2_flg22 <- table(as.data.frame(unlist(Rals2_flg22)))

# ---------------------------------- Dic, Eam, pseudomonas ----------------------------------
Dic_flg22 <- list()
Atu_variant <- subset(temp_hold_filtered_seq, temp_hold_filtered_seq$MAMP_Sequence == "ERLSSGLAINSAKDNAAGSGIV")
for (i in 1:nrow(Atu_variant)){Dic_flg22[[i]] <- strsplit(Atu_variant$File_Name[i],"_")[[1]][2]}
Dic_flg22 <- table(as.data.frame(unlist(Dic_flg22)))

Eam_flg22 <- list()
Atu_variant <- subset(temp_hold_filtered_seq, temp_hold_filtered_seq$MAMP_Sequence == "ERLSSGSRITSSKEDAAGQAIS")
for (i in 1:nrow(Atu_variant)){Eam_flg22[[i]] <- strsplit(Atu_variant$File_Name[i],"_")[[1]][2]}
Eam_flg22 <- table(as.data.frame(unlist(Eam_flg22)))

Pal_flg22 <- list()
Atu_variant <- subset(temp_hold_filtered_seq, temp_hold_filtered_seq$MAMP_Sequence == "SRLSSGLKVQNARDNVGVLSTI")
for (i in 1:nrow(Atu_variant)){Pal_flg22[[i]] <- strsplit(Atu_variant$File_Name[i],"_")[[1]][2]}
Pal_flg22 <- table(as.data.frame(unlist(Pal_flg22)))

# ---------------------------------- xanthomonas ----------------------------------
Xan1_flg22 <- list()
Atu_variant <- subset(temp_hold_filtered_seq, temp_hold_filtered_seq$MAMP_Sequence == "QQLSSGKRITSFSVDAAGGAIA")
for (i in 1:nrow(Atu_variant)){Xan1_flg22[[i]] <- strsplit(Atu_variant$File_Name[i],"_")[[1]][2]}
Xan1_flg22 <- table(as.data.frame(unlist(Xan1_flg22)))

Xan2_flg22 <- list()
Atu_variant <- subset(temp_hold_filtered_seq, temp_hold_filtered_seq$MAMP_Sequence == "QQLSSGKRITSFAVDAAGGAIA")
for (i in 1:nrow(Atu_variant)){Xan2_flg22[[i]] <- strsplit(Atu_variant$File_Name[i],"_")[[1]][2]}
Xan2_flg22 <- table(as.data.frame(unlist(Xan2_flg22)))



# manually input the frequency of eptiope variants and their respective species into prism to plot








