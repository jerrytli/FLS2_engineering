#-----------------------------------------------------------------------------------------------
# Coaker Lab - Plant Pathology Department UC Davis
# Author: Danielle M. Stevens
# Last Updated: 03/26/2024
# Script Purpose: Compare specific FLS2-flg22 pairs to understand why some receptors have expanded ligand receptors
# Inputs Necessary: FLS2 sequences (AtFLS2_LRR.fasta) found in FLS2_flg22_coordinates folder
# Outputs: 
#-----------------------------------------------------------------------------------------------


#########################################################
# funciton - convert DNAstringset attribute to dataframe
#########################################################

calculate_protein_properties <- function(data_frame, type = c("charge", "aIndex", "boman", "hydrophobicity", "bulk"), hydro_scale){
  if(type == "charge"){
    for (i in 1:nrow(data_frame)){
      data_frame$charge[[i]] <- Peptides::charge(data_frame$seq[[i]], pH = 5.4)}
    data_frame$charge <- unlist(data_frame$charge)
    }
  if(type == "aIndex"){
    for (i in 1:nrow(data_frame)){
      data_frame$aIndex[[i]] <- Peptides::aIndex(data_frame$seq[[i]])}
    data_frame$aIndex <- unlist(data_frame$aIndex)
    }
  if(type == "boman"){
    for (i in 1:nrow(data_frame)){
      data_frame$boman[[i]] <- Peptides::boman(data_frame$seq[[i]])}
    data_frame$boman <- unlist(data_frame$boman)
    }
  if(type == "bulk"){
    for (i in 1:nrow(data_frame)){
      data_frame$bulk[[i]] <- alakazam::bulk(data_frame$seq[[i]])}
    data_frame$bulk <- unlist(data_frame$bulk)
    }
  if(type == "hydrophobicity"){
    for (i in 1:nrow(data_frame)){
      data_frame$hydro_scale[[i]] <- Peptides::hydrophobicity(data_frame$seq[[i]], scale = hydro_scale)}
    data_frame$hydro_scale <- unlist(data_frame$hydro_scale)
    }
  return(data_frame)
}

# a small number were filtered out as they are likely not relativent for this work
possible_hydrophobicity_scales <- c("Aboderin", "AbrahamLeo", "Argos", "BlackMould", "BullBreese", "Casari",
                                    "Chothia","Cid", "Cowan3.4", "Cowan7.5", "Eisenberg", "Engelman", "Fasman",
                                    "Fauchere", "Goldsack", "Guy", "HoppWoods", "Janin", "Jones", "Juretic", "Kidera",
                                    "Kuhn","KyteDoolittle","Levitt","Manavalan", "Miyazawa", "Parker", "Ponnuswamy",
                                    "Prabhakaran", "Rao", "Rose", "Roseman", "Sweet", "Tanford", "Welling", "Wilson",
                                    "Wolfenden", "Zimmerman", "interfaceScale_pH8", "octanolScale_pH8",
                                    "oiScale_pH8")



#########################################################
# import fasta files of each FLS2 homolog studied
#########################################################

AtFLS2_seq <- aa2df_v2(Biostrings::readAAStringSet("../FLS2_flg22_coordinates/AtFLS2_LRR.fasta"))
SlFLS2_seq <- aa2df_v2(Biostrings::readAAStringSet("../FLS2_flg22_coordinates/SlFLS2_LRR.fasta"))
FcFLS2_seq <- aa2df_v2(Biostrings::readAAStringSet("../FLS2_flg22_coordinates/FcFLS2_LRR.fasta"))
VrFLS2_seq <- aa2df_v2(Biostrings::readAAStringSet("../FLS2_flg22_coordinates/VrFLS2_LRR.fasta"))
NbFLS2_seq <- aa2df_v2(Biostrings::readAAStringSet("../FLS2_flg22_coordinates/NbFLS2_1_LRR.fasta"))

FLS2XL_seq <- aa2df_v2(Biostrings::readAAStringSet("../FLS2_flg22_coordinates/FLS2xl_LRR.fasta"))
QvFLS2_seq <- aa2df_v2(Biostrings::readAAStringSet("../FLS2_flg22_coordinates/QvFLS2_LRR.fasta"))
TjFLS2_seq <- aa2df_v2(Biostrings::readAAStringSet("../FLS2_flg22_coordinates/TjFLS2_LRR.fasta"))


#########################################################
# import coordinates maps, these were based on 6 exposed residues for each LRR repeat in likely contact with the ligand peptide
#########################################################

# narrow range receptors 
At_FLS2_LRR_residues <- xlsx::read.xlsx("../FLS2_flg22_coordinates/AtFLS2_LRR_exposed_residue_coordinates.xlsx", sheetName = "brassicales_LRRexposed_residue_")
At_FLS2_LRR_residues <- At_FLS2_LRR_residues[1:28,1:8]
At_FLS2_LRR_residues <- cbind("LRR" = At_FLS2_LRR_residues$LRR, "X0" = c(At_FLS2_LRR_residues$First.leucine-1), At_FLS2_LRR_residues[,2:ncol(At_FLS2_LRR_residues)])


Sl_FLS2_LRR_residues <- xlsx::read.xlsx("../FLS2_flg22_coordinates/SlFLS2_LRR_exposed_residue_coordinates.xlsx", sheetName = "brassicales_LRRexposed_residue_")
Sl_FLS2_LRR_residues <- Sl_FLS2_LRR_residues[1:28,1:8]
Sl_FLS2_LRR_residues <- cbind("LRR" = Sl_FLS2_LRR_residues$LRR, "X0" = c(Sl_FLS2_LRR_residues$First.leucine-1), Sl_FLS2_LRR_residues[,2:ncol(Sl_FLS2_LRR_residues)])

Fc_FLS2_LRR_residues <- xlsx::read.xlsx("../FLS2_flg22_coordinates/FcFLS2_LRR_exposed_residue_coordinates.xlsx", sheetName = "brassicales_LRRexposed_residue_")
Fc_FLS2_LRR_residues <- Fc_FLS2_LRR_residues[1:28,1:8]
Fc_FLS2_LRR_residues <- cbind("LRR" = Fc_FLS2_LRR_residues$LRR, "X0" = c(Fc_FLS2_LRR_residues$First.leucine-1), Fc_FLS2_LRR_residues[,2:ncol(Fc_FLS2_LRR_residues)])

Vr_FLS2_LRR_residues <- xlsx::read.xlsx("../FLS2_flg22_coordinates/VrFLS2_LRR_exposed_residue_coordinates.xlsx", sheetName = "LRR_exposed_residue")
Vr_FLS2_LRR_residues <- Vr_FLS2_LRR_residues[1:28,1:8]
Vr_FLS2_LRR_residues <- cbind("LRR" = Vr_FLS2_LRR_residues$LRR, "X0" = c(Vr_FLS2_LRR_residues$First.leucine-1), Vr_FLS2_LRR_residues[,2:ncol(Vr_FLS2_LRR_residues)])

Nb_FLS2_LRR_residues <- xlsx::read.xlsx("../FLS2_flg22_coordinates/NbFLS2_1_LRR_exposed_residue_coordinates.xlsx", sheetName = "LRR_exposed_residue")
Nb_FLS2_LRR_residues <- Nb_FLS2_LRR_residues[1:28,1:8]
Nb_FLS2_LRR_residues <- cbind("LRR" = Nb_FLS2_LRR_residues$LRR, "X0" = c(Nb_FLS2_LRR_residues$First.leucine-1), Nb_FLS2_LRR_residues[,2:ncol(Nb_FLS2_LRR_residues)])

# broad range receptors
FLS2XL_LRR_residues <- xlsx::read.xlsx("../FLS2_flg22_coordinates/FLS2XL_LRR_exposed_residue_coordinates.xlsx", sheetName = "brassicales_LRRexposed_residue_")
FLS2XL_LRR_residues <- FLS2XL_LRR_residues[1:28,1:8]
FLS2XL_LRR_residues <- cbind("LRR" = FLS2XL_LRR_residues$LRR, "X0" = c(FLS2XL_LRR_residues$First.leucine-1), FLS2XL_LRR_residues[,2:ncol(FLS2XL_LRR_residues)])

QvFLS2_LRR_residues <- xlsx::read.xlsx("../FLS2_flg22_coordinates/QvFLS2_LRR_exposed_residue_coordinates.xlsx", sheetName = "brassicales_LRRexposed_residue_")
QvFLS2_LRR_residues <- QvFLS2_LRR_residues[1:28,1:8]
QvFLS2_LRR_residues <- cbind("LRR" = QvFLS2_LRR_residues$LRR, "X0" = c(QvFLS2_LRR_residues$First.leucine-1), QvFLS2_LRR_residues[,2:ncol(QvFLS2_LRR_residues)])

TjFLS2_LRR_residues <- xlsx::read.xlsx("../FLS2_flg22_coordinates/TjFLS2_LRR_exposed_residue_coordinates.xlsx", sheetName = "brassicales_LRRexposed_residue_")
TjFLS2_LRR_residues <- TjFLS2_LRR_residues[1:28,1:8]
TjFLS2_LRR_residues <- cbind("LRR" = TjFLS2_LRR_residues$LRR, "X0" = c(TjFLS2_LRR_residues$First.leucine-1), TjFLS2_LRR_residues[,2:ncol(TjFLS2_LRR_residues)])


#########################################################
# 
#########################################################

# narrow range receptors 
AtFLS2_database <- do.call(rbind, Map(data.frame, coordinates = 1:nchar(AtFLS2_seq$seq),
                                      seq = strsplit(AtFLS2_seq$seq, "")[[1]]))

SlFLS2_database <- do.call(rbind, Map(data.frame, coordinates = 1:nchar(SlFLS2_seq$seq),
                                      seq = strsplit(SlFLS2_seq$seq, "")[[1]]))

FcFLS2_database <- do.call(rbind, Map(data.frame, coordinates = 1:nchar(FcFLS2_seq$seq),
                                      seq = strsplit(FcFLS2_seq$seq, "")[[1]]))

VrFLS2_database <- do.call(rbind, Map(data.frame, coordinates = 1:nchar(VrFLS2_seq$seq),
                                      seq = strsplit(VrFLS2_seq$seq, "")[[1]]))

NbFLS2_database <- do.call(rbind, Map(data.frame, coordinates = 1:nchar(NbFLS2_seq$seq),
                                      seq = strsplit(NbFLS2_seq$seq, "")[[1]]))

# broad range receptors
FLS2XL_database <- do.call(rbind, Map(data.frame, coordinates = 1:nchar(FLS2XL_seq$seq),
                                      seq = strsplit(FLS2XL_seq$seq, "")[[1]]))

QvFLS2_database <- do.call(rbind, Map(data.frame, coordinates = 1:nchar(QvFLS2_seq$seq),
                                      seq = strsplit(QvFLS2_seq$seq, "")[[1]]))

TjFLS2_database <- do.call(rbind, Map(data.frame, coordinates = 1:nchar(TjFLS2_seq$seq),
                                      seq = strsplit(TjFLS2_seq$seq, "")[[1]]))



#########################################################
# filter for exposed residues
#########################################################

Athold_residues_list <- as.vector(matrix(t(At_FLS2_LRR_residues[,2:9])))
Slhold_residues_list <- as.vector(matrix(t(Sl_FLS2_LRR_residues[,2:9])))
Fchold_residues_list <- as.vector(matrix(t(Fc_FLS2_LRR_residues[,2:9])))
Vrhold_residues_list <- as.vector(matrix(t(Vr_FLS2_LRR_residues[,2:9])))
Nbhold_residues_list <- as.vector(matrix(t(Nb_FLS2_LRR_residues[,2:9])))

Xlhold_residues_list <- as.vector(matrix(t(FLS2XL_LRR_residues[,2:9])))
Qvhold_residues_list <- as.vector(matrix(t(QvFLS2_LRR_residues[,2:9])))
Tjhold_residues_list <- as.vector(matrix(t(TjFLS2_LRR_residues[,2:9])))


#########################################################
# Determine which property infulences FLS2 property data
#########################################################

# determine how large dataframe needs to be set based on total number of parameters
determine_chemistry_properties_rownames <- c("aliphatic", "bulk", "aIndex", "charge","boman", possible_hydrophobicity_scales)

determine_chemistry_properties <- data.frame("AtFLS2" = numeric(length(determine_chemistry_properties_rownames)), 
                                             "SlFLS2" = numeric(length(determine_chemistry_properties_rownames)), 
                                             "FcFLS2" = numeric(length(determine_chemistry_properties_rownames)), 
                                             "VrFLS2" = numeric(length(determine_chemistry_properties_rownames)), 
                                             "NbFLS2" = numeric(length(determine_chemistry_properties_rownames)), 
                                             
                                             "FLS2Xl" = numeric(length(determine_chemistry_properties_rownames)), 
                                             "QvFLS2" = numeric(length(determine_chemistry_properties_rownames)), 
                                             "TjFLS2" = numeric(length(determine_chemistry_properties_rownames)))

rownames(determine_chemistry_properties) <- determine_chemistry_properties_rownames


# seondary alipharic index overall score
determine_chemistry_properties$AtFLS2[1] <- alakazam::aliphatic(paste(unlist(AtFLS2_database[AtFLS2_database$coordinates %in% Athold_residues_list,2]), collapse = ""))
determine_chemistry_properties$SlFLS2[1] <- alakazam::aliphatic(paste(unlist(SlFLS2_database[SlFLS2_database$coordinates %in% Slhold_residues_list,2]), collapse = ""))
determine_chemistry_properties$FcFLS2[1] <- alakazam::aliphatic(paste(unlist(FcFLS2_database[FcFLS2_database$coordinates %in% Fchold_residues_list,2]), collapse = ""))
determine_chemistry_properties$VrFLS2[1] <- alakazam::aliphatic(paste(unlist(VrFLS2_database[VrFLS2_database$coordinates %in% Vrhold_residues_list,2]), collapse = ""))
determine_chemistry_properties$NbFLS2[1] <- alakazam::aliphatic(paste(unlist(NbFLS2_database[NbFLS2_database$coordinates %in% Nbhold_residues_list,2]), collapse = ""))

determine_chemistry_properties$FLS2Xl[1] <- alakazam::aliphatic(paste(unlist(FLS2XL_database[FLS2XL_database$coordinates %in% Xlhold_residues_list,2]), collapse = ""))
determine_chemistry_properties$QvFLS2[1] <- alakazam::aliphatic(paste(unlist(QvFLS2_database[QvFLS2_database$coordinates %in% Qvhold_residues_list,2]), collapse = ""))
determine_chemistry_properties$TjFLS2[1] <- alakazam::aliphatic(paste(unlist(TjFLS2_database[TjFLS2_database$coordinates %in% Tjhold_residues_list,2]), collapse = ""))


# bulkiness overall score
determine_chemistry_properties$AtFLS2[2] <- alakazam::bulk(paste(unlist(AtFLS2_database[AtFLS2_database$coordinates %in% Athold_residues_list,2]), collapse = ""))
determine_chemistry_properties$SlFLS2[2] <- alakazam::bulk(paste(unlist(SlFLS2_database[SlFLS2_database$coordinates %in% Slhold_residues_list,2]), collapse = ""))
determine_chemistry_properties$FcFLS2[2] <- alakazam::bulk(paste(unlist(FcFLS2_database[FcFLS2_database$coordinates %in% Fchold_residues_list,2]), collapse = ""))
determine_chemistry_properties$VrFLS2[2] <- alakazam::bulk(paste(unlist(VrFLS2_database[VrFLS2_database$coordinates %in% Vrhold_residues_list,2]), collapse = ""))
determine_chemistry_properties$NbFLS2[2] <- alakazam::bulk(paste(unlist(NbFLS2_database[NbFLS2_database$coordinates %in% Nbhold_residues_list,2]), collapse = ""))

determine_chemistry_properties$FLS2Xl[2] <- alakazam::bulk(paste(unlist(FLS2XL_database[FLS2XL_database$coordinates %in% Xlhold_residues_list,2]), collapse = ""))
determine_chemistry_properties$QvFLS2[2] <- alakazam::bulk(paste(unlist(QvFLS2_database[QvFLS2_database$coordinates %in% Qvhold_residues_list,2]), collapse = ""))
determine_chemistry_properties$TjFLS2[2] <- alakazam::bulk(paste(unlist(TjFLS2_database[TjFLS2_database$coordinates %in% Tjhold_residues_list,2]), collapse = ""))


# aIndex overall score
determine_chemistry_properties$AtFLS2[3] <- Peptides::aIndex(paste(unlist(AtFLS2_database[AtFLS2_database$coordinates %in% Athold_residues_list,2]), collapse = ""))
determine_chemistry_properties$SlFLS2[3] <- Peptides::aIndex(paste(unlist(SlFLS2_database[SlFLS2_database$coordinates %in% Slhold_residues_list,2]), collapse = ""))
determine_chemistry_properties$FcFLS2[3] <- Peptides::aIndex(paste(unlist(FcFLS2_database[FcFLS2_database$coordinates %in% Fchold_residues_list,2]), collapse = ""))
determine_chemistry_properties$VrFLS2[3] <- Peptides::aIndex(paste(unlist(VrFLS2_database[VrFLS2_database$coordinates %in% Vrhold_residues_list,2]), collapse = ""))
determine_chemistry_properties$NbFLS2[3] <- Peptides::aIndex(paste(unlist(NbFLS2_database[NbFLS2_database$coordinates %in% Nbhold_residues_list,2]), collapse = ""))

determine_chemistry_properties$FLS2Xl[3] <- Peptides::aIndex(paste(unlist(FLS2XL_database[FLS2XL_database$coordinates %in% Xlhold_residues_list,2]), collapse = ""))
determine_chemistry_properties$QvFLS2[3] <- Peptides::aIndex(paste(unlist(QvFLS2_database[QvFLS2_database$coordinates %in% Qvhold_residues_list,2]), collapse = ""))
determine_chemistry_properties$TjFLS2[3] <- Peptides::aIndex(paste(unlist(TjFLS2_database[TjFLS2_database$coordinates %in% Tjhold_residues_list,2]), collapse = ""))

# charge overall score
determine_chemistry_properties$AtFLS2[4] <- Peptides::charge(paste(unlist(AtFLS2_database[AtFLS2_database$coordinates %in% Athold_residues_list,2]), collapse = ""), pH = 5.4)
determine_chemistry_properties$SlFLS2[4] <- Peptides::charge(paste(unlist(SlFLS2_database[SlFLS2_database$coordinates %in% Slhold_residues_list,2]), collapse = ""), pH = 5.4)
determine_chemistry_properties$FcFLS2[4] <- Peptides::charge(paste(unlist(FcFLS2_database[FcFLS2_database$coordinates %in% Fchold_residues_list,2]), collapse = ""), pH = 5.4)
determine_chemistry_properties$VrFLS2[4] <- Peptides::charge(paste(unlist(VrFLS2_database[VrFLS2_database$coordinates %in% Vrhold_residues_list,2]), collapse = ""), pH = 5.4)
determine_chemistry_properties$NbFLS2[4] <- Peptides::charge(paste(unlist(NbFLS2_database[NbFLS2_database$coordinates %in% Nbhold_residues_list,2]), collapse = ""), pH = 5.4)

determine_chemistry_properties$FLS2Xl[4] <- Peptides::charge(paste(unlist(FLS2XL_database[FLS2XL_database$coordinates %in% Xlhold_residues_list,2]), collapse = ""), pH = 5.4)
determine_chemistry_properties$QvFLS2[4] <- Peptides::charge(paste(unlist(QvFLS2_database[QvFLS2_database$coordinates %in% Qvhold_residues_list,2]), collapse = ""), pH = 5.4)
determine_chemistry_properties$TjFLS2[4] <- Peptides::charge(paste(unlist(TjFLS2_database[TjFLS2_database$coordinates %in% Tjhold_residues_list,2]), collapse = ""), pH = 5.4)

# boman overall score
determine_chemistry_properties$AtFLS2[5] <- Peptides::boman(paste(unlist(AtFLS2_database[AtFLS2_database$coordinates %in% Athold_residues_list,2]), collapse = ""))
determine_chemistry_properties$SlFLS2[5] <- Peptides::boman(paste(unlist(SlFLS2_database[SlFLS2_database$coordinates %in% Slhold_residues_list,2]), collapse = ""))
determine_chemistry_properties$FcFLS2[5] <- Peptides::boman(paste(unlist(FcFLS2_database[FcFLS2_database$coordinates %in% Fchold_residues_list,2]), collapse = ""))
determine_chemistry_properties$VrFLS2[5] <- Peptides::boman(paste(unlist(VrFLS2_database[VrFLS2_database$coordinates %in% Vrhold_residues_list,2]), collapse = ""))
determine_chemistry_properties$NbFLS2[5] <- Peptides::boman(paste(unlist(NbFLS2_database[NbFLS2_database$coordinates %in% Nbhold_residues_list,2]), collapse = ""))

determine_chemistry_properties$FLS2Xl[5] <- Peptides::boman(paste(unlist(FLS2XL_database[FLS2XL_database$coordinates %in% Xlhold_residues_list,2]), collapse = ""))
determine_chemistry_properties$QvFLS2[5] <- Peptides::boman(paste(unlist(QvFLS2_database[QvFLS2_database$coordinates %in% Qvhold_residues_list,2]), collapse = ""))
determine_chemistry_properties$TjFLS2[5] <- Peptides::boman(paste(unlist(TjFLS2_database[TjFLS2_database$coordinates %in% Tjhold_residues_list,2]), collapse = ""))

# boman hydrophobicity scores
for (i in 6:nrow(determine_chemistry_properties)){
  determine_chemistry_properties$AtFLS2[i] <- Peptides::hydrophobicity(paste(unlist(AtFLS2_database[AtFLS2_database$coordinates %in% Athold_residues_list,2]), collapse = ""), scale = possible_hydrophobicity_scales[i-5])
  determine_chemistry_properties$SlFLS2[i] <- Peptides::hydrophobicity(paste(unlist(SlFLS2_database[SlFLS2_database$coordinates %in% Slhold_residues_list,2]), collapse = ""), scale = possible_hydrophobicity_scales[i-5])
  determine_chemistry_properties$FcFLS2[i] <- Peptides::hydrophobicity(paste(unlist(FcFLS2_database[FcFLS2_database$coordinates %in% Fchold_residues_list,2]), collapse = ""), scale = possible_hydrophobicity_scales[i-5])
  determine_chemistry_properties$VrFLS2[i] <- Peptides::hydrophobicity(paste(unlist(VrFLS2_database[VrFLS2_database$coordinates %in% Vrhold_residues_list,2]), collapse = ""), scale = possible_hydrophobicity_scales[i-5])
  determine_chemistry_properties$NbFLS2[i] <- Peptides::hydrophobicity(paste(unlist(NbFLS2_database[NbFLS2_database$coordinates %in% Nbhold_residues_list,2]), collapse = ""), scale = possible_hydrophobicity_scales[i-5])
  
  determine_chemistry_properties$FLS2Xl[i] <- Peptides::hydrophobicity(paste(unlist(FLS2XL_database[FLS2XL_database$coordinates %in% Xlhold_residues_list,2]), collapse = ""), scale = possible_hydrophobicity_scales[i-5])
  determine_chemistry_properties$QvFLS2[i] <- Peptides::hydrophobicity(paste(unlist(QvFLS2_database[QvFLS2_database$coordinates %in% Qvhold_residues_list,2]), collapse = ""), scale = possible_hydrophobicity_scales[i-5])
  determine_chemistry_properties$TjFLS2[i] <- Peptides::hydrophobicity(paste(unlist(TjFLS2_database[TjFLS2_database$coordinates %in% Tjhold_residues_list,2]), collapse = ""), scale = possible_hydrophobicity_scales[i-5])
}



#########################################################
# calculate PCA on property data to determine which may contribute to perception range
#########################################################

# removed aIndex as not a comparable metrix to other chemiysrt scales
calculate_property_PCA <- FactoMineR::PCA(determine_chemistry_properties[!rownames(determine_chemistry_properties) %in% "aIndex",],  graph = FALSE, scale.unit = TRUE)
individual_variables <- factoextra::get_pca_var(calculate_property_PCA)

# ----------- correlation plot to show which dims contribute to which individuals --------------------
corrplot::corrplot(individual_variables$contrib, is.corr=FALSE, tl.col = "black", col = COL1('Purples'), cl.ratio = 0.4, cl.align = 'l') 


# ----------- Contributions of variables to each PC dimensions --------------------

# dimensions 1 - size exported 1.5x2
factoextra::fviz_contrib(calculate_property_PCA, choice = "ind", axes = 1, top = 5, fill = "black", color = "black") + theme_linedraw() + 
  xlab('Chemistry Parameter') + 
  theme(axis.text = element_text(color = 'black', size = 6), axis.text.x = element_text(angle = 45, hjust = 1),
        axis.title = element_text(color = 'black', size = 7), 
        title = element_text(color = 'black', size = 7, face = 'bold'))

# dimensions 2 - size exported 1.5x2
factoextra::fviz_contrib(calculate_property_PCA, choice = "ind", axes = 2, top = 5, fill = "black", color = "black") + theme_linedraw() + 
  xlab('Chemistry Parameter') + 
  theme(axis.text = element_text(color = 'black', size = 6), axis.text.x = element_text(angle = 45, hjust = 1),
        axis.title = element_text(color = 'black', size = 7), 
        title = element_text(color = 'black', size = 7, face = 'bold'))

# dimensions 3 -- size exported 1.5x2
factoextra::fviz_contrib(calculate_property_PCA, choice = "ind", axes = 3, top = 5, fill = "black", color = "black") + theme_linedraw() + 
  xlab('Chemistry Parameter') + 
  theme(axis.text = element_text(color = 'black', size = 6), axis.text.x = element_text(angle = 45, hjust = 1),
        axis.title = element_text(color = 'black', size = 7), 
        title = element_text(color = 'black', size = 7, face = 'bold'))

# dimensions 4 - size exported
factoextra::fviz_contrib(calculate_property_PCA, choice = "ind", axes = 4, top = 5, fill = "black", color = "black") + theme_linedraw() + 
  xlab('Chemistry Parameter') + 
  theme(axis.text = element_text(color = 'black', size = 6), axis.text.x = element_text(angle = 45, hjust = 1),
        axis.title = element_text(color = 'black', size = 7), 
        title = element_text(color = 'black', size = 7, face = 'bold'))

# dimensions 5 - size exported
factoextra::fviz_contrib(calculate_property_PCA, choice = "ind", axes = 5, top = 5, fill = "black", color = "black") + theme_linedraw() + 
  xlab('Chemistry Parameter') + 
  theme(axis.text = element_text(color = 'black', size = 6), axis.text.x = element_text(angle = 45, hjust = 1),
        axis.title = element_text(color = 'black', size = 7), 
        title = element_text(color = 'black', size = 7, face = 'bold'))


# ----------- Contributions of variables to each PC dimensions --------------------
factoextra::fviz_pca_biplot(calculate_property_PCA, invisible ="var") + theme_classic2() + xlim(-9,14) + ylim(-0.5,4) +
  theme(axis.text = element_text(color = 'black', size = 15), 
        axis.title = element_text(color = 'black', size = 17))

# zoom into lower contributing factors
factoextra::fviz_pca_biplot(calculate_property_PCA, invisible ="var") + theme_classic2() + xlim(-1.85,1.02) + ylim(-0.14,0.01) +
  theme(axis.text = element_text(color = 'black', size = 15), 
        axis.title = element_text(color = 'black', size = 17))



######### ---------------------old code------------------------------
# need to z-score standardize each parameter
#for (i in 1:nrow(determine_chemistry_properties)){
#   a <- as.vector(unlist(determine_chemistry_properties[i,]))
#   m <- mean(a)
#   s <- sd(a)
#   print(paste(a,s, collapse = "_"))
#   z <- (a - m)/s
#   determine_chemistry_properties[i,] <- z
#}
#[!rownames(determine_chemistry_properties) %in% "aIndex",]
