#-----------------------------------------------------------------------------------------------
# Coaker Lab - Plant Pathology Department UC Davis
# Author: Danielle M. Stevens
# Last Updated: 08/21/2024
# Script Purpose: Compare specific FLS2-flg22 pairs to understand why some receptors have expanded ligand receptors
# Inputs Necessary: FLS2 protein sequences and residue coordinate map 
# Outputs: 
#-----------------------------------------------------------------------------------------------


#########################################################
# calculate residue propoeries - dim 1
#########################################################

AtFLS2_database <- calculate_protein_properties(AtFLS2_database, type = "bulk")
SlFLS2_database <- calculate_protein_properties(SlFLS2_database, type = "bulk")
FcFLS2_database <- calculate_protein_properties(FcFLS2_database, type = "bulk")
VrFLS2_database <- calculate_protein_properties(VrFLS2_database, type = "bulk")
NbFLS2_database <- calculate_protein_properties(NbFLS2_database, type = "bulk")

FLS2XL_database <- calculate_protein_properties(FLS2XL_database, type = "bulk")
QvFLS2_database <- calculate_protein_properties(QvFLS2_database, type = "bulk")
TjFLS2_database <- calculate_protein_properties(TjFLS2_database, type = "bulk")



hold_charge_matrix <- data.frame("AtFLS2" = AtFLS2_database[AtFLS2_database$coordinates %in% Athold_residues_list,3],
                                 "SlFLS2" = SlFLS2_database[SlFLS2_database$coordinates %in% Slhold_residues_list,3],
                                 "NbFLS2" = NbFLS2_database[NbFLS2_database$coordinates %in% Nbhold_residues_list,3],
                                 
                                 "VrFLS2" = VrFLS2_database[VrFLS2_database$coordinates %in% Vrhold_residues_list,3],
                                 "FcFLS2" = FcFLS2_database[FcFLS2_database$coordinates %in% Fchold_residues_list,3],
                                 
                                 "FLS2Xl" = FLS2XL_database[FLS2XL_database$coordinates %in% Xlhold_residues_list,3],
                                 "QvFLS2" = QvFLS2_database[QvFLS2_database$coordinates %in% Qvhold_residues_list,3],
                                 "TjFLS2" = TjFLS2_database[TjFLS2_database$coordinates %in% Tjhold_residues_list,3])


split = rep(1:28, each = 8)
column_ha = HeatmapAnnotation("NR" = anno_boxplot(t(hold_charge_matrix[,1:5]), height = unit(2, "cm")))
column_ha2 = HeatmapAnnotation("BR" = anno_boxplot(t(hold_charge_matrix[,6:8]), height = unit(2, "cm")))

ComplexHeatmap::Heatmap(t(as.matrix(hold_charge_matrix)),
                        cluster_columns = FALSE,
                        show_column_names = FALSE,
                        column_split = split,
                        #row_km = 2, 
                        cluster_rows = TRUE,
                        row_split = c(1,1,1,1,1,2,2,2),
                        
                        top_annotation = column_ha,
                        bottom_annotation = column_ha2,
                        
                        name = "Bulkiness", border = TRUE,
                        col = circlize::colorRamp2(c(0, 15, 25), c("#3690c0", "#ece2f0", "#016450"))) 



#########################################################
# calculate residue propoeries - dim 2
#########################################################

AtFLS2_database <- calculate_protein_properties(AtFLS2_database, type = "charge")
SlFLS2_database <- calculate_protein_properties(SlFLS2_database, type = "charge")
FcFLS2_database <- calculate_protein_properties(FcFLS2_database, type = "charge")
VrFLS2_database <- calculate_protein_properties(VrFLS2_database, type = "charge")
NbFLS2_database <- calculate_protein_properties(NbFLS2_database, type = "charge")

FLS2XL_database <- calculate_protein_properties(FLS2XL_database, type = "charge")
QvFLS2_database <- calculate_protein_properties(QvFLS2_database, type = "charge")
TjFLS2_database <- calculate_protein_properties(TjFLS2_database, type = "charge")



hold_charge_matrix <- data.frame("AtFLS2" = AtFLS2_database[AtFLS2_database$coordinates %in% Athold_residues_list,4],
                                 "SlFLS2" = SlFLS2_database[SlFLS2_database$coordinates %in% Slhold_residues_list,4],
                                 "NbFLS2" = NbFLS2_database[NbFLS2_database$coordinates %in% Nbhold_residues_list,4],
                                 "FcFLS2" = FcFLS2_database[FcFLS2_database$coordinates %in% Fchold_residues_list,4],
                                 
                                 "VrFLS2" = VrFLS2_database[VrFLS2_database$coordinates %in% Vrhold_residues_list,4],
                                 "FLS2Xl" = FLS2XL_database[FLS2XL_database$coordinates %in% Xlhold_residues_list,4],
                                 "QvFLS2" = QvFLS2_database[QvFLS2_database$coordinates %in% Qvhold_residues_list,4],
                                 "TjFLS2" = TjFLS2_database[TjFLS2_database$coordinates %in% Tjhold_residues_list,4])


split = rep(1:28, each = 8)
column_ha = HeatmapAnnotation("NR" = anno_boxplot(t(hold_charge_matrix[,1:5]), height = unit(2, "cm")))
column_ha2 = HeatmapAnnotation("BR" = anno_boxplot(t(hold_charge_matrix[,6:8]), height = unit(2, "cm")))

ComplexHeatmap::Heatmap(t(as.matrix(hold_charge_matrix)),
                        cluster_columns = FALSE,
                        show_column_names = FALSE,
                        column_split = split,
                        #row_km = 2, 
                        cluster_rows = TRUE,
                        row_split = c(1,1,1,1,2,2,2,2),
                        
                        top_annotation = column_ha,
                        bottom_annotation = column_ha2,
                        
                        name = "Charge", border = TRUE,
                        col = circlize::colorRamp2(c(-1, 0, 1), c("#D92121", "#ece2f0", "#3690c0"))) 




#########################################################
# calculate residue propoeries - lower rank parameter - Welling
#########################################################

AtFLS2_database <- calculate_protein_properties(AtFLS2_database, type = "hydrophobicity", "Welling")
SlFLS2_database <- calculate_protein_properties(SlFLS2_database, type = "hydrophobicity", "Welling")
FcFLS2_database <- calculate_protein_properties(FcFLS2_database, type = "hydrophobicity", "Welling")
VrFLS2_database <- calculate_protein_properties(VrFLS2_database, type = "hydrophobicity", "Welling")
NbFLS2_database <- calculate_protein_properties(NbFLS2_database, type = "hydrophobicity", "Welling")

FLS2XL_database <- calculate_protein_properties(FLS2XL_database, type = "hydrophobicity", "Welling")
QvFLS2_database <- calculate_protein_properties(QvFLS2_database, type = "hydrophobicity", "Welling")
TjFLS2_database <- calculate_protein_properties(TjFLS2_database, type = "hydrophobicity", "Welling")



hold_charge_matrix <- data.frame("AtFLS2" = AtFLS2_database[AtFLS2_database$coordinates %in% Athold_residues_list,5],
                                 "SlFLS2" = SlFLS2_database[SlFLS2_database$coordinates %in% Slhold_residues_list,5],
                                 "NbFLS2" = NbFLS2_database[NbFLS2_database$coordinates %in% Nbhold_residues_list,5],
                                 "FcFLS2" = FcFLS2_database[FcFLS2_database$coordinates %in% Fchold_residues_list,5],
                                 "VrFLS2" = VrFLS2_database[VrFLS2_database$coordinates %in% Vrhold_residues_list,5],
                          
                                 "FLS2Xl" = FLS2XL_database[FLS2XL_database$coordinates %in% Xlhold_residues_list,5],
                                 "QvFLS2" = QvFLS2_database[QvFLS2_database$coordinates %in% Qvhold_residues_list,5],
                                 "TjFLS2" = TjFLS2_database[TjFLS2_database$coordinates %in% Tjhold_residues_list,5])


split = rep(1:28, each = 8)
#column_ha = HeatmapAnnotation("NR" = anno_boxplot(t(hold_charge_matrix[,1:5]), height = unit(2, "cm")))
#column_ha2 = HeatmapAnnotation("BR" = anno_boxplot(t(hold_charge_matrix[,6:8]), height = unit(2, "cm")))

ComplexHeatmap::Heatmap(t(as.matrix(hold_charge_matrix)),
                        cluster_columns = FALSE,
                        show_column_names = FALSE,
                        column_split = split,
                        #row_km = 2, 
                        cluster_rows = TRUE,
                        row_split = c(1,1,1,1,1,2,2,2),
                        
                        #top_annotation = column_ha,
                        #bottom_annotation = column_ha2,
                        
                        name = "Welling", border = TRUE,
                        col = circlize::colorRamp2(c(0,0.5, 1), c("white", "grey" ,"black"))) 








