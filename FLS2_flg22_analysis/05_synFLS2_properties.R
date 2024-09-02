#-----------------------------------------------------------------------------------------------
# Coaker Lab - Plant Pathology Department UC Davis
# Author: Danielle M. Stevens
# Last Updated: 08/21/2024
# Script Purpose: Compare specific FLS2-flg22 pairs to understand why some receptors have expanded ligand receptors
# Inputs Necessary: flg22 sequences, FLS2 protein sequences and residue coordinate map 
# Outputs: 
#-----------------------------------------------------------------------------------------------



# ------------------------------------------------FLS2Xl vs. Vr FLS2 ---------------------------------------------------------

#########################################################
# load and prep certain chemical properties of synthetic FLS variants
#########################################################

# load Vr synthetic variants
synVrFLS2_15 <- aa2df_v2(Biostrings::readAAStringSet("../FLS2_flg22_coordinates/synVrFLS2_15.fasta"))
synVrFLS2_22 <- aa2df_v2(Biostrings::readAAStringSet("../FLS2_flg22_coordinates/synVrFLS2_22.fasta"))
synVrFLS2_30 <- aa2df_v2(Biostrings::readAAStringSet("../FLS2_flg22_coordinates/synVrFLS2_30.fasta"))


# load Vr database
synVrFLS2_database_15 <- do.call(rbind, Map(data.frame, coordinates = 1:nchar(synVrFLS2_15$seq), seq = strsplit(synVrFLS2_15$seq, "")[[1]]))
synVrFLS2_database_22 <- do.call(rbind, Map(data.frame, coordinates = 1:nchar(synVrFLS2_22$seq), seq = strsplit(synVrFLS2_22$seq, "")[[1]]))
synVrFLS2_database_30 <- do.call(rbind, Map(data.frame, coordinates = 1:nchar(synVrFLS2_30$seq), seq = strsplit(synVrFLS2_30$seq, "")[[1]]))


#########################################################
# calculate residue propoeries - dim 1
#########################################################


synVrFLS2_database_15 <- calculate_protein_properties(synVrFLS2_database_15, type = "bulk")
synVrFLS2_database_22 <- calculate_protein_properties(synVrFLS2_database_22, type = "bulk")
synVrFLS2_database_30 <- calculate_protein_properties(synVrFLS2_database_30, type = "bulk")


hold_charge_matrix <- data.frame("FLS2Xl" = FLS2XL_database[FLS2XL_database$coordinates %in% Xlhold_residues_list,3],
                                 "VrFLS2" = VrFLS2_database[VrFLS2_database$coordinates %in% Vrhold_residues_list,3],
                                 "synVrFLS2_22" = synVrFLS2_database_22[synVrFLS2_database_22$coordinates %in% Vrhold_residues_list,3]
                                 )

# hold values for property comparisons with flg22
#flg22_FLS2_property_comparisons <- hold_charge_matrix

split = rep(1:28, each = 8)
column_ha = HeatmapAnnotation("Average\n Bulk" = anno_boxplot(t(hold_charge_matrix), 
                                                              annotation_name_gp= gpar(fontsize = 8), height = unit(1.3, "cm")))

ComplexHeatmap::Heatmap(t(as.matrix(hold_charge_matrix)),
                        cluster_columns = FALSE,
                        show_column_names = FALSE,
                        column_split = split,
                        cluster_rows = FALSE,
                        column_title_gp = gpar(fontsize = 8), 
                        row_names_gp = grid::gpar(fontsize = 8),
                        #top_annotation = column_ha,
                        
                        name = "Bulkiness", border = TRUE, 
                        col = circlize::colorRamp2(c(0, 12, 25), c("#3690c0", "#ece2f0", "#016450")))


# export at about 9 x 1 in

#########################################################
# calculate residue propoeries - dim 2
#########################################################


synVrFLS2_database_15 <- calculate_protein_properties(synVrFLS2_database_15, type = "charge")
synVrFLS2_database_22 <- calculate_protein_properties(synVrFLS2_database_22, type = "charge")
synVrFLS2_database_30 <- calculate_protein_properties(synVrFLS2_database_30, type = "charge")

hold_charge_matrix <- data.frame("FLS2Xl" = FLS2XL_database[FLS2XL_database$coordinates %in% Xlhold_residues_list,4],
                                 "VrFLS2" = VrFLS2_database[VrFLS2_database$coordinates %in% Vrhold_residues_list,4],
                                 "synVrFLS2_22" = synVrFLS2_database_22[synVrFLS2_database_22$coordinates %in% Vrhold_residues_list,4]

)

# hold values for property comparisons with flg22
#flg22_FLS2_property_comparisons <- hold_charge_matrix

column_ha = HeatmapAnnotation("Average\n Charge" = anno_boxplot(t(hold_charge_matrix), 
                                                                annotation_name_gp= gpar(fontsize = 8), height = unit(1.3, "cm")))

ComplexHeatmap::Heatmap(t(as.matrix(hold_charge_matrix)),
                        cluster_columns = FALSE,
                        show_column_names = FALSE,
                        column_split = split,
                        cluster_rows = FALSE,
                        column_title_gp = gpar(fontsize = 8), 
                        row_names_gp = grid::gpar(fontsize = 8),
                        #top_annotation = column_ha,
                        
                        name = "Charge", border = TRUE,
                        col = circlize::colorRamp2(c(-1, 0, 1), c("#D92121", "#ece2f0", "#3690c0"))) 




# export at about 9 x 1 in


#########################################################
# calculate residue propoeries - lower rank parameter - Welling
#########################################################


synVrFLS2_database_15 <- calculate_protein_properties(synVrFLS2_database_15, type = "hydrophobicity", "Welling")
synVrFLS2_database_22 <- calculate_protein_properties(synVrFLS2_database_22, type = "hydrophobicity", "Welling")
synVrFLS2_database_30 <- calculate_protein_properties(synVrFLS2_database_30, type = "hydrophobicity", "Welling")

hold_charge_matrix <- data.frame("FLS2Xl" = FLS2XL_database[FLS2XL_database$coordinates %in% Xlhold_residues_list,5],
                                 "VrFLS2" = VrFLS2_database[VrFLS2_database$coordinates %in% Vrhold_residues_list,5],
                                 "synVrFLS2_22" = synVrFLS2_database_22[synVrFLS2_database_22$coordinates %in% Vrhold_residues_list,5]
                                 
)

# hold values for property comparisons with flg22
#flg22_FLS2_property_comparisons <- hold_charge_matrix

#column_ha = HeatmapAnnotation("Average\n Charge" = anno_boxplot(t(hold_charge_matrix), 
#                                                                annotation_name_gp= gpar(fontsize = 8), height = unit(1.3, "cm")))

ComplexHeatmap::Heatmap(t(as.matrix(hold_charge_matrix)),
                        cluster_columns = FALSE,
                        show_column_names = FALSE,
                        column_split = split,
                        cluster_rows = FALSE,
                        column_title_gp = gpar(fontsize = 8), 
                        row_names_gp = grid::gpar(fontsize = 8),
                        #top_annotation = column_ha,
                        
                        name = "Welling", border = TRUE,
                        col = circlize::colorRamp2(c(-3.5,0.5, 3.5), c("white", "grey" ,"black"))) 




# export at about 9 x 1 in



# ------------------------------------------------Qv vs. Fc FLS2 ---------------------------------------------------------


#########################################################
# load and prep certain chemical properties of synthetic FLS variants
#########################################################


# load Vr synthetic variants
synFcFLS2_31 <- aa2df_v2(Biostrings::readAAStringSet("../FLS2_flg22_coordinates/synFcFLS2_31.fasta"))
synFcFLS2_13 <- aa2df_v2(Biostrings::readAAStringSet("../FLS2_flg22_coordinates/synFcFLS2_13.fasta"))
synFcFLS2_15 <- aa2df_v2(Biostrings::readAAStringSet("../FLS2_flg22_coordinates/synFcFLS2_15.fasta"))


# load Vr database
synFcFLS2_database_31 <- do.call(rbind, Map(data.frame, coordinates = 1:nchar(synFcFLS2_31$seq), seq = strsplit(synFcFLS2_31$seq, "")[[1]]))
synFcFLS2_database_13 <- do.call(rbind, Map(data.frame, coordinates = 1:nchar(synFcFLS2_13$seq), seq = strsplit(synFcFLS2_13$seq, "")[[1]]))
synFcFLS2_database_15 <- do.call(rbind, Map(data.frame, coordinates = 1:nchar(synFcFLS2_15$seq), seq = strsplit(synFcFLS2_15$seq, "")[[1]]))


#########################################################
# calculate residue propoeries - dim 1
#########################################################


synFcFLS2_database_31 <- calculate_protein_properties(synFcFLS2_database_31, type = "bulk")
synFcFLS2_database_13 <- calculate_protein_properties(synFcFLS2_database_13, type = "bulk")
synFcFLS2_database_15 <- calculate_protein_properties(synFcFLS2_database_15, type = "bulk")


hold_charge_matrix <- data.frame("QvFLS2" = QvFLS2_database[QvFLS2_database$coordinates %in% Qvhold_residues_list,3],
                                 "FcFLS2" = FcFLS2_database[FcFLS2_database$coordinates %in% Fchold_residues_list,3],
                                 "synFcFLS2_13" = synFcFLS2_database_13[synFcFLS2_database_13$coordinates %in% Fchold_residues_list,3]

                                 
)

# hold values for property comparisons with flg22
#flg22_FLS2_property_comparisons <- cbind(flg22_FLS2_property_comparisons, hold_charge_matrix)

split = rep(1:28, each = 8)
column_ha = HeatmapAnnotation("Average\n Bulk" = anno_boxplot(t(hold_charge_matrix), 
                                                                annotation_name_gp= gpar(fontsize = 8), height = unit(1.3, "cm")))

ComplexHeatmap::Heatmap(t(as.matrix(hold_charge_matrix)),
                        cluster_columns = FALSE,
                        show_column_names = FALSE,
                        column_split = split,
                        cluster_rows = FALSE,
                        column_title_gp = gpar(fontsize = 8), 
                        row_names_gp = grid::gpar(fontsize = 8),
                        #top_annotation = column_ha,
                        
                        name = "Bulkiness", border = TRUE,
                        col = circlize::colorRamp2(c(0, 12, 25), c("#3690c0", "#ece2f0", "#016450")))



# export at about 9 x 1 in


#########################################################
# calculate residue propoeries - dim 2
#########################################################


synFcFLS2_database_31 <- calculate_protein_properties(synFcFLS2_database_31, type = "charge")
synFcFLS2_database_13 <- calculate_protein_properties(synFcFLS2_database_13, type = "charge")
synFcFLS2_database_15 <- calculate_protein_properties(synFcFLS2_database_15, type = "charge")


hold_charge_matrix <- data.frame("QvFLS2" = QvFLS2_database[QvFLS2_database$coordinates %in% Qvhold_residues_list,4],
                                 "FcFLS2" = FcFLS2_database[FcFLS2_database$coordinates %in% Fchold_residues_list,4],
                                 "synFcFLS2_13" = synFcFLS2_database_13[synFcFLS2_database_13$coordinates %in% Fchold_residues_list,4]

                                 
)

ComplexHeatmap::Heatmap(t(as.matrix(hold_charge_matrix)),
                        cluster_columns = FALSE,
                        show_column_names = FALSE,
                        column_split = split,
                        cluster_rows = FALSE,
                        column_title_gp = gpar(fontsize = 8), 
                        row_names_gp = grid::gpar(fontsize = 8),
                        #top_annotation = column_ha,
                        
                        name = "Charge", border = TRUE,
                        col = circlize::colorRamp2(c(-1, 0, 1), c("#D92121", "#ece2f0", "#3690c0"))) 





# export at about 9 x 1 in



