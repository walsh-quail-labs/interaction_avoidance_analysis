%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% Replace occurences of
%%% Neighborhood_Individual_Image_Ourversion_Histocat_withNon function in
%%% run_randomized_control_perm_histocat_SQCC.m with this file
%%% to save intermidiated computed results 

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function save_nbhd_matrix(nuc_path,regInds_path,celltype_path,dist,matFilePath)

[Neighbor_Matrix,allLabels,allLabelsInd,nOfCells,cellTypes]=compute_neighborhood_matrix(nuc_path,regInds_path,celltype_path,dist);
save(matFilePath,'Neighbor_Matrix','allLabels','allLabelsInd','nOfCells','cellTypes');
end 