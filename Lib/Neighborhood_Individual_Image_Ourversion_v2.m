function Neighborhood_Individual_Image_Ourversion_v2(nuc_path,regInds_path,celltype_path,dist,permutations,matFilePath)

[Neighbor_Matrix,allLabels,allLabelsInd,nOfCells,cellTypes]=compute_neighborhood_matrix(nuc_path,regInds_path,celltype_path,dist);

Neighbor_Matrix_index = Neighbor_Matrix+1;



patch_det = 0;



combos_oneside = nchoosek(allLabelsInd,2);
combos_all = [combos_oneside;fliplr(combos_oneside);[allLabelsInd,allLabelsInd]];
combos_all_histcount = [combos_all,zeros(size(combos_all,1),1)];
Phenograph_Vector = zeros(nOfCells,1);
for i =1 : nOfCells
    curType = cellTypes{i};
    for j =1 : length(allLabels)
        if strcmp(curType,allLabels{j}) == 1
            Phenograph_Vector(i) = allLabelsInd(j);
        end
    end
end
Phenograph_Vector_index =[0;Phenograph_Vector];


% Replace all neighbors with corresponding cell type
Phenograph_Neighor_Matrix = Phenograph_Vector_index(Neighbor_Matrix_index);


% Run through all combos_all_histcount
[combos_all_histcount_real] = Calculate_STDandMean(combos_all,Phenograph_Neighor_Matrix,...
    Phenograph_Vector,Neighbor_Matrix,patch_det);

for p=1:permutations
    combos_all_histcount_Perm_single = [];Phenograph_Vector_perm = [];
    Phenograph_Vector_index_perm = [];
    % Generate matrix for permutation
    Phenograph_Vector_perm = Phenograph_Vector(randperm(length(Phenograph_Vector)));
    Phenograph_Vector_index_perm = [0;Phenograph_Vector_perm];
    % Replace all neighbors with corresponding cell type
    Phenograph_Neighor_Matrix_perm = Phenograph_Vector_index_perm(Neighbor_Matrix_index);
    % Run through all combos_all_histcount
    [combos_all_histcount_Perm_single] = Calculate_STDandMean(combos_all,Phenograph_Neighor_Matrix_perm,...
        Phenograph_Vector_perm,Neighbor_Matrix,patch_det);
    combos_all_histcount_Perm(:,p+2) = combos_all_histcount_Perm_single(:,3);
end


% Calculate p-values
% Get real data and permutated data
real_data_mean = combos_all_histcount_real(:,3);
perm_data_mean = combos_all_histcount_Perm(:,3:end);
% Calculate amount higher or lower than mean (logic matrix)
% What is the likelihood realdata is higher than random?
Higher_perm_test=repmat(real_data_mean,1,size(perm_data_mean,2))<=perm_data_mean;
% What is the likelihood realdata is lower than random?
Lower_perm_test=repmat(real_data_mean,1,size(perm_data_mean,2))>=perm_data_mean;
% Calculate sum of lower and higherstat
Amount_higher = sum(Higher_perm_test,2);
Amount_lower = sum(Lower_perm_test,2);
% Calculate actuall pValues
pValue_higher = (Amount_higher+1)/(permutations+1);
pValue_lower = (Amount_lower+1)/(permutations+1);

save(matFilePath,'pValue_higher','pValue_lower','real_data_mean','combos_all');


