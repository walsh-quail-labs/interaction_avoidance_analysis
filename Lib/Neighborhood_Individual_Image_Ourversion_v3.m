function Neighborhood_Individual_Image_Ourversion_v3(nuc_path,regInds_path,celltype_path,dist,permutations,matFilePath)

[Neighbor_Matrix,allLabels,allLabelsInd,nOfCells,cellTypes]=compute_neighborhood_matrix(nuc_path,regInds_path,celltype_path,dist);


Neighbor_Matrix_index = Neighbor_Matrix+1;
Phenograph_Vector = zeros(nOfCells,1);
for i =1 : nOfCells
    curType = cellTypes{i};
    for j =1 : length(allLabels)
        if strcmp(curType,allLabels{j}) == 1
            Phenograph_Vector(i) = allLabelsInd(j);
        end
        if isempty(curType)
            Phenograph_Vector(i) = allLabelsInd(end);
        end
    end
end
Phenograph_Vector_index =[0;Phenograph_Vector];
nOfTypes = length(allLabels);

% Replace all neighbors with corresponding cell type
Phenograph_Neighor_Matrix = Phenograph_Vector_index(Neighbor_Matrix_index);

result_matrix_realdata = count_neighborhood_matrix(Phenograph_Neighor_Matrix,Phenograph_Vector_index,nOfTypes);
result_matrix_permdata = zeros(nOfTypes,nOfTypes,permutations);
for p=1:permutations
    combos_all_histcount_Perm_single = [];
    Phenograph_Vector_perm = Phenograph_Vector(randperm(length(Phenograph_Vector)));
    Phenograph_Vector_index_perm = [0;Phenograph_Vector_perm];
    Phenograph_Neighor_Matrix_perm = Phenograph_Vector_index_perm(Neighbor_Matrix_index);
    result_matrix_permdata(:,:,p) = count_neighborhood_matrix(Phenograph_Neighor_Matrix_perm,Phenograph_Vector_index_perm,nOfTypes);
end


% 
% % Calculate p-values
% % Get real data and permutated data
% real_data_mean = combos_all_histcount_real(:,3);
% perm_data_mean = combos_all_histcount_Perm(:,3:end);
% % Calculate amount higher or lower than mean (logic matrix)
% % What is the likelihood realdata is higher than random?


Higher_perm_test = zeros(size(result_matrix_permdata));
Lower_perm_test = zeros(size(result_matrix_permdata));
for p = 1 : permutations
    Higher_perm_test(:,:,p) = (result_matrix_permdata(:,:,p) >= result_matrix_realdata);
    Lower_perm_test(:,:,p) = (result_matrix_permdata(:,:,p) <= result_matrix_realdata);
end


% Higher_perm_test=repmat(real_data_mean,1,size(perm_data_mean,2))<=perm_data_mean;
% % What is the likelihood realdata is lower than random?
% Lower_perm_test=repmat(real_data_mean,1,size(perm_data_mean,2))>=perm_data_mean;




% % Calculate sum of lower and higherstat
Amount_higher = sum(Higher_perm_test,3);
Amount_lower = sum(Lower_perm_test,3);
% % Calculate actuall pValues
pValue_higher = 1 - (Amount_higher+1)/(permutations+1);
pValue_lower = 1 - (Amount_lower+1)/(permutations+1);


ML = pValue_lower;
ML(pValue_lower<0.01) = 1;
ML(pValue_lower>=0.01) = 0;
% imagesc(image_matrix_lower);
imagesc(ML);
colormap redblue
colorbar
caxis([0,1])
set(gca,'Xtick',1:length(allLabels),'Ytick',1:length(allLabels),'XtickLabel',allLabels,'YtickLabel',allLabels);

figure;
%avoidance
HL = pValue_higher;
HL(pValue_higher<0.01) = 1;
HL(pValue_higher>=0.01) = 0;
% imagesc(image_matrix_higher);
imagesc(HL);
colormap redblue
colorbar
caxis([0,1])
set(gca,'Xtick',1:length(allLabels),'Ytick',1:length(allLabels),'XtickLabel',allLabels,'YtickLabel',allLabels);

% 
% figure;
% imagesc(result_matrix_realdata)
% colormap redblue
% colorbar
% % caxis([0,1])
% set(gca,'Xtick',1:length(allLabels),'Ytick',1:length(allLabels),'XtickLabel',allLabels,'YtickLabel',allLabels);


save(matFilePath,'pValue_higher','pValue_lower','result_matrix_realdata','result_matrix_permdata');


