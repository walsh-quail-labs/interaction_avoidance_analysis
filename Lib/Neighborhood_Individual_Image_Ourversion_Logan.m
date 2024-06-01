function Neighborhood_Individual_Image_Ourversion_Logan(nuc_path,regInds_path,celltype_path,dist,permutations,matFilePath)

[Neighbor_Matrix,allLabels,allLabelsInd,nOfCells,cellTypes]=compute_neighborhood_matrix(nuc_path,regInds_path,celltype_path,dist);


Neighbor_Matrix_index = Neighbor_Matrix+1;
% 
% 
% 
% patch_det = 0;
% 
% 
% 
% combos_oneside = nchoosek(allLabelsInd,2);
% combos_all = [combos_oneside;fliplr(combos_oneside);[allLabelsInd,allLabelsInd]];
% combos_all_histcount = [combos_all,zeros(size(combos_all,1),1)];


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

% 
% 
% Replace all neighbors with corresponding cell type
Phenograph_Neighor_Matrix = Phenograph_Vector_index(Neighbor_Matrix_index);
% 
% 
% % Run through all combos_all_histcount
% [combos_all_histcount_real] = Calculate_STDandMean(combos_all,Phenograph_Neighor_Matrix,...
%     Phenograph_Vector,Neighbor_Matrix,patch_det);
% 

result_matrix_realdata = count_neighborhood_matrix_Logan(Phenograph_Neighor_Matrix,Phenograph_Vector_index,nOfTypes);
result_matrix_permdata = zeros(nOfTypes,nOfTypes,permutations);
% permutations = 1000;
for p=1:permutations
    combos_all_histcount_Perm_single = [];
    % Generate matrix for permutation
    Phenograph_Vector_perm = Phenograph_Vector(randperm(length(Phenograph_Vector)));
    Phenograph_Vector_index_perm = [0;Phenograph_Vector_perm];
%     % Replace all neighbors with corresponding cell type
    Phenograph_Neighor_Matrix_perm = Phenograph_Vector_index_perm(Neighbor_Matrix_index);
    result_matrix_permdata(:,:,p) = count_neighborhood_matrix_Logan(Phenograph_Neighor_Matrix_perm,Phenograph_Vector_index_perm,nOfTypes);
%     % Run through all combos_all_histcount
%     [combos_all_histcount_Perm_single] = Calculate_STDandMean(combos_all,Phenograph_Neighor_Matrix_perm,...
%         Phenograph_Vector_perm,Neighbor_Matrix,patch_det);
%     combos_all_histcount_Perm(:,p+2) = combos_all_histcount_Perm_single(:,3);
end

mean_permute = mean(result_matrix_permdata,3);

%ti = log(result_matrix_realdata./mean_permute);
ti = result_matrix_realdata./mean_permute;


% 
% st = find(ti < 1);
% ti(st) = -1./ti(st);


% figure;
% %avoidance
% imagesc(ti);
% colormap redblue
% colorbar
% caxis([-2,2])
% set(gca,'Xtick',1:length(allLabels),'Ytick',1:length(allLabels),'XtickLabel',allLabels,'YtickLabel',allLabels);

save(matFilePath,'ti','result_matrix_realdata','result_matrix_permdata');


