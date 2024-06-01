function Neighborhood_Individual_Image_Ourversion_Histocat_withNon(nuc_path,regInds_path,celltype_path,dist,permutations,matFilePath)
tic 
[Neighbor_Matrix,allLabels,allLabelsInd,nOfCells,cellTypes]=compute_neighborhood_matrix(nuc_path,regInds_path,celltype_path,dist);

%save intermidiate results for perm in python 
% save('tmp/20220810_Lung 5Slides_SQCC5_ROI_020_AS-18-020058 C5_20.mat','Neighbor_Matrix','allLabels','allLabelsInd','nOfCells','cellTypes');


% test neighborhood matrix 
test_neighbor_matrix = 0;
if test_neighbor_matrix
    allRegionInds = importdata(regInds_path);
    margin = dist + 1;
    nuc_data = importdata(nuc_path);
    celltype_data = importdata(celltype_path);

    Boundaries = nuc_data.Boundaries;
    nucleiImage = nuc_data.nucleiImage;
    cellTypes = celltype_data.cellTypes;
    allLabels = celltype_data.allLabels;
    figure;
    imshow(nucleiImage)
    hold on;
    for i =1 : length(Boundaries)
        cur_boundary = Boundaries{i};
        [bx,by] = ind2sub(size(nucleiImage),cur_boundary);
        plot(by,bx,'-');
    end
    hold on;
    for i = 1 : size(Neighbor_Matrix,1)
        cur_row = Neighbor_Matrix(i,:);
        cur_row = nonzeros(cur_row);
        cur_boundary = Boundaries{i};
        [bx,by] = ind2sub(size(nucleiImage),cur_boundary);
        mbx = mean(bx);
        mby = mean(by);
        for j = 1 : length(cur_row)
            neighbor = cur_row(j);
            cur_boundary = Boundaries{neighbor};
            [nx,ny] = ind2sub(size(nucleiImage),cur_boundary);
            mnx = mean(nx);
            mny = mean(ny);
            plot([mny;mby],[mnx;mbx],'-');
        end
    end
end


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
        if isempty(curType)
            %NONE  
            indexNone = -1;
            for k = 1 : length(allLabels)
                if strcmp(allLabels{k},'NONE')
                    indexNone = k;
                end
            end
            Phenograph_Vector(i) = allLabelsInd(indexNone);
        end
    end
end
Phenograph_Vector_index =[0;Phenograph_Vector];


% Replace all neighbors with corresponding cell type
Phenograph_Neighor_Matrix = Phenograph_Vector_index(Neighbor_Matrix_index);


% Run through all combos_all_histcount
% [combos_all_histcount_real] = Calculate_STDandMean_mod(combos_all,Phenograph_Neighor_Matrix,...
%     Phenograph_Vector,Neighbor_Matrix,patch_det);
[combos_all_histcount_real] = Calculate_STDandMean(combos_all,Phenograph_Neighor_Matrix,...
    Phenograph_Vector,Neighbor_Matrix,patch_det);
toc 
disp('comp finished')
tic 
for p=1:permutations
    combos_all_histcount_Perm_single = [];Phenograph_Vector_perm = [];
    Phenograph_Vector_index_perm = [];


    % Generate matrix for permutation
%     seed_ = randperm(length(Phenograph_Vector));
%     Phenograph_Vector_perm = Phenograph_Vector(seed_);


    Phenograph_Vector_perm = Phenograph_Vector(randperm(length(Phenograph_Vector)));

    Phenograph_Vector_index_perm = [0;Phenograph_Vector_perm];
    % Replace all neighbors with corresponding cell type
    Phenograph_Neighor_Matrix_perm = Phenograph_Vector_index_perm(Neighbor_Matrix_index);
    % Run through all combos_all_histcount
%     [combos_all_histcount_Perm_single] = Calculate_STDandMean_mod(combos_all,Phenograph_Neighor_Matrix_perm,...
%         Phenograph_Vector_perm,Neighbor_Matrix,patch_det);
    [combos_all_histcount_Perm_single] = Calculate_STDandMean(combos_all,Phenograph_Neighor_Matrix_perm,...
        Phenograph_Vector_perm,Neighbor_Matrix,patch_det);
    combos_all_histcount_Perm(:,p+2) = combos_all_histcount_Perm_single(:,3);
     
end

toc
disp('perm finished')
tic 
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
toc 
% visualize_p_values_for_one_image
save(matFilePath,'pValue_higher','pValue_lower','real_data_mean','combos_all');


