function Neighborhood_Individual_Image_Ourversion(nuc_path,celltype_path,dist,permutations,matFilePath)


nuc_data = importdata(nuc_path);
celltype_data = importdata(celltype_path);

cellTypes = celltype_data.cellTypes;
allLabels = celltype_data.allLabels;

for i = 1 : length(allLabels)
    if strcmp(allLabels{i},'Macrophages other') || strcmp(allLabels{i},'NONE')
        allLabels{i} = [];
    end
end
allLabels = allLabels(~cellfun(@isempty, allLabels));

Boundaries = nuc_data.Boundaries;
nucleiImage = nuc_data.nucleiImage;

nOfCells = length(Boundaries);  
N = length(allLabels);

ref_matrix = zeros(size(nucleiImage));
for i = 1: nOfCells  
    [xi,yi] = ind2sub(size(nucleiImage),Boundaries{i});
    mxi = round(mean(xi));
    myi = round(mean(yi));
    mxi = max(mxi,1);
    mxi = min(mxi,size(nucleiImage,1));
    myi = max(myi,1);
    myi = min(myi,size(nucleiImage,2));
    ref_matrix(mxi,myi) = i;    
end
list_neighbers = cell(nOfCells,1);
for i = 1: nOfCells
    binary_matrix = zeros(size(nucleiImage));
    binary_matrix(ref_matrix == i) = 1;
    binary_matrix = imdilate(binary_matrix,strel('disk',dist));
    cur_matrix = binary_matrix.*ref_matrix;   
    list_neighbers{i} = cur_matrix(cur_matrix ~= 0);
end
allLabelsInd = [1:N]';


Neighbor_Matrix = zeros(nOfCells,nOfCells);
for i= 1:nOfCells
    [xi,yi] = ind2sub(size(nucleiImage),Boundaries{i});
    mxi = mean(xi);
    myi = mean(yi);
    type_i = cellTypes{i};
    my_neigh = list_neighbers{i};
    for jj = 1:length(my_neigh)
        j = my_neigh(jj);
        if i~=j
            
            [xj,yj] = ind2sub(size(nucleiImage),Boundaries{j});
            mxj = mean(xj);
            myj = mean(yj);
            type_j = cellTypes{j};
            
            
            for k =1 : length(allLabels)
                if strcmp(type_j,allLabels{k}) == 1
                    type_j_index = allLabelsInd(k);
                end
            end
            
            
            
            if ~isempty(type_i) && ~isempty(type_j)
                if norm([mxi-mxj,myi-myj]) < dist
                    for itr = 1: length(allLabels)
                        if strcmp(allLabels{itr},type_i) == 1
                            cli = itr;
                            
                        end
                    end
                    for itr = 1: length(allLabels)
                        if strcmp(allLabels{itr},type_j) == 1
                            clj = itr;
                        end
                    end
                    Neighbor_Matrix(cli,clj) = type_j_index;
                    
                end
            end
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


