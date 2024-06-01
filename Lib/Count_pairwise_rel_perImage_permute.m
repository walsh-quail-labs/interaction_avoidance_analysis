function HM_perm = Count_pairwise_rel_perImage_permute(nuc_path,celltype_path,allLabels,dist,NofPerm)

% nuc_path = '/Users/elham/OneDrive - McGill University/SegmentationGBM/Data_BMR/BrM IMC run Dec 02/Slide6-Breast/Pano 07_7_ROI_17-236-C1_20/nuclei_multiscale.mat';
% celltype_path = '/Users/elham/OneDrive - McGill University/AxonNuclei_Brm/Main Files and Folders/CellTypes_March_Final/MATData/IMCData/BrainProject/Data_BMR/BrM IMC run Dec 02/Slide6-Breast/Pano 07_7_ROI_17-236-C1_20.mat';
% dist = 20;
% 
% allLabels = {'Cancer', 'B cell', 'Neutrophils', 'NK cell' , 'DCs cell', 'Endothelial cell', 'Mast cell', 'Astrocytes', 'Tc', 'Th', 'Treg', 'T other', 'Cl BMDM', 'Alt BMDM', 'Cl MG', 'Alt MG', 'Cl Mo', 'Non-Cl Mo', 'Int Mo','NONE'};

 

nuc_data = importdata(nuc_path);
celltype_data = importdata(celltype_path);

cellTypes = celltype_data.cellTypes;
Boundaries = nuc_data.Boundaries;
nucleiImage = nuc_data.nucleiImage;

nOfCells = length(Boundaries);  
allLabels(length(allLabels)) = [];
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

HM = zeros(N,N);
for perm= 1 : NofPerm
    HM_sub = zeros(N,N);
    new_cellTypes = cellTypes(randperm(nOfCells));
    
    for i= 1:nOfCells
        [xi,yi] = ind2sub(size(nucleiImage),Boundaries{i});
        mxi = mean(xi);
        myi = mean(yi);
        type_i = new_cellTypes{i};
        my_neigh = list_neighbers{i};
        for jj = 1:length(my_neigh)
            j = my_neigh(jj);
            if i~=j
                
                [xj,yj] = ind2sub(size(nucleiImage),Boundaries{j});
                mxj = mean(xj);
                myj = mean(yj);
                type_j = new_cellTypes{j};
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
                        HM(cli,clj) = HM(cli,clj)+1;
                        
                    end
                end
            end
            
        end
    end
end

HM_perm = HM/NofPerm;





