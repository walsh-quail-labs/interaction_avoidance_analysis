clear all
close all
nuc_path = '/Users/elham/OneDrive - McGill University/SegmentationGBM/Data_BMR/BrM IMC run Dec 02/Slide6-Breast/Pano 07_7_ROI_17-236-C1_20/nuclei_multiscale.mat';
celltype_path = '/Users/elham/OneDrive - McGill University/AxonNuclei_Brm/Main Files and Folders/CellTypes_March_Final/MATData/IMCData/BrainProject/Data_BMR/BrM IMC run Dec 02/Slide6-Breast/Pano 07_7_ROI_17-236-C1_20.mat';

nuc_data = importdata(nuc_path);
celltype_data = importdata(celltype_path);

cellTypes = celltype_data.cellTypes;
Boundaries = nuc_data.Boundaries;
nucleiImage = nuc_data.nucleiImage;
dist = 20;

%allLabels = {'Cancer','B Cell','T Cell','NONE'};
%N = length(allLabels)-1;
%HM = zeros(N,N);
%type_i = 'Cancer';
%type_j = 'T Cell';

%for itr = 1 : length(allLabels)
%    if strcmp(allLabels{itr},type_i) == 1
%        cli = itr;
%    end
%end
%for itr = 1 : length(allLabels)
%    if strcmp(allLabels{itr},type_j) == 1
%        clj = itr;
%    end
%end
%HM(cli,clj) = HM(cli,clj)+1;



allLabels = {'Cancer', 'B cell', 'Neutrophils', 'NK cell' , 'DCs cell', 'Endothelial cell', 'Mast cell', 'Astrocytes', 'Tc', 'Th', 'Treg', 'T other', 'Cl BMDM', 'Alt BMDM', 'Cl MG', 'Alt MG', 'Cl Mo', 'Non-Cl Mo', 'Int Mo'};
N = length(allLabels);
HM = zeros(N,N);

figure;
imshow(nucleiImage)


for i= 1:length(Boundaries)
    [xi,yi] = ind2sub(size(nucleiImage),Boundaries{i});
    mxi = mean(xi);
    myi = mean(yi);
    circle(myi,mxi,20)
%     type_i = cellTypes{i};
%     for j = 1:length(Boundaries)
%         if i~=j
%             
%             [xj,yj] = ind2sub(size(nucleiImage),Boundaries{j});            
%             mxj = mean(xj);
%             myj = mean(yj);            
%             type_j = cellTypes{j};
%             if ~isempty(type_i) && ~isempty(type_j)
%                 if norm([mxi-mxj,myi-myj]) < dist
%                    for itr = 1: length(allLabels)
%                       if strcmp(allLabels{itr},type_i) == 1
%                           cli = itr;
% 
%                       end
%                    end
%                    for itr = 1: length(allLabels)
%                       if strcmp(allLabels{itr},type_j) == 1
%                           clj = itr;
%                       end
%                    end
%                    HM(cli,clj) = HM(cli,clj)+1;
% 
%                 end
%             end
%         end
%         
%     end
end
% 
% figure('units','normalized','outerposition',[0 0 1 1])
% heatmap(allLabels,allLabels,HM)
% 

% Boundaries{i} => indices to [x,y]
% [xi,yi] = ind2sub(size(nucleiImage),Boundaries{i});
% [xi,yi] = ind2sub(size(nucleiImage),Boundaries{i}); % repeat j
% mxi = mean(xi);
% myi = mean(yi)

% sqrt((x1-x2)*(x1-x2)+(y1-y2)*(y1-y2)) === norm([x1-x2,y1-y2])

