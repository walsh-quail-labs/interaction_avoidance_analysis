function [Neighbor_Matrix,allLabels,allLabelsInd,nOfCells,cellTypes]=compute_neighborhood_matrix(nuc_path,regInds_path,celltype_path,dist)
allRegionInds = importdata(regInds_path);
margin = dist + 1;
nuc_data = importdata(nuc_path);
celltype_data = importdata(celltype_path);

Boundaries = nuc_data.Boundaries;
nucleiImage = nuc_data.nucleiImage;
if min(string(class(celltype_data))=="struct")~=1
    % load list of labels if they are not saved along the celltype
    disp("loading nada labels")
    cellTypes = celltype_data;
    allLabels = load("ConfigFiles/nada/allLabels.mat");
    allLabels = allLabels.allLabels;
else
    cellTypes = celltype_data.cellTypes;
    allLabels = celltype_data.allLabels;
end 

for i = 1 : length(allLabels)
    if strcmp(allLabels{i},'Macrophages other') %|| strcmp(allLabels{i},'NONE')
        allLabels{i} = [];
    end
end
allLabels = allLabels(~cellfun(@isempty, allLabels));
allLabels
nOfCells = length(Boundaries);  

N = length(allLabels);

imsize = size(nucleiImage);
allRegionIndsExpanded = cell(nOfCells,1);

ifCheck = 0;
if ifCheck
    BW1 = zeros(imsize);
    BW2 = zeros(imsize);
end 

for i = 1 : nOfCells
    [rx,ry] = ind2sub(imsize,allRegionInds{i});
    
    minx = min(rx);
    maxx = max(rx);
    miny = min(ry);
    maxy = max(ry);
    
    rx = rx-minx+margin;
    ry = ry-miny+margin;
    w = maxx-minx+2*margin;
    h = maxy-miny+2*margin;
    BW = zeros(w,h);
    BW(sub2ind(size(BW),rx,ry)) = 1;
%     BW(rx,ry) = 1;
%     se = strel('sphere',25);

    BW = imdilate(BW,strel('disk',dist));
    [nrx,nry] = ind2sub(size(BW),find(BW~=0));
    nrx = nrx+minx-margin;
    nry = nry+miny-margin;
    nrx = min(imsize(1),max(nrx,1));
    nry = min(imsize(2),max(nry,1));
    allRegionIndsExpanded{i} = sub2ind(imsize,nrx,nry);
    if ifCheck
        BW1(allRegionInds{i}) = 1;
        BW2(allRegionIndsExpanded{i}) = 1;
    end 
end
if ifCheck
    figure
    imshow(BW1)
    figure
    imshow(BW2)
end

intersection_matrix = zeros(nOfCells,nOfCells);

for i = 1 : nOfCells
    [xi,yi] = ind2sub(imsize,Boundaries{i});
    mxi = round(mean(xi));
    myi = round(mean(yi));
    mxi = max(mxi,1);
    mxi = min(mxi,size(nucleiImage,1));
    myi = max(myi,1);
    myi = min(myi,size(nucleiImage,2));
    for j = i+1 : nOfCells
        [xj,yj] = ind2sub(imsize,Boundaries{j});
        mxj = round(mean(xj));
        myj = round(mean(yj));
        mxj = max(mxj,1);
        mxj = min(mxj,size(nucleiImage,1));
        myj = max(myj,1);
        myj = min(myj,size(nucleiImage,2));
        if norm([mxj,myj]-[mxi,myi]) <= 8*dist
            A = allRegionIndsExpanded{i}';
            B = allRegionIndsExpanded{j}';
            CNum = length(intersect(A,B));
            intersection_matrix(i,j) = CNum;
            intersection_matrix(j,i) = CNum;
        end        
    end
end
intersection_cell = cell(nOfCells,1);
maxLength = 0;
for i = 1 : nOfCells
%     intersection_cell{i} = nonzeros(intersection_matrix(i,:));
    intersection_cell{i} = find(intersection_matrix(i,:)~=0);
    maxLength = max(length(intersection_cell{i}),maxLength);
end

Neighbor_Matrix = zeros(nOfCells,maxLength);

for i = 1 : nOfCells
    neigh_i = intersection_cell{i};
    for j = 1 : length(neigh_i)
        Neighbor_Matrix(i,j) = neigh_i(j);
    end
end


% 
% % neghiborhood_matrix = imbinarize(intersection_matrix);
% % average_connectivity = mean(sum(neghiborhood_matrix));
% 
% Neighbor_Matrix = zeros(nOfCells,length(allLabels));
allLabelsInd = [1:N]';
% 
% for i= 1:nOfCells
%     nOfCurrentInds = length(allRegionIndsExpanded{i});
%     
%     for j = 1:nOfCells
%         type_j = cellTypes{j};
%         for k =1 : length(allLabels)
%             if strcmp(type_j,allLabels{k}) == 1
%                 type_j_index = allLabelsInd(k);
%             end
%         end
%         Neighbor_Matrix(i,type_j_index) = Neighbor_Matrix(i,type_j_index)+(intersection_matrix(i,j)/nOfCurrentInds);
%     end
% end
end
