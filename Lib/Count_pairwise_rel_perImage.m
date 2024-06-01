function HM = Count_pairwise_rel_perImage(nuc_path,celltype_path,allLabels,dist)

nuc_data = importdata(nuc_path);
celltype_data = importdata(celltype_path);

cellTypes = celltype_data.cellTypes;
Boundaries = nuc_data.Boundaries;
nucleiImage = nuc_data.nucleiImage;

allLabels(length(allLabels)) = [];
N = length(allLabels);
HM = zeros(N,N);


for i= 1:length(Boundaries)
    [xi,yi] = ind2sub(size(nucleiImage),Boundaries{i});
    mxi = mean(xi);
    myi = mean(yi);
    type_i = cellTypes{i};
    for j = 1:length(Boundaries)
        if i~=j
            
            [xj,yj] = ind2sub(size(nucleiImage),Boundaries{j});            
            mxj = mean(xj);
            myj = mean(yj);            
            type_j = cellTypes{j};
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




