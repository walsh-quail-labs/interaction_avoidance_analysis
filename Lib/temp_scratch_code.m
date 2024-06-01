
threshold = 0.95;
maxIndFile = floor(nOfFiles*threshold);
heatMapMatrixNormalize = zeros(length(allLabels),length(IndicatorsList));
for i = 1 : length(allLabels)
    for j = 1 :length(IndicatorsList)
        V = zeros(nOfFiles,1);
        B = zeros(nOfFiles,1);
        for indexFile = 1:nOfFiles
            A = allHeatmapMatrices{indexFile};
            B = allNumCellTypes{indexFile};
            V(indexFile) = A(i,j);
            B(indexFile) = B(i);

        end
        [VSorted,Inds] = sort(V);
        BSorted = B(Inds);
        heatMapMatrixNormalize(i,j) = sum(VSorted(1:maxIndFile));
    end
end

overallNumCellTypes = zeros(length(allLabels),1);
for indexFile = 1:nOfFiles
    if indexFile~=57 && indexFile~=110
        overallNumCellTypes = overallNumCellTypes + allNumCellTypes{indexFile};
    end
end

for ind = 1 : size(overallHeatMapMatrix,2)
    heatMapMatrixNormalize(:,ind) = heatMapMatrixNormalize(:,ind)./overallNumCellTypes;
end

for ind = 1 : size(overallHeatMapMatrix,2)
    heatMapMatrixNormalize(:,ind) = zscore(heatMapMatrixNormalize(:,ind));
    %heatMapMatrixNormalize(:,ind) = (overallHeatMapMatrix(:,ind));
end
heatMapMatrixNormalize = fix(heatMapMatrixNormalize*1000)/1000;
figure('units','normalized','outerposition',[0 0 1 1])
heatmap(IndicatorsList,allLabels,heatMapMatrixNormalize)
%     heatmap(IndicatorsList,allLabels,overallHeatMapMatrix)
colormap redblue
caxis([-5,5])