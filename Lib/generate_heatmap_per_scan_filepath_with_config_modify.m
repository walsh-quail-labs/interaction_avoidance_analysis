% function heatMapMatrix = generate_heatmap_per_scan_filepath_with_config(files,i,panoDataSource,IndicatorsList,working_dir,celltype_date_path,allCellTypes)
function heatMapMatrix = generate_heatmap_per_scan_filepath_with_config_modify(CombinedFlag,panoFilePath,cellTypeFilePath,allRegIndsPath,IndicatorsList,allCellTypes,rad)


if CombinedFlag == 0

    % read channels (Pano)
    [~,Pano] = imreadText(panoFilePath);

    for lt = 1:length(IndicatorsList) 
        lookupQuery = IndicatorsList{lt};
        headerIndex = find(contains(Pano.textdata,lookupQuery));
        allChanelImages{lt} = imreadTextChannelRaw(Pano,headerIndex);  
    end


    % read nuclei
    allRegInds = importdata(allRegIndsPath);

    nOfCells = length(allRegInds);
    nOfInds = length(IndicatorsList);
    cellMeanIntensityMatrix = zeros(nOfCells,nOfInds);
    for i = 1 : nOfCells
        for j = 1 : nOfInds
            Ch = allChanelImages{j};  
            if ~isempty(allRegInds{i})
%                 im_rad = zeros(size(Ch));
%                 im_rad(allRegInds{i})=1;
%                 if rad > 0
%                     im_rad = imdilate(im_rad,strel('sphere',rad));
%                 else
%                     im_rad = imerode(im_rad,strel('sphere',-rad));
%                 end
%                 finalRadInds = im_rad==1;
                finalRadInds = dilate_erode_patch(allRegInds{i},size(Ch),rad);
                if  ~isempty(finalRadInds)
                    cellMeanIntensityMatrix(i,j) = mean(Ch(finalRadInds));
                end
%                 cellMeanIntensityMatrix(i,j) = mean(Ch(allRegInds{i}));
            end
        end
    end

    % read cell types
    data_cellTypes = importdata(cellTypeFilePath);
    cellTypes = data_cellTypes.cellTypes;
    cellTypes(cellfun(@isempty,cellTypes)) = {'NONE'};


    nOfCellTypes = length(allCellTypes);
    heatMapMatrix = zeros(nOfCellTypes,nOfInds);
    for i =1:length(allCellTypes)
        ccT = allCellTypes{i};
        headerIndex = find(ismember(cellTypes,{ccT}));
        if ~isempty(headerIndex)
            heatMapMatrix(i,:) = mean(cellMeanIntensityMatrix(headerIndex,:),1);
        end

    end
    
else
    
    nOfInds = length(IndicatorsList);
    
    cellMeanIntensityMatrix = zeros(1,nOfInds);
    
    for fileITR = 1 : length(panoFilePath)
        
        
        % read channels (Pano)
        [~,Pano] = imreadText(panoFilePath{fileITR});

        for lt = 1:length(IndicatorsList) 
            lookupQuery = IndicatorsList{lt};
            headerIndex = find(contains(Pano.textdata,lookupQuery));
            allChanelImages{lt} = imreadTextChannelRaw(Pano,headerIndex);  
        end


        % read nuclei
        allRegInds = importdata(allRegIndsPath{fileITR});

        nOfCells = length(allRegInds);
        
        cur_cellMeanIntensityMatrix = zeros(nOfCells,nOfInds);
        for i = 1 : nOfCells
            for j = 1 : nOfInds
                Ch = allChanelImages{j};
%                 im_rad = zeros(size(Ch));
%                 im_rad(allRegInds{i})=1;
%                 if rad > 0
%                     im_rad = imdilate(im_rad,strel('sphere',rad));
%                 else
%                     im_rad = imerode(im_rad,strel('sphere',-rad));
%                 end
%                 finalRadInds = im_rad==1;
                finalRadInds = dilate_erode_patch(allRegInds{i},size(Ch),rad);
                if ~isempty(finalRadInds)
                    cellMeanIntensityMatrix(i,j) = mean(Ch(finalRadInds));
%                     cur_cellMeanIntensityMatrix(i,j) = mean(Ch(allRegInds{i}));
                end
            end
        end
        cellMeanIntensityMatrix = [cellMeanIntensityMatrix;cur_cellMeanIntensityMatrix];
        
        
    end
    cellMeanIntensityMatrix(1,:) = [];
    
    
    % read cell types
    cellTypes = {};
    for fileITR = 1 : length(cellTypeFilePath)
        data_cellTypes = importdata(cellTypeFilePath{fileITR});
        cur_cellTypes = data_cellTypes.cellTypes;
        cur_cellTypes(cellfun(@isempty,cur_cellTypes)) = {'NONE'};
        cellTypes = [cellTypes;cur_cellTypes];
    end   
    
    
    
    nOfCellTypes = length(allCellTypes);
    heatMapMatrix = zeros(nOfCellTypes,nOfInds);
    for i =1:length(allCellTypes)
        ccT = allCellTypes{i};
        headerIndex = find(ismember(cellTypes,{ccT}));
        if ~isempty(headerIndex)
            heatMapMatrix(i,:) = mean(cellMeanIntensityMatrix(headerIndex,:),1);
        end

    end
    
    

end


end