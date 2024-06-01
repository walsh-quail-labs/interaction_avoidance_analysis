% function heatMapMatrix = generate_heatmap_per_scan_filepath_with_config(files,i,panoDataSource,IndicatorsList,working_dir,celltype_date_path,allCellTypes)
function cellMeanIntensityMatrix = generate_cell_intensity_per_scan_filepath_with_config_mask_2(CombinedFlag,panoFilePath,allRegIndsPath,IndicatorsList,TypeList,JMF_folderpath,input_pano_folder)


if CombinedFlag == 0

    % read channels (Pano)
    [~,Pano] = imreadText(panoFilePath);

    
    for lt = 1:length(IndicatorsList) 

        lookupQuery = IndicatorsList{lt};
        headerIndex = find(contains(Pano.textdata,lookupQuery));        
        allChanelImages{lt} = imreadTextChannelRaw(Pano,headerIndex);  
        JMF_filePath = extractBetween(panoFilePath,length(input_pano_folder)+1,length(panoFilePath)-4);
        JMF_filePath = JMF_filePath{1};
        JMF_filePath = fullfile(JMF_folderpath,JMF_filePath,strcat(TypeList{lt},'_2.png'));
        allJMFilePaths{lt} = JMF_filePath;
    end


    % read nuclei
    allRegInds = importdata(allRegIndsPath);

    nOfCells = length(allRegInds);
    nOfInds = length(IndicatorsList);
    cellMeanIntensityMatrix = zeros(nOfCells,nOfInds);
    for i = 1 : nOfCells
        for j = 1 : nOfInds

            Ch = double(imread(allJMFilePaths{lt})).*allChanelImages{j};  
            if ~isempty(allRegInds{i})
                cellMeanIntensityMatrix(i,j) = mean(Ch(allRegInds{i}));
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
        
        
%         % read channels (Pano)
        [~,Pano] = imreadText(panoFilePath{fileITR});
% 
        for lt = 1:length(IndicatorsList) 
            lookupQuery = IndicatorsList{lt};
            headerIndex = find(contains(Pano.textdata,lookupQuery));
            allChanelImages{lt} = imreadTextChannelRaw(Pano,headerIndex);  
            JMF_filePath = extractBetween(panoFilePath,length(input_pano_folder)+1,length(panoFilePath)-4);
            JMF_filePath = JMF_filePath{1};
            JMF_filePath = fullfile(JMF_folderpath,JMF_filePath,strcat(TypeList{lt},'_2.png'));
            allJMFilePaths{lt} = JMF_filePath;
        end


        % read nuclei
        allRegInds = importdata(allRegIndsPath{fileITR});

        nOfCells = length(allRegInds);
        
        cur_cellMeanIntensityMatrix = zeros(nOfCells,nOfInds);
        for i = 1 : nOfCells
            for j = 1 : nOfInds
                Ch = double(imread(JMF_filePath{fileITR})).*allChanelImages{j};  
                if ~isempty(allRegInds{i})
                    cur_cellMeanIntensityMatrix(i,j) = mean(Ch(allRegInds{i}));
                end
            end
        end
        cellMeanIntensityMatrix = [cellMeanIntensityMatrix;cur_cellMeanIntensityMatrix];
        
        
    end
    cellMeanIntensityMatrix(1,:) = [];
    
%     
%     % read cell types
%     cellTypes = {};
%     for fileITR = 1 : length(cellTypeFilePath)
%         data_cellTypes = importdata(cellTypeFilePath{fileITR});
%         cur_cellTypes = data_cellTypes.cellTypes;
%         cur_cellTypes(cellfun(@isempty,cur_cellTypes)) = {'NONE'};
%         cellTypes = [cellTypes;cur_cellTypes];
%     end   
%     
%     
%     
%     nOfCellTypes = length(allCellTypes);
%     heatMapMatrix = zeros(nOfCellTypes,nOfInds);
%     for i =1:length(allCellTypes)
%         ccT = allCellTypes{i};
%         headerIndex = find(ismember(cellTypes,{ccT}));
%         if ~isempty(headerIndex)
%             heatMapMatrix(i,:) = mean(cellMeanIntensityMatrix(headerIndex,:),1);
%         end
% 
%     end
    
    

end


end