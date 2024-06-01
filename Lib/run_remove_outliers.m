clear all;
clc;
addpath(genpath('Lib'))
config_file_path = 'template_ICI3.csv';
[~,configFileName,~] = fileparts(config_file_path);

CSVData = importdata(fullfile('ConfigFiles',config_file_path));


k = strfind(CSVData{1},'input_pano_folder');
input_pano_folder = extractBetween(CSVData{1},k+length('input_pano_folder')+1,length(CSVData{1}));
input_pano_folder = input_pano_folder{1};

k = strfind(CSVData{2},'input_cell_segmentation_folder');
input_cell_segmentation_folder = extractBetween(CSVData{2},k+length('input_cell_segmentation_folder')+1,length(CSVData{2}));
input_cell_segmentation_folder = input_cell_segmentation_folder{1};


k = strfind(CSVData{3},'input_cell_type_folder');
input_celltype_folder = extractBetween(CSVData{3},k+length('input_cell_type_folder')+1,length(CSVData{3}));
input_celltype_folder = input_celltype_folder{1};


k = strfind(CSVData{4},'input_keyPatientID');
input_keyPatientID = extractBetween(CSVData{4},k+length('input_keyPatientID')+1,length(CSVData{4}));
input_keyPatientID = input_keyPatientID{1};

k = strfind(CSVData{5},'group_PatientID');
group_PatientID = extractBetween(CSVData{5},k+length('group_PatientID')+1,length(CSVData{5}));
group_PatientID = group_PatientID{1};

k = strfind(CSVData{6},'allCellTypeLocation');
allCellTypesAddress = extractBetween(CSVData{6},k+length('allCellTypeLocation')+1,length(CSVData{6}));
allCellTypesAddress = allCellTypesAddress{1};
allLabelsTable = readtable(allCellTypesAddress);
allLabels = allLabelsTable.Variables;



k = strfind(CSVData{7},'output_folder');
output_folder = extractBetween(CSVData{7},k+length('output_folder')+1,length(CSVData{7}));
output_folder = output_folder{1};
mkdir(output_folder);


k = strfind(CSVData{8},'IndicatorsListPath');
IndicatorsListPath = extractBetween(CSVData{8},k+length('IndicatorsListPath')+1,length(CSVData{8}));
IndicatorsListPath = IndicatorsListPath{1};

IndicatorsListTable = readtable(IndicatorsListPath);
IndicatorsList = IndicatorsListTable.Variables; % the same as keep up table






k = strfind(CSVData{9},'csv_source_folder');
csv_source_folder = extractBetween(CSVData{9},k+length('csv_source_folder')+1,length(CSVData{9}));
csv_source_folder = csv_source_folder{1};




fprintf('The code is reading from pano data folder %s\n',input_pano_folder);
fprintf('The code is reading from cell segmentation data folder %s\n',input_cell_segmentation_folder);
fprintf('The code is reading from cell type data folder %s\n',input_celltype_folder);
fprintf('The code is reading key patient id file from %s\n',input_keyPatientID);
fprintf('The code is reading group patient ids file from %s\n',group_PatientID);
fprintf('The code will write its result to the folder %s\n',output_folder);

% Data = readKeyTable(input_keyPatientID);
groupPatientIDTable = readtable(group_PatientID);
nOfGroups = width(groupPatientIDTable);
groupPatientIDTableData = groupPatientIDTable.Variables;


inputKeyPatientIDTable = readtable(input_keyPatientID);

headerKeyTable = inputKeyPatientIDTable.Properties.VariableNames;
headerGroupPatientID = groupPatientIDTable.Properties.VariableNames;
%columnIndexFileName
%columnIndexImageID

columnIndexFileName = find(contains(headerKeyTable,'FileName'));
columnIndexImageID = find(contains(headerKeyTable,'UniqueImageID'));
inputKeyPatientIDData = inputKeyPatientIDTable.Variables;
fileNameColumn = inputKeyPatientIDData(:,columnIndexFileName);
imageIDColumn = inputKeyPatientIDData(:,columnIndexImageID);







% find all addresses here for all file names
middleString = [];
celltype_files =[];
for i = 1 : 100
    celltype_files = [celltype_files;rdir(fullfile(input_celltype_folder,middleString,'*.mat'))];
    middleString = fullfile(middleString,'*');
end



middleString = [];
regIndices_files =[];
for i = 1 : 100
    regIndices_files = [regIndices_files;rdir(fullfile(input_cell_segmentation_folder,middleString,'allRegionIndices.mat'))];
    middleString = fullfile(middleString,'*');
end


middleString = [];
pano_files =[];
for i = 1 : 100
    pano_files = [pano_files;rdir(fullfile(input_pano_folder,middleString,'*.txt'))];
    middleString = fullfile(middleString,'*');
end






for indexGroup = 1 : nOfGroups
    
    
    
    
    
    

    middleString = [];
    csv_files =[];
    for i = 1 : 100
        csv_files = [csv_files;rdir(fullfile(csv_source_folder,middleString,'*.csv'))];
        middleString = fullfile(middleString,'*');
    end



    
    
    
    
    
    groupFolderPath = fullfile(output_folder,headerGroupPatientID{indexGroup});
    mkdir(groupFolderPath);
    close all;
    cur_group_keys = groupPatientIDTableData(:,indexGroup);
    cur_group_keys = cur_group_keys(~cellfun(@isempty, cur_group_keys));
    fileNamesForThisGroup = cell(size(cur_group_keys));
    
    for i = 1 : length(cur_group_keys)
        cur_key = cur_group_keys{i};
        %rowIndex = find(contains(imageIDColumn,cur_key));
        exact_match_mask = strcmp(imageIDColumn, cur_key);
        rowIndex = find(exact_match_mask);

        
        
        if length(rowIndex)==1
            fileNamesForThisGroup{i} = fileNameColumn{rowIndex};
        else
            
            curCell = cell(length(rowIndex),1);
            for tt = 1 : length(rowIndex)
                curCell{tt} = fileNameColumn{rowIndex(tt)};
            end
            fileNamesForThisGroup{i} = curCell;
        end
        
    end
    
    
    nOfFiles = length(fileNamesForThisGroup);
    
    allHeatmapMatrices = cell(nOfFiles,1);
    % read channels
    for indexFile = 1  : nOfFiles
        
        
        
        
        %%% find pano (.txt) file addresses
        CombinedFlag = 0;
        if ~iscellstr(fileNamesForThisGroup{indexFile})
            for j = 1:length(pano_files)
                fullFilePath = pano_files(j).name;
                [folderPath,scanName,ext]=fileparts(fullFilePath);
                if strcmp(scanName,fileNamesForThisGroup{indexFile}) ==1
                    panoFilePath = fullFilePath;
                end
            end
        else
            CombinedFlag = 1;
            curCell = fileNamesForThisGroup{indexFile};
            scanName = '';
            panoFilePath = cell(length(curCell),1);
            for tt = 1 : length(curCell)
                queryStr = curCell{tt};
                for j = 1:length(pano_files)
                    fullFilePath = pano_files(j).name;
                    [folderPath,scanName,ext]=fileparts(fullFilePath);
                    if strcmp(scanName,queryStr) ==1
                        panoFilePath{tt} = fullFilePath;
                    end
                end
            end
        end
        
        %%% find cell type file addresses
        if ~iscellstr(fileNamesForThisGroup{indexFile})
            for j = 1:length(celltype_files)
                fullFilePath = celltype_files(j).name;
                [folderPath,scanName,ext]=fileparts(fullFilePath);
                if strcmp(scanName,fileNamesForThisGroup{indexFile}) ==1
                    cellTypeFilePath = fullFilePath;
                end
            end
            [~,finalScanName,~]=fileparts(cellTypeFilePath);
            
        else
            curCell = fileNamesForThisGroup{indexFile};
            finalScanName = '';
            cellTypeFilePath = cell(length(curCell),1);
            for tt = 1 : length(curCell)
                queryStr = curCell{tt};
                for j = 1:length(celltype_files)
                    fullFilePath = celltype_files(j).name;
                    [folderPath,scanName,ext]=fileparts(fullFilePath);
                    if strcmp(scanName,queryStr) ==1
                        cellTypeFilePath{tt} = fullFilePath;
                    end
                end
                [~,cur_scanName,~]=fileparts(cellTypeFilePath{tt});
                if tt > 1
                    if length(cur_scanName) < length(finalScanName)
                        finalScanName = cur_scanName;
                    end
                else
                    finalScanName = cur_scanName;
                end
            end
            
        end
        
        
        %%% find all region indices file addresses
        if ~iscellstr(fileNamesForThisGroup{indexFile})
            for j = 1:length(regIndices_files)
                fullFilePath = regIndices_files(j).name;
                [folderPath,~,ext]=fileparts(fullFilePath);
                endout=regexp(folderPath,filesep,'split');
                scanName = endout{end};
                if strcmp(scanName,fileNamesForThisGroup{indexFile}) ==1
                    allRegIndsPath = fullFilePath;
                end
            end
        else
            curCell = fileNamesForThisGroup{indexFile};
            scanName = '';
            allRegIndsPath = cell(length(curCell),1);
            for tt = 1 : length(curCell)
                queryStr = curCell{tt};
                for j = 1:length(regIndices_files)
                    fullFilePath = regIndices_files(j).name;
                    [folderPath,~,ext]=fileparts(fullFilePath);
                    endout=regexp(folderPath,filesep,'split');
                    scanName = endout{end};
                    if strcmp(scanName,queryStr) ==1
                        allRegIndsPath{tt} = fullFilePath;
                    end
                end
            end
        end
        
        %%% find csv (.csv) file addresses
        CombinedFlag = 0;
        csvFilePathSourceFound = 0;
        if ~iscellstr(fileNamesForThisGroup{indexFile})
            for j = 1:length(csv_files)
                fullFilePath = csv_files(j).name;
                [folderPath,scanName,ext]=fileparts(fullFilePath);
                if strcmp(scanName,fileNamesForThisGroup{indexFile}) ==1
                    csvFilePathSourceFound = 1;
                end
            end
        else
            CombinedFlag = 1;
            curCell = fileNamesForThisGroup{indexFile};
            scanName = '';
            csvFilePathSource = cell(length(curCell),1);
            for tt = 1 : length(curCell)
                queryStr = curCell{tt};
                for j = 1:length(csv_files)
                    fullFilePath = csv_files(j).name;
                    [folderPath,scanName,ext]=fileparts(fullFilePath);
                    if strcmp(scanName,queryStr) ==1
                        csvFilePathSourceFound = 1;
                    end
                end
            end
        end
        
        
        
        
        
        
        if csvFilePathSourceFound==0        
            heatMapMatrix = generate_heatmap_per_scan_filepath_with_config_outliers(CombinedFlag,panoFilePath,cellTypeFilePath,allRegIndsPath,IndicatorsList,allLabels);

            allHeatmapMatrices{indexFile} = heatMapMatrix;        
            csvFilePath = fullfile(csv_source_folder,strcat(finalScanName,'.csv'));
            mkdir(csv_source_folder)
            fileID = fopen(csvFilePath,'w');

            fprintf(fileID,'%s,',finalScanName);
            for ii = 1 : length(IndicatorsList)
                if ii~=length(IndicatorsList)
                    fprintf(fileID,'%s,',IndicatorsList{ii});
                else

                    fprintf(fileID,'%s',IndicatorsList{ii});
                end
            end
            fprintf(fileID,'\n');
            for ii = 1 : size(heatMapMatrix,1)
                fprintf(fileID,'%s,',allLabels{ii});
                for jj = 1 : size(heatMapMatrix,2)
                    if jj~=size(heatMapMatrix,2)
                        fprintf(fileID,'%f,',heatMapMatrix(ii,jj));
                    else
                        fprintf(fileID,'%f',heatMapMatrix(ii,jj));
                    end
                end
                fprintf(fileID,'\n');
            end
            fclose(fileID);
        else
            disp('reading from source');
            csvFilePath = fullfile(csv_source_folder,strcat(finalScanName,'.csv'));
            T = readtable(csvFilePath);
            T(:,1) = [];
            heatMapMatrix= T.Variables;
            allHeatmapMatrices{indexFile} = heatMapMatrix;
        end
        
        
    end
    
    overallHeatMapMatrix = zeros(length(allLabels),length(IndicatorsList));

    for indexFile = 1:nOfFiles
        overallHeatMapMatrix = overallHeatMapMatrix + allHeatmapMatrices{indexFile};
    end
    heatMapMatrixNormalize = overallHeatMapMatrix;
    for ind = 1 : size(overallHeatMapMatrix,2)
         heatMapMatrixNormalize(:,ind) = zscore(overallHeatMapMatrix(:,ind));
         %heatMapMatrixNormalize(:,ind) = (overallHeatMapMatrix(:,ind));
    end
    heatMapMatrixNormalize = fix(heatMapMatrixNormalize*1000)/1000;
    figure('units','normalized','outerposition',[0 0 1 1])
    heatmap(IndicatorsList,allLabels,heatMapMatrixNormalize)
    colormap redblue
    caxis([-5,5])
    pdfFilePath = fullfile(groupFolderPath,strcat(headerGroupPatientID{indexGroup},'.pdf'));
    export_fig(strcat(pdfFilePath));
    
    
    
    
    
    
    
    
    csvFilePath = fullfile(groupFolderPath,strcat(headerGroupPatientID{indexGroup},'.csv'));
    fileID = fopen(csvFilePath,'w');

    %%%%%%
    overallHeatMapMatrix = heatMapMatrixNormalize;
    fprintf(fileID,'%s,',headerGroupPatientID{indexGroup});
    for ii = 1 : length(IndicatorsList)
        if ii~=length(IndicatorsList)
            fprintf(fileID,'%s,',IndicatorsList{ii});
        else
            
            fprintf(fileID,'%s',IndicatorsList{ii});
        end
    end
    fprintf(fileID,'\n');
    
    
    
    for ii = 1 : size(overallHeatMapMatrix,1)
        fprintf(fileID,'%s,',allLabels{ii});
        for jj = 1 : size(overallHeatMapMatrix,2)
            if jj~=size(overallHeatMapMatrix,2)
                fprintf(fileID,'%f,',overallHeatMapMatrix(ii,jj));
            else
                fprintf(fileID,'%f',overallHeatMapMatrix(ii,jj));
            end
        end
        fprintf(fileID,'\n');
    end
    fclose(fileID);
    
    
end




