clear all;
clc;
addpath(genpath('Lib'))
config_file_path = 'config_templateFoxp3.csv';
[~,configFileName,~] = fileparts(config_file_path);

CSVData = importdata(fullfile('ConfigFiles',config_file_path));


lineCounter = 1;
k = strfind(CSVData{lineCounter},'input_pano_folder');
input_pano_folder = extractBetween(CSVData{lineCounter},k+length('input_pano_folder')+1,length(CSVData{1}));
input_pano_folder = input_pano_folder{lineCounter};


lineCounter = lineCounter +1;
k = strfind(CSVData{lineCounter},'input_cell_segmentation_folder');
input_cell_segmentation_folder = extractBetween(CSVData{lineCounter},k+length('input_cell_segmentation_folder')+1,length(CSVData{lineCounter}));
input_cell_segmentation_folder = input_cell_segmentation_folder{1};


lineCounter = lineCounter +1;
k = strfind(CSVData{lineCounter},'input_keyPatientID');
input_keyPatientID = extractBetween(CSVData{lineCounter},k+length('input_keyPatientID')+1,length(CSVData{lineCounter}));
input_keyPatientID = input_keyPatientID{1};

lineCounter = lineCounter +1;
k = strfind(CSVData{lineCounter},'group_PatientID');
group_PatientID = extractBetween(CSVData{lineCounter},k+length('group_PatientID')+1,length(CSVData{lineCounter}));
group_PatientID = group_PatientID{1};

% 
% lineCounter = lineCounter +1;
% k = strfind(CSVData{lineCounter},'output_folder');
% output_folder = extractBetween(CSVData{lineCounter},k+length('output_folder')+1,length(CSVData{lineCounter}));
% output_folder = output_folder{1};
% mkdir(output_folder);

lineCounter = lineCounter +1;
k = strfind(CSVData{lineCounter},'IndicatorsListPath');
IndicatorsListPath = extractBetween(CSVData{lineCounter},k+length('IndicatorsListPath')+1,length(CSVData{lineCounter}));
IndicatorsListPath = IndicatorsListPath{1};
IndicatorsListTable = readtable(IndicatorsListPath);
IndicatorsList = IndicatorsListTable.Indicators; % the same as keep up table
TypeList = IndicatorsListTable.Type;


lineCounter = lineCounter +1;
k = strfind(CSVData{lineCounter},'csv_source_folder');
csv_source_folder = extractBetween(CSVData{lineCounter},k+length('csv_source_folder')+1,length(CSVData{lineCounter}));
csv_source_folder = csv_source_folder{1};
mkdir(csv_source_folder);


lineCounter = lineCounter +1;
k = strfind(CSVData{lineCounter},'JMF_folderpath');
JMF_folderpath = extractBetween(CSVData{lineCounter},k+length('JMF_folderpath')+1,length(CSVData{lineCounter}));
JMF_folderpath = JMF_folderpath{1};



fprintf('The code is reading from pano data folder %s\n',input_pano_folder);
fprintf('The code is reading from cell segmentation data folder %s\n',input_cell_segmentation_folder);
fprintf('The code is reading key patient id file from %s\n',input_keyPatientID);
fprintf('The code is reading group patient ids file from %s\n',group_PatientID);
% fprintf('The code will write its result to the folder %s\n',output_folder);

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
regIndices_files =[];
for i = 1 : 10
    regIndices_files = [regIndices_files;rdir(fullfile(input_cell_segmentation_folder,middleString,'allRegionIndices.mat'))];
    middleString = fullfile(middleString,'*');
end


middleString = [];
pano_files =[];
for i = 1 : 10
    pano_files = [pano_files;rdir(fullfile(input_pano_folder,middleString,'*.txt'))];
    middleString = fullfile(middleString,'*');
end






for indexGroup = 1 %: nOfGroups
    
    
    
    
    
    

    middleString = [];
    csv_files =[];
    for i = 1 : 100
        csv_files = [csv_files;rdir(fullfile(csv_source_folder,middleString,'*.csv'))];
        middleString = fullfile(middleString,'*');
    end



    
    
    
    
    
%     groupFolderPath = fullfile(output_folder,headerGroupPatientID{indexGroup});
%     mkdir(groupFolderPath);
    close all;
    cur_group_keys = groupPatientIDTableData(:,indexGroup);
    cur_group_keys = cur_group_keys(~cellfun(@isempty, cur_group_keys));
    fileNamesForThisGroup = cell(size(cur_group_keys));
    
    for i = 1 : length(cur_group_keys)
        cur_key = cur_group_keys{i};
        rowIndex = find(contains(imageIDColumn,cur_key));
        
        
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
    

    parfor indexFile = 1  : nOfFiles
        
        
        
        
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
            [~,finalScanName,~]=fileparts(fileparts(allRegIndsPath));
        else
            curCell = fileNamesForThisGroup{indexFile};
            scanName = '';
            finalScanName = '';
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
                [~,cur_scanName,~]=fileparts(fileparts(allRegIndsPath{tt}));
                if tt > 1
                    if length(cur_scanName) < length(finalScanName)
                        finalScanName = cur_scanName;
                    end
                else
                    finalScanName = cur_scanName;
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

            if CombinedFlag==0


                JMF_filePath = extractBetween(panoFilePath,length(input_pano_folder)+1,length(panoFilePath)-4);
                JMF_filePath = JMF_filePath{1};
                JMF_filePath = fullfile(JMF_folderpath,JMF_filePath,'FoxP3_2.png');
            else

                for ikt = 1:length(panoFilePath)
                    JMF_filePath_cur = extractBetween(panoFilePath{ikt},length(input_pano_folder)+1,length(panoFilePath)-4);
                    JMF_filePath_cur = JMF_filePath_cur{1};
                    JMF_filePath{ikt} = fullfile(JMF_folderpath,JMF_filePath_cur,'FoxP3_2.png');
                    JMF_filePath
                    indexFile
                end
            end

            cellMeanIntensityMatrix = generate_cell_intensity_per_scan_filepath_with_config_mask_2(CombinedFlag,panoFilePath,allRegIndsPath,IndicatorsList,TypeList,JMF_folderpath,input_pano_folder);
            csvFilePath = fullfile(csv_source_folder,strcat(finalScanName,'.csv'));
            fileID = fopen(csvFilePath,'w');

            for ii = 1 : length(IndicatorsList)
                if ii~=length(IndicatorsList)
                    fprintf(fileID,'%s,',IndicatorsList{ii});
                else

                    fprintf(fileID,'%s',IndicatorsList{ii});
                end
            end
            fprintf(fileID,'\n');
            for ii = 1 : size(cellMeanIntensityMatrix,1)
                 fprintf(fileID,'%s,',allLabels{ii});
                for jj = 1 : size(cellMeanIntensityMatrix,2)
                    if jj~=size(cellMeanIntensityMatrix,2)
                        fprintf(fileID,'%f,',cellMeanIntensityMatrix(ii,jj));
                    else
                        fprintf(fileID,'%f',cellMeanIntensityMatrix(ii,jj));
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
            cellMeanIntensityMatrix= T.Variables;
        end
        
        
    end
    allLabels = TypeList;
    
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




