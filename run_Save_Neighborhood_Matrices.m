delete(gcp('nocreate'));
clear all;
clc;
addpath(genpath('Lib'))
output_folder = 'inter_nada_final';
mkdir(output_folder); 
config_file_path = 'template_nada_final.csv';


mkdir(output_folder);
[~,configFileName,~] = fileparts(config_file_path);
CSVData = importdata(fullfile('ConfigFiles/nada_final',config_file_path));

dist = 6;

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

% find all addresses here for all file names
middleString = [];
celltype_files =[];
for i = 1 : 10
    celltype_files = [celltype_files;rdir(fullfile(input_celltype_folder,middleString,'*.mat'))];
    middleString = fullfile(middleString,'*');
end

fprintf('The code is reading from pano data folder %s\n',input_pano_folder);
fprintf('The code is reading from cell segmentation data folder %s\n',input_cell_segmentation_folder);
fprintf('The code is reading from cell type data folder %s\n',input_celltype_folder);
fprintf('The code is reading key patient id file from %s\n',input_keyPatientID);
fprintf('The code is reading group patient ids file from %s\n',group_PatientID);

groupPatientIDTable = readtable(group_PatientID);
nOfGroups = width(groupPatientIDTable);
groupPatientIDTableData = groupPatientIDTable.Variables;

inputKeyPatientIDTable = readtable(input_keyPatientID);

headerKeyTable = inputKeyPatientIDTable.Properties.VariableNames;
headerGroupPatientID = groupPatientIDTable.Properties.VariableNames;

columnIndexFileName = find(contains(headerKeyTable,'FileName'));
columnIndexImageID = find(contains(headerKeyTable,'UniqueImageID'));
inputKeyPatientIDData = inputKeyPatientIDTable.Variables;
fileNameColumn = inputKeyPatientIDData(:,columnIndexFileName);
imageIDColumn = inputKeyPatientIDData(:,columnIndexImageID);

% find all addresses here for all file names
middleString = [];
celltype_files =[];
for i = 1 : 10
    celltype_files = [celltype_files;rdir(fullfile(input_celltype_folder,middleString,'*.mat'))];
    middleString = fullfile(middleString,'*');
end

middleString = [];
regIndices_files =[];
for i = 1 : 10
    regIndices_files = [regIndices_files;rdir(fullfile(input_cell_segmentation_folder,middleString,'nuclei_multiscale.mat'))];
    middleString = fullfile(middleString,'*');
end

middleString = [];
pano_files =[];
for i = 1 : 10
    pano_files = [pano_files;rdir(fullfile(input_pano_folder,middleString,'*.txt'))];
    middleString = fullfile(middleString,'*');
end


for indexGroup = 1: nOfGroups
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
    allHeatmapMatrices = cell(nOfFiles,1);
    list_files = cell(nOfFiles,2);

    for indexFile = 1 : nOfFiles
        disp(fileNamesForThisGroup{indexFile})

        if ~iscellstr(fileNamesForThisGroup{indexFile})
            for j = 1:length(celltype_files)
                fullFilePath = celltype_files(j).name;
                [folderPath,scanName,ext]=fileparts(fullFilePath);
                if strcmp(scanName,fileNamesForThisGroup{indexFile}) ==1
                    disp(scanName)
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
        
        if ~iscell(allRegIndsPath)
            nuc_path = allRegIndsPath;
            celltype_path = cellTypeFilePath;

            [folder_path,finalScanName,~] = fileparts(celltype_path);
            matFilePath = fullfile(output_folder,strcat(finalScanName,'.mat'));

            tmp = split(nuc_path, '/');
            tmp2 = join(tmp(1:length(tmp)-1),'/');
            regInds_path = strcat(tmp2{1}, '/allRegionindices.mat');
         
            
            if ~exist(matFilePath,'file')
                  save_nbhd_matrix(nuc_path,regInds_path,celltype_path,dist,matFilePath);

            else
                fprintf('%s exists\n',matFilePath);
            end
        else
            for iCellITR = 1 : length(allRegIndsPath)
                nuc_path = allRegIndsPath{iCellITR};
                celltype_path = cellTypeFilePath{iCellITR};
                
                
                
                [folder_path,finalScanName,~] = fileparts(celltype_path);
                matFilePath = fullfile(output_folder,strcat(finalScanName,'.mat'));
                if length(input_celltype_folder)+2 < length(folder_path)
                    folder_structure = extractBetween(folder_path,length(input_celltype_folder)+2,length(folder_path));
                    regInds_path = fullfile(input_cell_segmentation_folder,folder_structure,finalScanName,'allRegionIndices.mat');
                    regInds_path = regInds_path{1};
                else
                    folder_structure = '';
                    regInds_path = fullfile(input_cell_segmentation_folder,folder_structure,finalScanName,'allRegionIndices.mat');
                end
                if ~exist(matFilePath,'file')
                    save_nbhd_matrix(nuc_path,regInds_path,celltype_path,dist,matFilePath);
               
                else
                    fprintf('%s exists\n',matFilePath);
                end
            end
        end



        
  
    end
    
  
end







