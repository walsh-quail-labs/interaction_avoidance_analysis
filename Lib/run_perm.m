clear all;
clc;
addpath(genpath('Lib'))
config_file_path = 'config_template3.csv';
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






for indexGroup = 1 : nOfGroups
    
    
    
    
    groupFolderPath = fullfile(output_folder,headerGroupPatientID{indexGroup});
    mkdir(groupFolderPath);
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
    % read channels
    parfor indexFile = 1  : nOfFiles
        
        
                dist = 20;
        NofPerm = 1000;
        run_parfor_perm_image(fileNamesForThisGroup,indexFile,pano_files,celltype_files,regIndices_files,dist,NofPerm,output_folder,allLabels);
        
        
        
        
    end
    
    
    
end




