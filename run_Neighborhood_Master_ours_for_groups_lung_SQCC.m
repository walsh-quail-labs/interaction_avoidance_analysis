delete(gcp('nocreate'));
clear all;
clc;
addpath(genpath('Lib'))

%!! check alllabels and celltype location
%if they are included in the celltype files or not 
%if not comment 62 and add label file 63-64
config_file_path = 'template_nada_final.csv';

[~,configFileName,~] = fileparts(config_file_path);
CSVData = importdata(fullfile('ConfigFiles/nada_final',config_file_path));


%-- loop through all the permutation folders

for perm_folder_idx =1:length(filelist)
    
    output_folder        = fullfile(perm_folder,filelist(perm_folder_idx).name);
    custom_gatesfolder   = fullfile(save_folder,filelist(perm_folder_idx).name);
    mkdir(custom_gatesfolder);
    
    t = strsplit(output_folder,"-");
    t = char(t{length(t)});
    
    permutations = str2num(t(1:length(t)-1))*1000;
    disp(filelist(perm_folder_idx).name)
    disp(permutations)
    
    pVal_sig = 0.01;
    Special_clusters_name = '1';
    Extra_information = '';
    pixelexpansion = 6;
    cut_off_percent = 0.3;
   
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
%     allLabels = allLabelsTable.Variables;
    allLabels = load('ConfigFiles/nada_final/allLabels.mat');
    allLabels = allLabels.allLabels;
    
    % k = strfind(CSVData{7},'output_folder');
    % output_folder = extractBetween(CSVData{7},k+length('output_folder')+1,length(CSVData{7}));
    % output_folder = output_folder{1};
    % mkdir(output_folder);
    
    fprintf('The code is reading from pano data folder %s\n',input_pano_folder);
    fprintf('The code is reading from cell segmentation data folder %s\n',input_cell_segmentation_folder);
    fprintf('The code is reading from cell type data folder %s\n',input_celltype_folder);
    fprintf('The code is reading key patient id file from %s\n',input_keyPatientID);
    fprintf('The code is reading group patient ids file from %s\n',group_PatientID);
    
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
    
    
    
    
    
    %merge them to one sheet
    %first col group 1 and second col for group 2
    %compare the two groups 
    %heatmaps plotted separately
    for indexGroup = 1: nOfGroups
        close all;
        custom_gatesfolder_for_group = fullfile(custom_gatesfolder,headerGroupPatientID{indexGroup});
        mkdir(custom_gatesfolder_for_group);
        GroupName = headerGroupPatientID{indexGroup};
        
        
        
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
                 
    
        
        allHeatmapMatrices = cell(nOfFiles,1);
        % read channels
        list_files = cell(nOfFiles,2);
        for indexFile = 1 : nOfFiles
            
            
            
            
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
            
            
            nuc_path = allRegIndsPath;
            celltype_path = cellTypeFilePath;
            
            %matFilePath = fullfile(csv_source_folder,strcat(finalScanName,'.mat'));
            matFilePath = fullfile(output_folder,strcat(finalScanName,'_perm.mat'));
            
            
            
            
            
    %         index_row = str2double(selected_rows{indexFile});
            list_files{indexFile,1} = finalScanName;
            list_files{indexFile,2}= cur_group_keys{indexFile};
            
                
            
            
            input_data = importdata(matFilePath);
    
    
    
    
            pValue_higher = input_data.pValue_higher;
            pValue_lower = input_data.pValue_lower;
            combos_all = input_data.combos_all;
            real_data_mean = input_data.real_data_mean;
    
    
            real_data_mean = real_data_mean';
            pValue_higher = pValue_higher';
            pValue_lower = pValue_lower';
            combos_all = double(combos_all);
    
            Higher_logic = pValue_higher<pVal_sig;
            Lower_logic = pValue_lower<pVal_sig;
    
            % Extract higher and lower logic after pVal_sig for each image
            parfor_gates_high(indexFile) = {[combos_all,Higher_logic]};
            parfor_gates_low(indexFile) = {[combos_all,Lower_logic]};
        end
        
        [Matrix_Delta,Matrix_low,Unique_all,Unique_low_all,Matrix_high,Unique_high_all,pheno_name,Matrix_Delta_cut_noNaN,Unique_all_string_names_cut]...
        = Heatmap_individual_images_ours(parfor_gates_high,parfor_gates_low,permutations,custom_gatesfolder_for_group,Special_clusters_name,Extra_information,pVal_sig);
    
    
        % Generate an assymmetric heatmap
        Delta_allvsall = Asymmetric_heatmap_ours_pdf(parfor_gates_high,parfor_gates_low, Matrix_high,...
            Matrix_low,Unique_high_all,Unique_low_all,pheno_name,pixelexpansion,permutations,custom_gatesfolder_for_group,Extra_information,pVal_sig,cut_off_percent,allLabels,GroupName);
        
        
    %     writetable(array2table(Delta_allvsall),strcat(custom_gatesfolder,'/',GroupName,'/Delta_allvsall.xlsx'));
        
        csvFilePath = fullfile(custom_gatesfolder_for_group,strcat(GroupName,'.csv'));
        fileID = fopen(csvFilePath,'w');
    
    %     
        fprintf(fileID,'%s,','Delta');
        for ii = 1 : length(allLabels)
            if ii~=length(allLabels)
                fprintf(fileID,'%s,',allLabels{ii});
            else
    
                fprintf(fileID,'%s',allLabels{ii});
            end
        end
        fprintf(fileID,'\n');
    
        for ii = 1 : size(Delta_allvsall,1)
            fprintf(fileID,'%s,',allLabels{ii});
            for jj = 1 : size(Delta_allvsall,2)
                if jj~=size(Delta_allvsall,2)
                    fprintf(fileID,'%f,',Delta_allvsall(ii,jj));
                else
                    fprintf(fileID,'%f',Delta_allvsall(ii,jj));
                end
            end
            fprintf(fileID,'\n');
        end
        %%%%%%%%%%
        for ii = 1 : length(allLabels)+1
            fprintf(fileID,'\t');
        end
        fprintf(fileID,'\n');
    
        for ii = 1 : length(allLabels)+1
            fprintf(fileID,'\t');
        end
        fprintf(fileID,'\n');
    
    
    
    end
    close all force 
end 









