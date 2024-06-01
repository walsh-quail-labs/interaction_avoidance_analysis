delete(gcp('nocreate'));
clear all;
clc;
addpath(genpath('Lib'))
config_file_path = 'template_nada_final.csv';
[~,configFileName,~] = fileparts(config_file_path);
CSVData = importdata(fullfile('ConfigFiles/nada_final',config_file_path));

for perm_folder_idx = 1:length(filelist)
output_folder        = fullfile(perm_folder,filelist(perm_folder_idx).name);
custom_gatesfolder   = fullfile(save_folder,filelist(perm_folder_idx).name);
mkdir(custom_gatesfolder);

t = strsplit(output_folder,"-");
% t = strsplit(output_folder,"/");
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
% allLabels = allLabelsTable.Variables;
allLabels = load('ConfigFiles/nada_final/allLabels.mat');
allLabels = allLabels.allLabels;

% k = strfind(CSVData{7},'output_folder');
% output_folder = extractBetween(CSVData{7},k+length('output_folder')+1,length(CSVData{7}));
% output_folder = output_folder{1};
% mkdir(output_folder);

% Data = readKeyTable(input_keyPatientID);



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







% find Label addresses here for all file names
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




allTables_interaction = cell(nOfGroups,1);
allTables_avoidance = cell(nOfGroups,1);

allTables_interaction_save = cell(nOfGroups,1);
allTables_avoidance_save = cell(nOfGroups,1);
for indexGroup = 1: nOfGroups
    close all;
    TI = table;
    TA = table;
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
        %Before
        %         t  = pValue_higher*1001;
        %         t = t-1;
        %
        %         T.(finalScanName) = t;
        
        %Interaction
        t1  = pValue_higher*(permutations+1);
        t1 = t1-1; 
        t2  = pValue_lower*(permutations+1);
        t2 = t2-1;

        numEqs = (t1+t2)-permutations;

        t_interactin = t2-numEqs;
        %Avoidance
        t1  = pValue_lower*(permutations+1);
        t1 = t1-1;
        t2  = pValue_higher*(permutations+1);
        t2 = t2-1;

        numEqs = (t1+t2)-permutations;

        t_avoidance = t2-numEqs;


        %%%%%%%%%%%%%%%%
        TI.(finalScanName) = t_interactin;
        TA.(finalScanName) = t_avoidance;

        %         parfor_gates_high(indexFile) = {[combos_all,Higher_logic]};
        %         parfor_gates_low(indexFile) = {[combos_all,Lower_logic]};

        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


        %     CurMatInt = reshape_values_permutations(t_interactin,combos_all);
        %
        %     TT = table();
        %     TT.Interaction_Avoidance = allLabels;
        %     for iii = 1 : length(allLabels)
        %         column_name = allLabels{iii};
        %         TT.(column_name) = CurMatInt(:,iii);
        %     end

    end

    allTables_interaction{indexGroup} = TI;
    allTables_avoidance{indexGroup} = TA;
    


%     colNames = {'Interaction'};
%     for i=1:length(combos_all)
%         colNames = [colNames,sprintf('%s-%s',allLabels{combos_all(i,1)},allLabels{combos_all(i,2)})];
%     end 
% 
%     TI_transpose = rows2vars(TI);
%     TI_transpose.Properties.VariableNames = colNames;
%     part1 = TI_transpose(:, 1);
%     part2 = TI_transpose(:, 2:end);
%     group_column = repmat({GroupName}, size(TI_transpose, 1), 1);
%     TI_updated = [part1, table(group_column, 'VariableNames', {'Group'}), part2];
%     
%     colNames = {'Avoidance'};
%     for i=1:length(combos_all)
%         colNames = [colNames,sprintf('%s-%s',allLabels{combos_all(i,1)},allLabels{combos_all(i,2)})];
%     end 
%     TA_transpose = rows2vars(TA);
%     TA_transpose.Properties.VariableNames = colNames;
%     part1 = TA_transpose(:, 1);
%     part2 = TA_transpose(:, 2:end);
%     group_column = repmat({GroupName}, size(TA_transpose, 1), 1);
%     TA_updated = [part1, table(group_column, 'VariableNames', {'Group'}), part2];
%     
%     allTables_interaction_save{indexGroup} = TI_updated;
%     allTables_avoidance_save{indexGroup} = TA_updated;

%     allTables_interaction_pvalues{indexGroup} = TI;
%     allTables_avoidance_p_values{indexGroup} = TA;
end

% savedir = fullfile(custom_gatesfolder,"perfile");
% if exist(savedir, 'dir')
%     rmdir(savedir,'s');
% end 
% mkdir(savedir);
% 
% interaction_table = [allTables_interaction_save{1};allTables_interaction_save{2}];
% writetable(interaction_table,fullfile(savedir,'Interaction_per_file.csv'));
% 
% avoidance_table = [allTables_avoidance_save{1};allTables_avoidance_save{2}];
% writetable(avoidance_table,fullfile(savedir,'Avoidance_per_file.csv'));


T_interaction1 = allTables_interaction{1};
T_interaction2 = allTables_interaction{2};

T_avoidance1 = allTables_avoidance{1};
T_avoidance2 = allTables_avoidance{2};

[pvaluesI,stats_values] = compute_pvalues_between_two_tables(T_interaction1,T_interaction2);
[muvaluesI,md1valuesI,md2valuesI] = compute_meanstd_between_two_tables(T_interaction1,T_interaction2);
%%%final results
pvalue_final_interaction = reshape_values_permutations(pvaluesI,combos_all);
muvalues_final_interaction = reshape_values_permutations(abs(muvaluesI),combos_all);
md1values_final_interaction = reshape_values_permutations(md1valuesI,combos_all);
md2values_final_interaction = reshape_values_permutations(md2valuesI,combos_all);


%%%%%%%
csvFilePath = fullfile(custom_gatesfolder,'Interaction.csv');
fileID = fopen(csvFilePath,'w');


fprintf(fileID,'%s,','PValues');
for ii = 1 : length(allLabels)
    if ii~=length(allLabels)
        fprintf(fileID,'%s,',allLabels{ii});
    else

        fprintf(fileID,'%s',allLabels{ii});
    end
end
fprintf(fileID,'\n');

for ii = 1 : size(pvalue_final_interaction,1)
    fprintf(fileID,'%s,',allLabels{ii});
    for jj = 1 : size(pvalue_final_interaction,2)
        if jj~=size(pvalue_final_interaction,2)
            fprintf(fileID,'%f,',pvalue_final_interaction(ii,jj));
        else
            fprintf(fileID,'%f',pvalue_final_interaction(ii,jj));
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


fprintf(fileID,'%s,','diffrence between means');
for ii = 1 : length(allLabels)
    if ii~=length(allLabels)
        fprintf(fileID,'%s,',allLabels{ii});
    else

        fprintf(fileID,'%s',allLabels{ii});
    end
end
fprintf(fileID,'\n');
for ii = 1 : size(muvalues_final_interaction,1)
    fprintf(fileID,'%s,',allLabels{ii});
    for jj = 1 : size(muvalues_final_interaction,2)
        if jj~=size(muvalues_final_interaction,2)
            fprintf(fileID,'%f,',muvalues_final_interaction(ii,jj));
        else
            fprintf(fileID,'%f',muvalues_final_interaction(ii,jj));
        end
    end
    fprintf(fileID,'\n');
end
fprintf(fileID,'\n');
%%%%%%%%%%%%%%%%%%%%
for ii = 1 : length(allLabels)+1
    fprintf(fileID,'\t');
end
fprintf(fileID,'\n');

for ii = 1 : length(allLabels)+1
    fprintf(fileID,'\t');
end
fprintf(fileID,'\n');


fprintf(fileID,'%s,',headerGroupPatientID{1});
for ii = 1 : length(allLabels)
    if ii~=length(allLabels)
        fprintf(fileID,'%s,',allLabels{ii});
    else

        fprintf(fileID,'%s',allLabels{ii});
    end
end
fprintf(fileID,'\n');
for ii = 1 : size(md1values_final_interaction,1)
    fprintf(fileID,'%s,',allLabels{ii});
    for jj = 1 : size(md1values_final_interaction,2)
        if jj~=size(md1values_final_interaction,2)
            fprintf(fileID,'%f,',md1values_final_interaction(ii,jj));
        else
            fprintf(fileID,'%f',md1values_final_interaction(ii,jj));
        end
    end
    fprintf(fileID,'\n');
end
fprintf(fileID,'\n');
%%%%%%%%%%%%%%%%%%%%
for ii = 1 : length(allLabels)+1
    fprintf(fileID,'\t');
end
fprintf(fileID,'\n');

for ii = 1 : length(allLabels)+1
    fprintf(fileID,'\t');
end
fprintf(fileID,'\n');


fprintf(fileID,'%s,',headerGroupPatientID{2});
for ii = 1 : length(allLabels)
    if ii~=length(allLabels)
        fprintf(fileID,'%s,',allLabels{ii});
    else

        fprintf(fileID,'%s',allLabels{ii});
    end
end
fprintf(fileID,'\n');
for ii = 1 : size(md2values_final_interaction,1)
    fprintf(fileID,'%s,',allLabels{ii});
    for jj = 1 : size(md2values_final_interaction,2)
        if jj~=size(md2values_final_interaction,2)
            fprintf(fileID,'%f,',md2values_final_interaction(ii,jj));
        else
            fprintf(fileID,'%f',md2values_final_interaction(ii,jj));
        end
    end
    fprintf(fileID,'\n');
end
fprintf(fileID,'\n');

fclose(fileID);

%%%%%%%%%



[pvaluesA,stats_values] = compute_pvalues_between_two_tables(T_avoidance1,T_avoidance2);
[muvaluesA,md1valuesA,md2valuesA] = compute_meanstd_between_two_tables(T_avoidance1,T_avoidance2);
%%%final results
pvalue_final_avoidance = reshape_values_permutations(pvaluesA,combos_all);
muvalues_final_avoidance = reshape_values_permutations(abs(muvaluesA),combos_all);
md1values_final_avoidance = reshape_values_permutations(md1valuesA,combos_all);
md2values_final_avoidance = reshape_values_permutations(md2valuesA,combos_all);




csvFilePath = fullfile(custom_gatesfolder,'Avoidance.csv');
fileID = fopen(csvFilePath,'w');


fprintf(fileID,'%s,','PValues');
for ii = 1 : length(allLabels)
    if ii~=length(allLabels)
        fprintf(fileID,'%s,',allLabels{ii});
    else

        fprintf(fileID,'%s',allLabels{ii});
    end
end
fprintf(fileID,'\n');

for ii = 1 : size(pvalue_final_avoidance,1)
    fprintf(fileID,'%s,',allLabels{ii});
    for jj = 1 : size(pvalue_final_avoidance,2)
        if jj~=size(pvalue_final_avoidance,2)
            fprintf(fileID,'%f,',pvalue_final_avoidance(ii,jj));
        else
            fprintf(fileID,'%f',pvalue_final_avoidance(ii,jj));
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


fprintf(fileID,'%s,','diffrence between means');
for ii = 1 : length(allLabels)
    if ii~=length(allLabels)
        fprintf(fileID,'%s,',allLabels{ii});
    else

        fprintf(fileID,'%s',allLabels{ii});
    end
end
fprintf(fileID,'\n');
for ii = 1 : size(muvalues_final_avoidance,1)
    fprintf(fileID,'%s,',allLabels{ii});
    for jj = 1 : size(muvalues_final_avoidance,2)
        if jj~=size(muvalues_final_avoidance,2)
            fprintf(fileID,'%f,',muvalues_final_avoidance(ii,jj));
        else
            fprintf(fileID,'%f',muvalues_final_avoidance(ii,jj));
        end
    end
    fprintf(fileID,'\n');
end
fprintf(fileID,'\n');
%%%%%%%%%%%%%%%%%%%%
for ii = 1 : length(allLabels)+1
    fprintf(fileID,'\t');
end
fprintf(fileID,'\n');

for ii = 1 : length(allLabels)+1
    fprintf(fileID,'\t');
end
fprintf(fileID,'\n');


fprintf(fileID,'%s,',headerGroupPatientID{1});
for ii = 1 : length(allLabels)
    if ii~=length(allLabels)
        fprintf(fileID,'%s,',allLabels{ii});
    else

        fprintf(fileID,'%s',allLabels{ii});
    end
end
fprintf(fileID,'\n');
for ii = 1 : size(md1values_final_avoidance,1)
    fprintf(fileID,'%s,',allLabels{ii});
    for jj = 1 : size(md1values_final_avoidance,2)
        if jj~=size(md1values_final_avoidance,2)
            fprintf(fileID,'%f,',md1values_final_avoidance(ii,jj));
        else
            fprintf(fileID,'%f',md1values_final_avoidance(ii,jj));
        end
    end
    fprintf(fileID,'\n');
end
fprintf(fileID,'\n');
%%%%%%%%%%%%%%%%%%%%
for ii = 1 : length(allLabels)+1
    fprintf(fileID,'\t');
end
fprintf(fileID,'\n');

for ii = 1 : length(allLabels)+1
    fprintf(fileID,'\t');
end
fprintf(fileID,'\n');


fprintf(fileID,'%s,',headerGroupPatientID{2});
for ii = 1 : length(allLabels)
    if ii~=length(allLabels)
        fprintf(fileID,'%s,',allLabels{ii});
    else

        fprintf(fileID,'%s',allLabels{ii});
    end
end
fprintf(fileID,'\n');
for ii = 1 : size(md2values_final_avoidance,1)
    fprintf(fileID,'%s,',allLabels{ii});
    for jj = 1 : size(md2values_final_avoidance,2)
        if jj~=size(md2values_final_avoidance,2)
            fprintf(fileID,'%f,',md2values_final_avoidance(ii,jj));
        else
            fprintf(fileID,'%f',md2values_final_avoidance(ii,jj));
        end
    end
    fprintf(fileID,'\n');
end
fprintf(fileID,'\n');

fclose(fileID);

end 