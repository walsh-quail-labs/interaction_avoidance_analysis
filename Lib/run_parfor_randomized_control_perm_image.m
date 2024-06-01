function run_parfor_randomized_control_perm_image(fileNamesForThisGroup,indexFile,pano_files,celltype_files,regIndices_files,dist,NofPerm,output_folder,allLabels)
        
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
        if CombinedFlag == 0
            Neighborhood_Individual_Image_Ourversion(nuc_path,celltype_path,allLabels,dist,NofPerm,matFilePath);
%             heatMapMatrix = Count_pairwise_rel_perImage_permute(nuc_path,celltype_path,allLabels,dist,NofPerm);
%             save(matFilePath,'heatMapMatrix');                
        else
            heatMapMatrix = 0;
            for fileITR = 1 : length(nuc_path)
                HM = Count_pairwise_rel_perImage_permute(nuc_path{fileITR},celltype_path{fileITR},allLabels,dist,NofPerm);
                heatMapMatrix = heatMapMatrix + HM;
            end
            save(matFilePath,'heatMapMatrix');
        end