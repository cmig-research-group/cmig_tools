function info = FEMA_save(outputType, dirOutput, varargin)
% Function to save different output from FEMA in different format

% Start timer
tSaveOverall = tic;
logging(FEMA_info)

% Parse input arguments
p = inputParser;
p.KeepUnmatched = true;
allowEmpty = @(f) @(x) isempty(x) || f(x);
addRequired(p, 'outputType', @(x) iscell(x) || ischar(x));
addRequired(p, 'dirOutput', @ischar);
addParameter(p, 'saveDesignMatrix', true, @islogical);
addParameter(p, 'outPrefix', [], allowEmpty(@(x) ischar(x)));
addParameter(p, 'matType', 'uncompressed', @(x) ismember(lower(x), {'compressed', 'uncompressed'}));
parse(p, outputType, dirOutput, varargin{:});

unMatched = p.Unmatched;
if isfield(unMatched, 'nonDEAP')
    writeJSON = ~logical(unMatched.nonDEAP);
else
    writeJSON = false;
end
unMatched_flds = fieldnames(unMatched);
saveDesignMatrix = p.Results.saveDesignMatrix;
outPrefix = p.Results.outPrefix;
matType = p.Results.matType;
info = unMatched.info;

% Append mat type to info
FEMA_save.matType = matType;

if isfield(unMatched, 'RandomEffects')
    RandomEffects = unMatched.RandomEffects;

    % Clean up mapping names
    RandomEffects_name = RandomEffects;

    RandomEffects_name(ismember(RandomEffects_name, 'F')) = {'Family effect'};
    RandomEffects_name(ismember(RandomEffects_name, 'S')) = {'Subject effect'};
    RandomEffects_name(ismember(RandomEffects_name, 'E')) = {'Error'};
    RandomEffects_name(ismember(RandomEffects_name, 'A')) = {'Additive genetic effect'};
    RandomEffects_name(ismember(RandomEffects_name, 'D')) = {'Dominant genetic effect'};
    RandomEffects_name(ismember(RandomEffects_name, 'M')) = {'Maternal effect'};
    RandomEffects_name(ismember(RandomEffects_name, 'P')) = {'Paternal effect'};
    RandomEffects_name(ismember(RandomEffects_name, 'H')) = {'Home effect'};
    RandomEffects_name(ismember(RandomEffects_name, 'T')) = {'Twin effect'};
end

% Ensure outputType is specified as lower case
outputType = lower(outputType);

if ischar(outputType)
    outputType = cellstr(outputType);
end

if ~exist(dirOutput, 'dir')
    mkdir(dirOutput);
end

nii_list = {'nii', 'nii.gz', 'nifti', 'voxel', 'voxelwise'}; 
gii_list = {'gii', 'gifti', 'vertex', 'vertexwise'};
splitLR = true; 
% M_atl     = [0 -1 0 101; 0 0 1 -131; -1 0 0 101; 0 0 0 1];
M_atl_sub = [0 -2 0 102; 0 0 2 -132; -2 0 0 102; 0 0 0 1];

if any(ismember(outputType, nii_list))
    if isfield(unMatched, 'mask')
        mask = unMatched.mask;
    else
        error('mask is required if outputType is specified as nifti');
    end
end

% Check for unstructured covariance parameters
if isfield(unMatched.unstructParams, 'eidOrd')
    eidOrd_all = unMatched.unstructParams.eidOrd;
end

%% output summary json
if any(ismember(outputType, 'summary'))
    info_tmp = info;
    info_flds = fieldnames(info_tmp);
    for ii=1:length(info_flds)
        if isfield(info_tmp.(info_flds{ii}), 'missingness')
            idx = ~cellfun(@isempty, regexpi(fieldnames(info_tmp.(info_flds{ii}).missingness), '^id_')); 
            rm_fld = fieldnames(info_tmp.(info_flds{ii}).missingness);
            info_tmp.(info_flds{ii}).missingness = rmfield(info_tmp.(info_flds{ii}).missingness, rm_fld(idx));
        end
    end
    if ~isempty(info_tmp.FEMA_makeDesign)
        if isfield(info_tmp.FEMA_makeDesign.settings, 'splines')
            info_tmp.FEMA_makeDesign.settings = rmfield(info_tmp.FEMA_makeDesign.settings.splines, 'basisSubset');
        end 
    end
    if writeJSON
        info_json = save_jsonencode(info_tmp, PrettyPrint=true);
        fid = fopen(fullfile(dirOutput, 'FEMA_summary.json'), 'w');
        fwrite(fid, info_json);
        fclose(fid);
    end
end

%% Output mat files
if any(ismember(outputType, {'mat', 'corrmat', 'external'}))
    % All possible things to save for estimates
    tSaveMat = tic; 
    toSave_estimates = {'beta_hat', 'beta_se', 'zmat', 'logpmat', 'sig2tvec', 'sig2mat', ...
                        'beta_hat_perm', 'beta_se_perm', 'zmat_perm', 'sig2tvec_perm',   ...
                        'sig2mat_perm', 'Wald', 'logp_Wald', 'logLikvec',      ...
                        'logLikvec_perm', 'Hessmat', 'coeffCovar', 'binvec_save',        ...
                        'nvec_bins', 'tvec_bins', 'residuals_GLS', 'unstructParams',     ...
                        'contrasts', 'hypValues', 'info', 'colnames_model'};

    toSave_design    = {'designMatrix', 'colnames_model', 'X', 'fid', 'iid', 'eid', 'agevec', ...
                        'FamilyStruct', 'MotherID', 'FatherID', 'HomeID', 'PregID'};
    
    % First: save main (permuted and non-permuted) statistics
    if ~isempty(outPrefix)
        outName = [outPrefix, '_estimates.mat'];
    else 
        outName = 'FEMA_estimates.mat';
    end
    saveName = fullfile(dirOutput, outName);
    toDelete = unMatched_flds(~ismember(unMatched_flds, toSave_estimates));
    toSave_struct = rmfield(unMatched, toDelete);
    checkSize = whos('toSave_struct');
    if checkSize.bytes > 2^31
        if strcmpi(matType, 'compressed')
            save(saveName, '-struct', 'toSave_struct', '-v7.3');
        else
            save(saveName, '-struct', 'toSave_struct', '-v7.3', '-nocompression');
        end
    else
        save(saveName, '-struct', 'toSave_struct');
    end
    fname_mat = saveName; 
    FEMA_save.timing.tSaveMat = toc(tSaveMat);
    
    % Second: save filtered design matrix and related variables
    if saveDesignMatrix
        tSaveDesign = tic;
        if ~isempty(outPrefix)
            outName = [outPrefix, '_designMatrix.mat'];
        else
            outName = 'FEMA_designMatrix.mat';
        end
        saveName = fullfile(dirOutput, outName);
        toDelete = unMatched_flds(~ismember(unMatched_flds, toSave_design));
        toSave_struct = rmfield(unMatched, toDelete);
        checkSize = whos('toSave_struct');
        if checkSize.bytes > 2^31
            if strcmpi(matType, 'compressed')
                save(saveName, '-struct', 'toSave_struct', '-v7.3');
            else
                save(saveName, '-struct', 'toSave_struct', '-v7.3', '-nocompression');
            end
        else 
            save(saveName, '-struct', 'toSave_struct'); 
        end
        FEMA_save.timing.tSaveDesign = toc(tSaveDesign);
    end
end

%% Save id list
if any(ismember(outputType, {'ids'}))
    tSaveIDs = tic;
    % If designMatrix table variable exists, rename the columns to
    % 'participant_id', and 'session_id', and then save it as a parquet file
    if ismember('designMatrix', unMatched_flds)
        if istable(unMatched.designMatrix)
            % Output name
            if ~isempty(outPrefix)
                outName = [outPrefix, '_designMatrix.parquet'];
            else
                outName = 'FEMA_designMatrix.parquet';
            end
            saveName = fullfile(dirOutput, outName);

            % Rename the design matrix columns before saving as a parquet
            unMatched.designMatrix.Properties.VariableNames{strcmpi(unMatched.designMatrix.Properties.VariableNames, 'fid')} = 'family_id';
            unMatched.designMatrix.Properties.VariableNames{strcmpi(unMatched.designMatrix.Properties.VariableNames, 'iid')} = 'participant_id';
            unMatched.designMatrix.Properties.VariableNames{strcmpi(unMatched.designMatrix.Properties.VariableNames, 'eid')} = 'session_id';

            % Save as a parquet file
            parquetwrite(saveName, unMatched.designMatrix);

            % Extract and save IDs
            ids  = unMatched.designMatrix;
            ids  = ids(:, ismember(ids.Properties.VariableNames, {'participant_id', 'session_id'}));

            if ~isempty(outPrefix)
                outName = [outPrefix, '_idList.parquet'];
            else
                outName = 'FEMA_idList.parquet';
            end
            saveName = fullfile(dirOutput, outName);
            if isempty(ids)
                warning('ID list is empty; something went wrong');
            else
                parquetwrite(saveName, ids);
            end
        end
    else
        if all(ismember({'iid', 'eid'}, unMatched_flds))
            if isnumeric(unMatched.iid)
                unMatched.iid = num2cell(unMatched.iid);
            end
            if isnumeric(unMatched.eid)
                unMatched.eid = num2cell(unMatched.eid);
            end
            ids = cell2table([unMatched.iid, unMatched.eid], ...
                             'VariableNames', {'participant_id', 'session_id'});
            if isempty(ids)
                warning('ID list is empty; something went wrong');
            else
                if ~isempty(outPrefix)
                    outName = [outPrefix, '_idList.parquet'];
                else
                    outName = 'FEMA_idList.parquet';
                end
                saveName = fullfile(dirOutput, outName);
                parquetwrite(saveName, ids);
            end
        else
            error('Unable to extract id list');
        end
    end
    FEMA_save.timing.tSaveIDs = toc(tSaveIDs);
end

%% NIfTI and GIfTI
if any(ismember(outputType, [nii_list gii_list]))
    tSaveImages = tic;
    toSave_FFX_stats  = {'beta_hat', 'beta_se', 'zmat', 'logpmat'};
    toSave_RFX_stats  = {'sig2tvec', 'sig2mat', 'sig2mat_normalized'};
    toSave_Wald_stats = {'Wald', 'logp_Wald'};
    toSave_FFX        = toSave_FFX_stats(ismember(toSave_FFX_stats,   unMatched_flds));
    toSave_RFX        = toSave_RFX_stats(ismember(toSave_RFX_stats,   unMatched_flds));    
    toSave_Wald       = toSave_Wald_stats(ismember(toSave_Wald_stats, unMatched_flds) & ...
                        ~cellfun(@(x) isempty(unMatched.(x)), toSave_Wald_stats));

    
    % Loop over all variables that can be saved, then call writeNIfTI 
    % or writeGIfTI 
    for ff = 1:length(toSave_FFX)
        if ~isempty(outPrefix)
            outName = [outPrefix, '_FFX'];
        else
            outName = 'FFX';
        end
        statName = toSave_FFX{ff}; 
        dispName = dispNameStat(statName);
        workVar = unMatched.(toSave_FFX{ff});

        % Assign values 
        if ~isfield(unMatched, 'colsinterest')
            unMatched.colsinterest = 1:size(workVar,1);
        end 
        for j = 1:length(unMatched.colsinterest)
            jj = unMatched.colsinterest(j);
            
            if ~isfield(unMatched, 'vars_of_interest')
                % Generate vars_of_interest once
                % unMatched.vars_of_interest{j} = ['var', num2str(jj)];
                unMatched.vars_of_interest = strrep(strcat({'var'}, num2str(colvec(unMatched.colsinterest))), ' ', '');
            end 

            if any(ismember(outputType, nii_list))
                saveName_tmp = [outName, '_', toSave_FFX{ff}, '_col', num2str(jj, '%03d'), '_', unMatched.vars_of_interest{j}, '.nii.gz'];
                saveName = fullfile(dirOutput, saveName_tmp);
                saveData = fullvol(workVar(jj,:), mask);
                niftiwrite_amd(saveData, saveName, M_atl_sub);
                json.fixed(j).params.(toSave_FFX{ff}).file_name = saveName_tmp;
            else
                % save
                saveName_tmp = [outName, '_', toSave_FFX{ff}, '_col', num2str(jj, '%03d'), '_', unMatched.vars_of_interest{j}];
                saveName = fullfile(dirOutput, saveName_tmp);
                saveData = workVar(jj,:);
                writeGIfTI(saveData, [], saveName, splitLR);
                json.fixed(j).params.(toSave_FFX{ff}).file_name = {[saveName_tmp '_lh.gii'] ; [saveName_tmp '_rh.gii']};
            end 
            value_range = saveDataRange(saveData);
            json.fixed(j).name = unMatched.vars_of_interest{j}; 
            json.fixed(j).params.(toSave_FFX{ff}).display_name = dispName;
            json.fixed(j).params.(toSave_FFX{ff}).value_range  = value_range;
        end
    end

    for ff=1:length(toSave_RFX)
        if ~isempty(outPrefix)
            outName = [outPrefix, '_RFX'];
        else
            outName = 'RFX';
        end
        if strcmp(toSave_RFX{ff}, 'sig2tvec') 
            workVar = unMatched.(toSave_RFX{ff});

            if any(ismember(outputType, nii_list))
                saveData = fullvol(workVar, mask);
                saveName_tmp = [outName, '_TotVar.nii.gz'];
                saveName = fullfile(dirOutput, saveName_tmp);
                niftiwrite_amd(saveData, saveName, M_atl_sub);
                json.random(1).params.file_name = saveName_tmp;
            else 
                saveData = workVar;
                saveName_tmp = [outName, '_TotVar'];
                saveName = fullfile(dirOutput, saveName_tmp);
                writeGIfTI(saveData, [], saveName, splitLR);
                json.random(1).params.file_name = {[saveName_tmp '_lh.gii'] ; [saveName_tmp '_rh.gii']};
            end 
            % json.random.total_var.file_name = saveName_tmp;
            % json.random.total_var.display_name = 'Total Variance';
            value_range = saveDataRange(saveData);
            json.random(1).name = 'Total variance';
            json.random(1).params.display_name = 'Total Variance';
            json.random(1).params.value_range = value_range;
        elseif  strcmp(toSave_RFX{ff}, 'sig2mat') || strcmp(toSave_RFX{ff}, 'sig2mat_normalized')
            workVar = unMatched.(toSave_RFX{ff});

            dim_workVar = size(workVar);
            if length(dim_workVar) > 2
                for rr = 1:dim_workVar(3)
                    if ~isfield(unMatched, 'RandomEffects')
                        % json.random.sig2mat(rr).name = ['Random Effect ', num2str(rr)];
                        json.random(rr+1).name = ['Random Effect ', num2str(rr)];
                    else
                        % json.random.sig2mat(rr).name = unMatched.RandomEffects{rr};
                        json.random(rr+1).name = RandomEffects_name{rr};
                    end
                    count_var = 1;
                    % count_corr = 1;
                    for i1 = 1:dim_workVar(2)
                        for i2 = 1:i1                            
                            if i1 == i2
                                if isnumeric(eidOrd_all)
                                    eidOrd = num2str(eidOrd_all(i1));
                                else 
                                    eidOrd = eidOrd_all{i1};
                                end 
                                volName = ['Var_', RandomEffects{rr}, '_', eidOrd];   
                                if any(ismember(outputType, nii_list))
                                    saveName_tmp = [outName, '_', volName, '.nii.gz'];
                                    saveName = fullfile(dirOutput, saveName_tmp);
                                    vecData = reshape(squeeze(workVar(i1,i2,rr,:)), 1, []);
                                    saveData = fullvol(vecData, mask);
                                    niftiwrite_amd(saveData, saveName, M_atl_sub);
                                    json.random(rr+1).params.(matlab.lang.makeValidName(volName)).file_name = ...
                                    saveName_tmp;           
                                else
                                    saveName_tmp = [outName, '_', volName];
                                    saveName = fullfile(dirOutput, saveName_tmp);
                                    saveData = reshape(squeeze(workVar(i1,i2,rr,:)), 1, []);
                                    writeGIfTI(saveData, [], saveName, splitLR);
                                    json.random(rr+1).params.(matlab.lang.makeValidName(volName)).file_name = ...
                                                              {[saveName_tmp '_lh.gii'] ; [saveName_tmp '_rh.gii']};
                                end
                                % json.random.sig2mat(rr).variance(count_var).display_name = volName;
                                % json.random.sig2mat(rr).variance(count_var).file_name = saveName_tmp;
                                value_range = saveDataRange(saveData);
                                json.random(rr+1).params.(matlab.lang.makeValidName(volName)).display_name = strrep(volName, 'Var_', 'Variance ');
                                json.random(rr+1).params.(matlab.lang.makeValidName(volName)).value_range = value_range;
                                count_var = count_var + 1;
                            else
                                if isnumeric(eidOrd_all)
                                    eidOrd1 = num2str(eidOrd_all(i1));
                                    eidOrd2 = num2str(eidOrd_all(i2));
                                else
                                    eidOrd1 = eidOrd_all{i1};
                                    eidOrd2 = eidOrd_all{i2};
                                end
                                volName = ['Corr_', RandomEffects{rr}, '_', eidOrd2, '-', eidOrd1];  
                                if any(ismember(outputType, nii_list)) 
                                    saveName_tmp = [outName, '_', volName, '.nii.gz']; 
                                    saveName = fullfile(dirOutput, saveName_tmp);
                                    mask = unMatched.mask;
                                    vecData = reshape(squeeze(workVar(i1,i2,rr,:)), 1, []);
                                    saveData = fullvol(vecData, mask);
                                    niftiwrite_amd(saveData, saveName, M_atl_sub);     
                                    json.random(rr+1).params.(matlab.lang.makeValidName(volName)).file_name = ...
                                        saveName_tmp;                            
                                else
                                    saveName_tmp = [outName, '_', volName];
                                    saveName = fullfile(dirOutput, saveName_tmp);
                                    saveData = reshape(squeeze(workVar(i1,i2,rr,:)), 1, []);
                                    writeGIfTI(saveData, [], saveName, splitLR);
                                    json.random(rr+1).params.(matlab.lang.makeValidName(volName)).file_name = ...
                                        {[saveName_tmp '_lh.gii'] ; [saveName_tmp '_rh.gii']};
                                end
                                % json.random.sig2mat(rr).correlation(count_corr).display_name = volName;
                                % json.random.sig2mat(rr).correlation(count_corr).file_name = saveName_tmp;
                                value_range = saveDataRange(saveData);
                                json.random(rr+1).params.(matlab.lang.makeValidName(volName)).display_name = strrep(volName, 'Corr_', 'Correlation ');
                                
                                json.random(rr+1).params.(matlab.lang.makeValidName(volName)).value_range = value_range;
                                count_var = count_var + 1;
                            end  
                        end
                    end
                end
            else
                for rr = 1:dim_workVar(1)
                    if ~isfield(unMatched, 'RandomEffects')
                        % json.random.sig2mat(rr).name = ['RFX_', num2str(rr, '%02d')];
                        json.random(rr+1).name = ['RFX_', num2str(rr, '%02d')];
                    else
                        % json.random.sig2mat(rr).name = unMatched.RandomEffects{rr};
                        json.random(rr+1).name = RandomEffects_name{rr};
                    end
                    volName = ['Var_' RandomEffects{rr}];
                    % json.random(rr+1).params.sig2mat(rr).variance.display_name = volName;
                    if any(ismember(outputType, nii_list))
                        vecData = workVar(rr,:);
                        saveData = fullvol(vecData, mask);
                        saveName_tmp = [outName, '_', volName, '.nii.gz']; 
                        saveName = fullfile(dirOutput, saveName_tmp);
                        niftiwrite_amd(saveData, saveName, M_atl_sub);
                        json.random(rr+1).params.variance.file_name = saveName_tmp;
                    else 
                        saveName_tmp = [outName, '_', volName];
                        saveName = fullfile(dirOutput, saveName_tmp);  
                        saveData = workVar(rr,:);                        
                        writeGIfTI(saveData, [], saveName, splitLR);
                        json.random(rr+1).params.variance.file_name = {[saveName_tmp '_lh.gii'] ; [saveName_tmp '_rh.gii']};
                    end
                    % json.random.sig2mat(rr).variance.file_name = saveName_tmp;
                    % json.random.sig2mat(rr).file_name = saveName_tmp;
                    value_range = saveDataRange(saveData);
                    json.random(rr+1).params.variance.display_name = 'Variance';
                    json.random(rr+1).params.variance.value_range = value_range;
                end
            end
        end 
    end

    % Save Wald test
    for ff = 1:length(toSave_Wald)
        if ~isempty(outPrefix)
            outName = [outPrefix, '_FFX'];
        else
            outName = 'FFX';
        end
        statName = toSave_Wald{ff}; 
        dispName = dispNameStat(statName);
        workVar = unMatched.(toSave_Wald{ff});

        % Assign values
        for j = 1:size(workVar,1)

            if any(ismember(outputType, nii_list))
                saveName_tmp = [outName, '_', toSave_Wald{ff}, '_', unMatched.splines_of_interest{j,2}, '.nii.gz'];
                saveName = fullfile(dirOutput, saveName_tmp);
                saveData = fullvol(workVar(j,:), mask);
                niftiwrite_amd(saveData, saveName, M_atl_sub);
                json.fixed(length(unMatched.colsinterest)+j).params.(toSave_Wald{ff}).file_name = saveName_tmp;
            else
                saveName_tmp = [outName, '_', toSave_Wald{ff}, '_', unMatched.splines_of_interest{j,2}];
                saveName = fullfile(dirOutput, saveName_tmp);
                saveData = workVar(j,:);
                writeGIfTI(saveData, [], saveName, splitLR);
                json.fixed(length(unMatched.colsinterest)+j).params.(toSave_Wald{ff}).file_name = {[saveName_tmp '_lh.gii'] ; [saveName_tmp '_rh.gii']};
            end
            value_range = saveDataRange(saveData);
            json.fixed(length(unMatched.colsinterest)+j).name = unMatched.splines_of_interest{j,2}; 
            json.fixed(length(unMatched.colsinterest)+j).params.(toSave_Wald{ff}).display_name = dispName;
            json.fixed(length(unMatched.colsinterest)+j).params.(toSave_Wald{ff}).value_range = value_range;
        end
    end    

    if writeJSON
        fname_json = fullfile(dirOutput, 'FEMA_mapping.json');
        tmp = save_jsonencode(json, PrettyPrint=true); 
        fid = fopen(fname_json, 'w');
        fprintf(fid, tmp);
        fclose(fid);
    end

    FEMA_save.timing.tSaveImages = toc(tSaveImages);
end

%% Save corrmat results
if any(ismember(outputType, {'corrmat'}))
    tSaveCorr         = tic;
    toSave_FFX_stats  = {'beta_hat', 'beta_se', 'zmat', 'logpmat'};
    toSave_RFX_stats  = {'sig2tvec', 'sig2mat', 'sig2mat_normalized'};
    toSave_Wald_stats = {'Wald', 'logp_Wald'};
    toSave_FFX        = toSave_FFX_stats(ismember(toSave_FFX_stats,   unMatched_flds));
    toSave_RFX        = toSave_RFX_stats(ismember(toSave_RFX_stats,   unMatched_flds));    
    toSave_Wald       = toSave_Wald_stats(ismember(toSave_Wald_stats, unMatched_flds) & ...
                        ~cellfun(@(x) isempty(unMatched.(x)), toSave_Wald_stats));

    % Loop over all variables that can be saved
    for ff = 1:length(toSave_FFX)
        if ~isempty(outPrefix)
            outName = [outPrefix, '_FFX'];
        else
            outName = 'FFX';
        end
        % statName = toSave_FFX{ff}; 
        % switch statName
        %     case 'beta_hat'
        %         dispName = 'Beta'; 
        %     case 'beta_se'
        %         dispName = 'Standard Error';
        %     case 'zmat'
        %         dispName = 'T statistics';
        %     case 'logpmat'
        %         dispName = 'Signed -log10(p)';
        % end
        workVar = unMatched.(toSave_FFX{ff});

        % Assign values
        if ~isfield(unMatched, 'colsinterest')
            unMatched.colsinterest = 1:size(workVar,1);
        end

        % Loop over every fixed effect of interest
        for j = 1:length(unMatched.colsinterest)
            jj = unMatched.colsinterest(j);
            
            if ~isfield(unMatched, 'vars_of_interest')
                % Generate vars_of_interest once
                unMatched.vars_of_interest = strrep(strcat({'var'}, num2str(colvec(unMatched.colsinterest))), ' ', '');
            end

            % Create JSON object and write
            % tmp = struct;
            % tmp.upperTriangle = workVar(jj,:);
            fname_json = [outName, '_', toSave_FFX{ff}, '_col', num2str(jj, '%03d'), '_', unMatched.vars_of_interest{j}, '.json'];
            if writeJSON
                FEMA_vec2json(workVar(jj,:), fullfile(dirOutput, fname_json));
                % tmp = FEMA_save_jsonencode(tmp, PrettyPrint=true);
                % fid = fopen(fullfile(dirOutput, fname_json), 'w');
                % fprintf(fid, tmp);
                % fclose(fid);
                % Keep track of all files being created
                json.fixed(j).params.(toSave_FFX{ff}).file_name = fname_json;
            end
            json.fixed(j).name = unMatched.vars_of_interest{j}; 
            % json.fixed(j).params.(toSave_FFX{ff}).display_name = dispName;
        end
    end

    % Random effects
    for ff=1:length(toSave_RFX)
        if ~isempty(outPrefix)
            outName = [outPrefix, '_RFX'];
        else
            outName = 'RFX';
        end

        if strcmpi(toSave_RFX{ff}, 'sig2tvec') 

             % Create JSON object and write
            % tmp = struct;
            % tmp.upperTriangle = unMatched.(toSave_RFX{ff});
            fname_json = [outName, '_TotVar.json'];
            if writeJSON
                FEMA_vec2json(unMatched.(toSave_RFX{ff}), fullfile(dirOutput, fname_json));
                % tmp = FEMA_save_jsonencode(tmp, PrettyPrint=true);
                % fid = fopen(fullfile(dirOutput, fname_json), 'w');
                % fprintf(fid, tmp);
                % fclose(fid);
                % Keep track of all files being created
                json.random(1).params.file_name = fname_json;
            end
            json.random(1).name = 'Total variance';
            % json.random(1).params.display_name = 'Total Variance';
        elseif strcmp(toSave_RFX{ff}, 'sig2mat') || strcmp(toSave_RFX{ff}, 'sig2mat_normalized')
            workVar     = unMatched.(toSave_RFX{ff});
            dim_workVar = size(workVar);
            if length(dim_workVar) > 2
                for rr = 1:dim_workVar(3)
                    if ~isfield(unMatched, 'RandomEffects')
                        json.random(rr+1).name = ['Random Effect ', num2str(rr)];
                    else
                        json.random(rr+1).name = RandomEffects_name{rr};
                    end
                    count_var = 1;
                    for i1 = 1:dim_workVar(2)
                        for i2 = 1:i1                            
                            if i1 == i2
                                if isnumeric(eidOrd_all)
                                    eidOrd = num2str(eidOrd_all(i1));
                                else 
                                    eidOrd = eidOrd_all{i1};
                                end

                                % Create JSON object and write
                                % tmp = struct;
                                % tmp.upperTriangle = reshape(squeeze(workVar(i1,i2,rr,:)), 1, []);
                                % tmp = FEMA_save_jsonencode(tmp, PrettyPrint=true);
                                volName = ['Var_', RandomEffects{rr}, '_', eidOrd];
                                fname_json = [outName, '_', volName, '.json'];
                                if writeJSON
                                    FEMA_vec2json(reshape(squeeze(workVar(i1,i2,rr,:)), 1, []), fullfile(dirOutput, fname_json));
                                    % fid = fopen(fullfile(dirOutput, fname_json), 'w');
                                    % fprintf(fid, tmp);
                                    % fclose(fid);
                                    % Keep track of all files being created
                                    json.random(rr+1).params.(matlab.lang.makeValidName(volName)).file_name = fname_json;
                                end
                                % json.random(rr+1).params.(matlab.lang.makeValidName(volName)).display_name = strrep(volName, 'Var_', 'Variance ');
                                count_var = count_var + 1;
                            else
                                if isnumeric(eidOrd_all)
                                    eidOrd1 = num2str(eidOrd_all(i1));
                                    eidOrd2 = num2str(eidOrd_all(i2));
                                else
                                    eidOrd1 = eidOrd_all{i1};
                                    eidOrd2 = eidOrd_all{i2};
                                end

                                % Create JSON object and write
                                % tmp = struct;
                                % tmp.upperTriangle = reshape(squeeze(workVar(i1,i2,rr,:)), 1, []);
                                % tmp = FEMA_save_jsonencode(tmp, PrettyPrint=true);
                                volName = ['Corr_', RandomEffects{rr}, '_', eidOrd2, '-', eidOrd1];
                                fname_json = [outName, '_', volName, '.json'];
                                if writeJSON
                                    FEMA_vec2json(reshape(squeeze(workVar(i1,i2,rr,:)), 1, []), fullfile(dirOutput, fname_json));
                                    % fid = fopen(fullfile(dirOutput, fname_json), 'w');
                                    % fprintf(fid, tmp);
                                    % fclose(fid);
                                    % Keep track of all files being created
                                    json.random(rr+1).params.(matlab.lang.makeValidName(volName)).file_name = fname_json;
                                end
                                % json.random(rr+1).params.(matlab.lang.makeValidName(volName)).display_name = strrep(volName, 'Corr_', 'Correlation ');
                                count_var = count_var + 1;
                            end
                        end
                    end
                end
            else
                for rr = 1:dim_workVar(1)
                    if ~isfield(unMatched, 'RandomEffects')
                        json.random(rr+1).name = ['RFX_', num2str(rr, '%02d')];
                    else
                        json.random(rr+1).name = RandomEffects_name{rr};
                    end

                    % Create JSON object and write
                    % tmp = struct;
                    % tmp.upperTriangle = workVar(rr,:);
                    % tmp = FEMA_save_jsonencode(tmp, PrettyPrint=true);
                    volName = ['Var_' RandomEffects{rr}];
                    fname_json = [outName, '_', volName, '.json'];
                    if writeJSON
                        FEMA_vec2json(workVar(rr,:), fullfile(dirOutput, fname_json));
                        % fid = fopen(fullfile(dirOutput, fname_json), 'w');
                        % fprintf(fid, tmp);
                        % fclose(fid);
                        % Keep track of all files being created
                        json.random(rr+1).params.variance.file_name = fname_json;
                    end
                    % json.random(rr+1).params.variance.display_name = 'Variance';
                end
            end
        end 
    end

    % Save Wald test
    for ff = 1:length(toSave_Wald)
        if ~isempty(outPrefix)
            outName = [outPrefix, '_FFX'];
        else
            outName = 'FFX';
        end
        % statName = toSave_Wald{ff}; 
        % switch statName
        %     case 'Wald'
        %         dispName = 'Wald statistics'; 
        %     case 'p_Wald'
        %         dispName = 'Wald p value';
        %     case 'logp_Wald'
        %         dispName = 'Wald -log10(p)';
        % end
        workVar = unMatched.(toSave_Wald{ff});

        % Assign values
        for j = 1:size(workVar,1)

            % Create JSON object and write
            % tmp = struct;
            % tmp.upperTriangle = workVar(j,:);
            fname_json = [outName, '_', toSave_Wald{ff}, '_', unMatched.splines_of_interest{j,2}, '.json'];
            if writeJSON
                FEMA_vec2json(workVar(j,:), fullfile(dirOutput, fname_json));
                % tmp = FEMA_save_jsonencode(tmp, PrettyPrint=true);
                % fid = fopen(fullfile(dirOutput, fname_json), 'w');
                % fprintf(fid, tmp);
                % fclose(fid);
                % Keep track of all files being created
                json.fixed(length(unMatched.colsinterest)+j).params.(toSave_Wald{ff}).file_name = fname_json;
            end
            json.fixed(length(unMatched.colsinterest)+j).name = unMatched.splines_of_interest{j,2};
            % json.fixed(length(unMatched.colsinterest)+j).params.(toSave_Wald{ff}).display_name = dispName;
        end
    end

    % Write out mapping file
    if writeJSON
        fname_json = fullfile(dirOutput, 'FEMA_mapping.json');
        tmp = save_jsonencode(json, PrettyPrint=true); 
        fid = fopen(fname_json, 'w');
        fprintf(fid, tmp);
        fclose(fid);
    end

    % Write out label file: single file with ROI names
    if writeJSON && isfield(unMatched, 'ymat_names')
        fname_json = fullfile(dirOutput, 'Corrmat_labels.json');
        tmp = save_jsonencode(unMatched.ymat_names, PrettyPrint=true); 
        fid = fopen(fname_json, 'w');
        fprintf(fid, tmp);
        fclose(fid);
    
        % fname_json = fullfile(dirOutput, 'Corrmat_label_source.json');
        % tmp = FEMA_save_jsonencode(unMatched.ymat_names(:,1), PrettyPrint=true); 
        % fid = fopen(fname_json, 'w');
        % fprintf(fid, tmp);
        % fclose(fid);
        % 
        % fname_json = fullfile(dirOutput, 'Corrmat_label_target.json');
        % tmp = FEMA_save_jsonencode(unMatched.ymat_names(:,2), PrettyPrint=true); 
        % fid = fopen(fname_json, 'w');
        % fprintf(fid, tmp);
        % fclose(fid);
    end
    FEMA_save.timing.tSaveCorr = toc(tSaveCorr);
end

%% Save tabulated results
if any(ismember(outputType, {'tables', 'external'}))
    tSaveTables = tic;

    % Object to keep track of filenames
    json_fileTracker = struct;
    track = 1;

    % Part I: grouping by dependent variable - every dependent variable
    %         gets its own separate JSON file

    % Create a structure for display names: fixed effects
    display_fixed.beta_hat  = 'Beta';
    display_fixed.beta_se   = 'Standard Error';
    display_fixed.zmat      = 'T statistics';
    display_fixed.logpmat   = 'Signed -log10(p)';
    display_fixed.Wald      = 'Wald statistics';
    display_fixed.logp_Wald = 'Wald -log10(p)';

    % Create a structure for display names: random effects
    if ndims(unMatched.sig2mat) == 2
        display_random.variance = 'Variance';
    else
        display_random.variance = 'Variance';
        display_random.varcorr  = 'Variance-Correlation';
    end
    % display_random.totVar = 'Total variance';

    % Prepare fixed effects part of the formula
    unMatched.colnames_model = rowvec(unMatched.colnames_model);
    %ffx_formula = [strcat(unMatched.colnames_model(1:length(unMatched.colnames_model)-1), {' + '}), ...
    %                      unMatched.colnames_model(end)];
    if writeJSON
        ffx_formula = [strcat(unMatched.FFX_conceptMapping(1:size(unMatched.FFX_conceptMapping, 1)-1,1)', {' + '}), ...
                              unMatched.FFX_conceptMapping(end,1)];
        ffx_formula = horzcat(ffx_formula{:});

        % Prepare random effects part of the formula
        if ndims(unMatched.sig2mat) == 2
            rfx_formula = [{' + '}, strcat({'(1|'}, strrep(rowvec(RandomEffects_name(1:length(RandomEffects_name)-1)), ' ', ''), {') + '}), strrep(RandomEffects_name(end), ' ', '')];
            rfx_formula = horzcat(rfx_formula{:});
        else
            % unstructured covariance matrix, variances and covariances
            if isnumeric(eidOrd_all)
                eidOrd_all = cellstr(num2str(eidOrd_all'));
            end

            % Prepare random effects part of the formula
            rfx_formula = [{' + '}, strcat({'us(eid + 1|'}, strrep(rowvec(RandomEffects_name(1:length(RandomEffects_name)-1)), ' ', ''), {') + '}), strrep(RandomEffects_name(end), ' ', '')];
            rfx_formula = horzcat(rfx_formula{:});
        end
    end 

    % How many coefficients?
    howMany = size(unMatched.beta_hat,1);

    % Name of the fields: source and target
    sourceFieldNames = {'beta_hat', 'beta_se', 'zmat', 'logpmat'};
    destFieldNames   = {'beta_hat', 'beta_se', 'zmat', 'logpmat'};

    % Work out variable source
    if ~isempty(unMatched.FFX_conceptMapping)
        varSourceInfo = cell(howMany, 1);
        for v = 1:howMany
            tmpX = unMatched.colnames_model{v};
            if strcmpi(tmpX, 'Intercept')
                varSourceInfo{v} = 'Intercept';
            else
                varSourceInfo{v} = unMatched.FFX_conceptMapping{cellfun(@(x) ismember(tmpX, x), unMatched.FFX_conceptMapping(:,2)), 1};
            end
        end
        % Where do these variables come from?
        [tab_fixed(1:howMany).variableSource] = varSourceInfo{:};
    end

    % Which fixed effects are of interest?
    tmp_ofInterest = false(howMany, 1);
    tmp_ofInterest(unMatched.colsinterest) = true;
    tmp_ofInterest = num2cell(tmp_ofInterest);

    for yy = 1:length(unMatched.ymat_names)
        % Initialize structure
        json_tabularY = struct;
        tab_fixed     = struct;

        % Dependent variable name
        depName = unMatched.ymat_names{yy};

        % Prepare JSON data
        json_tabularY.dependent = depName;
        
        if writeJSON
            formula = strcat(depName, {' ~ '}, ffx_formula, rfx_formula);
            json_tabularY.formula = formula;
        end
        
        % Prepare fixed effects statistics
        % Assign names
        [tab_fixed(1:howMany).name] = unMatched.colnames_model{:};
        
        % Assign statistics
        for ff = 1:length(sourceFieldNames)
            tmp = num2cell(unMatched.(sourceFieldNames{ff})(:,yy));
            [tab_fixed(1:size(tmp,1)).(destFieldNames{ff})] = tmp{:};        
        end

        % Which of these are of interest?
        [tab_fixed(1:howMany).of_interest] = tmp_ofInterest{:};

        % Now work on Wald test
        count = howMany+1;
        for ff = 1:size(unMatched.splines_of_interest,1)
            tab_fixed(count).name           = unMatched.splines_of_interest{ff,2};
            tab_fixed(count).Wald           = unMatched.Wald(ff,yy);
            tab_fixed(count).logp_Wald      = unMatched.logp_Wald(ff,yy);
            tab_fixed(count).variableSource = strrep(unMatched.splines_of_interest{ff,2}, '_splines', ''); % unMatched.FFX_conceptMapping{cellfun(@(x) ismember(strrep(unMatched.splines_of_interest{ff,2}, '_splines', ''), x), unMatched.FFX_conceptMapping(:,2)), 1};
            tab_fixed(count).of_interest    = true;
            count = count + 1;
        end

        % Variables that are of interest should be sorted to the top
        [~, b] = sort([(tab_fixed(:).of_interest)], 'descend');
        tab_fixed = tab_fixed(1,b);

        % Assign id
        ids = num2cell(0:size(tab_fixed,2)-1);
        [tab_fixed(1:size(tab_fixed,2)).id] = ids{:};

        % Assign fixed effects structure to json_tabularY
        json_tabularY.fixed       = tab_fixed;
        json_tabularY.display_name = display_fixed;

        % Save as JSON file
        if writeJSON
            table_json = save_jsonencode(json_tabularY, PrettyPrint=true);
            if ~isempty(outPrefix)
                outName = [outPrefix, '_FFX_byY_', unMatched.ymat_names{yy}, '.json'];
            else
                outName = ['FFX_byY_', unMatched.ymat_names{yy}, '.json'];
            end
            fid = fopen(fullfile(dirOutput, outName), 'w');
            fwrite(fid, table_json);
            fclose(fid);
            json_fileTracker.fixed.byY{track} = outName;
            track = track + 1;
        end

        % Now work on random effects
        % Initialize structure
        json_tabularY = struct;
        tab_random    = struct;

        % Prepare JSON data
        json_tabularY.dependent = depName;
        if isfield(unMatched.unstructParams, 'sig2mat_normalized')
            if isnumeric(eidOrd_all)
                eidOrd_all = cellstr(num2str(eidOrd_all'));
            end
            [a, b]                      = ndgrid(eidOrd_all);
            json_tabularX.rowNames      = a(:,1);
            json_tabularX.colNames      = b(:,1);
            % json_tabularY.varcorr_names = strcat(b, {'-'}, a);
        end

        % Regular compound symmetry / diagonal covariance
        if ndims(unMatched.sig2mat) == 2
            for rr = 1:length(unMatched.RandomEffects)
                tab_random(rr).id       = rr-1;
                tab_random(rr).name     = RandomEffects_name{rr};
                tab_random(rr).variance = unMatched.sig2mat(rr,yy);
            end
            tab_random(rr+1).id       = rr;
            tab_random(rr+1).name     = 'Total variance';
            tab_random(rr+1).variance = unMatched.sig2tvec(yy);
        else
            if isfield(unMatched.unstructParams, 'sig2mat_normalized')
                for rr=1:size(unMatched.unstructParams.sig2mat_normalized, 3)
                    tab_random(rr).id      = rr-1;
                    tmp_varcorr            = squeeze(unMatched.unstructParams.sig2mat_normalized(:,:,rr,yy));
                    tab_random(rr).varcorr = tmp_varcorr;
                    if strcmpi(RandomEffects_name{rr}, 'Subject effect')
                         tab_random(rr).name = 'Subject effect and Error';
                    else
                        tab_random(rr).name = RandomEffects_name{rr};
                    end
                end
                tab_random(rr+1).id       = rr;
                tab_random(rr+1).name     = 'Total variance';
                tab_random(rr+1).variance = unMatched.sig2tvec(yy);
                json_tabularY.rowNames    = eidOrd_all;
                json_tabularY.colNames    = eidOrd_all;
            else
                error('sig2mat_normalized field missing');
            end
        end

        % Assign random effects structure to json_tabularY
        json_tabularY.random       = tab_random;
        json_tabularY.display_name = display_random;

        % Save as JSON file
        if writeJSON
            table_json = save_jsonencode(json_tabularY, PrettyPrint=true);
            if ~isempty(outPrefix)
                outName = [outPrefix, '_RFX_byY_', unMatched.ymat_names{yy}, '.json'];
            else
                outName = ['RFX_byY_', unMatched.ymat_names{yy}, '.json'];
            end
            fid = fopen(fullfile(dirOutput, outName), 'w');
            fwrite(fid, table_json);
            fclose(fid);
            json_fileTracker.random.byY{track} = outName;
            track = track + 1;
        end
    end

    % Part II: each fixed effect gets its own separate JSON file
    % Prepare fixed effects part of the formula
    % Create a structure for display names: fixed effects without Wald
    display_fixed           = struct;
    display_fixed.beta_hat  = 'Beta';
    display_fixed.beta_se   = 'Standard Error';
    display_fixed.zmat      = 'T statistics';
    display_fixed.logpmat   = 'Signed -log10(p)';

    % Create a structure for display names: fixed effects with Wald
    display_fixed_Wald.Wald      = 'Wald statistics';
    display_fixed_Wald.logp_Wald = 'Wald -log10(p)';

    % Create a structure for display names: random effects
    if ndims(unMatched.sig2mat) == 2
        display_random.variance = 'Variance';
    else
        display_random.varcorr  = 'Variance-Correlation';
    end
    % display_random.totVar = 'Total variance';

    unMatched.colnames_model = rowvec(unMatched.colnames_model);
    if writeJSON
        %ffx_formula = [strcat(unMatched.colnames_model(1:length(unMatched.colnames_model)-1), {' + '}), ...
        %                      unMatched.colnames_model(end)];
        ffx_formula = [strcat(unMatched.FFX_conceptMapping(1:size(unMatched.FFX_conceptMapping, 1)-1,1)', {' + '}), ...
                              unMatched.FFX_conceptMapping(end,1)];
        ffx_formula = horzcat(ffx_formula{:});

        % Prepare random effects part of the formula
        if ndims(unMatched.sig2mat) == 2
            rfx_formula = [{' + '}, strcat({'(1|'}, strrep(rowvec(RandomEffects_name(1:length(RandomEffects_name)-1)), ' ', ''), {') + '}), strrep(RandomEffects_name(end), ' ', '')];
            rfx_formula = horzcat(rfx_formula{:});
        else
            % unstructured covariance matrix, variances and covariances
            if isnumeric(eidOrd_all)
                eidOrd_all = cellstr(num2str(eidOrd_all'));
            end

            % Prepare random effects part of the formula
            rfx_formula = [{' + '}, strcat({'us(eid + 1|'}, strrep(rowvec(RandomEffects_name(1:length(RandomEffects_name)-1)), ' ', ''), {') + '}), strrep(RandomEffects_name(end), ' ', '')];
            rfx_formula = horzcat(rfx_formula{:});
        end

        formula = strcat({'y ~ '}, ffx_formula, rfx_formula);
    end
   
    % How many dependent variables?
    howMany = size(unMatched.beta_hat,2);

    % Name of the fields: source and target
    sourceFieldNames = {'beta_hat', 'beta_se', 'zmat', 'logpmat'};
    destFieldNames   = {'beta_hat', 'beta_se', 'zmat', 'logpmat'};

    % Go over fixed effects
    for xx = 1:size(unMatched.beta_hat,1)
        % Is this covariate of interest?
        if any(xx == unMatched.colsinterest)
            tmp_ofInterest = true;
        else
            tmp_ofInterest = false;
        end

        % Initialize structure
        json_tabularX = struct;
        tab_dependent = struct;

        % Fixed effect name
        fixName = unMatched.colnames_model{xx};

        % Prepare JSON data
        json_tabularX.fixed       = fixName;
        if writeJSON
            json_tabularX.formula = formula;
        end
        json_tabularX.of_interest = tmp_ofInterest;

        % Prepare fixed effects statistics
        % Assign names
        [tab_dependent(1:howMany).name] = unMatched.ymat_names{:};

        % Assign statistics
        for ff = 1:length(sourceFieldNames)
            tmp = num2cell(unMatched.(sourceFieldNames{ff})(xx,:));
            [tab_dependent(1:size(tmp,2)).(destFieldNames{ff})] = tmp{:};
        end

        % Assign id
        ids = num2cell(0:size(tab_dependent,2)-1);
        [tab_dependent(1:size(tab_dependent,2)).id] = ids{:};

        % Assign fixed effects structure to json_tabularY
        json_tabularX.dependent    = tab_dependent;
        json_tabularX.display_name = display_fixed;

        % Save as JSON file
        if writeJSON
            table_json = save_jsonencode(json_tabularX, PrettyPrint=true);
            if ~isempty(outPrefix)
                outName = [outPrefix, '_FFX_byX_', fixName, '.json'];
            else
                outName = ['FFX_byX_', fixName, '.json'];
            end
            fid = fopen(fullfile(dirOutput, outName), 'w');
            fwrite(fid, table_json);
            fclose(fid);
            json_fileTracker.fixed.byX{track} = outName;
            track = track + 1;
        end
    end

    % Now work on Wald test
    for xx = 1:size(unMatched.splines_of_interest,1)
        % Initialize structure
        json_tabularX = struct;
        tab_dependent = struct;

        % Fixed effect name
        fixName = unMatched.splines_of_interest{xx,2};

        % Prepare JSON data
        json_tabularX.fixed       = fixName;
        if writeJSON
            json_tabularX.formula = formula;
        end
        json_tabularX.of_interest = true;

        % Assign names
        [tab_dependent(1:howMany).name] = unMatched.ymat_names{:};

        % Assign statistics
        tmp = num2cell(unMatched.Wald(xx,:));
        [tab_dependent(1:howMany).Wald] = tmp{:};

        tmp = num2cell(unMatched.logp_Wald(xx,:));
        [tab_dependent(1:howMany).logp_Wald] = tmp{:};

        % Assign id
        ids = num2cell(0:size(tab_dependent,2)-1);
        [tab_dependent(1:size(tab_dependent,2)).id] = ids{:};
        
        % Assign fixed effects structure to json_tabularY
        json_tabularX.dependent    = tab_dependent;
        json_tabularX.display_name = display_fixed_Wald;

        % Save as JSON file
        if writeJSON
            table_json = save_jsonencode(json_tabularX, PrettyPrint=true);
            if ~isempty(outPrefix)
                outName = [outPrefix, '_FFX_byX_', fixName, '.json'];
            else
                outName = ['FFX_byX_', fixName, '.json'];
            end
            fid = fopen(fullfile(dirOutput, outName), 'w');
            fwrite(fid, table_json);
            fclose(fid);
            json_fileTracker.fixed.byX{track} = outName;
            track = track + 1;
        end
    end

    % Now work on random effects
    % Regular compound symmetry / diagonal covariance
    if ndims(unMatched.sig2mat) == 2
        for rr = 1:length(unMatched.RandomEffects)
            % Initialize structure
            json_tabularX = struct;
            tab_dependent = struct;
    
            % Prepare JSON data
            json_tabularX.random = RandomEffects_name{rr};

            % Assign names
            [tab_dependent(1:howMany).name] = unMatched.ymat_names{:};

            % Assign statistics
            tmp = num2cell(unMatched.sig2mat(rr,:));
            [tab_dependent(1:size(tmp,2)).variance] = tmp{:};

            % Assign id
            ids = num2cell(0:size(tab_dependent,2)-1);
            [tab_dependent(1:size(tab_dependent,2)).id] = ids{:};

            % Assign fixed effects structure to json_tabularY
            json_tabularX.dependent    = tab_dependent;
            json_tabularX.display_name = display_random;

            % Save as JSON file
            if writeJSON
                table_json = save_jsonencode(json_tabularX, PrettyPrint=true);
                if ~isempty(outPrefix)
                    outName = [outPrefix, '_RFX_byX_', unMatched.RandomEffects{rr}, '.json'];
                else
                    outName = ['RFX_byX_', unMatched.RandomEffects{rr}, '.json'];
                end
                fid = fopen(fullfile(dirOutput, outName), 'w');
                fwrite(fid, table_json);
                fclose(fid);
                json_fileTracker.random.byX{track} = outName;
                track = track + 1;
            end
        end
    else
        for rr = 1:length(unMatched.RandomEffects)-1
            % Initialize structure
            json_tabularX = struct;
            tab_dependent = struct;

            if isnumeric(eidOrd_all)
                eidOrd_all = cellstr(num2str(eidOrd_all'));
            end
            % [a, b]                      = ndgrid(eidOrd_all);
            json_tabularX.rowNames      = eidOrd_all;
            json_tabularX.colNames      = eidOrd_all;
            % json_tabularX.varcorr_names = strcat(b, {'-'}, a);

            % Assign names
            % [tab_dependent(1:howMany).name] = unMatched.ymat_names{:};
    
            % Unstructured covariance
            if isfield(unMatched.unstructParams, 'sig2mat_normalized')
                for yy = 1:length(unMatched.ymat_names)
                    tmp                       = squeeze(unMatched.unstructParams.sig2mat_normalized(:,:,rr,yy));
                    tab_dependent(yy).id      = yy-1;
                    tab_dependent(yy).name    = unMatched.ymat_names{yy};
                    tab_dependent(yy).varcorr = tmp;
                end
            else
                error('sig2mat_normalized field missing');
            end

            % Assign id
            ids = num2cell(0:size(tab_dependent,2)-1);
            [tab_dependent(1:size(tab_dependent,2)).id] = ids{:};

            % Assign fixed effects structure to json_tabularY
            json_tabularX.dependent    = tab_dependent;
            json_tabularX.display_name = struct('varcorr', 'Variance-correlation');

            % Save as JSON file
            if writeJSON
                table_json = save_jsonencode(json_tabularX, PrettyPrint=true);
                if ~isempty(outPrefix)
                    outName = [outPrefix, '_RFX_byX_', unMatched.RandomEffects{rr}, '.json'];
                else
                    outName = ['RFX_byX_', unMatched.RandomEffects{rr}, '.json'];
                end
                fid = fopen(fullfile(dirOutput, outName), 'w');
                fwrite(fid, table_json);
                fclose(fid);
                json_fileTracker.random.byX{track} = outName;
                track = track + 1;
            end

        end
    end

    % Now do the total variance
    % Initialize structure
    json_tabularX = struct;
    tab_dependent = struct;

    % Prepare JSON data
    json_tabularX.random = 'Total variance';

    % Assign names
    [tab_dependent(1:howMany).name] = unMatched.ymat_names{:};

    % Assign variances
    tmp = num2cell(unMatched.sig2tvec);
    [tab_dependent(1:size(tmp,2)).variance] = tmp{:};

    % Assign id
    ids = num2cell(0:size(tab_dependent,2)-1);
    [tab_dependent(1:size(tab_dependent,2)).id] = ids{:};    

    % Assign fixed effects structure to json_tabularY
    json_tabularX.dependent    = tab_dependent;
    json_tabularX.display_name = struct('variance', 'Variance');

    % Save as JSON file
    if writeJSON
        table_json = save_jsonencode(json_tabularX, PrettyPrint=true);
        if ~isempty(outPrefix)
            outName = [outPrefix, '_RFX_byX_TotVar.json'];
        else
            outName = 'RFX_byX_TotVar.json';
        end
        fid = fopen(fullfile(dirOutput, outName), 'w');
        fwrite(fid, table_json);
        fclose(fid);
        json_fileTracker.random.byX{track} = outName;

        % Any empty entries in json_fileTracker should be deleted
        json_fileTracker.fixed.byX(cellfun(@isempty, json_fileTracker.fixed.byX)) = [];
        json_fileTracker.fixed.byY(cellfun(@isempty, json_fileTracker.fixed.byY)) = [];

        json_fileTracker.random.byX(cellfun(@isempty, json_fileTracker.random.byX)) = [];
        json_fileTracker.random.byY(cellfun(@isempty, json_fileTracker.random.byY)) = [];

        % Write out json_fileTracker
        json_fileTracker = save_jsonencode(json_fileTracker, PrettyPrint=true);
        fid = fopen(fullfile(dirOutput, 'FEMA_mapping.json'), 'w');
        fwrite(fid, json_fileTracker);
        fclose(fid);
    end
    
    FEMA_save.timing.tSaveTables = toc(tSaveTables);
end 

% Prepare fixed effects concept mapping and write it as a JSON file
if ismember('FFX_conceptMapping', unMatched_flds)
    tConceptMapping = tic;

    if writeJSON
        if ~isempty(outPrefix)
            outName = [outPrefix, '_FEMA_conceptMapping.json'];
        else
            outName = 'FEMA_conceptMapping.json';
        end

        fid = fopen(fullfile(dirOutput, outName), 'w');
        fwrite(fid, save_jsonencode(unMatched.FFX_conceptMapping', PrettyPrint=true));
        fclose(fid);
    end

    FEMA_save.timing.tSaveConceptMapping = toc(tConceptMapping);
end

if any(ismember(outputType, {'roi'}))
    roi_atlas = unMatched.roi_atlas;
    if any(ismember({'at', 'aseg', 'dsk', 'dst', 'aparc', 'aparc_a2009s'}, roi_atlas))
        tSaveROI = tic;
        
        json_roi = struct('fixed', [], 'random', [], 'atlas', atlas2tabNames(roi_atlas, 'tab2atlas'));
        % map ymat_names to tabulated names
        switch roi_atlas
            case {'aseg'}
                [ymat_names ymat_codes] = aseg2tabroinames(unMatched.ymat_names, 'tab2atlas');
            case {'at'}
                [ymat_names ymat_codes] = fiber2tabroinames(unMatched.ymat_names, 'tab2atlas');
            case {'dsk', 'aparc'}
                [ymat_names ymat_codes] = aparc2tabroinames(unMatched.ymat_names, 'tab2atlas');
            case {'dst', 'aparc_a2009s'}
                [ymat_names ymat_codes] = aparc2009s2tabroinames(unMatched.ymat_names, 'tab2atlas');
            otherwise
                error('Unknown ROI atlas: %s. No mapping to tabulated names.', roi_atlas);
        end
    
        % remove "lh", "rh"  and "all" from ymat_names
        idx_rm = find(ismember(ymat_names, {'lh', 'rh', 'all'}));
        ymat_names(idx_rm) = [];
        ymat_codes(idx_rm) = [];
    
        % get indices of left and right hemisphere (should this be a function?)
        isLeft = ~cellfun(@isempty, regexp(lower(ymat_names), '(^|[_\.-])(lh|left|l)([_\.-]|$)', 'once'));
        isRight = ~cellfun(@isempty, regexp(lower(ymat_names), '(^|[_\.-])(rh|right|r)([_\.-]|$)', 'once'));
        idxLeft = find(isLeft);
        idxRight = find(isRight);
    
        % Prepare statistics to save
        toSave_FFX_stats  = {'beta_hat', 'beta_se', 'zmat', 'logpmat'};
        toSave_RFX_stats  = {'sig2tvec', 'sig2mat', 'sig2mat_normalized'};
        toSave_Wald_stats = {'Wald', 'logp_Wald'};
        toSave_FFX        = toSave_FFX_stats(ismember(toSave_FFX_stats,   unMatched_flds));
        toSave_RFX        = toSave_RFX_stats(ismember(toSave_RFX_stats,   unMatched_flds));    
        toSave_Wald       = toSave_Wald_stats(ismember(toSave_Wald_stats, unMatched_flds) & ...
                            ~cellfun(@(x) isempty(unMatched.(x)), toSave_Wald_stats));
    
        % Fixed effects
        if ~isempty(toSave_FFX)
            if ~isfield(unMatched, 'colsinterest')
                unMatched.colsinterest = 1:size(unMatched.(toSave_FFX{1}),1);
            end 
            if ~isfield(unMatched, 'vars_of_interest')
                % Generate vars_of_interest once
                % unMatched.vars_of_interest{j} = ['var', num2str(jj)];
                unMatched.vars_of_interest = strrep(strcat({'var'}, num2str(colvec(unMatched.colsinterest))), ' ', '');
            end       
            % save out one json for each fixed effect
            for j = 1:length(unMatched.colsinterest)
                jj = unMatched.colsinterest(j);
                saveName_tmp = [outName, '_col', num2str(jj, '%03d'), '_', unMatched.vars_of_interest{j}, '.json'];
                saveName = fullfile(dirOutput, saveName_tmp);
    
                for ff = 1:length(toSave_FFX)
                    workVar = unMatched.(toSave_FFX{ff});
                    statName = toSave_FFX{ff};
                    dispName = dispNameStat(statName);
                    saveData = workVar(jj,:);
                    value_range = saveDataRange(saveData);
    
                    % data to save
                    json_roi.fixed(j).name = unMatched.vars_of_interest{j};
                    json_roi.fixed(j).params.(toSave_FFX{ff}).display_name = dispName;
                    json_roi.fixed(j).params.(toSave_FFX{ff}).value_range = value_range;
                    json_roi.fixed(j).params.(toSave_FFX{ff}).values.left = roiDataStruct(saveData, ymat_codes, idxLeft);
                    json_roi.fixed(j).params.(toSave_FFX{ff}).values.right = roiDataStruct(saveData, ymat_codes, idxRight);
                end 
            end 
        end
    
        % Save Wald test
        if ~isempty(toSave_Wald)
            nFFX = length(json_roi.fixed);
            for j = 1:size(unMatched.splines_of_interest, 1)
                saveName_tmp = [outName, '_', unMatched.splines_of_interest{j,2}, '.json'];
                saveName = fullfile(dirOutput, saveName_tmp);
    
                for ff = 1:length(toSave_Wald)
                    workVar = unMatched.(toSave_Wald{ff});
                    statName = toSave_Wald{ff};
                    dispName = dispNameStat(statName);
                    saveData = workVar(j,:);
                    value_range = saveDataRange(saveData);
                    % data to save
                    %json_fileTracker.fixed.roi(j).params.(statName).display_name = dispName;
                    %json_fileTracker.fixed.roi(j).params.(statName).value_range = value_range;
                    json_roi.fixed(nFFX+j).name = unMatched.splines_of_interest{j,2};
                    json_roi.fixed(nFFX+j).params.(toSave_Wald{ff}).display_name = dispName;
                    json_roi.fixed(nFFX+j).params.(toSave_Wald{ff}).value_range = value_range;
                    json_roi.fixed(nFFX+j).params.(toSave_Wald{ff}).values.left = roiDataStruct(saveData, ymat_codes, idxLeft);
                    json_roi.fixed(nFFX+j).params.(toSave_Wald{ff}).values.right = roiDataStruct(saveData, ymat_codes, idxRight);
                end 
            end 
        end
        
        if ~isempty(toSave_RFX)
            if ~isempty(outPrefix)
                outName = [outPrefix, '_RFX'];
            else
                outName = 'RFX';
            end
            saveName_tmp = [outName, '_all.json'];
            saveName = fullfile(dirOutput, saveName_tmp);
    
            for ff = 1:length(toSave_RFX)
                if strcmp(toSave_RFX{ff}, 'sig2tvec') 
                    statName = toSave_RFX{ff};
                    dispName = dispNameStat(statName);
                    workVar = unMatched.(toSave_RFX{ff});
                    value_range = saveDataRange(workVar);
    
                    json_roi.random(1).name = dispName;
                    json_roi.random(1).params.display_name = dispName;
                    json_roi.random(1).params.value_range = value_range;
    
                    json_roi.random(1).params.values.left = roiDataStruct(workVar, ymat_codes, idxLeft);
                    json_roi.random(1).params.values.right = roiDataStruct(workVar, ymat_codes, idxRight);
    
                elseif  strcmp(toSave_RFX{ff}, 'sig2mat') || strcmp(toSave_RFX{ff}, 'sig2mat_normalized')
                    statName = toSave_RFX{ff};
                    if isempty(unMatched.unstructParams) 
                        dispName = dispNameStat(statName);                    
                        workVar = unMatched.(toSave_RFX{ff});
                        dim_workVar = size(workVar);
                        for rr = 1:dim_workVar(1)
                            vecData = workVar(rr, :);
                            value_range = saveDataRange(vecData);
                            json_roi.random(rr+1).name = RandomEffects_name{rr};
                            json_roi.random(rr+1).params.variance.display_name = dispName;
                            json_roi.random(rr+1).params.variance.value_range = value_range;
                            json_roi.random(rr+1).params.variance.values.left = roiDataStruct(vecData, ymat_codes, idxLeft);
                            json_roi.random(rr+1).params.variance.values.right = roiDataStruct(vecData, ymat_codes, idxRight);
                        end
                    else 
                        workVar = unMatched.(toSave_RFX{ff});
                        dim_workVar = size(workVar);
                        for rr = 1:dim_workVar(3)
                            for i1 = 1:dim_workVar(2)
                                for i2 = 1:i1   
                                    vecData = reshape(squeeze(workVar(i1,i2,rr,:)), 1, []);     
                                    if i1 == i2
                                        if isnumeric(eidOrd_all)
                                            eidOrd = num2str(eidOrd_all(i1));
                                        else 
                                            eidOrd = eidOrd_all{i1};
                                        end 
                                        dispName = ['Variance ', RandomEffects{rr}, '_', eidOrd];
                                        volName = ['Var_', RandomEffects{rr}, '_', eidOrd];
                                        volName = matlab.lang.makeValidName(volName);
                                        value_range = saveDataRange(vecData);
                                        json_roi.random(rr+1).name =  RandomEffects_name{rr}; 
                                        json_roi.random(rr+1).params.(volName).display_name = dispName;
                                        json_roi.random(rr+1).params.(volName).value_range = value_range;
                                        json_roi.random(rr+1).params.(volName).values.left = roiDataStruct(vecData, ymat_codes, idxLeft);
                                        json_roi.random(rr+1).params.(volName).values.right = roiDataStruct(vecData, ymat_codes, idxRight);
                                    else                                
                                        if isnumeric(eidOrd_all)
                                            eidOrd1 = num2str(eidOrd_all(i1));
                                            eidOrd2 = num2str(eidOrd_all(i2));
                                        else
                                            eidOrd1 = eidOrd_all{i1};
                                            eidOrd2 = eidOrd_all{i2};
                                        end
                                        dispName = ['Correlation ', RandomEffects{rr}, '_', eidOrd2, '-', eidOrd1];
                                        volName = ['Corr_', RandomEffects{rr}, '_', eidOrd2, '-', eidOrd1];
                                        volName = matlab.lang.makeValidName(volName);
                                        value_range = saveDataRange(vecData);
                                        json_roi.random(rr+1).name = RandomEffects_name{rr};
                                        json_roi.random(rr+1).params.(volName).display_name = dispName;
                                        json_roi.random(rr+1).params.(volName).value_range = value_range;
                                        json_roi.random(rr+1).params.(volName).values.left = roiDataStruct(vecData, ymat_codes, idxLeft);
                                        json_roi.random(rr+1).params.(volName).values.right = roiDataStruct(vecData, ymat_codes, idxRight);
                                    end  
                                end
                            end
                        end
                    end 
                end 
            end
        end
        % save json_roi
        if ~isempty(outPrefix)
            outName = [outPrefix, 'FFX_RFX_ROI_all.json'];
        else
            outName = 'FFX_RFX_ROI_all.json';
        end
        saveName = fullfile(dirOutput, outName);
        if writeJSON
            roi_json = save_jsonencode(json_roi, PrettyPrint=true);
            roi_json = regexprep(roi_json, '"x([0-9]+)"\s*:', '"$1":');
            fid = fopen(saveName, 'w');
            fwrite(fid, roi_json);
            fclose(fid);
        end
    
        FEMA_save.timing.tSaveROI = toc(tSaveROI);
        % save FEMA_mapping.json (json_fileTracker is the encoded string from tables block when present)
        if writeJSON
            if exist('json_fileTracker', 'var') && (ischar(json_fileTracker) || isstring(json_fileTracker))
                json_fileTracker = jsondecode(json_fileTracker);
            else
                json_fileTracker = struct(struct('roi', {{}}), 'random', struct('roi', {{}}));
            end
            json_fileTracker.roi = saveName;
            json_fileTracker = save_jsonencode(json_fileTracker, PrettyPrint=true);
            fid = fopen(fullfile(dirOutput, 'FEMA_mapping.json'), 'w');
            fwrite(fid, json_fileTracker);
            fclose(fid);
        end
    else
        warning('Unknown ROI atlas: %s. No mapping outputs.', roi_atlas);
    end
end

FEMA_save.timing.tSaveOverall = toc(tSaveOverall);
if exist('fname_mat', 'var') && exist(fname_mat, 'file')
    info.FEMA_save = FEMA_save;
    save(fname_mat, 'info', '-append');
end 
end

function value_range = saveDataRange(saveData)
    cal_min = min(saveData, [], 'all', 'omitmissing');
    cal_max = max(saveData, [], 'all', 'omitmissing');
    if isinf(cal_max)
        cal_max = max(saveData(~isinf(saveData)));
    end
    if isinf(cal_min)
        cal_min = min(saveData(~isinf(saveData)));
    end
    value_range = [cal_min, cal_max];
end

function S = roiDataStruct(saveData, ymat_codes, idx)
    vals = saveData(idx);
    codes = ymat_codes(idx);
    S = struct;
    for k = 1:numel(codes)
        S.(['x', num2str(codes(k))]) = vals(k);
    end
end

function FEMA_vec2json(vector, filename)
    % "Impoverished" function for converting a 1xv vector to JSON
    % Specifically designed for corrmat
    
    % Open file for writing
    fid = fopen(filename, 'w');
    
    % Opening content and brackets
    fprintf(fid, '%s\n%s\n', '{', '"upperTriangle":[');
    
    % Write all lines but the last entry, comma ending
    fprintf(fid, '%g,\n', vector(1,1:end-1));
    
    % Write the last line, no comma
    fprintf(fid, '%g\n', vector(1,end));
    
    % Closing brackets
    fprintf(fid, '%s\n%s', ']', '}');
    
    % Close the file
    fclose(fid);
end
 