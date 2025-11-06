function FEMA_save(outputType, dirOutput, outPrefix, varargin)
% Function to save different output from FEMA in different format

p = inputParser;
p.KeepUnmatched = true;
addParameter(p, 'saveDesignMatrix', true, @islogical);

parse(p, varargin{:});

toSave_struct = unMatched;
toSave_struct_flds = fieldnames(toSave_struct);
saveDesignMatrix = p.Results.saveDesignMatrix;
unMatched = p.Unmatched;
% Ensure outputType is specified as lower case
outputType = lower(outputType);

nii_list = {'.nii', '.nii.gz', 'nifti', 'voxel'}; 
gii_list = {'.gii', 'gifti', 'vertex'};
splitLR = true; 
M_atl     = [0    -1     0   101; 0     0     1  -131; -1     0     0   101; 0     0     0     1];
M_atl_sub = [0    -2     0   102; 0     0     2  -132; -2     0     0   102; 0     0     0     1];
    
% Maybe allow cell type outputType and use if matching to write multiple
% things in a single call?
if ismember('.mat', outputType)
    % All possible things to save for estimates
    toSave_estimates = {'beta_hat', 'beta_se', 'zmat', 'logpmat', 'sig2tvec', 'sig2mat', ...
                        'beta_hat_perm', 'beta_se_perm', 'zmat_perm', 'sig2tvec_perm',   ...
                        'sig2mat_perm', 'logLikvec', 'logLikvec_perm', 'Hessmat',        ...
                        'coeffCovar', 'binvec_save', 'nvec_bins', 'tvec_bins',           ...
                        'residuals_GLS', 'unstructParams', 'contrasts', 'hypValues', 'info', ...
                        'colnames_model'};

    toSave_design    = {'designMatrix', 'colnames_model', 'X', 'fid', 'iid', 'eid', 'agevec', ...
                        'FamilyStruct', 'MotherID', 'FatherID', 'HomeID', 'PregID'};
    
    % First: save main (permuted and non-permuted) statistics
    saveName = fullfile(dirOutput, [outPrefix, '_estimates.mat']);
    %toSave   = toSave_estimates(ismember(toSave_estimates, fieldnames(toSave_struct)));
    toDelete = toSave_struct_flds(~ismember(toSave_struct_flds, toSave_estimates));
    toSave_struct = rmfield(unMatched, toDelete);
    checkSize = whos('toSave_struct');
    if checkSize.bytes > 2^31
        save(saveName,  '-struct', 'toSave_struct', '-v7.3');
    else
        save(saveName,  '-struct', 'toSave_struct');
    end
    
    % Second: save filtered design matrix and related variables
    if saveDesignMatrix
        saveName = fullfile(dirOutput, [outPrefix, '_designMatrix.mat']);
        toSave   = toSave_design(ismember(toSave_design, fieldnames(toSave_struct)));
        toDelete = toSave_struct_flds(~ismember(toSave_struct_flds, toSave_design));
        toSave_struct = rmfield(unMatched, toDelete);
        checkSize = whos('toSave_struct');
        if checkSize.bytes > 2^31
            save(saveName,  '-struct', 'toSave_struct', '-v7.3');
        else
            save(saveName,  '-struct', 'toSave_struct');
        end
    end
end 

if any(ismember(outputType, [nii_list gii_list]))
    toSave_FFX_NIfTI = {'beta_hat', 'beta_se', 'zmat', 'logpmat'};
    toSave_RFX_NIfTI = {'sig2tvec', 'sig2mat', 'sig2mat_normalised'};
    toSave_FFX      = toSave_FFX_NIfTI(ismember(toSave_FFX_NIfTI, toSave_struct_flds));
    toSave_RFX      = toSave_RFX_NIfTI(ismember(toSave_RFX_NIfTI, toSave_struct_flds));
    
    % Loop over all variables that can be saved, then call writeNIfTI 
    % or writeGIfTI 
    for ff = 1:length(toSave_FFX)
        outPrefix = [outName '_FFX'];
        statName = toSave_FFX{ff}; 
        switch statName
            case 'beta_hat'
                dispName = 'Beta'; 
            case 'beta_se'
                dispName = 'Standard Error';
            case 'zmat'
                dispName = 'T statistics';
            case 'logpmat'
                dispName = 'Signed -log10(p)';
        end
        workVar = unMatched.(toSave_FFX{ff});

        % Assign values 
        if ~isfield(unMatched, 'colsinterest')
            unMatched.colsinterest = 1:size(workVar,1);
        end 
        for j = 1:length(unMatched.colsinterest)
            jj = unMatched.colsinterest(j);
            
            if ~isfield(unMatched, 'vars_of_interest')
                unMatched.vars_of_interest{j} = ['var', num2str(jj)];
            end 

            if any(ismember(outputType, nii_list))
                mask = unMatched.mask;
                saveName_tmp = [outPrefix, '_', toSave_FFX{ff}, '_col', num2str(jj, '%03d'), '_', unMatched.    vars_of_interest{j}, '.nii.gz'];
                saveName = fullfile(dirOutput, saveName_tmp);
                saveData = fullvol(workVar(jj,:), mask);
                niftiwrite_amd(saveData, saveName, M_atl_sub, true);
            else
                saveName_tmp = [outPrefix, '_', toSave_FFX{ff}, '_col', num2str(jj, '%03d'), '_', unMatched.    vars_of_interest{j}, '.gii.gz'];
                saveName = fullfile(dirOutput, saveName_tmp);
                saveData = workVar(jj,:);
                writeGIfTI(saveData, [], saveName, [], [], splitLR);
            end 
            json.fixed(j).name = unMatched.vars_of_interest{j}; 
            json.fixed(j).params.(toSave_FFX{ff}).file_name = saveName_tmp; 
            json.fixed(j).params.(toSave_FFX{ff}).display_name = dispName;
        end
    end    
    for ff=1:length(toSave_RFX)
        outPrefix = 'RFX';
        if strcmp(toSave_RFX{ff}, 'sig2tvec') 
            workVar = unMatched.(toSave_RFX{ff});
            if any(ismember(outputType, nii_list))
                mask = unMatched.mask;
                saveData = fullvol(workVar, mask);
                saveName_tmp = [outPrefix, '_TotVar.nii.gz'];
                saveName = fullfile(dirOutput, saveName_tmp);
                niftiwrite_amd(saveData, saveName, M_atl_sub, true);
            else 
                saveData = workVar;
                saveName_tmp = [outPrefix, '_TotVar.gii.gz'];
                saveName = fullfile(dirOutput, saveName_tmp);
                writeGIfTI(saveData, [], saveName, [], [], splitLR);
            end 
            json.random.total_var.file_name = saveName_tmp;
            json.random.total_var.display_name = 'Total Variance';
        elseif  strcmp(toSave_RFX{ff}, 'sig2mat')|| strcmp(toSave_RFX{ff}, 'sig2mat_normalised')
            workVar = unMatched.(toSave_RFX{ff});
            dim_workVar = size(workVar);
            if length(dim_workVar) > 2
                for rr = 1:dim_workVar(3)
                    if ~isfield(unMatched, 'RandomEffects')
                        json.random.sig2mat(rr).name = ['Random Effect ', num2str(rr)];
                    else
                        json.random.sig2mat(rr).name = unMatched.RandomEffects{rr};
                    end
                    count_var = 1; 
                    count_corr = 1;
                    for i1 = 1:dim_workVar(2)
                        for i2 = 1:i1                            
                            if i1 == i2
                                if isnumeric(unMatched.eidOrd)
                                    eidOrd = num2str(unMatched.eidOrd(i1));
                                else 
                                    eidOrd = unMatched.eidOrd{i1};
                                end 
                                volName = ['Var_', eidOrd];   
                                if any(ismember(outputType, nii_list))
                                    saveName_tmp = [outPrefix, '_', volName, '.nii.gz'];
                                    saveName = fullfile(dirOutput, saveName_tmp);
                                    mask = unMatched.mask;
                                    vecData = squeeze(workVar(i1,i2,rr,:));
                                    saveData = fullvol(vecData, mask);
                                    niftiwrite_amd(saveData, saveName, M_atl_sub, true);           
                                else
                                    saveName_tmp = [outPrefix, '_', volName, '.gii.gz'];
                                    saveName = fullfile(dirOutput, saveName_tmp);
                                    saveData = squeeze(workVar(i1,i2,rr,:));                        
                                    writeGIfTI(saveData, [], saveName, [], [], splitLR);
                                end
                                json.random.sig2mat(rr).variance(count_var).display_name = volName;
                                json.random.sig2mat(rr).variance(count_var).file_name = saveName_tmp;
                                count_var = count_var + 1;
                            else
                                if isnumeric(unMatched.eidOrd)
                                    eidOrd1 = num2str(unMatched.eidOrd(i1));
                                    eidOrd2 = num2str(unMatched.eidOrd(i2));
                                else
                                    eidOrd1 = unMatched.eidOrd{i1};
                                    eidOrd2 = unMatched.eidOrd{i2};
                                end
                                volName = ['Corr_', eidOrd2, '-', eidOrd1];  
                                if any(ismember(outputType, nii_list)) 
                                    saveName_tmp = [outPrefix, '_', volName, '.nii.gz']; 
                                    saveName = fullfile(dirOutput, saveName_tmp);
                                    mask = unMatched.mask;
                                    vecData = squeeze(workVar(i1,i2,rr,:));
                                    saveData = fullvol(vecData, mask);
                                    niftiwrite_amd(saveData, saveName, M_atl_sub, true);                                 
                                else
                                    saveName_tmp = [outPrefix, '_', volName, '.gii.gz'];
                                    saveName = fullfile(dirOutput, saveName_tmp);  
                                    saveData = squeeze(workVar(i1,i2,rr,:));                        
                                    writeGIfTI(saveData, [], saveName, [], [], splitLR);
                                end
                                json.random.sig2mat(rr).correlation(count_corr).display_name = volName;
                                json.random.sig2mat(rr).correlation(count_corr).file_name = saveName_tmp;
                                count_corr = count_corr + 1;
                            end  
                        end
                    end
                end
            else
                for rr = 1:dim_workVar(1)
                    if ~isfield(unMatched, 'RandomEffects')
                        json.random.sig2mat(rr).name = ['RFX_', num2str(rr, '%02d')];
                    else
                        json.random.sig2mat(rr).name = unMatched.RandomEffects{rr};
                    end
                    volName = ['Var'];
                    json.random.sig2mat(rr).variance.display_name = volName;
                    if any(ismember(outputType, nii_list))
                        vecData = workVar(rr,:);
                        saveData = fullvol(vecData, mask);
                        saveName_tmp = [outPrefix, '_', volName, '.nii.gz']; 
                        saveName = fullfile(dirOutput, saveName_tmp);
                        niftiwrite_amd(saveData, saveName, M_atl_sub, true);
                    else 
                        saveName_tmp = [outPrefix, '_', volName, '.gii.gz'];
                        saveName = fullfile(dirOutput, saveName_tmp);  
                        saveData = workVar(rr,:);                        
                        writeGIfTI(saveData, [], saveName, [], [], splitLR);
                    end
                    json.random.sig2mat(rr).variance.file_name = saveName_tmp;
                end
            end
        end 
    end
    fname_json = fullfile(dirOutput, 'FEMA_NIfTI_mapping.json');
    tmp = jsonencode(json, PrettyPrint=true); 
    fid = fopen(fname_json, 'w');
    fprintf(fid, tmp);
    fclose(fid);
end 



    case {'tables'}
        % Need to dump regression tables into text/csv files

    case {'summary'}
        % This needs to be a JSON file
end  
