function FEMA_save(outputType, dirOutput, outPrefix, saveDesignMatrix, varargin)
% Function to save different output from FEMA in different format

% Get a sense for all variables in the workspace
tmpInfo = whos;

% Ensure outputType is specified as lower case
outputType = lower(outputType);

% Maybe allow cell type outputType and use if matching to write multiple
% things in a single call?
switch(outputType)

    case '.mat'

    % All possible things to save for estimates
    toSave_estimates = {'beta_hat', 'beta_se', 'zmat', 'logpmat', 'sig2tvec', 'sig2mat', ...
                        'beta_hat_perm', 'beta_se_perm', 'zmat_perm', 'sig2tvec_perm',   ...
                        'sig2mat_perm', 'logLikvec', 'logLikvec_perm', 'Hessmat',        ...
                        'coeffCovar', 'binvec_save', 'nvec_bins', 'tvec_bins',           ...
                        'reusableVars', 'contrasts', 'hypValues'};

    toSave_design    = {'X', 'fid', 'iid', 'eid', 'agevec', 'FamilyStruct', ...
                        'MotherID', 'FatherID', 'HomeID', 'PregID'};

    toSave_settings  = {'niter', 'contrasts', 'nbins', 'RandomEffects', 'nperms',      ...
                        'CovType', 'FixedEstType', 'RandomEstType', 'GroupByFamType',  ...
                        'NonnegFlag', 'SingleOrDouble', 'logLikflag', 'PermType',      ...
                        'returnReusable', 'doPar', 'numWorkers', 'numThreads', 'hypValues'};
    
    % First: save main (permuted and non-permuted) statistics
    saveName = fullfile(dirOutput, [outPrefix, '_estimates.mat']);
    toSave   = toSave_estimates(ismember(toSave_estimates, {tmpInfo(:).name}));
    
    if sum([tmpInfo(ismember({tmpInfo(:).name}', toSave)).bytes]) > 2^31
        save(saveName, toSave{:}, '-v7.3');
        % save(saveName, 'beta_hat', 'beta_se', 'zmat', 'logpmat', 'sig2tvec', 'sig2mat', ...
        %         'beta_hat_perm', 'beta_se_perm', 'zmat_perm', 'sig2tvec_perm',          ...
        %         'sig2mat_perm', 'logLikvec', 'logLikvec_perm', 'Hessmat',               ...
        %         'coeffCovar', 'binvec_save', 'nvec_bins', 'tvec_bins',                  ...
        %         'reusableVars', 'contrasts', 'hypValues', '-v7.3');
    else
        save(saveName, toSave{:});
        % save(saveName, 'beta_hat', 'beta_se', 'zmat', 'logpmat', 'sig2tvec', 'sig2mat', ...
        %         'beta_hat_perm', 'beta_se_perm', 'zmat_perm', 'sig2tvec_perm',          ...
        %         'sig2mat_perm', 'logLikvec', 'logLikvec_perm', 'Hessmat',               ...
        %         'coeffCovar', 'binvec_save', 'nvec_bins', 'tvec_bins',                  ...
        %         'reusableVars', 'contrasts', 'hypValues');
    end
    
    % Second: save filtered design matrix and related variables
    if saveDesignMatrix
        saveName = fullfile(dirOutput, [outPrefix, '_designMatrix.mat']);
        toSave   = toSave_design(ismember(toSave_design, {tmpInfo(:).name}));
        if sum([tmpInfo(ismember({tmpInfo(:).name}', toSave)).bytes]) > 2^31
            save(saveName, toSave{:}, '-v7.3');
            % save(saveName, 'X', 'fid', 'iid', 'eid', 'agevec', 'FamilyStruct', ...
            %                'MotherID', 'FatherID', 'HomeID', 'PregID', '-v7.3');
        else
            % save(saveName, 'X', 'fid', 'iid', 'eid', 'agevec', 'FamilyStruct', ...
            save(saveName, toSave{:});
            %                'MotherID', 'FatherID', 'HomeID', 'PregID');
        end
    end
    
    % Third: save settings (this does not require v7.3 check)
    saveName = fullfile(dirOutput, [outPrefix, '_settings.mat']);
    toSave   = toSave_settings(ismember(toSave_settings, {tmpInfo(:).name}));
    save(saveName, toSave{:});
    % save(saveName, 'niter', 'contrasts', 'nbins', 'RandomEffects', 'nperms',      ...
    %                'CovType', 'FixedEstType', 'RandomEstType', 'GroupByFamType',  ...
    %                'NonnegFlag', 'SingleOrDouble', 'logLikflag', 'PermType',      ...
    %                'returnReusable', 'doPar', 'numWorkers', 'numThreads', 'hypValues');

    case {'.nii', '.nii.gz', 'nifti', 'voxel'}
        % Need to additionally save JSON colnames
        % Bring in the updated file naming convention
        toSave_NIfTI = {'beta_hat', 'beta_se', 'zmat', 'logpmat', 'sig2tvec', 'sig2mat', ...
                        'beta_hat_perm', 'beta_se_perm', 'zmat_perm', 'sig2tvec_perm',   ...
                        'sig2mat_perm'};
        toSave       = toSave_NIfTI(ismember(toSave_NIfTI, {tmpInfo(:).name}));

        % Loop over all variables that can be saved, call them a generic
        % variable, expand to mask (need to ensure this is passed in), and
        % then call writeNIfTI (need to allow non structure character
        % input to writeNIfTI): maybe make a different writeNIfTI that only
        % takes data, output directory and name as inputs?
        for ff = 1:length(toSave)
            % Assign to generic variable
            workVar = eval(toSave{ff});

            % Initialize: 3D or 4D
            tmp_dim = size(workVar,2);
            volData = zeros([size(mask) tmp_dim]);
            
            % Assign values and call writeNIfTI
            for j = 1:tmp
                volData(:,:,:,j) = single(fullvol(workVar(j,:),mask));
                writeNIFTI(volData, dirOutput, outPrefix);
            end
        end

    case {'.gii', 'gifti', 'vertex'}
        % Need to call writeGIfTI

    case {'tables'}
        % Need to dump regression tables into text/csv files

    case {'summary'}
        % This needs to be a JSON file


        %%


    %write column names to json for DEAP
    fname_col = sprintf('%s/FEMA_results_colnames.json',dirname_out{des});
    out = struct('colnames_model',{colnames_model},'RandomEffects',{RandomEffects});
    jsonStr = jsonencode(out);
    fid = fopen(fname_col,'w');
    fprintf(fid,'%s\n',jsonStr);
    fclose(fid);


        % =========================================================================
        % Write VERTEX results
        % =========================================================================
    elseif strcmpi(datatype, 'vertex')

        if contains(outputFormat,'mat')

            if nperms>0 & tfce==0
                save(fpath_out,base_variables_to_save{:},'zmat_perm','beta_hat_perm','colnames_interest','colsinterest','-v7.3');
            elseif nperms>0 & tfce==1
                save(fpath_out,base_variables_to_save{:},'zmat_perm','beta_hat_perm','tfce_perm','colnames_interest','colsinterest','-v7.3');
            elseif nperms==0
                save(fpath_out,base_variables_to_save{:},'-v7.3');
            end

            logging('Results written to %s',fpath_out);
        end

        if contains(outputFormat, 'nifti') %FIXME: these are much smaller, so haven't added the same optimization as for voxelwise

            randomFields = {'sig2tvec', 'sig2mat'};

            results = struct('beta_hat',beta_hat,'beta_se',beta_se,'zmat',zmat,'logpmat',logpmat,'sig2tvec',sig2tvec,'sig2mat',sig2mat);
            % Write out in FreeSurfer curv format, with naming consistent with volumes
            fieldnamelist = setdiff(fieldnames(results),randomFields);
            icnum = ico+1;
            load(fullfile(fileparts(fileparts(which('FEMA_wrapper'))), 'showSurf', 'SurfView_surfs.mat'), 'icsurfs');
            % load('~/matlab/cache/SurfView_surfs.mat'); % this does not include white
            S = struct;
            S.nverts = 2*size(icsurfs{icnum}.vertices,1);
            S.nfaces = 2*size(icsurfs{icnum}.faces,1);
            S.faces = cat(1,icsurfs{icnum}.faces,icsurfs{icnum}.faces+size(icsurfs{icnum}.vertices,1));
            % parse IVs
            if isempty(ivnames)
                excludeCol = strmatch('mri_info_',colnames_model);
                nCol = length(colnames_model);
                ivCol = setdiff(1:nCol, excludeCol);
            else
                [~,ivCol,~] = intersect(colnames_model,ivnames);
            end
            if length(ivCol) < 1, error('No IVs found! Not writing nifti.'), end
            % write out the fixed effects
            for fi = 1:length(fieldnamelist)
                fieldname = fieldnamelist{fi};
                vol_nifti = results.(fieldname);
                for iv = ivCol(:)'
                    colname = sprintf('FE%02d',iv);
                    valvec = vol_nifti(iv,:);
                    fname_out = sprintf('%s/FEMA_results_vertexwise_%s_%s_%s.fsvals',dirname_out{1},fstem_imaging,fieldname,colname);
                    fs_write_curv(fname_out,valvec,S.nfaces);
                    fprintf(1,'file %s written\n',fname_out);
                end
            end
            % write out the random effects
            fieldnamelist = randomFields;
            for fi = 1:length(fieldnamelist)
                fieldname = fieldnamelist{fi};
                vol_nifti = results.(fieldname);
                for iv = 1:size(vol_nifti,4)
                    colname = sprintf('RE%02d',iv);
                    valvec = vol_nifti(iv,:);
                    fname_out = sprintf('%s/FEMA_results_vertexwise_%s_%s_%s.fsvals',dirname_out{1},fstem_imaging,fieldname,colname);
                    fs_write_curv(fname_out,valvec,S.nfaces);
                    fprintf(1,'file %s written\n',fname_out);
                end
            end

        end

    elseif strcmpi(datatype, 'corrmat')

        if 0
            figure; im = reshape(beta_hat(1,:),dims(2:end)); imagesc(im,max(abs(im(:)))*[-1 1]); colormap(blueblackred); axis equal tight;
            figure; im = reshape(beta_hat(2,:),dims(2:end)); imagesc(im,max(abs(im(:)))*[-1 1]); colormap(blueblackred); axis equal tight;
            figure; im = reshape(zmat(1,:),dims(2:end)); imagesc(im,max(abs(im(:)))*[-1 1]); colormap(blueblackred); axis equal tight;
            figure; im = reshape(zmat(2,:),dims(2:end)); imagesc(im,max(abs(im(:)))*[-1 1]); colormap(blueblackred); axis equal tight;
        end
        %beta_hat = reshape(beta_hat,[size(beta_hat,1) dims(2:end)]);
        %beta_se = reshape(beta_se,[size(beta_hat,1) dims(2:end)]);
        %zmat = reshape(zmat,[size(beta_hat,1) dims(2:end)]);
        %logpmat = reshape(logpmat,[size(beta_hat,1) dims(2:end)]);
        %sig2tvec = reshape(sig2tvec,[size(sig2tvec,1) dims(2:end)]);
        %sig2mat = reshape(sig2mat,[size(sig2mat,1) dims(2:end)]);

        if nperms>0
            save(fpath_out,base_variables_to_save{:},'colnames_imaging','zmat_perm','beta_hat_perm','colnames_interest','colsinterest','-v7.3');
        else
            save(fpath_out,base_variables_to_save{:},'colnames_imaging','-v7.3');
        end
        logging('Results written to %s',fpath_out);

    end

    if ~isempty(fpath_out)
        fpaths_out = cat(2,fpaths_out,fpath_out);
    end
end