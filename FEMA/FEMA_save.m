function FEMA_save(outputType, dirOutput, outPrefix, saveDesignMatrix, varargin)
% Function to save different output from FEMA in different format

p = inputParser;
p.KeepUnmatched = true;

parse(p, varargin{:});

toSave_struct = p.Unmatched;
toSave_struct_flds = fieldnames(toSave_struct);

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

    toSave_info  = {'niter', 'contrasts', 'nbins', 'RandomEffects', 'nperms',      ...
                        'CovType', 'FixedEstType', 'RandomEstType', 'GroupByFamType',  ...
                        'NonnegFlag', 'SingleOrDouble', 'logLikflag', 'PermType',      ...
                        'returnReusable', 'doPar', 'numWorkers', 'numThreads', 'hypValues'};
    
    % First: save main (permuted and non-permuted) statistics
    saveName = fullfile(dirOutput, [outPrefix, '_estimates.mat']);
    %toSave   = toSave_estimates(ismember(toSave_estimates, fieldnames(toSave_struct)));
    toDelete = toSave_struct_flds(~ismember(toSave_struct_flds, toSave_estimates));
    toSave_struct = rmfield(p.Unmatched, toDelete);
    checkSize = whos('toSave_struct');
    if checkSize.bytes > 2^31
        save(saveName,  '-struct', 'toSave_struct', '-v7.3');
    else
        save(saveName,  '-struct', 'toSave_struct');
    end
    
    % Second: save filtered design matrix and related variables
    if saveDesignMatrix
        saveName = fullfile(dirOutput, [outPrefix, '_designMatrix.mat']);
        %toSave   = toSave_design(ismember(toSave_design, fieldnames(toSave_struct)));
        toDelete = toSave_struct_flds(~ismember(toSave_struct_flds, toSave_design));
        toSave_struct = rmfield(p.Unmatched, toDelete);
        checkSize = whos('toSave_struct');
        if checkSize.bytes > 2^31
            save(saveName,  '-struct', 'toSave_struct', '-v7.3');
        else
            save(saveName,  '-struct', 'toSave_struct');
        end
    end
    
    % Third: save info (this does not require v7.3 check)
    saveName = fullfile(dirOutput, [outPrefix, '_info.mat']);
    %toSave   = toSave_info(ismember(toSave_info, {tmpInfo(:).name}));
    toDelete = toSave_struct_flds(~ismember(toSave_struct_flds, toSave_info));
    toSave_struct = rmfield(p.Unmatched, toDelete);
    save(saveName,  '-struct', 'toSave_struct');


    case {'.nii', '.nii.gz', 'nifti', 'voxel'}
        % Need to additionally save JSON colnames
        % Bring in the updated file naming convention
        toSave_NIfTI = {'beta_hat', 'beta_se', 'zmat', 'logpmat', 'sig2tvec', 'sig2mat', ...
                        'beta_hat_perm', 'beta_se_perm', 'zmat_perm', 'sig2tvec_perm',   ...
                        'sig2mat_perm'};
        toSave       = toSave_NIfTI(ismember(toSave_NIfTI, toSave_struct_flds));
        mask = p.Unmatched.mask;
        % Loop over all variables that can be saved, call them a generic
        % variable, expand to mask (need to ensure this is passed in), and
        % then call writeNIfTI (need to allow non structure character
        % input to writeNIfTI): maybe make a different writeNIfTI that only
        % takes data, output directory and name as inputs?
        for ff = 1:length(toSave)
            % Assign to generic variable
            workVar = p.Unmatched.(toSave{ff});

            M_atl     = [0    -1     0   101; 0     0     1  -131; -1     0     0   101; 0     0     0     1];
            M_atl_sub = [0    -2     0   102; 0     0     2  -132; -2     0     0   102; 0     0     0     1];
            % Assign values and call writeNIfTI
            for j = 1:length(p.Unmatched.colsinterest)
                jj = p.Unmatched.colsinterest(j);
                volData = single(fullvol(workVar(jj,:), mask));
                saveName = fullfile(dirOutput, [outPrefix, '_', toSave{ff}, '_', num2str(jj, '%03d'), '_', p.Unmatched.vars_of_interest{j}, '.nii.gz']);
                niftiwrite_amd(volData, saveName, M_atl_sub, true);
            end
        end

    case {'.gii', 'gifti', 'vertex'}
        % Need to call writeGIfTI
            % Need to additionally save JSON colnames
        % Bring in the updated file naming convention
        toSave_GIfTI = {'beta_hat', 'beta_se', 'zmat', 'logpmat', 'sig2tvec', 'sig2mat', ...
                        'beta_hat_perm', 'beta_se_perm', 'zmat_perm', 'sig2tvec_perm',   ...
                        'sig2mat_perm'};
        toSave       = toSave_GIfTI(ismember(toSave_GIfTI, toSave_struct_flds));
        % Loop over all variables that can be saved, call them a generic
        % variable, expand to mask (need to ensure this is passed in), and
        % then call writeNIfTI (need to allow non structure character
        % input to writeNIfTI): maybe make a different writeNIfTI that only
        % takes data, output directory and name as inputs?
        for ff = 1:length(toSave)
            % Assign to generic variable
            workVar = p.Unmatched.(toSave{ff});
            surfData = single(workVar(p.Unmatched.colsinterest,:));
            basename = [outPrefix, '_', toSave{ff}];
            colnames_gii = strcat(num2str(p.Unmatched.colsinterest, '%03d'), {'_'}, p.Unmatched.vars_of_interest);
            % Assign values and call writeGIfTI
            for j = 1:length(p.Unmatched.colsinterest)
                jj = p.Unmatched.colsinterest(j);
                saveName = fullfile(dirOutput, [outPrefix, '_', toSave{ff}, '_', num2str(jj), '_', p.Unmatched.vars_of_interest{j}, '.gii.gz']);
                writeGIfTI(surfData, [], outDir, basename, colnames_gii);
            end

    case {'tables'}
        % Need to dump regression tables into text/csv files

    case {'summary'}
        % This needs to be a JSON file


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
%    elseif strcmpi(datatype, 'vertex')
%
%        if contains(outputFormat,'mat')
%
%            if nperms>0 & tfce==0
%                save(fpath_out,base_variables_to_save{:},'zmat_perm','beta_hat_perm','colnames_interest','colsinterest',%'-v7.3');
%            elseif nperms>0 & tfce==1
%                save(fpath_out,base_variables_to_save{:},'zmat_perm','beta_hat_perm','tfce_perm','colnames_interest',%'colsinterest','-v7.3');
%            elseif nperms==0
%                save(fpath_out,base_variables_to_save{:},'-v7.3');
%            end
%
%            logging('Results written to %s',fpath_out);
%        end
%
%        if contains(outputFormat, 'nifti') %FIXME: these are much smaller, so haven't added the same optimization as %for voxelwise
%
%            randomFields = {'sig2tvec', 'sig2mat'};
%
%            results = struct('beta_hat',beta_hat,'beta_se',beta_se,'zmat',zmat,'logpmat',logpmat,'sig2tvec',sig2tvec,%'sig2mat',sig2mat);
%            % Write out in FreeSurfer curv format, with naming consistent with volumes
%            fieldnamelist = setdiff(fieldnames(results),randomFields);
%            icnum = ico+1;
%            load(fullfile(fileparts(fileparts(which('FEMA_wrapper'))), 'showSurf', 'SurfView_surfs.mat'), 'icsurfs');
%            % load('~/matlab/cache/SurfView_surfs.mat'); % this does not include white
%            S = struct;
%            S.nverts = 2*size(icsurfs{icnum}.vertices,1);
%            S.nfaces = 2*size(icsurfs{icnum}.faces,1);
%            S.faces = cat(1,icsurfs{icnum}.faces,icsurfs{icnum}.faces+size(icsurfs{icnum}.vertices,1));
%            % parse IVs
%            if isempty(ivnames)
%                excludeCol = strmatch('mri_info_',colnames_model);
%                nCol = length(colnames_model);
%                ivCol = setdiff(1:nCol, excludeCol);
%            else
%                [~,ivCol,~] = intersect(colnames_model,ivnames);
%            end
%            if length(ivCol) < 1, error('No IVs found! Not writing nifti.'), end
%            % write out the fixed effects
%            for fi = 1:length(fieldnamelist)
%                fieldname = fieldnamelist{fi};
%                vol_nifti = results.(fieldname);
%                for iv = ivCol(:)'
%                    colname = sprintf('FE%02d',iv);
%                    valvec = vol_nifti(iv,:);
%                    fname_out = sprintf('%s/FEMA_results_vertexwise_%s_%s_%s.fsvals',dirname_out{1},fstem_imaging,%fieldname,colname);
%                    fs_write_curv(fname_out,valvec,S.nfaces);
%                    fprintf(1,'file %s written\n',fname_out);
%                end
%            end
%            % write out the random effects
%            fieldnamelist = randomFields;
%            for fi = 1:length(fieldnamelist)
%                fieldname = fieldnamelist{fi};
%                vol_nifti = results.(fieldname);
%                for iv = 1:size(vol_nifti,4)
%                    colname = sprintf('RE%02d',iv);
%                    valvec = vol_nifti(iv,:);
%                    fname_out = sprintf('%s/FEMA_results_vertexwise_%s_%s_%s.fsvals',dirname_out{1},fstem_imaging,%fieldname,colname);
%                    fs_write_curv(fname_out,valvec,S.nfaces);
%                    fprintf(1,'file %s written\n',fname_out);
%                end
%            end
%
%        end
%
%    elseif strcmpi(datatype, 'corrmat')
%
%        if 0
%            figure; im = reshape(beta_hat(1,:),dims(2:end)); imagesc(im,max(abs(im(:)))*[-1 1]); colormap%(blueblackred); axis equal tight;
%            figure; im = reshape(beta_hat(2,:),dims(2:end)); imagesc(im,max(abs(im(:)))*[-1 1]); colormap%(blueblackred); axis equal tight;
%            figure; im = reshape(zmat(1,:),dims(2:end)); imagesc(im,max(abs(im(:)))*[-1 1]); colormap(blueblackred); %axis equal tight;
%            figure; im = reshape(zmat(2,:),dims(2:end)); imagesc(im,max(abs(im(:)))*[-1 1]); colormap(blueblackred); %axis equal tight;
%        end
%        %beta_hat = reshape(beta_hat,[size(beta_hat,1) dims(2:end)]);
%        %beta_se = reshape(beta_se,[size(beta_hat,1) dims(2:end)]);
%        %zmat = reshape(zmat,[size(beta_hat,1) dims(2:end)]);
%        %logpmat = reshape(logpmat,[size(beta_hat,1) dims(2:end)]);
%        %sig2tvec = reshape(sig2tvec,[size(sig2tvec,1) dims(2:end)]);
%        %sig2mat = reshape(sig2mat,[size(sig2mat,1) dims(2:end)]);
%
%        if nperms>0
%            save(fpath_out,base_variables_to_save{:},'colnames_imaging','zmat_perm','beta_hat_perm','colnames_interest',%'colsinterest','-v7.3');
%        else
%            save(fpath_out,base_variables_to_save{:},'colnames_imaging','-v7.3');
%        end
%        logging('Results written to %s',fpath_out);
%
%    end
%
%    if ~isempty(fpath_out)
%        fpaths_out = cat(2,fpaths_out,fpath_out);
%    end
end