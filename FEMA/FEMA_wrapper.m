function [fpaths_out beta_hat beta_se zmat logpmat sig2tvec sig2mat beta_hat_perm beta_se_perm zmat_perm sig2tvec_perm sig2mat_perm inputs mask tfce_perm colnames_interest save_params logLikvec Hessmat coeffCovar] = FEMA_wrapper(varargin)
%
% FEMA_wrapper can be called either as:
%   FEMA_wrapper(fstem_imaging, fname_design, dirname_out, dirname_imaging, datatype, ...)
%   FEMA_wrapper(fname_json)
% Wrapper function to run whole FEMA pipeline:
%     1) To load and process imaging data (FEMA_process_data)
%     2) To intersect with design matrices (FEMA_intersect_design)
%     3) To run linear mixed effects models (FEMA_fit)
%     4) To save outputs in specified format
%
% USAGE: FEMA_wrapper(fstem_imaging,fname_design,dirname_out,dirname_imaging,datatype,varargin)
%
% INPUTS
%   fstem_imaging <char>       :  name of vertex/voxel-mapped phenotype (e.g., 'thickness-sm16', 'FA')
%   fname_design <cell>        :  cell array with path to file with design matrix file: 
%                                 this can either be a csv, txt or parquet file of the complete design matrix
%                                 created to the FEMA specifications or it can a JSON config file to be passed to 
%                                 FEMA_makeDesign to create a design matrix. 
%   dirname_out <cell>         :  cell array with path to output directory --> if batching design matrices must include separate output directory as rows within dirname_out as a cell array
%   dirname_imaging <char>     :  path to imaging data directory
%   datatype <char>            :  'voxel','vertex','external', 'corrmat'
%                                   NB: Other than 'external' all code is written to expect ABCD data
%
% Optional input arguments:
%   contrasts <num> OR <path>  :  contrast matrix, or path to file containing contrast matrix (readable by readtable)
%   ico <num>                  :  ico-number for vertexwise analyses (0-based, default 5)
%   transformY <string>        :  transformation to apply to Y (default 'none')
%   output <string>            :  'mat' (default) or 'nifti' or 'deap' or concatenations to write multiple formats.
%   vars_of_interest <string>           :  comma-separated list of IVs to write [this is used only for DEAP]
%   RandomEffects <cell>       :  list of random effects to estimate (default {'F','S','E'}):
%                                   family relatedness (F)
%                                   subject - required for longitudinal analyses (S)
%                                   error (E) - always required
%                                   additive genetic relatedness (A) - must include file path to genetic relatedness data (GRM) for this option
%   GRM_file <char>            :  path to genetic relatedness data (GRM) - default [] - only required if A random effect specified
%   preg_file <char>           :  path to pregnancy data - default [] - only required if T random effect specified
%   address_file <char>        :  path to address data - default [] - only required if H random effect specified
%   nperms <num>               :  default 0 --> if >0 will run and output permuted effects
%   mediation <num>            :  default 0 --> if 1 will ensure same seed used for resampling of the two models used for a mediation analysis
%   niter <num>                :  input for FEMA_fit - default 1
%   nbins <num>                :  input for FEMA_fit - default 20 - number of bins across Y for estimating random effects
%   CovType <char>             :  input for FEMA_fit - default 'analytic' --> no other options currently available
%   FixedEstType <char>        :  input for FEMA_fit - default 'GLS' --> other option: 'OLS'
%   GroupByFamType <boolean>   :  input for FEMA_fit - default true
%   NonnegFlag <blooean>       :  input for FEMA_fit - default true - non-negativity constraint on random effects estimation
%   precision <char>      :  input for FEMA_fit - default 'double' --> other option: 'single' - for precision
%   logLikflag <boolean>       :  input for FEMA_fit - default 0
%   permtype <char>            :  input for FEMA_fit - options:
%                                   'wildbootstrap' - residual boostrap --> creates null distribution by randomly flipping the sign of each observation
%                                   'wildbootstrap-nn' - non-null boostrap --> estimates distribution around effect of interest using sign flipping (used for sobel test)
%   tfce <num>                 :  default 0 --> if 1 will run TFCE

%
% OUTPUTS
%   fpaths_out                 :  results will be saved here in the specified format
%   All other outputs are optional if user wants results to be output in the MATLAB workspace
%   All permutation outputs will be empty is nperms==0
%


% This software is Copyright (c) 2021 The Regents of the University of California. All Rights Reserved.
% See LICENSE.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% TO DOs
% - FEMA_get_defaults 
% - contrasts
% - add study and release variable for non-DEAP mode. is there a way to get this from the data rather than user input?
% - DEAP mode: do we need to create the makeDEsign config file? 
% - DEAP mode: how do we set dfirname_out? 
%
%
%
%
%
%
%
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% logging('***Start FEMA v2.3, 11/15/2023)');
logging('***Start***');
disp(FEMA_info);
starttime = now();
rng shuffle %Set random number generator so different every time

%PARSING INPUTS
% check for required params
requiredSet1 = {'config', 'data', 'output'};
requiredSet2 = {{'fstem_imaging', 'fname_design', 'dirname_out', 'dirname_imaging', 'datatype'}, ...
                {'fstem_imaging', 'config_design', 'dirname_out', 'dirname_imaging', 'datatype'}};

existSet1 = all(ismember(requiredSet1, varargin));
existSet2 = any(cellfun(@(x) all(ismember(x, varargin)), requiredSet2));

if ~existSet1 & ~existSet2
    error('Missing required parameters.\nMust provide %s or %s.', ...
    strjoin(requiredSet1, ', '), strjoin(requiredSet2, ', '));
end

if existSet1
    % check that it's a json file and it exists
    [fstem_imaging, config_design, dirname_out, dirname_imaging, datatype, dataFile, extraArgs] = ...
        FEMA_parseInputs(varargin{:});
    varargin = extraArgs;
end

% if isdeployed
%   logging('***FEMA_wrapper_app Compiled 11/8/2021, FEMA v2.0***\n\n'); %TODO: remember to change this before compiling
% end


inputs = inputParser;
% required inputs
addRequired(inputs, 'fstem_imaging', @ischar);
addRequired(inputs, 'dirname_out', @ischar);
addRequired(inputs, 'dirname_imaging', @ischar);
addRequired(inputs, 'datatype', @(x) ischar(x) && ...
                     ismember(x, {'voxel', 'vertex', 'corrmat', ...
                                  'roi', 'external'}));
% optional inputs
addParameter(inputs, 'fname_design', @(x) iscell(x) || ischar(x));
addParameter(inputs, 'config_design', @(x) iscell(x) || ischar(x));
addParameter(inputs, 'iid_filter', [], @(x) iscell(x) || ischar(x));
addParameter(inputs, 'eid_filter', [], @(x) iscell(x) || ischar(x));
addParameter(inputs, 'fname_qc', [], @(x) ischar(x));
addParameter(inputs, 'qc_var', [], @(x) ischar(x));
addParameter(inputs, 'transformY', 'none', ...
                      @(x) ischar(x) && ...
                      ismember(x, {'center' 'centre' 'demean' ...
                                   'std' 'standardize' 'normalize' ...
                                   'logn' 'log10' 'inverseranknorm' ...
                                   'ranknorm' 'int'}));
addParameter(inputs, 'study', 'abcd', @(x) ischar(x) && ismember(x, {'abcd', 'hbcd'}));
addParameter(inputs, 'release', '6.0', @(x) ischar(x));
addParameter(inputs, 'outPrefix', ['FEMA_', char(datetime('now', 'Format', 'yyyyMMdd_HHmmSS'))], @(x) ischar(x));
addParameter(inputs, 'saveDesignMatrix', true, @(x) isscalar(x) && islogical(x) || isnumeric(x));
addParameter(inputs, 'returnResiduals', false, @(x) islogical(x) || isnumeric(x));


addParameter(inputs, 'ico', 5, @(x) isscalar(x) && isnumeric(x));
addParameter(inputs, 'fname_contrast', [], @(x) ischar(x));
addParameter(inputs, 'outputType', '.mat', @(x) ischar(x) && ...
                     ismember(x, {'.mat' '.nii' '.nii.gz' ... 
                                  'nifti' 'voxel' '.gii' ...
                                  'gifti' 'vertex' 'tables' ...
                                  'summary'}));
addParameter(inputs, 'synth',  false, @(x) isscalar(x) && islogical(x) || isnumeric(x)); 
addParameter(inputs, 'vars_of_interest', '', @(x) ischar(x) || iscell(x));
addParameter(inputs, 'RandomEffects', {'F' 'S' 'E'}, @(x) iscell(x) || ischar(x));
addParameter(inputs, 'GRM_file', [], @(x) ischar(x));
addParameter(inputs, 'preg_file', [], @(x) ischar(x));
addParameter(inputs, 'address_file', [], @(x) ischar(x));
addParameter(inputs, 'nperms', 0, @(x) isscalar(x) && isnumeric(x));
addParameter(inputs, 'mediation', false, @(x) isscalar(x) && islogical(x) || isnumeric(x));
addParameter(inputs, 'tfce', false, @(x) isscalar(x) && islogical(x) || isnumeric(x));
addParameter(inputs, 'corrvec_thresh', 0.8, @(x) isscalar(x) && isnumeric(x));
%FEMA_fit variable inputs
addParameter(inputs, 'niter', 1, @(x) isscalar(x) && isnumeric(x));
addParameter(inputs, 'nbins', 20, @(x) isscalar(x) && isnumeric(x));
addParameter(inputs, 'CovType', 'analytic', @(x) ischar(x) && ...
                      ismember(x, {'analytic', 'unstructured'}));
addParameter(inputs, 'FixedEstType', 'GLS', @(x) ischar(x) && ismember(x, {'OLS', 'GLS'}));
addParameter(inputs, 'RandomEstType', 'MoM', @(x) ischar(x) && ismember(x, {'MoM', 'ML'}));
addParameter(inputs, 'GroupByFamType', true, @(x) isscalar(x) && islogical(x) || isnumeric(x));
addParameter(inputs, 'NonnegFlag', true, @(x) isscalar(x) && islogical(x) || isnumeric(x)); % Perform lsqnonneg on random effects estimation
addParameter(inputs, 'precision', 'double', @(x) ischar(x) && ...
                      ismember(x, {'single', 'double'}));
addParameter(inputs, 'logLikflag', false, @(x) isscalar(x) && islogical(x) || isnumeric(x));
addParameter(inputs, 'Hessflag', false, @(x) isscalar(x) && islogical(x) || isnumeric(x));
addParameter(inputs, 'ciflag', false, @(x) isscalar(x) && islogical(x) || isnumeric(x));
addParameter(inputs, 'permtype', 'wildbootstrap', @(x) ischar(x) && ...
                      ismember(x, {'wildbootstrap', 'wildbootstrap-nn', 'none'}));


parse(inputs, fstem_imaging, dirname_out, dirname_imaging, datatype, varargin{:})

% Display input arguments for log
%disp(inputs.Results)

niter = inputs.Results.niter;
ico = inputs.Results.ico;
fname_contrast = inputs.Results.fname_contrast;
transformY = inputs.Results.transformY;
outputType = inputs.Results.outputType;
if strcmp(datatype,'external') 
        outputType = 'tables'; 
end 
nbins = inputs.Results.nbins;
if ~isempty(inputs.Results.vars_of_interest)
    vars_of_interest = split(inputs.Results.vars_of_interest,','); 
else
    vars_of_interest = {};
end

if ~isempty(vars_of_interest) || isdeployed
    logging('%d Variables of interest specified: %s',length(vars_of_interest), strjoin(vars_of_interest, ', '));
end
outPrefix = inputs.Results.outPrefix;
saveDesignMatrix = inputs.Results.saveDesignMatrix;
iid_filter = inputs.Results.iid_filter;
eid_filter = inputs.Results.eid_filter;
RandomEffects = inputs.Results.RandomEffects;
fname_GRM = inputs.Results.GRM_file;
fname_address = inputs.Results.address_file;
fname_pregnancy = inputs.Results.preg_file;
CovType = inputs.Results.CovType;
FixedEstType = inputs.Results.FixedEstType;
RandomEstType = inputs.Results.RandomEstType;
GroupByFamType = inputs.Results.GroupByFamType;
NonnegFlag = inputs.Results.NonnegFlag;
precision = inputs.Results.precision;
logLikflag = inputs.Results.logLikflag;
Hessflag = inputs.Results.Hessflag;
ciflag = inputs.Results.ciflag;
nperms = inputs.Results.nperms;
permtype = inputs.Results.permtype;
mediation = inputs.Results.mediation;
synth = inputs.Results.synth;
tfce = inputs.Results.tfce;
corrvec_thresh = inputs.Results.corrvec_thresh;

if ~exist('config_design', 'var') & ~ exist('fname_design', 'var')
    error('Either config_design or fname_design must be specified.');
end

if exist('fname_design', 'var') && ~isempty(fname_design)
    if ~iscell(fname_design)
        fname_design = {fname_design};
    end 
    designExists = true;
    n_desmat = length(fname_design);
end

if exist('config_design', 'var') && ~isempty(config_design)
    if ~iscell(config_design)
        config_design = {config_design};
    end
    designExists = false;
    n_desmat = length(config_design);
end

if ~iscell(dirname_out)
    dirname_out = {dirname_out};
end

if n_desmat ~= length(dirname_out) 
    error('if using cell array specifying multiple designs, fname_design and dirname_out must both have an equal number of items.')
end

fprintf('Random effects: %s\n', strjoin(RandomEffects, ', '));

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% LOAD AND PROCESS IMAGING DATA FOR ANALYSIS - ABCD specific function unless datatype='external'
[ymat, iid_concat, eid_concat, ivec_mask, mask, colnames_imaging, GRM, preg, address, missing_process_data] = ...
    FEMA_process_data(fstem_imaging, dirname_imaging, datatype, 'ico', ico, ...
                      'GRM_file', fname_GRM, 'preg_file', fname_pregnancy, ...
                      'address_file', fname_address, 'corrvec_thresh', corrvec_thresh, ...
                      'iid', iid_filter, 'eid', eid_filter);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% back up data 
GRM_bak=GRM;
ymat_bak=ymat;
            
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% IF RUNNING MEDIATION ANALYSIS
if mediation==1

    seed=rng; % Save rng in order to ensure resampling is the same for both models

    if ~ismember(lower(permtype),{'wildbootstrap-nn'})
        error('Change permtype input to wildboostrap-nn to perform non-null WB needed for mediation analysis OR set mediation=0.')
    end

    if length(fname_design)~=2
        error('fname_design must be a cell array containing 2 nested design matrices to test for mediation.')
    end

    if nperms==0
        error('Change nperms>0 to run mediation analysis OR set mediation=0.')
    end

    if tfce==1
        error('Cannot run mediation and TFCE as different resampling scheme required. Set tfce=0.')
    end

end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% CREATE DESIGN MATRIX and  INTERSECT YMAT
%Loops over multiple design matrices (rows in fname_design cell array) to run several models with the same imaging data
fpaths_out = {};

for des=1:n_desmat
    % read in or create design matrix
    if designExists
        designMatrix = readtable(fname_design{des});
    else 
        logging('Creating design matrix.');
        [designMatrix, vars_of_interest] = FEMA_makeDesign(config_design{des},'dataFile', dataFile, ...
                                       'iid', iid_concat, 'eid', eid_concat); 
    end
    % get column names
    exp_colnames = {'iid' 'eid' 'fid' 'agevec'};
    exp_colnames_exist = all(ismember(exp_colnames, designMatrix.Properties.VariableNames));
    if exp_colnames_exist
        colnames_model = designMatrix.Properties.VariableNames(~ismember(designMatrix.Properties.VariableNames, exp_colnames));
    else 
        warning('Assuming first 4 columns of design matrix are participant ID, session ID, family ID and age, in that order.')
        colnames_model = designMatrix.Properties.VariableNames(5:end);
    end
    if ~isempty(vars_of_interest)
        colsinterest = find(ismember(colnames_model, vars_of_interest)); 
    end 

    % contrasts 
    if ~isempty(fname_contrast)
        [contrasts, hypValues, contrast_names] = FEMA_parse_contrastFile(fname_contrast, colnames_model); 
        colnames_model = cat(2, contrast_names, colnames_model);
    else
        contrasts = [];
        hypValues = [];
    end

    % intersect ymat with design matrix
    [X, iid, eid, fid, agevec, ymat, GRM, PregID, HomeID, missing_intersect_design] = ...
        FEMA_intersect_design(designMatrix, ymat_bak, iid_concat, eid_concat, ...
                             'GRM', GRM_bak, 'preg', preg, 'address', address);

    if synth==1 % Make synthesized data
        [ymat sig2tvec_true sig2mat_true] = FEMA_synthesize(X, iid, eid, fid, agevec, ymat, GRM, 'nbins', nbins, 'RandomEffects', RandomEffects); % Make GRM and zygmat optional arguments? % Need to update SSE_synthesize_dev to accept list of random effects to include, and range of values

        % sig2mat_true(length(RandomEffects),:) = 1-sum(sig2mat_true(1:length(RandomEffects)-1,:),1); % This shouldn't be needed, if RandomEffects include 'E'
    else
        sig2tvec_true = []; sig2mat_true = [];
    end
    synthstruct = struct('sig2tvec_true',sig2tvec_true,'sig2mat_true',sig2mat_true);

    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % transform ymat
    % always winsorize tfmri data
    if contains(fstem_imaging, 'beta', 'IgnoreCase', true) 
        [ymat settingsTransform] = doTransformation(ymat, 'winsorize', ...
                                                    'lower_bound', 2, 'upper_bound', 98);
    end 
    % transform according to user input 
    if ~strcmpi(transformY, 'none')
        [ymat settingsTransform] = doTransformation(ymat, transformY);
    end

    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % IF RUNNING MODELS FOR MEDIATION ANALYSIS
    if mediation==1

        rng(seed); %set same seed for both runs of FEMA_fit

        if des==1
            X_bak=X;
        end

        if des==2
            if size(X_bak,1)~=size(X,1)
                error('Design matrices must have the same number of observations for mediation analysis.') % Check after intersection
            end
        end

    end

    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % FIT MODEL
    [beta_hat, beta_se, zmat, logpmat, sig2tvec, sig2mat, Hessmat, logLikvec,                            ...
     beta_hat_perm, beta_se_perm, zmat_perm, sig2tvec_perm, sig2mat_perm, logLikvec_perm,                ...
     binvec_save, nvec_bins, tvec_bins, FamilyStruct, coeffCovar, unstructParams, residuals_GLS, info] = ...
     FEMA_fit(X, iid, eid, fid, agevec, ymat, contrasts, nbins, GRM.GRM,                                 ...
             'niter', niter, 'RandomEffects', RandomEffects, 'nperms', nperms, 'CovType', CovType,       ...
             'FixedEstType', FixedEstType, 'RandomEstType', RandomEstType, ...
             'GroupByFamType', GroupByFamType, 'NonnegFlag', NonnegFlag, ...
             'SingleOrDouble', precision, 'logLikflag', logLikflag, ...
             'Hessflag', Hessflag, 'ciflag', ciflag, 'permtype', permtype, ...
             'PregID', PregID, 'HomeID', HomeID, 'synthstruct', synthstruct, 'returnResiduals', returnResiduals);

    if sum(~mask)>0

        z_tmp=zeros(size(zmat,1),size(mask,2));
        p_tmp=zeros(size(logpmat,1),size(mask,2));
        beta_tmp=zeros(size(beta_hat,1),size(mask,2));
        betase_tmp=zeros(size(beta_se,1),size(mask,2));
        sig2mat_tmp=zeros(size(sig2mat,1),size(mask,2));
        sig2tvec_tmp=zeros(size(sig2tvec,1),size(mask,2));
        coeffCovar_tmp=zeros(size(coeffCovar, 1), size(coeffCovar, 2), size(mask, 2));

        z_tmp(:,ivec_mask)=zmat;
        p_tmp(:,ivec_mask)=logpmat;
        beta_tmp(:,ivec_mask)=beta_hat;
        betase_tmp(:,ivec_mask)=beta_se;
        sig2mat_tmp(:,ivec_mask)=sig2mat;
        sig2tvec_tmp(:,ivec_mask)=sig2tvec;
        coeffCovar_tmp(:, :, ivec_mask)=coeffCovar;

        zmat=z_tmp;
        logpmat=p_tmp;
        beta_hat=beta_tmp;
        beta_se=betase_tmp;
        sig2mat=sig2mat_tmp;
        sig2tvec=sig2tvec_tmp;
        coeffCovar=coeffCovar_tmp;

        if nperms>0

            zperm_tmp=zeros(size(zmat_perm,1),size(mask,2),size(zmat_perm,3));
            betaperm_tmp=zeros(size(beta_hat_perm,1),size(mask,2),size(beta_hat_perm,3));
            betaseperm_tmp=zeros(size(beta_se_perm,1),size(mask,2),size(beta_se_perm,3));
            sig2matperm_tmp=zeros(size(sig2mat_perm,1),size(mask,2),size(sig2mat_perm,3));
            sig2tvecperm_tmp=zeros(size(sig2tvec_perm,1),size(mask,2),size(sig2tvec_perm,3));
            coeffCovar_perm_tmp=zeros(size(coeffCovar, 1), size(coeffCovar, 2), size(mask, 2), size(coeffCovar, 4));

            zperm_tmp(:,ivec_mask,:)=zmat_perm;
            betaperm_tmp(:,ivec_mask,:)=beta_hat_perm;
            betaseperm_tmp(:,ivec_mask,:)=beta_se_perm;
            sig2matperm_tmp(:,ivec_mask,:)=sig2mat_perm;
            sig2tvecperm_tmp(:,ivec_mask,:)=sig2tvec_perm;
            coeffCovar_perm_tmp(:, :, ivec_mask, :)=coeffCovar_perm;

            zmat_perm=zperm_tmp;
            beta_hat_perm=betaperm_tmp;
            beta_se_perm=betaseperm_tmp;
            sig2mat_perm=sig2matperm_tmp;
            sig2tvec_perm=sig2tvecperm_tmp;
            coeffCovar_perm=coeffCovar_perm_tmp;

        end
    end

    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    % IF RUNNING TFCE

    if tfce==1 && nperms>0
        tfce_perm = FEMA_tfce(zmat_perm,colsinterest,datatype,mask);
    elseif tfce==1 && nperms==0
        error('Change nperms>0 to run TFCE OR set tfce=0.')
    else
        tfce_perm=[];
    end

    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    % IF RUNNING RESAMPLING

    % Only save permutations for IVs of interest to reduce size of output
    % If wanting to save permutations for all columns of X make colsinterest=[1:size(X,2)]

    if nperms>0
        beta_hat_perm=beta_hat_perm(colsinterest,:,:);
        beta_se_perm=beta_se_perm(colsinterest,:,:);
        zmat_perm=zmat_perm(colsinterest,:,:);
        colnames_interest=colnames_model(colsinterest);
    elseif nperms==0
        colnames_interest=colnames_model;
    end

    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % SAVE OUTPUT

    % always save out a .mat file 
    outputType = '.mat';
    FEMA_save(outputType, dirname_out{des}, outPrefix, ...
              'beta_hat', beta_hat, 'beta_se', beta_se, 'zmat', zmat, 'logpmat', logpmat, ...
              'sig2tvec', sig2tvec, 'sig2mat', sig2mat, ... % 'sig2mat_normalized', sig2mat_normalized, 
              'beta_hat_perm', beta_hat_perm, 'beta_se_perm', beta_se_perm, 'zmat_perm', zmat_perm, ...
              'sig2tvec_perm', sig2tvec_perm, 'sig2mat_perm', sig2mat_perm, 'logLikvec', logLikvec, ...
              'logLikvec_perm', logLikvec_perm, 'Hessmat', Hessmat, 'coeffCovar', coeffCovar, ...
              'binvec_save', binvec_save, 'nvec_bins', nvec_bins, 'tvec_bins', tvec_bins, ...
              'reusableVars', reusableVars, 'contrasts', contrasts, 'hypValues', hypValues, ...
              'designMatrix', designMatrix, ... % toSave_design
               'niter', niter, 'nbins', nbins, ... % start of toSave_info
              'RandomEffects', RandomEffects, 'nperms', nperms, 'CovType', CovType, ...
              'FixedEstType', FixedEstType, 'RandomEstType', RandomEstType, 'GroupByFamType', GroupByFamType, ...
              'NonnegFlag', NonnegFlag, 'precision', precision, 'logLikflag', logLikflag, ...
              'permtype', permtype); 
              %'FamilyStruct', FamilyStruct, 'MotherID', MotherID, 'FatherID', ...
              %FatherID, 'HomeID', HomeID, 'PregID', PregID, ... \
              %'numWorkers', numWorkers, 'numThreads', numThreads

    % different save options depending on datatype
    % NB corrmat is always saved as .mat only -> do we need another file type for DEAP mode? 
    switch datatype
        case {'voxel', 'vertex', 'roi'}
            % need to add how to handle roi 
            % save nifti/gifti of fixed effects variables of itnerest
            outPrefix = 'FFX'; 
            FEMA_save(datatype, dirname_out{des}, outPrefix, 'beta_hat', beta_hat, 'beta_se', beta_se, ...
                      'zmat', zmat,'logpmat', logpmat, 'colsinterest', colsinterest, ...
                      'vars_of_interest', vars_of_interest); 
            % save nifti/gifti of random effects 
            outPrefix = 'RFX'; 
            if strcmp(CovType, 'analytic')
                FEMA_save(datatype, dirname_out{des}, outPrefix, ...
                          'sig2mat', sig2mat, 'sig2tvec', sig2tvec, 'RandomEffects', RandomEffects);      
            else
                FEMA_save(datatype, dirname_out{des}, outPrefix, ...
                          'sig2mat_normalized', unstructuredParams.sig2mat_normalized, ...
                          'sig2tvec', sig2tvec, 'RandomEffects', RandomEffects, ...
                          'eidOrd', unstructuredParams.eidOrd);
            end
        case 'external'
            outPrefix = 'regression_tables';
            outputType = 'tables'; 
            FEMA_save(outputType, dirname_out{des}, outPrefix, 'beta_hat', beta_hat, 'beta_se', beta_se, 'zmat', zmat,'logpmat', logpmat,'sig2tvec',sig2tvec,'sig2mat',sig2mat, 'sig2mat_normalized', sig2mat_normalized, 'info', FEMA_wrapper_info, 'colnames_model', colnames_model, 'ymat_names', fstem_imaging) % check fstem_imaging makes sense for tabulated 
    end 
    
end  % LOOP over design matrices

endtime = now();
logging('***Done*** (%0.2f seconds)',(endtime-starttime)*3600*24);

end 

