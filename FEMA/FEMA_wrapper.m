function [beta_hat, beta_se, zmat, logpmat, sig2tvec, sig2mat, Hessmat, logLikvec,              ...
          beta_hat_perm, beta_se_perm, zmat_perm, sig2tvec_perm, sig2mat_perm, logLikvec_perm,  ...
          binvec_save, nvec_bins, tvec_bins, FamilyStruct, coeffCovar, unstructParams,          ...
          residuals_GLS, info_fit, inputs, mask, tfce_perm, info] = FEMA_wrapper(varargin)
%
% FEMA_wrapper can be called either as:
%   FEMA_wrapper(fstem_imaging, fname_design, dirname_out, dirname_imaging, datatype, ...)
%   FEMA_wrapper('config', fname_json, 'data', fname_data, 'output', dirname_out)
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
%   dirname_out <cell>         :  cell array with path to output directory --> if batching design matrices must include 
%                                 separate output directory as rows within dirname_out as a cell array
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
%   precision <char>           :  input for FEMA_fit - default 'double' --> other option: 'single' - for precision
%   logLikflag <boolean>       :  input for FEMA_fit - default 0
%   PermType <char>            :  input for FEMA_fit - options:
%                                   'wildbootstrap' - residual boostrap --> creates null distribution by randomly 
%                                                                           flipping the sign of each observation
%                                   'wildbootstrap-nn' - non-null boostrap --> estimates distribution around effect 
%                                                                              of interest using sign flipping (used for sobel test)
%   tfce <num>                 :  default 0 --> if 1 will run TFCE

%
% OUTPUTS
%   All  outputs are optional if user wants results to be output in the MATLAB workspace
%   All permutation outputs will be empty is nperms==0
%
% See LICENSE.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% TO DOs
% - FEMA_get_defaults 
% - contrasts
% - add study and release variable for non-DEAP mode. is there a way to get this from the data rather than user input?
% - DEAP mode: do we need to create the makeDesign config file? 
% - DEAP mode: how do we set dirname_out? 
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

logging('***Start***');
disp(FEMA_info);
tStart = tic;
rng shuffle %Set random number generator so different every time

%PARSING INPUTS
% check for required params
requiredSet1 = {'config', 'data', 'output'};
requiredSet2 = {{'fstem_imaging', 'fname_design',  'dirname_out', 'dirname_imaging', 'datatype'}, ...
                {'fstem_imaging', 'config_design', 'dirname_out', 'dirname_imaging', 'datatype'}};

% if any(cellfun(@(x) iscell(x), varargin))
%     existSet1 = all(ismember(requiredSet1, horzcat(varargin{:})));
%     existSet2 = any(cellfun(@(x) all(ismember(x, horzcat(varargin{:}))), requiredSet2));
% else
    % existSet1 = all(cellfun(@(x) any(ismember(requiredSet1, x)), varargin, 'UniformOutput', true));
    existSet1 = all(ismember(requiredSet1, varargin(1:2:end)));
    existSet2 = any(cellfun(@(x) all(ismember(x, varargin(1:2:end))), requiredSet2));
% end

if ~existSet1 && ~existSet2
    error('Missing required parameters.\nMust provide one of \n%s \nOR\n%s\nOR\n%s', ...
    strjoin(requiredSet1, ', '), strjoin(requiredSet2{1}, ', '), strjoin(requiredSet2{2}, ', '));
end

if existSet1
    % check that it's a json file and it exists
    [fstem_imaging, config_design, dirname_out, dirname_imaging, datatype, dataFile, extraArgs] = ...
     FEMA_parseInputs(varargin{:});

    % Any additional arguments?
    leftOverArgs = setdiff(setdiff(varargin(1:2:end), extraArgs(1:2:end)), {'config', 'data', 'output'});
    for ii = 1:length(leftOverArgs)
        loc = ismember(varargin(1:2:end), leftOverArgs{ii});
        extraArgs{end+1} = leftOverArgs{ii}; %#ok<AGROW>
        extraArgs{end+1} = varargin{find(loc)*2}; %#ok<AGROW>
        % extraArgs{end+1} = varargin{find(ismember(varargin, leftOverArgs{ii}))+1}; %#ok<AGROW>
    end
    varargin = extraArgs;
else
    fstem_imaging   = varargin{find(strcmpi(varargin, {'fstem_imaging'})) + 1};
    dirname_out     = varargin{find(strcmpi(varargin, {'dirname_out'})) + 1};
    dirname_imaging = varargin{find(strcmpi(varargin, {'dirname_imaging'})) + 1};
    datatype        = varargin{find(strcmpi(varargin, {'datatype'})) + 1};
end

inputs = inputParser;

% required inputs
addRequired(inputs, 'fstem_imaging', @(x) iscell(x) || ischar(x));
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
                      ismember(x, {'none' 'center' 'centre' 'demean' ...
                                   'std' 'standardize' 'normalize' ...
                                   'logn' 'log10' 'inverseranknorm' ...
                                   'ranknorm' 'int'}));
addParameter(inputs, 'study', 'abcd', @(x) ischar(x) && ismember(x, {'abcd', 'hbcd'}));
addParameter(inputs, 'release', '6.0', @(x) ischar(x));
addParameter(inputs, 'outPrefix', [], @(x) ischar(x));
addParameter(inputs, 'saveDesignMatrix', true, @(x) isscalar(x) && islogical(x) || isnumeric(x));
addParameter(inputs, 'returnResiduals', false, @(x) islogical(x) || isnumeric(x));
addParameter(inputs, 'matType', 'uncompressed', @(x) ischar(x) && ismember(lower(x), {'compressed', 'uncompressed'}));

addParameter(inputs, 'ico', 5, @(x) isscalar(x) && isnumeric(x));
addParameter(inputs, 'fname_contrast', [], @(x) ischar(x));
addParameter(inputs, 'outputType', [], @(x) (iscell(x) || ischar(x)) && ...
                     all(ismember(x, {'mat' 'nii' 'nii.gz' ... 
                                  'nifti' 'voxel' 'gii' ...
                                  'gifti' 'vertex' 'tables' ...
                                  'corrmat' 'external' ...
                                  'summary' 'none' 'cache'})));
addParameter(inputs, 'synth',  false, @(x) isscalar(x) && islogical(x) || isnumeric(x)); 
addParameter(inputs, 'vars_of_interest', '', @(x) ischar(x) || iscell(x));
addParameter(inputs, 'RandomEffects', {'F' 'S' 'E'}, @(x) iscell(x) || ischar(x));
addParameter(inputs, 'GRM_file', [], @(x) ischar(x) || isempty(x));
addParameter(inputs, 'preg_file', [], @(x) ischar(x));
addParameter(inputs, 'address_file', [], @(x) ischar(x));
addParameter(inputs, 'nperms', 0, @(x) isscalar(x) && isnumeric(x));
addParameter(inputs, 'mediation', false, @(x) isscalar(x) && islogical(x) || isnumeric(x));
addParameter(inputs, 'tfce', false, @(x) isscalar(x) && islogical(x) || isnumeric(x));
addParameter(inputs, 'corrvec_thresh', 0.8, @(x) isscalar(x) && isnumeric(x));

% FEMA_fit variable inputs
addParameter(inputs, 'niter', 1, @(x) isscalar(x) && isnumeric(x));
addParameter(inputs, 'nbins', 20, @(x) isscalar(x) && isnumeric(x));
addParameter(inputs, 'CovType', 'analytic', @(x) ischar(x) && ...
                      ismember(x, {'analytical', 'analytic', 'unstructured'}));
addParameter(inputs, 'FixedEstType', 'GLS', @(x) ischar(x) && ismember(x, {'OLS', 'GLS'}));
addParameter(inputs, 'RandomEstType', 'MoM', @(x) ischar(x) && ismember(x, {'MoM', 'ML'}));
addParameter(inputs, 'GroupByFamType', true, @(x) isscalar(x) && islogical(x) || isnumeric(x));
addParameter(inputs, 'NonnegFlag', true, @(x) isscalar(x) && islogical(x) || isnumeric(x)); % Perform lsqnonneg on random effects estimation
addParameter(inputs, 'precision', 'double', @(x) ischar(x) && ...
                      ismember(x, {'single', 'double'}));
addParameter(inputs, 'logLikflag', false, @(x) isscalar(x) && islogical(x) || isnumeric(x));
addParameter(inputs, 'Hessflag', false, @(x) isscalar(x) && islogical(x) || isnumeric(x));
addParameter(inputs, 'ciflag', false, @(x) isscalar(x) && islogical(x) || isnumeric(x));
addParameter(inputs, 'PermType', 'wildbootstrap', @(x) ischar(x) && ...
                      ismember(x, {'wildbootstrap', 'wildbootstrap-nn', 'none'}));
addParameter(inputs, 'doPar', false, @islogical);
addParameter(inputs, 'numWorkers', 2, @isnumeric);

parse(inputs, fstem_imaging, dirname_out, dirname_imaging, datatype, varargin{:})

% Display input arguments for log
%disp(inputs.Results)

niter = inputs.Results.niter;
ico = inputs.Results.ico;
fname_contrast = inputs.Results.fname_contrast;
transformY = inputs.Results.transformY;
outputType = inputs.Results.outputType;
nbins = inputs.Results.nbins;
if ~isempty(inputs.Results.vars_of_interest)
    vars_of_interest = split(inputs.Results.vars_of_interest,','); 
else
    vars_of_interest = {};
end
study = inputs.Results.study;
release = inputs.Results.release;
returnResiduals = inputs.Results.returnResiduals;
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
PermType = inputs.Results.PermType;
mediation = inputs.Results.mediation;
synth = inputs.Results.synth;
tfce = inputs.Results.tfce;
corrvec_thresh = inputs.Results.corrvec_thresh;
fname_design = inputs.Results.fname_design;
matType = inputs.Results.matType;
doPar = inputs.Results.doPar;
numWorkers = inputs.Results.numWorkers;

if ~exist('config_design', 'var') && ~ exist('fname_design', 'var')
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
    error(['If using cell array specifying multiple design matrices, ', ...
           'fname_design and dirname_out must both have an equal number of items.']);
end

% if hbcd and deap mode 
if strcmpi(study, 'hbcd') && existSet1
    fname_fam = '/data/hbcd/support_files/hbcd_family_id.parquet';
else 
    fname_fam = [];
end 


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
[ymat, iid_concat, eid_concat, ivec_mask, mask, GRM, preg, address, info_process_data] = ...
FEMA_process_data(fstem_imaging, dirname_imaging, datatype, 'ico', ico, 'GRM_file', fname_GRM, 'preg_file', ...
                  fname_pregnancy, 'address_file', fname_address, 'corrvec_thresh', corrvec_thresh, ...
                  'iid', iid_filter, 'eid', eid_filter, 'study', study, 'release', release);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% back up data 
GRM_bak=GRM;
ymat_bak=ymat;
            
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% check permutations if mediation or tfce is requested
if mediation==1

    seed=rng; % Save rng in order to ensure resampling is the same for both models

    if ~ismember(lower(PermType),{'wildbootstrap-nn'})
        error('Change PermType input to wildboostrap-nn to perform non-null WB needed for mediation analysis OR set mediation=0.')
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
% create or load design matrix and intersect with imaging data


for des=1:n_desmat % loop if multiple design matrices
    % read in or create design matrix
    if designExists
        tRead = tic;
        designMatrix = readtable(fname_design{des});
        tmp = struct('settings', [], 'timing', []);
        tmp.settings.fname = fname_design{des};
        tmp.timing.tRead   = toc(tRead);
        info_makeDesign    = tmp;
    else 
        disp(['In Wrapper: about to call FEMA_makeDesign; full path to dataFile is: ', dataFile]);
        [designMatrix, vars_of_interest, splines_of_interest, FFX_conceptMapping, info_makeDesign] = ...
         FEMA_makeDesign(config_design{des},'dataFile', dataFile, ...
                         'iid', iid_concat, 'eid', eid_concat,    ...
                         'fname_fam', fname_fam, 'study', study);
    end
    if ~isempty(vars_of_interest) 
        logging('%d Variables of interest specified: %s.',length(vars_of_interest), strjoin(vars_of_interest, ', '));
    else 
        logging('No variables of interest specified.');
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
        [~, colsinterest] = ismember(vars_of_interest, colnames_model);
    else
        colsinterest = [];
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
    [X, iid, eid, fid, agevec, ymat, GRM, PregID, HomeID, info_intersect_design] = ...
        FEMA_intersect_design(designMatrix, ymat_bak, iid_concat, eid_concat, ...
                              'GRM', GRM_bak, 'preg', preg, 'address', address);

    if synth==1 % Make synthesized data
        % Make GRM and zygmat optional arguments? 
        % Need to update SSE_synthesize_dev to accept list of random
        % effects to include, and range of values
        [ymat, sig2tvec_true, sig2mat_true] = ...
         FEMA_synthesize(X, iid, eid, fid, agevec, ymat, GRM, ...
                        'nbins', nbins, 'RandomEffects', RandomEffects);

        % This shouldn't be needed, if RandomEffects include 'E'
        % sig2mat_true(length(RandomEffects),:) = 1-sum(sig2mat_true(1:length(RandomEffects)-1,:),1); 
    else
        sig2tvec_true = []; 
        sig2mat_true = [];
    end
    synthstruct = struct('sig2tvec_true',sig2tvec_true,'sig2mat_true',sig2mat_true);

    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % transform ymat according to user input 
    if ~strcmpi(transformY, 'none')
        logging('Applying %s to Y matrix.', transformY);
        [ymat, info_process_data.transformY] = doTransformation(ymat, transformY);
    end

    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % mediation 
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
    % fit model 
    [beta_hat, beta_se, zmat, logpmat, sig2tvec, sig2mat, Hessmat, logLikvec,                            ...
     beta_hat_perm, beta_se_perm, zmat_perm, sig2tvec_perm, sig2mat_perm, logLikvec_perm,                ...
     binvec_save, nvec_bins, tvec_bins, FamilyStruct, coeffCovar, unstructParams, residuals_GLS, info_fit] = ...
     FEMA_fit(X, iid, eid, fid, agevec, ymat, contrasts, nbins, GRM.GRM,                                 ...
             'niter', niter, 'RandomEffects', RandomEffects, 'nperms', nperms, 'CovType', CovType,       ...
             'FixedEstType', FixedEstType, 'RandomEstType', RandomEstType, ...
             'GroupByFamType', GroupByFamType, 'NonnegFlag', NonnegFlag, ...
             'precision', precision, 'logLikflag', logLikflag, ...
             'Hessflag', Hessflag, 'ciflag', ciflag, 'PermType', PermType, ...
             'PregID', PregID, 'HomeID', HomeID, 'synthstruct', synthstruct, ...
             'returnResiduals', returnResiduals, 'doPar', doPar, 'numWorkers', numWorkers);

    %% Do Wald test on splines_of_interest
    % Initialize results from Wald test
    info_FEMA_WaldTest = struct('name', {{}}, 'L', {{}}, 'hypValue', {{}}, 'doF', {{}}, 'numObs', {{}}, 'timing', {{}});
    if ~isempty(splines_of_interest)
        [Wald, logp_Wald] = deal(zeros(size(splines_of_interest,1), size(beta_hat,2)));
        
        % ToDo: hypValue needs to be picked from contrast file
        %       handle Wald test for permutations
        hypValue      = 0;
        doF           = false;
        numObs        = NaN;
        tInit_overall = tic;

        for ii = 1:size(splines_of_interest,1)

            tInit = tic;

            % Which variables to work with?
            [~, toWork] = ismember(splines_of_interest{ii,1}, colnames_model);

            % Contrast matrix
            L = zeros(size(beta_hat, 1));
            L(toWork, toWork) = eye(length(toWork));
            
            % Do Wald test
            [Wald(ii,:), logp_Wald(ii,:)] = ...
            FEMA_WaldTest(L, beta_hat, coeffCovar, hypValue, doF, numObs);

            % Record some information
            info_FEMA_WaldTest.name{ii}     = splines_of_interest(ii,2);
            info_FEMA_WaldTest.L{ii}        = L;
            info_FEMA_WaldTest.hypValue{ii} = hypValue;
            info_FEMA_WaldTest.doF{ii}      = doF;
            info_FEMA_WaldTest.numObs{ii}   = numObs;
            info_FEMA_WaldTest.timing{ii}   = toc(tInit);
        end

        % Record overall timing
        info_FEMA_WaldTest.tOverall = toc(tInit_overall);
    else
        [Wald, logp_Wald] = deal([]);
    end

    if sum(~mask)>0
        z_tmp           = zeros(size(zmat,1), size(mask,2));
        p_tmp           = zeros(size(logpmat,1), size(mask,2));
        beta_tmp        = zeros(size(beta_hat,1), size(mask,2));
        betase_tmp      = zeros(size(beta_se,1), size(mask,2));
        sig2tvec_tmp    = zeros(size(sig2tvec,1), size(mask,2));
        coeffCovar_tmp  = zeros(size(coeffCovar, 1), size(coeffCovar, 2), size(mask, 2));
        if strcmpi(CovType, 'unstructured')
            sig2mat_tmp = zeros([size(sig2mat,1:3), size(mask,2)]);
        else
            sig2mat_tmp = zeros(size(sig2mat,1), size(mask,2));
        end

        if ~isempty(splines_of_interest)
            W_tmp           = zeros(size(Wald,1), size(mask,2));
            logp_Wald_tmp   = zeros(size(logp_Wald,1), size(mask,2));
        end

        z_tmp(:, ivec_mask)             = zmat;
        p_tmp(:, ivec_mask)             = logpmat;
        beta_tmp(:, ivec_mask)          = beta_hat;
        betase_tmp(:, ivec_mask)        = beta_se;
        sig2tvec_tmp(:, ivec_mask)      = sig2tvec;
        coeffCovar_tmp(:, :, ivec_mask) = coeffCovar;
        if strcmpi(CovType, 'unstructured')
            sig2mat_tmp(:, :, :, ivec_mask) = sig2mat;
        else
            sig2mat_tmp(:, ivec_mask) = sig2mat;
        end

        if ~isempty(splines_of_interest)
            W_tmp(:,ivec_mask)         = Wald;
            logp_Wald_tmp(:,ivec_mask) = logp_Wald;
        end

        zmat       = z_tmp;
        logpmat    = p_tmp;
        beta_hat   = beta_tmp;
        beta_se    = betase_tmp;
        sig2mat    = sig2mat_tmp;
        sig2tvec   = sig2tvec_tmp;
        coeffCovar = coeffCovar_tmp;

        if ~isempty(splines_of_interest)
            Wald      = W_tmp;
            logp_Wald = logp_Wald_tmp;
        end

        if nperms>0
            % Permutations are not currently handled for unstructured
            % covariance; this part needs to be updated when permutations
            % are incorporated for unstructured covariance
            zperm_tmp           = zeros(size(zmat_perm,1), size(mask,2), size(zmat_perm,3));
            betaperm_tmp        = zeros(size(beta_hat_perm,1), size(mask,2), size(beta_hat_perm,3));
            betaseperm_tmp      = zeros(size(beta_se_perm,1), size(mask,2), size(beta_se_perm,3));
            sig2matperm_tmp     = zeros(size(sig2mat_perm,1), size(mask,2), size(sig2mat_perm,3));
            sig2tvecperm_tmp    = zeros(size(sig2tvec_perm,1), size(mask,2), size(sig2tvec_perm,3));
            coeffCovar_perm_tmp = zeros(size(coeffCovar, 1), size(coeffCovar, 2), size(mask, 2), size(coeffCovar, 4));

            zperm_tmp(:, ivec_mask,:)               = zmat_perm;
            betaperm_tmp(:, ivec_mask,:)            = beta_hat_perm;
            betaseperm_tmp(:, ivec_mask,:)          = beta_se_perm;
            sig2matperm_tmp(:, ivec_mask,:)         = sig2mat_perm;
            sig2tvecperm_tmp(:, ivec_mask,:)        = sig2tvec_perm;
            coeffCovar_perm_tmp(:, :, ivec_mask, :) = coeffCovar_perm;

            zmat_perm       = zperm_tmp;
            beta_hat_perm   = betaperm_tmp;
            beta_se_perm    = betaseperm_tmp;
            sig2mat_perm    = sig2matperm_tmp;
            sig2tvec_perm   = sig2tvecperm_tmp;
            coeffCovar_perm = coeffCovar_perm_tmp;
        end
    end

    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % tfce
    if tfce==1 && nperms>0
        tfce_perm = FEMA_tfce(zmat_perm,colsinterest,datatype,mask);
    elseif tfce==1 && nperms==0
        error('Change nperms>0 to run TFCE OR set tfce=0.')
    else
        tfce_perm=[];
    end

    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % resampling 
    % only save permutations for IVs of interest to reduce size of output
    % if wanting to save permutations for all columns of X make colsinterest=[1:size(X,2)]
    if nperms>0
        beta_hat_perm=beta_hat_perm(colsinterest,:,:);
        beta_se_perm=beta_se_perm(colsinterest,:,:);
        zmat_perm=zmat_perm(colsinterest,:,:);
        colnames_interest=colnames_model(colsinterest);
    end

    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % collate info data 
    info.FEMA_process_data = info_process_data; 
    info.FEMA_makeDesign = info_makeDesign;
    info.FEMA_intersect_design = info_intersect_design;
    info.FEMA_fit = info_fit;
    info.FEMA_WaldTest = info_FEMA_WaldTest;
    info.summary = struct('version', info.FEMA_process_data.FEMA_version, ...
                          'nObservations', info.FEMA_fit.nObservations, ...
                          'nUqSubjects', info.FEMA_fit.nUqSubjects, ...
                          'nFamilies', info.FEMA_fit.nFamilies, ... 
                          'df', info.FEMA_fit.df, ... 
                          'modelSingularity', info.FEMA_fit.modelSingularity, ...
                          'n_cols_X', info.FEMA_fit.num_X, ...
                          'n_cols_ymat', info.FEMA_fit.num_ymat, ...
                          'RandomEffects', strjoin(info.FEMA_fit.settings.RandomEffects, ', '), ...
                          'CovType', info.FEMA_fit.settings.CovType, ...
                          'nbins', info.FEMA_fit.settings.nbins);

    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % save output 
    % always save .mat different save options depending on datatype
    % NB corrmat is always saved as .mat only 
    
    % check roi atlas type if datatype is 'roi' - should proably do this earlier
    if strcmpi(datatype, 'roi')
        roi_atlas_match = regexp(fstem_imaging, '(?<=__)(at|aseg|dsk|dst)(?=$|\.)', 'match');
        if ~isempty(roi_atlas_match)
            roi_atlas = roi_atlas_match{1};
        else
            error('Cannot find the atlas type for ROI datatype.')
        end
    else 
        roi_atlas = []; 
    end
    
    if isempty(outputType)
        % Determine output format based on datatype and roi_atlas
        if strcmpi(datatype, 'roi')
            %if ismember(roi_atlas, {'at', 'aseg'})
            %    outputType = {'roi', 'nifti', 'tables', 'mat', 'summary'};
            %elseif ismember(roi_atlas, {'dsk', 'dst'})
            %    outputType = {'roi', 'gifti', 'tables', 'mat', 'summary'};
            %end
            outputType = {'roi', 'tables', 'mat', 'summary'};
        else
            outputType = {datatype, 'mat', 'summary'};
        end
        % output ids for DEAP 
        if existSet1
            outputType = [outputType, {'ids'}]; %#ok<AGROW>
        end
    end
    
    inputs = inputs.Results;
    if any(ismember(outputType, 'cache'))
        % In this case, additional FEMA_fit parameters can be saved into
        % inputs_FEMA, including ymat
        inputs.X = X;
        inputs.iid = iid;
        inputs.eid = eid;
        inputs.fid = fid;
        inputs.agevec = agevec;
        inputs.ymat = ymat;
        inputs.GRM = GRM.GRM;
        inputs.PregID = PregID;
        inputs.HomeID = HomeID;
        inputs.contrasts = contrasts;
        inputs.synthstruct = synthstruct;
    end

    if any(~ismember(outputType, 'none'))
        info = FEMA_save(outputType, dirname_out{des}, 'outPrefix', outPrefix, ...
                        'beta_hat', beta_hat, 'beta_se', beta_se, 'zmat', zmat, 'logpmat', logpmat, ...
                        'sig2tvec', sig2tvec, 'sig2mat', sig2mat, 'unstructParams', unstructParams, ...
                        'beta_hat_perm', beta_hat_perm, 'beta_se_perm', beta_se_perm, 'zmat_perm', zmat_perm, ...
                        'sig2tvec_perm', sig2tvec_perm, 'sig2mat_perm', sig2mat_perm, 'logLikvec', logLikvec, ...
                        'logLikvec_perm', logLikvec_perm, 'Hessmat', Hessmat, 'coeffCovar', coeffCovar, ...
                        'Wald', Wald, 'logp_Wald', logp_Wald, ...
                        'splines_of_interest', splines_of_interest, 'FFX_conceptMapping', FFX_conceptMapping, ...
                        'binvec_save', binvec_save, 'nvec_bins', nvec_bins, 'tvec_bins', tvec_bins, ...
                        'residuals_GLS', residuals_GLS, 'contrasts', contrasts, 'hypValues', hypValues, ...
                        'designMatrix', designMatrix, 'colnames_model', colnames_model, 'mask', mask, ...
                        'niter', niter, 'nbins', nbins, 'FamilyStruct', FamilyStruct, ...
                        'RandomEffects', RandomEffects, 'nperms', nperms, 'CovType', CovType, ...
                        'FixedEstType', FixedEstType, 'RandomEstType', RandomEstType, 'GroupByFamType', GroupByFamType, ...
                        'NonnegFlag', NonnegFlag, 'precision', precision, 'logLikflag', logLikflag, ...
                        'PermType', PermType, 'colsinterest', colsinterest, 'vars_of_interest', vars_of_interest,...
                        'ymat_names', info.FEMA_process_data.ymat_names, 'roi_atlas', roi_atlas, ...
                        'fstem_imaging', fstem_imaging, 'saveDesignMatrix', saveDesignMatrix, 'info', info, 'matType', matType);
    end

end % looping over desmat 
         
% endtime = now();
logging('***Done*** (%0.2f seconds)', toc(tStart));
% logging('***Done*** (%0.2f seconds)',(endtime-starttime)*3600*24);

end 


