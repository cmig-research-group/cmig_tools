function caller_FEMA_fit_GWAS(file_PLINK, file_FEMA_fit, GWASType, dirOutput, outPrefix, varargin)
% Caller function for FEMA_fit_GWAS (compiled version)
%% Mandatory inputs:
% ------------------
% file_PLINK:       full path to a PLINK file (no extension)
%
% file_FEMA_fit:    full path to output from caller_FEMA_fit (point to
%                   *_estimates.mat file)
%
% GWASType:         should be one of the following (default: long_interact):
%                       * cross:            standard cross-sectional GWAS
%                       * cross_interact:   cross-sectional GWAS with
%                                           interaction of SNPs with basis
%                                           function
%                       * long:             standard longitudinal GWAS
%                       * long_interact:    longitudinal GWAS with
%                                           interaction of SNPs with basis
%                                           function
% 
% dirOutput:        full path to where the results should be written out
% 
% outPrefix:        prefix to the output file names (without extension;
%                   default: FEMA_GWAS-GWASType-yyyyMMMdd-HHmmSS)
%
%% Optional inputs (comma-separated name-value pairs):
% ----------------------------------------------------
% file_basisFunc:   full path to the output file(s) from calling 
%                   caller_createBasisFunctions (either point to the .mat
%                   file or the *-basisFunction.csv file); required for
%                   cross_interact and long_interact analyses (default: [])
%                   If this is a csv file, we assume that the files have
%                   header rows
% 
% appendMainEffect: true or false indicating if the main effect of SNP (an
%                   intercept term) should be appended as the first column
%                   of the basis functions specified in file_BasisFunc; if
%                   yes, the first estimate is the main effect of the SNP,
%                   followed by the estimates of the interaction terms; if
%                   not, the main effect is not estimated - this can be a
%                   scenario where the user has already appended the main
%                   effec of the SNP as part of the basis function OR can
%                   be a scenario where the basis functions were not SVD
%                   modified, therefore, the main effect is captured as
%                   part of the overall basis functions (default: true)
%
% file_contrast:    full path to a contrast file relevant for GWAS: note
%                   that this file is different from the contrast file for
%                   fixed effects contrasts - this file should have column
%                   names corresponding to the basis functions (relevant
%                   for cross_interact and long_intereact GWAS types): by
%                   default, an omnibus eye(k) test is run, where k is the
%                   number of basis functions after accounting for
%                   appendMainEffect being true (default: [])
% 
% splitBy:          how should the genome be split-up for processing;
%                   should be one of the following (default: snp) 
%                       * 'snp'
%                       * 'chromosome'
% 
% chunkSize         if the genome is being split by SNPs, this controls the
%                   number of SNPs in each chunk (default: 1000)
% 
% meanImpute:       should missing SNPs be replaced with mean imputed
%                   values? If not, will resule in NaN (default: true)
% 
% roundOff:         after mean imputation, should the values be converted
%                   to 0/1/2 or left as fractions? (default: false)
% 
% transform:        should be one of the following:
%                       * 'center'
%                       * 'centre'
%                       * 'std'
%                       * 'none' (default)
% 
% stdType:          if transform is 'std', then should be one of the
%                   following (controls how the standardization happens): 
%                       * 'emperical' (default)
%                       * 'gcta'
% 
% pValType:         distribution for calculating p values; should be one of:
%                       * 'z' (default)
%                       * 't'
%                       * 'chi'
% 
% df:               degrees of freedom (required for computing the
%                   appropriately penalised mean squared error); the
%                   default is: size(ymat, 1) - (size(Xvars, 2) + numBasisFunc)
% 
% OLSflag:          true or false indicating if OLS solution should be used
%                   for estimating SNP effect (default: false)
% 
% doF:              true or false indicating if F test should be conducted
%                   instead of the default Wald test, and then use F
%                   distribution for calculating the p values (default: false)
% 
% doCoeffCovar:     if coefficient covariance for SNPs should be saved
%                   (default: false); note that this is false by default as
%                   the variable size can become quite massive; coefficient
%                   covariance is required for any post fitting contrast
%                   estimation
% 
% SingleOrDouble:   control the numerical precision; should be one of the
%                   following (default: setting that was used for FEMA_fit): 
%                       * 'double'
%                       * 'single'
% 
% doPar:            should parallel processing be used to process different
%                   chunks of the genome? (default: false)
% 
% numWorkers:       if doPar is true, how many parallel workers to start
%                   (default: 2)
% 
% numThreads:       if doPar is true, how many threads per parallel worker
%                   (default: 2)
%
%% To Do:
% * Handle scneario where *_designMatrix file is not saved as part of
% caller_FEMA_fit - this will require inputs of file_X, file_ymat,
% re-intersecting the data, calling FEMA_parse_family to generate
% clusterinfo, and then proceeding with analysis
% 
% * incorporate interaction by sex - what should be the input file?

%% Start
updateString = [char(datetime('now')), ': caller_FEMA_fit_GWAS: job started'];
disp(updateString);

%% Check mandatory inputs
if ~exist('file_PLINK', 'var') || isempty(file_PLINK)
    error('Please provide a full path to a PLINK file (without extension)');
end

if ~exist('file_FEMA_fit', 'var') || isempty(file_FEMA_fit)
    error('Please provide a full path to the output from caller_FEMA_fit');
else
    % Check if the estimates, the designMatrix, and the settings files exist
    if ~exist(file_FEMA_fit, 'file')
        error(['Unable to find file: ', file_FEMA_fit]);
    else
        file_FEMA_designMatrix = strrep(file_FEMA_fit, '_estimates.mat', '_designMatrix.mat');
        if ~exist(file_FEMA_designMatrix, 'file')
            error(['Unable to find file: ', file_FEMA_designMatrix]);
        else
            file_FEMA_settings = strrep(file_FEMA_fit, '_estimates.mat', '_settings.mat');
            if ~exist(file_FEMA_settings, 'file')
                error(['Unable to find file: ', file_FEMA_settings]);
            end
        end
    end
end

if ~exist('GWASType', 'var') || isempty(GWASType)
    error('Please specify the type of GWAS to be performed');
else
    GWASType = lower(GWASType);
    chkGWAS  = {'cross', 'cross_interact', 'long', 'long_interact'};
    if ~ismember(GWASType, chkGWAS)
        error(['Unknown GWASType specified: ', GWASType]);
    else
        if strcmpi(GWASType, 'cross_interact') || strcmpi(GWASType, 'long_interact')
            doInteraction = true;
        else
            doInteraction = false;
        end
    end
end

if ~exist('dirOutput', 'var') || isempty(dirOutput)
    error('Please provide a full path to where results should be saved');
else
    if ~exist(dirOutput, 'dir')
        mkdir(dirOutput);
    end
end

if ~exist('outPrefix', 'var') || isempty(outPrefix)
    outPrefix = ['FEMA_GWAS-', GWASType, '-', char(datetime('now', 'Format', 'yyyyMMMdd-HHmmSS'))];
end

%% Assign default optional inputs
p = inputParser;
addParameter(p, 'file_basisFunc', '');
addParameter(p, 'appendMainEffect', true);
addParameter(p, 'file_contrast', '');
addParameter(p, 'splitBy', 'snp');
addParameter(p, 'chunkSize', 1000);
addParameter(p, 'meanImpute', true);
addParameter(p, 'roundOff', false);
addParameter(p, 'transform', 'none');
addParameter(p, 'stdType', 'emperical');
addParameter(p, 'pValType', 'z');
addParameter(p, 'df', []);
addParameter(p, 'OLSflag', false);
addParameter(p, 'doF', false);
addParameter(p, 'doCoeffCovar', false);
addParameter(p, 'SingleOrDouble', 'double');
addParameter(p, 'doPar', false);
addParameter(p, 'numWorkers', 2);
addParameter(p, 'numThreads', 2);

%% Parse optional inputs
parse(p, varargin{:})
file_basisFunc      = p.Results.file_basisFunc;
appendMainEffect    = p.Results.appendMainEffect;
file_contrast       = p.Results.file_contrast;
splitBy             = p.Results.splitBy;
meanImpute          = p.Results.meanImpute;
roundOff            = p.Results.roundOff;
transform           = p.Results.transform;
stdType             = p.Results.stdType;
pValType            = p.Results.pValType;
df                  = p.Results.df;
OLSflag             = p.Results.OLSflag;
doF                 = p.Results.doF;
doCoeffCovar        = p.Results.doCoeffCovar;
SingleOrDouble      = p.Results.SingleOrDouble;
doPar               = p.Results.doPar;

if ~isnumeric(p.Results.chunkSize)
    chunkSize = str2double(p.Results.chunkSize);
else
    chunkSize = p.Results.chunkSize;
end

if ~isnumeric(p.Results.numWorkers)
    numWorkers = str2double(p.Results.numWorkers);
else
    numWorkers = p.Results.numWorkers;
end

if ~isnumeric(p.Results.numThreads)
    numThreads = str2double(p.Results.numThreads);
else
    numThreads = p.Results.numThreads;
end

% Not quite sure why the logical values need to be converted from character
% here but not in caller_FEMA_fit
if isdeployed
    chkFun = @(x) (ischar(x) | isstring(x)) && ismember(lower(x), {'true'});

    if chkFun(appendMainEffect)
        appendMainEffect = true;
    else
        appendMainEffect = false;
    end

    if chkFun(meanImpute)
        meanImpute = true;
    else
        meanImpute = false;
    end
    
    if chkFun(roundOff)
        roundOff = true;
    else
        roundOff = false;
    end

    if chkFun(OLSflag)
        OLSflag = true;
    else
        OLSflag = false;
    end

    if chkFun(doF)
        doF = true;
    else
        doF = false;
    end

    if chkFun(doCoeffCovar)
        doCoeffCovar = true;
    else
        doCoeffCovar = false;
    end

    if chkFun(doPar)
        doPar = true;
    else
        doPar = false;
    end
end
disp(p.Results);

%% Make some decisions based on optional input
% If interaction is true, is the basis function file provided?
if doInteraction
    if isempty(file_basisFunc)
        error(['Please provide full path to a file containing the ', ...
               'basis functions which will be used for analyses']);
    else
        % If this is a mat file, read the basis function
        [~, ~, ext] = fileparts(file_basisFunc);
        if isempty(ext)
            error('Please provide full path to the basis function file');
        else
            if strcmpi(ext, '.mat')
                temp_bf  = load(file_basisFunc, 'basisFunction', 'varNames_basisFunction');
                toChk_bf = {'basisFunction', 'varNames_basisFunction'};
                chk_bf   = ismember(toChk_bf, fieldnames(temp_bf));
                if sum(chk_bf) ~= length(toChk_bf)
                    missMsg = sprintf('%s, ', toChk_bf{~chk_bf});
                    error(['The following required variables are missing from ', ...
                           file_basisFunc, ': ', missMsg(1:end-2)]);
                else
                    bfSNP                  = temp_bf.basisFunction;
                    varNames_basisFunction = temp_bf.varNames_basisFunction;
                    clear temp_bf;
                end
            else
                if strcmpi(ext, '.csv')
                    temp_bf                = readtable(file_basisFunc);
                    bfSNP                  = temp_bf{:,:};
                    varNames_basisFunction = temp.Properties.VariableNames;
                    clear temp_bf;
                else
                    error(['Unknown extension provided for file_basisFunc: ', ext]);
                end
            end
        end
    end
else
    bfSNP = [];
end

% Do basis functions need to be appended with intercept?
if exist('bfSNP', 'var') && ~isempty(bfSNP)
    if appendMainEffect
        bfSNP = [ones(size(bfSNP, 1), 1), bfSNP];
    end
end

% Do the basis functions have associated contrast file(s)?
if exist('file_contrast', 'var') && ~isempty(file_contrast)
    if ~exist(file_contrast, 'file')
        error(['Unable to find: ', file_contrast]);
    else
        [contrasts, hypValues] = FEMA_parse_contrastFile(file_contrast, varNames_basisFunction);
    end
else
    contrasts = [];
    hypValues = [];
end

updateString = [char(datetime('now')), ': caller_FEMA_fit_GWAS: finished input parsing'];
disp(updateString);

%% Get relevant variables from FEMA_fit results
updateString = [char(datetime('now')), ': caller_FEMA_fit_GWAS: loading FEMA_fit results'];
disp(updateString);

% We need the following:
% ymat residuals:   saved in estimates file (reusableVars)
% sig2mat:          saved in estimates file
% binvec_save:      saved in estimates file
% FamilyStruct:     saved in design matrix
% iid:              saved in design matrix
% X:                saved in design matrix
% SingleOrDouble:   saved in settings file
% RandomEffects:    saved in settings file
% CovType:          saved in settings file

temp_estimates  = load(file_FEMA_fit, 'reusableVars', 'sig2mat', 'binvec_save');
temp_designvars = load(file_FEMA_designMatrix, 'FamilyStruct', 'iid', 'X');
temp_settings   = load(file_FEMA_settings, 'RandomEffects', 'SingleOrDouble', 'CovType');

% Make sure all fields are present
toChk_estimates = {'reusableVars', 'sig2mat', 'binvec_save'};
chk_estimates   = ismember(toChk_estimates, fieldnames(temp_estimates));
if sum(chk_estimates) ~= length(toChk_estimates)
    missMsg = sprintf('%s, ', toChk_estimates{~chk_estimates});
    error(['The following required variables are missing from ', ...
           file_FEMA_fit, ': ', missMsg(1:end-2)]);
else
    reusableVars = temp_estimates.reusableVars;
    sig2mat      = temp_estimates.sig2mat;
    binvec_save  = temp_estimates.binvec_save;
    clear temp_estimates;
end

toChk_design = {'FamilyStruct', 'iid', 'X'};
chk_design   = ismember(toChk_design, fieldnames(temp_designvars));
if sum(chk_design) ~= length(toChk_design)
    missMsg = sprintf('%s, ', toChk_design{~chk_design});
    error(['The following required variables are missing from ', ...
           file_FEMA_designMatrix, ': ', missMsg(1:end-2)]);
else
    FamilyStruct = temp_designvars.FamilyStruct;
    iid          = temp_designvars.iid;
    X            = temp_designvars.X;
    clear temp_designvars;
end

toChk_settings = {'RandomEffects', 'SingleOrDouble', 'CovType'};
chk_settings   = ismember(toChk_settings, fieldnames(temp_settings));
if sum(chk_settings) ~= length(toChk_settings)
    missMsg = sprintf('%s, ', toChk_settings{~chk_settings});
    error(['The following required variables are missing from ', ...
           file_FEMA_settings, ': ', missMsg(1:end-2)]);
else
    RandomEffects = temp_settings.RandomEffects;
    CovType       = temp_settings.CovType;
    if isempty(SingleOrDouble)
        SingleOrDouble = temp_settings.SingleOrDouble;
    else
        if ~strcmpi(temp_settings.SingleOrDouble, SingleOrDouble)
            warning(['Different settings of SingleOrDouble between ', ...
                     'FEMA_fit and current call to FEMA_fit_GWAS; using: ', SingleOrDouble]);
        end
    end
    clear temp_settings
end

%% Get some information about genetics
updateString = [char(datetime('now')), ': caller_FEMA_fit_GWAS: extracting basic genetics information'];
disp(updateString);

onlyCheck = true;
lowMem    = false;

[~, Chr, SNPID, BP, check, errMsg, genInfo] =             ...
 FEMA_parse_PLINK(file_PLINK, iid, [], onlyCheck, lowMem, ...
                  meanImpute, roundOff, transform, stdType);
if ~check
    error(errMsg);
end

%% Compile some variables
updateString = [char(datetime('now')), ': caller_FEMA_fit_GWAS: compiling V term'];
disp(updateString);

[allWsTerms, timing_compile] = FEMA_compileTerms(FamilyStruct.clusterinfo, binvec_save,      ...
                                           sig2mat, RandomEffects, FamilyStruct.famtypevec,  ...
                                           reusableVars.GroupByFamType, CovType,             ...
                                           SingleOrDouble, reusableVars.visitnum);

%% Divide genome into chunks
updateString = [char(datetime('now')), ': caller_FEMA_fit_GWAS: splitting genome into parts'];
disp(updateString);

[splitInfo, timing_divide] = divideSNPs(file_PLINK, splitBy, chunkSize, outPrefix, ...
                                        SNPID,      Chr,     BP,        genInfo);

%% Extract out ymat residuals and clear up some space
ymat_res_gls = reusableVars.ymat_res_gls;
clear reusableVars FamilyStruct iid

%% Initialize parallel cluster
if doPar
    updateString = [char(datetime('now')), ': caller_FEMA_fit_GWAS: initializing parallel cluster'];
    disp(updateString);
    
    local = parcluster('local');
    local.NumThreads = numThreads;
    pool = local.parpool(numWorkers);    
end

%% Do model fitting
t_GWAS = tic;
updateString = [char(datetime('now')), ': caller_FEMA_fit_GWAS: calling FEMA_fit_GWAS'];
disp(updateString);

if doPar
    parfor parts = 1:length(splitInfo)
        t1 = tic;
        FEMA_fit_GWAS(splitInfo{parts}, ymat_res_gls, binvec_save, X, allWsTerms, ...
                      'bfSNP', bfSNP, 'pValType', pValType, 'df', df,             ...
                      'OLSflag', OLSflag, 'SingleOrDouble', SingleOrDouble,       ...
                      'outDir', dirOutput, 'outName', splitInfo{parts}.outName,   ...
                      'L', contrasts, 'hypValue', hypValues, 'doF', doF,          ...
                      'doCoeffCovar', doCoeffCovar);
    
        % Also display progress
        disp(['Finished: ', num2str(parts, '%04d'), ' in ', num2str(toc(t1), '%.2f'), 's']);
    end
    delete(pool);
else
    for parts = 1:length(splitInfo)
        t1 = tic;
        FEMA_fit_GWAS(splitInfo{parts}, ymat_res_gls, binvec_save, X, allWsTerms, ...
                      'bfSNP', bfSNP, 'pValType', pValType, 'df', df,             ...
                      'OLSflag', OLSflag, 'SingleOrDouble', SingleOrDouble,       ...
                      'outDir', dirOutput, 'outName', splitInfo{parts}.outName,   ...
                      'L', contrasts, 'hypValue', hypValues, 'doF', doF,          ...
                      'doCoeffCovar', doCoeffCovar);
        
        % Also display progress
        disp(['Finished: ', num2str(parts, '%04d'), ' in ', num2str(toc(t1), '%.2f'), 's']);
    end
end
timing_GWAS = toc(t_GWAS);

%% Save miscellenous information as a mat file
% GWAS results are already saved
updateString = [char(datetime('now')), ': caller_FEMA_fit_GWAS: finished FEMA_fit_GWAS'];
disp(updateString);

% Get a sense for all variables in the workspace
tmpInfo = whos;

saveName = fullfile(dirOutput, [outPrefix, '_miscInfo.mat']);
toSave   = {'timing_compile', 'timing_divide', 'timing_GWAS',                    ...
            'file_PLINK', 'file_FEMA_fit', 'file_basisFunc',  'file_contrast',   ...
            'bfSNP', 'pValType', 'df', 'OLSflag', 'SingleOrDouble', 'contrasts', ...
            'doF', 'doCoeffCovar', 'GWASType',  'appendMainEffect', 'splitBy',   ...
            'chunkSize', 'meanImpute', 'roundOff', 'transform', 'stdType',       ...
            'doPar', 'numWorkers', 'numThreads', 'hypValues'};

if sum([tmpInfo(ismember({tmpInfo(:).name}', toSave)).bytes]) > 2^31
    save(saveName, 'timing_compile', 'timing_divide', 'timing_GWAS',                     ...
                    'file_PLINK', 'file_FEMA_fit', 'file_basisFunc',  'file_contrast',   ...
                    'bfSNP', 'pValType', 'df', 'OLSflag', 'SingleOrDouble', 'contrasts', ...
                    'doF', 'doCoeffCovar', 'GWASType',  'appendMainEffect', 'splitBy',   ...
                    'chunkSize', 'meanImpute', 'roundOff', 'transform', 'stdType',       ...
                    'doPar', 'numWorkers', 'numThreads', 'hypValues', '-v7.3');
else
    save(saveName, 'timing_compile', 'timing_divide', 'timing_GWAS',                     ...
                    'file_PLINK', 'file_FEMA_fit', 'file_basisFunc',  'file_contrast',   ...
                    'bfSNP', 'pValType', 'df', 'OLSflag', 'SingleOrDouble', 'contrasts', ...
                    'doF', 'doCoeffCovar', 'GWASType',  'appendMainEffect', 'splitBy',   ...
                    'chunkSize', 'meanImpute', 'roundOff', 'transform', 'stdType',       ...
                    'doPar', 'numWorkers', 'numThreads', 'hypValues');
end