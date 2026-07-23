function [beta_hat,      beta_se,        zmat,        logpmat,          ...
          sig2tvec,      sig2mat,        Hessmat,     logLikvec,        ...
          beta_hat_perm, beta_se_perm,   zmat_perm,   sig2tvec_perm,    ...
          sig2mat_perm,  logLikvec_perm, binvec_save, nvec_bins,        ...
          tvec_bins,     FamilyStruct,   coeffCovar,  unstructParams,   ...
          residuals_GLS, info] = FEMA_fit_binary(X, iid, eid, fid,      ...
                                                 agevec, ymat, maxIter, ...
                                                 contrasts, nbins, GRM, varargin)
% Function to fit fast and efficient linear mixed effects model
%
% For notation below:
% n = observations,
% p = predictors (fixed effects),
% v = number of outcome variables - currently set as 1
% c = number of contrasts to evaluate
% r = number of random effects
%

%% Inputs:
% X               <num>            [n x p]    design matrix, with intercept if needed
% iid             <cell>           [n x 1]    subject IDs to match imaging data
% eid             <cell>           [n x 1]    eventname
% fid             <num>            [n x 1]    family ID (members of the same family unit have same value)
% agevec          <num>            [n x 1]    participants age
% ymat            <num>            [n x v]    matrix of imaging data
% niter           <num>            [1 x 1]    maximal number of iterations (default 200) 
% contrasts       <num> OR <path>  [c x p]    contrast matrix, where c is number of contrasts to compute,
%                                             OR path to file containing contrast matrix (readable by readtable)
% pihatmat        <num>            [n x n]    matrix of genetic relatedness --> already intersected to match X and Y sample
%


%% Optional input arguments:
% RandomEffects   <cell>           list of random effects to estimate (default {'F','S'}):
%                                       * F:  family relatedness
%                                       * S:  subject - required for longitudinal analyses
%                                       * A:  additive genetic relatedness - must include file path to genetic relatedness data (pihat) for this option
%                                       * D:  dominant genetic relatedness - square of A
%                                       * M:  maternal effect - effect of having same mother
%                                       * P:  paternal effect  - effect of having same father
%                                       * H:  home effect - effect of living at the same address
%                                       * T:  twin effect - effect of having the same pregnancy ID
% nperms          <num>            deault 0 --> if >0 will run and output permuted effects
% CovType         <char>           default 'analytic' --> no other options currently available
% RandomEstType   <char>           default 'MoM' --> other option: 'ML' (much slower)
% GroupByFamType  <boolean>        default true
% NonnegFlag      <blooean>        default true - non-negativity constraint on random effects estimation
% precision  <char>           default 'double' --> other option: 'single' - for precision
% logLikflag      <boolean>        default true - compute log-likelihood
% PermType        <char>           permutation type:
%                                       * 'wildbootstrap':    residual boostrap --> creates null distribution by randomly flipping the sign of each observation
%                                       * 'wildbootstrap-nn': non-null boostrap --> estimates distribution around effect of interest using sign flipping (used for sobel test)
% returnReusable  <boolean>        default false - if true, additionally returns reusableVars as a structure with some variables that can be reused (primarily by FEMA-GWAS)
% maxIter         <num>            default 200
% tol             <num>            default 1e-4


%% Outputs:
% beta_hat                         [c+p x v]  estimated beta coefficients
% beta_se                          [c+p x v]  estimated beta standard errors
% zmat                             [c+p x v]  z statistics
% logpmat                          [c+p x v]  log10 p-values
% sig2tvec                         [1   x v]  total residual error of model at each vertex/voxel
% sig2mat                          [r   x v]  normalized random effect variances
% FamilyStruct                                structure type (can be passed as input to avoid re-parsing family structure etc.)

%% Parse inputs
tInit = tic;
logging(FEMA_info);
logging('***Start***');

% before = memory().MemUsedMATLAB / 1024^2;
% Extremely quick sanity check on X and y variables
% First, make sure that they are either single or double precision; if not,
% cast as double precision
if ~ismember(class(X), {'single', 'double'})
    X = cast(X, 'double');
end
if ~ismember(class(ymat), {'single', 'double'})
    ymat = cast(ymat, 'double');
end

% Now make sure that there are no NaNs or Infs
if logical(sum(any(isnan(X)))) || logical(sum(any(isnan(ymat)))) || ...
        logical(sum(any(isinf(X)))) || logical(sum(any(isinf(ymat))))
    error('X and/or ymat have NaN or Inf; please check your data');
else
    % Additional check for constant values in y variables
    if any(var(ymat) == 0)
        warning('One or more columns in ymat are constant');
    end
end

p = inputParser;

if ~exist('maxIter', 'var') || isempty(maxIter)
    maxIter = 200;
end

if ~exist('contrasts', 'var')
    contrasts = [];
end

if ~isfinite(contrasts)
    fname_contrasts = p.Results.contrasts;
    logging('Reading contrast matrix from %s', fname_contrasts);
    contrasts = readtable(fname_contrasts);
end

% Zeros-pad contrasts, if needed
if ~isempty(contrasts) && size(contrasts,2) < size(X,2)
    contrasts = cat(2, contrasts, zeros([size(contrasts, 1) size(X, 2) - size(contrasts, 2)]));
end

if ~exist('nbins', 'var') || isempty(nbins)
    nbins = 20;
end

if ~exist('GRM', 'var')
    GRM = [];
end

% Should change to allow p to be passed in, so as to avoid having to
% duplicate input argument parsing in FEMA_wrapper and FEMA_fit
p = inputParser;
addParamValue(p,'CovType', 'analytic'); %#ok<*NVREPLA>
addParameter(p, 'FixedEstType', 'GLS');
addParamValue(p,'RandomEstType', 'MoM');
addParamValue(p,'PermType', 'wildbootstrap');
addParamValue(p,'GroupByFamType', true);
addParamValue(p,'NonnegFlag', true); % Perform lsqnonneg on random effects estimation
addParamValue(p,'precision', 'double');
addParamValue(p,'RandomEffects', {'F' 'S' 'E'}); % Default to Family, Subject, and eps
addParamValue(p,'logLikflag', false);
addParamValue(p,'Hessflag', false);
addParamValue(p,'ciflag', false);
addParamValue(p,'nperms', 0);
addParamValue(p,'FatherID', {}); % Father ID, ordered same as GRM
addParamValue(p,'MotherID', {}); % Mother ID, ordered same as GRM
addParamValue(p,'PregID', {}); % Pregnancy effect (same ID means twins), ordered same as GRM
addParamValue(p,'HomeID', {}); % Home effect (defined as same address ID), ordered same as GRM
addParamValue(p,'FamilyStruct', {}); % Avoids recomputing family strucutre et al
addParameter(p, 'returnResiduals', false); % Additionally returns GLS residuals
addParamValue(p,'synthstruct', ''); % True / synthesized random effects
addParamValue(p,'maxIter', 200); % Maximum iteration numbers
addParamValue(p,'tol', 1e-4); % tolerance set for both fixed and random effects - Line 448
addParamValue(p,'AddIntercept', true); % default add intercept
addParameter(p, 'doPar', false);
addParameter(p, 'numWorkers', 2);
addParameter(p, 'numThreads', 2);

parse(p,varargin{:})

CovType              = p.Results.CovType; %#ok<*NASGU>
FixedEstType         = lower(p.Results.FixedEstType);
RandomEstType        = p.Results.RandomEstType;
GroupByFamType       = p.Results.GroupByFamType;
NonnegFlag           = p.Results.NonnegFlag;
precision            = p.Results.precision;
RandomEffects        = p.Results.RandomEffects;
MoMflag              = ismember(lower(RandomEstType), {'mom'});
MLflag               = ismember(lower(RandomEstType), {'ml'});
logLikflag           = p.Results.logLikflag;
Hessflag             = p.Results.Hessflag;
ciflag               = p.Results.ciflag;
nperms               = p.Results.nperms;
PermType             = p.Results.PermType;
FamilyStruct         = p.Results.FamilyStruct;
returnResiduals      = p.Results.returnResiduals;
synthstruct          = p.Results.synthstruct;
maxIter              = p.Results.maxIter;
tol                  = p.Results.tol;
AddIntercept         = p.Results.AddIntercept;
doPar                = p.Results.doPar;
numWorkers           = p.Results.numWorkers;
numThreads           = p.Results.numThreads;

% Assign some logical operators
OLSflag         = ismember(FixedEstType,  {'ols'});
MLflag          = ismember(RandomEstType, {'ml'});
unstructuredCov = ismember(CovType, {'unstructured'});

% Check if lsqminnorm can be used
if exist('lsqminnorm', 'file')
    useLSQ = true;
else
    useLSQ = false;
end

% Ensure CovType is valid
if ~ismember(CovType, {'analytic', 'analytical', 'unstructured'})
    warning(['Unknown CovType specified: ', CovType, '; setting CovType to analytic']);
    CovType = 'analytic';
else
    if strcmpi(CovType, 'analytical')
        CovType = 'analytic';
    end
end

% Ensure permType is valid
if ~isempty(PermType)
    if strcmpi(PermType, 'none')
        PermType = [];
        nperms   = 0;
    else
        if ~ismember(PermType, {'wildbootstrap', 'wildbootstrap-nn'})
            error(['Unknown resampling scheme specified: ', PermType, '; PermType should be either wildbootstrap or wildbootstrap-nn']);
        end
    end
end

% If permutation and unstructured covariance, warn the user
if unstructuredCov && nperms > 0
    warning('Permutations not yet implemented for unstructured covariance');
    nperms = 0;
end

% Examine RandomEffects and ensure E is always the last term - relevant for
% unstructured covariance
RandomEffects = rowvec(RandomEffects);
tmp           = strcmpi(RandomEffects, 'E');
if ~any(tmp)
    warning('RandomEffects did not include E term; appending E as the last random effect');
    RandomEffects = [RandomEffects, 'E'];
else
    if find(tmp) ~= length(RandomEffects)
        RandomEffects = [RandomEffects(~tmp), RandomEffects(tmp)];
        logging(['Re-arranging RandomEffects as: ', sprintf('%s ', RandomEffects{:})]);
    end
end

% Grouping by family type is only supported for RandomEffects 'F' 'S' 'E'
if ~isempty(setdiff(RandomEffects,{'F' 'S' 'E'}))
    GroupByFamType = false;
end

if ~unstructuredCov
    unstructParams = [];
end

if ~returnResiduals
    residuals_GLS = [];
end

% Check if lsqminnorm can be used
if exist('lsqminnorm', 'file')
    useLSQ = true;
else
    useLSQ = false;
end

% Get some basic info
[num_obs, num_y] = size(ymat);
num_RFX          = length(RandomEffects);

% Add intercept, if required
if AddIntercept
    if ~all(X(:,1)==1)
        X = [ones(num_obs,1), X];
    end
end

% Number of X variables (updated, if intercept was added)
num_X = size(X, 2);

% Check if X is rank deficient
if rank(double(X)) < num_X
    lowRank = true;
else
    lowRank = false;
end

% permutation initialization and ensuring that all outputs are initialized
if nperms>0
    beta_hat_perm = zeros(num_X, nperms);
    beta_se_perm  = zeros(num_X, nperms);
    zmat_perm     = zeros(num_X, nperms);
    sig2tvec_perm = zeros(1, nperms);
    sig2mat_perm  = zeros(num_RFX, nperms);
    if logLikflag
        logLikvec_perm = zeros(1, nperms);
    else
        logLikvec_perm = [];
    end
else
    [logLikvec, beta_hat_perm, beta_se_perm, zmat_perm, ...
     sig2tvec_perm, sig2mat_perm, logLikvec_perm] = deal([]);
end

[binvec_save, nvec_bins, tvec_bins] = deal([]);

%% Save all input parameters
info.FEMA_version                = FEMA_info;
info.provenance                  = 'FEMA_fit_binary';
info.settings.nbins              = nbins;
info.settings.GRM_input          = ~isempty(GRM);
info.settings.contrasts_input    = ~isempty(contrasts);
info.settings.GroupByFamType     = GroupByFamType;
info.settings.RandomEffects      = RandomEffects;
info.settings.precision          = precision;
info.settings.OLSflag            = OLSflag;
info.settings.useLSQ             = useLSQ;
info.settings.unstructuredCov    = unstructuredCov;
info.settings.lowRank            = lowRank;
info.settings.CovType            = CovType;
info.settings.maxIter            = maxIter;
info.settings.FixedEstType       = FixedEstType;
info.settings.RandomEstType      = RandomEstType;
info.settings.NonnegFlag         = NonnegFlag;
info.settings.logLikflag         = logLikflag;
info.settings.Hessflag           = Hessflag;
info.settings.ciflag             = ciflag;
info.settings.nperms             = nperms;
info.settings.PermType           = PermType;
info.settings.FamilyStruct_input = ~isempty(FamilyStruct);
info.settings.synthstruct        = synthstruct;
info.settings.returnResiduals    = returnResiduals;
info.settings.doPar              = doPar;
info.settings.numWorkers         = numWorkers;
info.settings.numThreads         = numThreads;

%% Report model singularity
% Should perhaps report a more standard measure of model singularity?
modelSingularity = cond(X'*X)/cond(diag(diag(X'*X)));
logging('Model singularity index = %g', modelSingularity);

%% Save some basic information
info.FEMA_version       = FEMA_info;
info.num_X              = num_X;
info.num_ymat           = num_y;
info.num_RFX            = num_RFX;
info.lowRank            = lowRank;
info.modelSingularity   = modelSingularity;
info.timing.parseInputs = toc(tInit);

%% Parse family structure, if necessary
if ~exist('FamilyStruct', 'var') || isempty(FamilyStruct)
    % Save some information
    info.nObservations = length(iid);
    info.nUqSubjects   = length(unique(iid));

    tInit_parseFamily = tic;
    [clusterinfo, Ss, iid, famtypevec, famtypelist, subj_famtypevec] =                  ...
     FEMA_parse_family(iid, eid, fid, agevec, GRM, 'RandomEffects', RandomEffects, ...
                       'FatherID', p.Results.FatherID,  'MotherID', p.Results.MotherID, ...
                       'PregID',   p.Results.PregID,    'HomeID',   p.Results.HomeID); %#ok<*ASGLU>
    
    numUqSubjs = info.nUqSubjects;
    [~, ~, IC_subj] = unique(iid,'stable'); % nsubj = length(iid_list); num_obs = length(iid);
    [~, ~, IC_fam]  = unique(fid,'stable'); % nfam  = length(fid_list);
    nfam      = length(unique(fid));
    nfamtypes = length(famtypelist);

    % Save some more information
    info.nFamilies = length(clusterinfo);
    info.nFamTypes = nfamtypes;
    
    % remove after the if-else
    % RandomVar{1} = sparse(1:num_obs, IC_fam, ones(num_obs,1), num_obs, nfam);
    % RandomVar{2} = sparse(1:num_obs, IC_subj, ones(num_obs,1), num_obs, nsubj);
    
    % Prepare generalized matrix version of MoM estimator
    % tic
    S_sum = Ss{1};
    for i = 2:length(Ss)
        S_sum = S_sum + Ss{i};
    end
    [subvec1, subvec2] = find(S_sum); % Use full matrix, to simplify IGLS -- should be possible to limit to tril
    %[subvec1 subvec2] = find(tril(S_sum)); % Should exclude diagonals: tril(S_sum,-1)
    indvec = sub2ind([num_obs num_obs],subvec1,subvec2);

    % Should delete as there's no need to do binning
    % F_num = S_sum;
    % for fi = 1:nfam
    %     F_num(clusterinfo{fi}.jvec_fam,clusterinfo{fi}.jvec_fam) = fi;
    % end
    % fnumvec = F_num(indvec);
    % 
    % for fi = 1:nfam
    %     jvec_tmp  = clusterinfo{fi}.jvec_fam;
    %     [sv, si]  = sort(jvec_tmp);
    %     I_tmp     = reshape(1:length(jvec_tmp)^2, length(jvec_tmp) * [1 1]);
    %     ivec_fam  = find(fnumvec==fi);
    %     ivec_fam  = ivec_fam(colvec(I_tmp(si, si)));
    %     %  ivec_fam = find(fnumvec==fi); ivec_fam(colvec(I_tmp(si,si))) = ivec_fam;
    %     clusterinfo{fi}.ivec_fam = ivec_fam;
    % end

    % Scale back to using tril on S_sum
    % [subvec1, subvec2] = find(tril(S_sum)); % Should exclude diagonals: tril(S_sum,-1)
    % indvec             = sub2ind([num_obs num_obs],subvec1,subvec2);

    M = zeros(length(indvec),length(Ss));
    for i = 1:length(Ss)
        M(:,i) = Ss{i}(indvec);
    end

    % Should delete as there's no need to do binning
    % Create grid of normalized random effects
    % binvals_edges       = linspace(0,1,nbins+1); 
    % binvals_edges(end)  = binvals_edges(end)+0.0001;
    % 
    % % New ND version
    % if length(RandomEffects) == 2
    %     sig2gridi = colvec(1:length(binvals_edges)-1);
    %     sig2gridl = colvec(binvals_edges(1:end-1));
    %     sig2gridu = colvec(binvals_edges(2:end));
    % else
    %     sig2gridi = ndgrid_amd(repmat({1:length(binvals_edges)-1}, [1 length(RandomEffects)-1]));
    %     sig2gridl = ndgrid_amd(repmat({binvals_edges(1:end-1)},    [1 length(RandomEffects)-1]));
    %     sig2gridu = ndgrid_amd(repmat({binvals_edges(2:end)},      [1 length(RandomEffects)-1]));
    % end
    % sig2grid_ivec = find(sum(sig2gridl,2)<=1); % Get rid of "impossible" bins
    % sig2gridl     = sig2gridl(sig2grid_ivec,:);
    % sig2gridu     = sig2gridu(sig2grid_ivec,:);
    % sig2gridi     = sig2gridi(sig2grid_ivec,:);
    % sig2grid      = (sig2gridl+sig2gridu)/2;
    % sig2gridind   = sub2ind_amd(nbins*ones(1,length(RandomEffects)-1),sig2gridi);
    % nsig2bins     = size(sig2gridl,1); % Should handle case of no binning

    % Prepare FamilyStruct
    FamilyStruct = struct('clusterinfo', {clusterinfo}, 'M', {M},                     ...
                          'famtypevec',  {famtypevec},  'famtypelist', {famtypelist}, ...
                          'nfamtypes',   nfamtypes,     'nfam', nfam,                 ...
                          'subvec1',     subvec1,       'subvec2', subvec2,           ...
                          'Ss',          {Ss});
    
    info.timing.tParseFamily = toc(tInit_parseFamily);
else
    clusterinfo = FamilyStruct.clusterinfo;
    M           = FamilyStruct.M;
    % nsig2bins   = FamilyStruct.nsig2bins;
    nfam        = FamilyStruct.nfam; % nfam defined by the fid? duplicate
    famtypevec  = FamilyStruct.famtypevec;
    nfamtypes   = FamilyStruct.nfamtypes;
    % sig2grid    = FamilyStruct.sig2grid;
    % sig2gridl   = FamilyStruct.sig2gridl;
    % sig2gridu   = FamilyStruct.sig2gridu;
    subvec1     = FamilyStruct.subvec1;
    subvec2     = FamilyStruct.subvec2;
    Ss          = FamilyStruct.Ss;

    % Remove outside the if-else - duplicated
    % [iid_list, IA, IC_subj] = unique(iid,'stable'); nsubj = length(iid_list); num_obs = length(iid);
    % [fid_list, IA, IC_fam]  = unique(fid,'stable'); nfam  = length(fid_list);
    % RandomVar{1} = sparse(1:num_obs, IC_fam, ones(num_obs,1), num_obs, nfam);
    % RandomVar{2} = sparse(1:num_obs, IC_subj, ones(num_obs,1), num_obs, nsubj);
end

% Should be generalized ------- later
RandomVar = struct();
RandomVar.("V_F") = sparse(1:num_obs, IC_fam, ones(num_obs,1),  num_obs, nfam);
RandomVar.("V_S") = sparse(1:num_obs, IC_subj, ones(num_obs,1), num_obs, numUqSubjs);

Mi = single(pinv(M));
Cov_MoM = Mi*Mi'; % Variance  / covariance of MoM estimates, per unit of residual error variance

logging('size(M) = [%d %d]',size(M));
logging('Cov_MoM:'); disp(Cov_MoM);
logging('Mi*M:'); disp(Mi*M);

if ~isempty(synthstruct) % do we need this anymore?
    sig2mat_true  = synthstruct.sig2mat_true;
    sig2tvec_true = synthstruct.sig2tvec_true;

    % nvec_bins_true = NaN(nsig2bins,1);
    % binvec_true    = NaN(1,size(ymat,2));
    % for sig2bini = 1:nsig2bins
    %     tmpvec = true;
    %     for ri = 1:size(sig2mat_true,1)-1
    %         tmpvec = tmpvec & sig2mat_true(ri,:) >= sig2gridl(sig2bini,ri) & ...
    %                           sig2mat_true(ri,:) <  sig2gridu(sig2bini,ri);
    %     end
    %     ivec_bin = find(tmpvec);
    %     nvec_bins_true(sig2bini) = length(ivec_bin);
    %     binvec_true(ivec_bin) = sig2bini;
    % end
end

% % Various initialization
% beta_hat                                = zeros(size(X,2), size(ymat,2), class(ymat));
% [beta_se, zmat, ymat_hat, ymat_res]     = deal(zeros(size(beta_hat), class(ymat)));
% [betacon_hat, betacon_se]               = deal(zeros(size(contrasts,1), size(ymat,2), class(ymat)));
% binvec                                  = NaN(1, size(ymat,2));
% Randombeta                              = cell(1, size(ymat,2));

if Hessflag
    Hessmat = NaN([num_RFX num_RFX num_y]);
else
    Hessmat = [];
end

for permi = 0:nperms
    if permi == 0
        ymat_current = ymat;
    else
        beta_hat_null = zeros(size(X,2), 1);
        if ~exist('sig2mat_null', 'var')
            X_null = ones(num_obs,1);
            [~, ~, ~, ~, ~, sig2mat_null] = FEMA_fit_binary(X_null, iid, eid, fid, ...
                                                            agevec, ymat, maxIter, ...
                                                            [], [], [], ...
                                                            'RandomEffects', {'F','S','E'}, ...
                                                            'RandomEstType','MoM');
        end

        ymat_current = generate_null_sample(X, iid, fid, agevec, beta_hat_null, ...
                                            sig2mat_null, clusterinfo, RandomEffects);
    end
    
    beta_hat_current                = zeros(num_X, 1, precision);
    [beta_se_current, zmat_current] = deal(zeros(size(beta_hat_current), precision));
    % [betacon_hat_current, betacon_se_current] = deal(zeros(size(contrasts,1), 1, class(ymat_current)));
    
    %% Initialization
    % record should be deleted after
    beta_record     = zeros(maxIter, num_X, num_y); % used for convergence check
    sigmat_record   = zeros(maxIter, num_RFX, num_y);
    deviance_record = zeros(maxIter, num_y);

    iter = 1;
    
    converged = false;

    while ~converged && (iter <= maxIter)
        if iter == 1
            % --- Initialization (iter = 1) ---
            % initialization of working response and working weights
            u   = log((ymat_current + 0.5) ./ (1.5 - ymat_current));
            W   = (ymat_current + 0.5) .* (1.5 - ymat_current) ./4;
            W_1 = 1 ./ W; % avoid redundant calculation
        
            % Initially use OLS estimate
            XtX = X' * X;
            if lowRank
                if useLSQ
                    iXtX = lsqminnorm(XtX, eye(size(XtX)));
                else
                    iXtX = pinv(X);
                end
            else
                iXtX     = XtX \ eye(size(XtX));
            end
            df               = (num_obs - num_X); 
            beta_hat_current = iXtX * (X' * u);
            u_res            = u - X * beta_hat_current;
            sig2tvec_current = sum(u_res.^2,1)/df;

            p_marginal       =  1 ./ (1 + exp(-X * beta_hat_current));
            r                = (ymat_current - p_marginal)./(p_marginal.*(1-p_marginal));
        
            % Compute random variances (only fid and iid)
            [~, sig2mat_current] =    ...
             FEMA_fit_simplified(X, iid, eid, fid, u_res, sig2tvec_current,  ...
                                 GRM, W_1, 'MLflag', MLflag, 'FamilyStruct', ...
                                 FamilyStruct, 'NonnegFlag', NonnegFlag);
        
            % GLS updating fixed effects estimates
            [allWsTerms, beta_hat_current, beta_cov_current] =   ...
             FEMA_GLS(u, X, W_1, sig2mat_current, RandomEffects, ...
                      clusterinfo, nfamtypes, famtypevec,        ...
                      'GroupByFamType', GroupByFamType, 'useLSQ', useLSQ);

            % calculate the probability according to the updated parameters
            [u_update, prob, deviance] = compute_BLUP(u, X, beta_hat_current, sig2mat_current, allWsTerms, ...
                                                   RandomEffects, RandomVar, ymat_current);
            deviance_record(iter,:) = deviance;
    
            % record the parameter estimation
            beta_record(iter,:,:)   = beta_hat_current;
            sigmat_record(iter,:,:) = sig2mat_current;
        
        else % for iter>0
            step = 1;
            previous_deviance = deviance_record(iter-1,:);

            % record old state
            u_old          = u;
            u_update_old   = u_update;
            prob_old       = prob;
            W_1_old        = W_1;
            beta_old       = beta_hat_current;
            beta_cov_old   = beta_cov_current;
            sig2mat_old    = sig2mat_current;
            sig2tvec_old   = sig2tvec_current;
            allWsTerms_old = allWsTerms;

            W   = prob .* (1 - prob);
            W_1 = 1 ./ W;

            % Firth penalization
            U        = (allWsTerms * X) * beta_cov_current;
            h_diag   = sum(U .* X, 2);
            bias_adj = h_diag .* (0.5 - prob);

            % Line searching
            while true

                u = u_update_old + step * (ymat_current - prob_old + bias_adj) ./ (prob_old .* (1-prob_old));
    
                % reupdate the beta_hat and sig2mat
                [~, beta_hat_current, beta_cov_current] =                     ...
                 FEMA_GLS(u, X, W_1, sig2mat_old, RandomEffects, clusterinfo, ...
                          nfamtypes, famtypevec, 'allWsTerms', allWsTerms,    ...
                          'GroupByFamType', GroupByFamType, 'useLSQ', useLSQ);

                u_res            = u - (X * beta_hat_current);
                sig2tvec_current = sum(u_res.^2,1)/df; 

                [~, sig2mat_current] = FEMA_fit_simplified(X, iid, eid, fid,             ...
                                                           u_res, sig2tvec_current,      ...
                                                           GRM, W_1, 'MLflag', MLflag,   ...
                                                           'FamilyStruct', FamilyStruct, ...
                                                           'NonnegFlag', NonnegFlag);
    
                [allWsTerms, ~, ~] = FEMA_GLS(u, X, W_1, sig2mat_current, RandomEffects, ...
                                              clusterinfo, nfamtypes, famtypevec, ...
                                              'GroupByFamType', GroupByFamType, 'useLSQ', useLSQ, ...
                                              'GLSflag',false);
    
                % calculate the probability according to the updated parameters
                [u_update, prob, deviance] = compute_BLUP(u, X, beta_hat_current, ...
                                                          sig2mat_current, allWsTerms, ...
                                                          RandomEffects, RandomVar, ymat_current);
    
                % if all(deviance < previous_deviance)
                if all(deviance < previous_deviance - 1e-3 * 2 * step * sum((ymat_current - prob_old).^2 ./ (prob_old .* (1-prob_old)),1))
                % if true
                    % accept current step
                    deviance_record(iter,:) = deviance;

                    break;
    
                else
    
                    step = step * 0.5;
    
                    if step < 1e-2
                        % take the (t-1) state
                        u                       = u_old;
                        u_update                = u_update_old;
                        prob                    = prob_old;
                        W_1                     = W_1_old;
                        beta_hat_current        = beta_old;
                        beta_cov_current        = beta_cov_old;
                        sig2mat_current         = sig2mat_old;
                        sig2tvec_current        = sig2tvec_old;
                        allWsTerms              = allWsTerms_old;
                        deviance_record(iter,:) = previous_deviance;
                        break;
                    end
                end
            end
     
            % Convergence Check
            beta_change = norm(beta_hat_current - beta_old) / (norm(beta_old) + tol);
            sig_change  = norm(sig2mat_current - sig2mat_old) / (norm(sig2mat_old) + tol);
            if beta_change <= tol && sig_change <= tol
                converged = true;
            end
    
            beta_record(iter,:,:) = beta_hat_current;
            sigmat_record(iter,:,:) = sig2mat_current;
        
        end
        
        iter = iter + 1;
    
    end

    % if converged
    %     fprintf('Converged at iteration %d.\n', iter-1);
    % else
    %     fprintf('Maximum iterations (%d) reached.\n', maxIter);
    % end


    % Laplace correction for all random effects
    % Currently resorting to loop because multiple things in
    % laplace_correct will need to be changed to make this compatible with
    % multiple y variables
    coeffCovar     = zeros(num_X, num_X, num_y, precision);
    beta_se_robust = zeros(num_X, num_y, precision);
    
    for yy = 1:num_y
        [sig2mat_lap, intercept_lap] = laplace_correct(ymat_current(:,yy), X, beta_hat_current(:,yy), ...
                                                       RandomEffects, RandomVar,          ...
                                                       sig2mat_current(:,yy), IC_fam, IC_subj);

        sig2mat_current(1:end-1,yy) = sig2mat_lap;
        beta_hat_current(1,yy)      = intercept_lap;
    
        % degree of freedom for random effects
        sigma2_total = sum(sig2mat_current(1:end-1,yy), 1);
        W            = prob(:,yy) .* (1-prob(:,yy)); 
        [G, ~]       = findgroups(fid);
        W_sum        = splitapply(@sum, W, G')';
        h_edf        = sum(sigma2_total * W_sum ./ (sigma2_total * W_sum + 1));
    
        % sandwich estimates of beta_se
        % S0 = X' * diag(1 ./ W) * X;

        % PP: Are we sure of the calculation here?
        % Here is an reference: https://stats.stackexchange.com/a/284939
        % solve(t(X) %*% X) %*% t(X) %*% diag(r^2) %*% X %*% solve(t(X) %*% X)
        % 
        % Part 1): Note that X has been modified! X_mod = X * sqrt(W)
        % 
        % Part 2): r are the working residuals: r = (y - yHat) * (d(eta)/d(mu))
        % https://www.rdocumentation.org/packages/binomTools/versions/1.0-1/topics/Residuals
        % eta and mu are calculated at the end of IRLS:
        % https://github.com/wch/r-source/blob/84961436a6409c44ef226c0fadb16f5f84bc2b40/src/library/stats/R/glm.R#L263
        % 
        % Now, Google Gemini insists that the d(mu)/d(eta) simplifies to
        % p(1-p), and therefore, the working residuals should be: 
        % (y - p)/(p * (1-p))
        % 
        % See, here, for a reference on how to compute these residuals:
        % https://stats.stackexchange.com/a/485734
        % 
        % Using the code from above, 
        % predictedValue = 1/(1 + e^-z), where z = X * beta
        % mu = exp(predictedValue)/(1+exp(predictedValue))
        % workingResiduals = (Y-mu) / (mu*(1-mu))
        % (Y-exp(predictedValue)/(1+exp(predictedValue))) / (exp(predictedValue)/(1+exp(predictedValue))*(1-exp(predictedValue)/(1+exp(predictedValue))))
        % (Y-exp(1/(1 + e^-z))/(1+exp(1/(1 + e^-z)))) / (exp(1/(1 + e^-z))/(1+exp(1/(1 + e^-z)))*(1-exp(1/(1 + e^-z))/(1+exp(1/(1 + e^-z)))))
        % (Y-exp(1/(1 + e^-(X * beta)))/(1+exp(1/(1 + e^-(X * beta))))) / (exp(1/(1 + e^-(X * beta)))/(1+exp(1/(1 + e^-(X * beta))))*(1-exp(1/(1 + e^-(X * beta)))/(1+exp(1/(1 + e^-(X * beta))))))
        % Asking Google Gemini to simplify the above gives me (y - p)/(p * (1-p))
        % 
        % So, r = (y - p)/(p * (1-p)), needs to be scaled by sqrt(weights)
        % r_mod = r * sqrt(W)
        % r_mod = ((y - p)/(p * (1-p))) * sqrt(W)
        % 
        % Therefore, practically
        % solve(t(X) %*% X) %*% t(X) %*% diag(r^2) %*% X %*% solve(t(X) %*% X)
        % becomes
        % solve(t(X_mod) %*% X_mod) %*% t(X_mod) %*% diag(r_mod^2) %*% X_mod %*% solve(t(X_mod) %*% X_mod)
        % Solve is just a way to compute the inverse, so
        % inv(t(X_mod) * X_mod) * t(X_mod) * diag(r_mod^2) * X_mod * inv(t(X_mod) * X_mod)
        % inv(X_mod' * X_mod) * X_mod' * diag(r_mod^2) * X_mod * inv(X_mod' * X_mod)
        % inv((X * sqrt(W))' * (X * sqrt(W))) * (X * sqrt(W))' * diag(r_mod^2) * (X * sqrt(W)) * inv((X * sqrt(W))' * (X * sqrt(W)))
        % inv((X * sqrt(W))' * (X * sqrt(W))) * (X * sqrt(W))' * diag((((y - p)/(p * (1-p))) * sqrt(W))^2) * (X * sqrt(W)) * inv((X * sqrt(W))' * (X * sqrt(W)))
        
        % Separate weights
        W_bread    = prob(:,yy) .* (1-prob(:,yy));
        W_meat     = (ymat_current(:,yy) - prob(:,yy)).^2;

        % Compute bread and meat
        bread = pinv(X' * (W_bread .* X));
        meat  = X' * (W_meat .* X);

        % V sandwich
        V_sandwich = bread * meat * bread;
        
        % Older solution: not sure if the above is entirely correct
        % Cannot use iXtX directly since it is not weighted
        % Do we need the sigma2_total weighting? Isn't that sort of mean
        % squared error which is not used in logistic regression?
        % S0                     = X' * (X ./ W);
        % V1                     = iXtX * S0 * iXtX;
        % V2                     = sigma2_total * iXtX;
        % [coeffCovar, converge] = nearestSPD_timeout(V1 + V2);
        [tmpCoeff, converge] = nearestSPD_timeout(V_sandwich);
        if ~converge
            warning('Could not convert coefficient covariance matrix into positive semidefinite; results might be inaccurate');
            coeffCovar(:,:,yy) = V_sandwich;
        else
            V_sandwich = tmpCoeff;
            coeffCovar(:,:,yy) = tmpCoeff;
        end
        beta_se_robust(:,yy) = sqrt(diag(V_sandwich) * num_obs / (num_obs-num_X-h_edf));
    end

    % z-statistics
    zmat_current    = double(beta_hat_current) ./ double(beta_se_robust);
    logpmat_current = -log10(normcdf(-abs(zmat_current))*2);

    if permi == 0
        beta_hat = beta_hat_current;
        beta_se  = beta_se_robust;
        zmat     = zmat_current;
        logpmat  = logpmat_current;
        sig2tvec = sig2tvec_current;
        sig2mat  = sig2mat_current;
        if logLikflag
            logLikvec = logLikvec_current;
        end
    else
        beta_hat_perm(:,:,permi) = beta_hat_current;
        beta_se_perm(:,:,permi)  = beta_se_robust;
        zmat_perm(:,:,permi)     = zmat_current;
        sig2tvec_perm(:,:,permi) = sig2tvec_current;
        sig2mat_perm(:,:,permi)  = sig2mat_current;
        if logLikflag
            logLikvec_perm(:,:,permi) = logLikvec_current;
        end
    end

    % after_mem = memory().MemUsedMATLAB / 1024^2;
    % memory_usage = after_mem - before;

    % fprintf("Memory used: %.2f MB\n", after_mem-before);
    
    % fprintf('Permutation at permi %d. Converged at iteration %d.\n', [permi, iter-1]);
end