function [beta_hat,      beta_se,        zmat,        logpmat,          ...
          sig2tvec,      sig2mat,        Hessmat,     logLikvec,        ...
          beta_hat_perm, beta_se_perm,   zmat_perm,   sig2tvec_perm,    ...
          sig2mat_perm,  logLikvec_perm, binvec_save, nvec_bins,        ...
          tvec_bins,     FamilyStruct,   coeffCovar,  unstructParams,   ...
          residuals_GLS, info] = FEMA_fit_binary(X, iid, eid, fid,      ...
                                                 agevec, ymat, contrasts, ...
                                                 nbins, GRM, varargin)
% Function to fit fast and efficient linear mixed effects model
%
% For notation below:
% n = observations,
% p = predictors (fixed effects),
% v = number of outcome variables
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
% maxIter         <num>            [1 x 1]    maximal number of IRLS iterations (default 200) 
% contrasts       <num> OR <path>  [c x p]    contrast matrix, where c is number of contrasts to compute,
%                                             OR path to file containing contrast matrix (readable by readtable)
% nbins           <num>            [1 x 1]    number of covariance bins (default: 20)
%


%% Optional input arguments:
% RandomEffects   <cell>           list of random effects to estimate (default {'F','S'}):
%                                       * F:  family relatedness
%                                       * S:  subject - required for longitudinal analyses
%                                       * A:  additive genetic relatedness - must include file path to genetic relatedness data (GRM) for this option
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
% precision       <char>           default 'double' --> other option: 'single' - for precision
% max_pages       <num>            default 100000 - maximum number of matrix pages processed together
% logLikflag      <boolean>        default true - compute log-likelihood
% PermType        <char>           permutation type:
%                                       * 'wildbootstrap':    residual boostrap --> creates null distribution by randomly flipping the sign of each observation
%                                       * 'wildbootstrap-nn': non-null boostrap --> estimates distribution around effect of interest using sign flipping (used for sobel test)
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
addParamValue(p,'max_pages', 100000);
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
max_pages            = p.Results.max_pages;
RandomEffects        = p.Results.RandomEffects;
MoMflag              = ismember(lower(RandomEstType), {'mom'});
MLflag               = ismember(lower(RandomEstType), {'ml'});
logLikflag           = p.Results.logLikflag;
Hessflag             = p.Results.Hessflag;
ciflag               = p.Results.ciflag;
nperms               = p.Results.nperms;
PermType             = lower(p.Results.PermType);
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
    beta_hat_perm = zeros(num_X, num_y, nperms, precision);
    beta_se_perm  = zeros(num_X, num_y, nperms, precision);
    zmat_perm     = zeros(num_X, num_y, nperms, precision);
    sig2tvec_perm = zeros(1, num_y, nperms, precision);
    sig2mat_perm  = zeros(num_RFX, num_y, nperms, precision);
    if logLikflag
        logLikvec_perm = zeros(1, num_y, nperms, precision);
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
    
    nfam      = length(unique(fid));
    nfamtypes = length(famtypelist);

    % Save some more information
    info.nFamilies = length(clusterinfo);
    info.nFamTypes = nfamtypes;
    
    % Prepare generalized matrix version of MoM estimator
    S_sum = Ss{1};
    for i = 2:length(Ss)
        S_sum = S_sum + Ss{i};
    end
    [subvec1, subvec2] = find(S_sum); % Use full matrix, to simplify IGLS -- should be possible to limit to tril
    %[subvec1 subvec2] = find(tril(S_sum)); % Should exclude diagonals: tril(S_sum,-1)
    indvec = sub2ind([num_obs num_obs],subvec1,subvec2);

    fam_of_obs = zeros(num_obs, 1);
    for fi = 1:nfam
        fam_of_obs(clusterinfo{fi}.jvec_fam) = fi;
    end
    fnumvec_d = fam_of_obs(subvec1);

    if ~all(fnumvec_d >= 1)
        F_num = S_sum;
        for fi = 1:nfam
            F_num(clusterinfo{fi}.jvec_fam,clusterinfo{fi}.jvec_fam) = fi;
        end
        fnumvec = F_num(indvec);
    
        for fi = 1:nfam
            jvec_tmp  = clusterinfo{fi}.jvec_fam;
            [sv, si]  = sort(jvec_tmp);
            I_tmp     = reshape(1:length(jvec_tmp)^2, length(jvec_tmp) * [1 1]);
            ivec_fam  = find(fnumvec==fi);
            ivec_fam  = ivec_fam(colvec(I_tmp(si, si)));
            %  ivec_fam = find(fnumvec==fi); ivec_fam(colvec(I_tmp(si,si))) = ivec_fam;
            clusterinfo{fi}.ivec_fam = ivec_fam;
        end
    else
        npos   = length(fnumvec_d);
        [~, isort] = sort(fnumvec_d * (npos + 1) + (1:npos)');
        counts = accumarray(fnumvec_d, 1, [nfam 1]);
        bounds = [0; cumsum(counts)];
        for fi = 1:nfam
            jvec_tmp  = clusterinfo{fi}.jvec_fam;
            [~, si]   = sort(jvec_tmp);
            I_tmp     = reshape(1:length(jvec_tmp)^2, length(jvec_tmp) * [1 1]);
            ivec_fam  = isort(bounds(fi)+1 : bounds(fi+1));
            clusterinfo{fi}.ivec_fam = ivec_fam(colvec(I_tmp(si, si)));
        end
    end

    M = zeros(length(indvec),length(Ss));
    for i = 1:length(Ss)
        M(:,i) = Ss{i}(indvec);
    end

    % Create grid of normalized random effects
    binvals_edges       = linspace(0,1,nbins+1); 

    sig2gridi = ndgrid_amd(repmat({1:nbins}, [1 num_RFX]));
    sig2gridl = ndgrid_amd(repmat({binvals_edges(1:end-1)},    [1 num_RFX]));
    sig2gridu = ndgrid_amd(repmat({binvals_edges(2:end)},      [1 num_RFX]));
    sig2grid      = (sig2gridl+sig2gridu)/2;
    sig2gridind   = sub2ind_amd(nbins*ones(1,num_RFX),sig2gridi);
    nsig2bins     = size(sig2gridl,1); % Should handle case of no binning

    % Prepare FamilyStruct
    FamilyStruct = struct('clusterinfo', {clusterinfo}, 'M', {M},                     ...
                          'famtypevec',  {famtypevec},  'famtypelist', {famtypelist}, ...
                          'nfamtypes',   nfamtypes,     'nfam', nfam,                 ...
                          'sig2grid',    sig2grid,      'sig2gridl', sig2gridl,       ...
                          'sig2gridu',   sig2gridu,     'sig2gridi', sig2gridi,       ...
                          'sig2gridind', sig2gridind,   'nsig2bins', nsig2bins,       ...
                          'subvec1',     subvec1,       'subvec2', subvec2,           ...
                          'Ss',          {Ss});
    
    info.timing.tParseFamily = toc(tInit_parseFamily);
else
    clusterinfo = FamilyStruct.clusterinfo;
    M           = FamilyStruct.M;
    nsig2bins   = FamilyStruct.nsig2bins;
    nfam        = FamilyStruct.nfam; % nfam defined by the fid? duplicate -- no duplicate
    famtypevec  = FamilyStruct.famtypevec;
    nfamtypes   = FamilyStruct.nfamtypes;
    sig2grid    = FamilyStruct.sig2grid;
    sig2gridl   = FamilyStruct.sig2gridl;
    sig2gridu   = FamilyStruct.sig2gridu;
    subvec1     = FamilyStruct.subvec1;
    subvec2     = FamilyStruct.subvec2;
    Ss          = FamilyStruct.Ss;

end

% Should be generalized ------- later
numUqSubjs = info.nUqSubjects;
[~, ~, IC_subj] = unique(iid,'stable');
[~, ~, IC_fam]  = unique(fid,'stable');

% precompute iteration-invariant GLS sparse-assembly coordinates once
GLSCache = FEMA_GLS_make_cache(clusterinfo, RandomEffects, nfamtypes, famtypevec, GroupByFamType);

RandomVar = struct();
RandomVar.("V_F") = sparse(1:num_obs, IC_fam, ones(num_obs,1),  num_obs, nfam);
RandomVar.("V_S") = sparse(1:num_obs, IC_subj, ones(num_obs,1), num_obs, numUqSubjs);

Mi = single(pinv(M));
Cov_MoM = Mi*Mi'; % Variance  / covariance of MoM estimates, per unit of residual error variance

% logging('size(M) = [%d %d]',size(M));
% logging('Cov_MoM:'); disp(Cov_MoM);
% logging('Mi*M:'); disp(Mi*M);

if ~isempty(synthstruct) % do we need this anymore?----------------------------------
    sig2mat_true  = synthstruct.sig2mat_true;
    sig2tvec_true = synthstruct.sig2tvec_true;
    nvec_bins_true = NaN(nsig2bins,1);
    binvec_true    = NaN(1,size(ymat,2));
    for sig2bini = 1:nsig2bins
        tmpvec = true;
        for ri = 1:size(sig2mat_true,1)-1
            tmpvec = tmpvec & sig2mat_true(ri,:) >= sig2gridl(sig2bini,ri) & ...
                              sig2mat_true(ri,:) <  sig2gridu(sig2bini,ri);
        end
        ivec_bin = find(tmpvec);
        nvec_bins_true(sig2bini) = length(ivec_bin);
        binvec_true(ivec_bin) = sig2bini;
    end
end

if Hessflag
    Hessmat = NaN([num_RFX num_RFX num_y]);
else
    Hessmat = [];
end

for permi = 0:nperms
    if permi == 0
        ymat_current = ymat;
    else
        beta_hat_null = zeros(num_X, num_y, precision);
        if ~exist('sig2mat_null', 'var')
            X_null = ones(num_obs,1);
            [~, ~, ~, ~, ~, sig2mat_null] = FEMA_fit_binary(X_null, iid, eid, fid, ...
                                                            agevec, ymat, [], 0, [], ...
                                                            'RandomEffects', RandomEffects, ...
                                                            'RandomEstType','MoM','AddIntercept', false);
        end

        ymat_current = zeros(size(ymat), class(ymat));
        for yy = 1:num_y
            ymat_current(:,yy) = generate_null_sample( ...
                X, iid, fid, agevec, beta_hat_null(:,yy), ...
                sig2mat_null(:,yy), clusterinfo, RandomEffects);
        end
    end
    
    beta_hat_current                = zeros(num_X, num_y, precision);
    betacon_hat_current = zeros(size(contrasts,1), num_y, precision);
    
    %% Initialization
    % record should be deleted after --------------------------------------
    beta_record     = zeros(maxIter, num_X, num_y); % used for convergence check
    sigmat_record   = zeros(maxIter, num_RFX, num_y);
    deviance_record = zeros(maxIter, num_y);

    % --- Initialization ---
    % initialization of working response and working weights
    u   = log((ymat_current + 0.5) ./ (1.5 - ymat_current));
    W   = (ymat_current + 0.5) .* (1.5 - ymat_current) ./ 4;
    W_1 = 1 ./ W;

    % Initially use OLS estimate
    XtX = X' * X;
    if lowRank
        if useLSQ
            iXtX = lsqminnorm(XtX, eye(size(XtX)));
        else
            iXtX = pinv(XtX);
        end
    else
        iXtX = XtX \ eye(size(XtX));
    end

    df               = num_obs - num_X; 
    beta_hat_current = iXtX * (X' * u);
    u_res            = u - X * beta_hat_current;
    sig2tvec_current = sum(u_res.^2, 1) / df;

    % Compute random variances
    [~, sig2mat_current] = FEMA_fit_simplified(X, iid, eid, fid, u_res, sig2tvec_current,  ...
                                               GRM, W_1, 'MLflag', MLflag, 'FamilyStruct', ...
                                               FamilyStruct, 'NonnegFlag', NonnegFlag);

    % Snap to random effects grid
    [binvec_current, nvec_bins_current] = FEMA_snap_to_sig2grid(sig2mat_current, nbins, ...
                                                                FamilyStruct, W_1);


    % GLS updating fixed effects estimates
    [allWsTerms, beta_hat_current, beta_cov_current, tvec_bins_current] =   ...
     FEMA_GLS(u, X, W_1, sig2mat_current, RandomEffects, clusterinfo, ...
              nfamtypes, famtypevec, 'GroupByFamType', GroupByFamType, ...
              'precision', precision, 'useLSQ', useLSQ, 'max_pages', max_pages, ...
              'nbins', nbins, 'binvec', binvec_current, 'GLSCache', GLSCache);


    %Restore the absolute working scale for the coefficient covariance
    beta_cov_current = beta_cov_current .* reshape(sig2tvec_current, [1 1 num_y]);

    % calculate the probability according to the updated parameters
    [u_update, prob, deviance] = compute_BLUP(u, X, beta_hat_current, sig2mat_current, allWsTerms, ...
                                           RandomEffects, RandomVar, ymat_current);

    % record the parameter estimation---------------------------------------
    beta_record(1,:,:)   = reshape(beta_hat_current, [1 num_X num_y]);
    sigmat_record(1,:,:) = reshape(sig2mat_current, [1 num_RFX num_y]);
    deviance_record(1,:) = deviance;
    converged            = false(1, num_y);
    
    %% Iterate PQL/IRLS updates
    for iter = 2:maxIter
        ivec_active = find(~converged);
        if isempty(ivec_active)
            break;
        end

        previous_deviance = deviance(ivec_active);

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

        W_1(:,ivec_active) = 1 ./ (prob_old(:,ivec_active) .* (1 - prob_old(:,ivec_active)));

        % Firth penalization
        bias_adj = zeros(num_obs, num_y, precision);
        for yy = ivec_active
            PX             = allWsTerms_old{yy} * X;
            U              = PX * (beta_cov_old(:,:,yy) / sig2tvec_old(yy));
            h_diag         = sum(U .* X, 2);
            bias_adj(:,yy) = h_diag .* (0.5 - prob_old(:,yy));
        end

        step         = ones(1, num_y);
        ivec_pending = ivec_active;

        % Line searching
        while ~isempty(ivec_pending)
            u_trial = u_update_old(:,ivec_pending) + ...
                (ymat_current(:,ivec_pending) - prob_old(:,ivec_pending) + bias_adj(:,ivec_pending)) ./ ...
                (prob_old(:,ivec_pending) .* (1 - prob_old(:,ivec_pending))) .* step(ivec_pending);

            % Snap to random effects grid before the GLS update
            binvec_old = FEMA_snap_to_sig2grid(sig2mat_old(:,ivec_pending), ...
                                               nbins, FamilyStruct, W_1(:,ivec_pending));

            % reupdate the beta_hat and sig2mat
            [~, beta_trial, beta_cov_trial] = FEMA_GLS(u_trial, X, W_1(:,ivec_pending), ...
                                                       sig2mat_old(:,ivec_pending), RandomEffects, ...
                                                       clusterinfo, nfamtypes, famtypevec, ...
                                                       'GroupByFamType', GroupByFamType, ...
                                                       'precision', precision, 'useLSQ', useLSQ, 'max_pages', max_pages, ...
                                                       'binvec', binvec_old, 'GLSCache', GLSCache, 'nbins', nbins);
            
            beta_cov_trial = beta_cov_trial .* reshape(sig2tvec_old(ivec_pending), [1 1 length(ivec_pending)]);

            u_res_trial    = u_trial - X * beta_trial;
            sig2tvec_trial = sum(u_res_trial.^2, 1) / df;

            [~, sig2mat_trial] = FEMA_fit_simplified(X, iid, eid, fid, u_res_trial, ...
                                                     sig2tvec_trial, GRM, W_1(:,ivec_pending), ...
                                                     'MLflag', MLflag, 'FamilyStruct', FamilyStruct, ...
                                                     'NonnegFlag', NonnegFlag);

            W_1_trial = W_1(:,ivec_pending);

            % Snap to random effects grid again with the updated random
            % effects before rebuilding the precision matrices
            binvec_trial = FEMA_snap_to_sig2grid(sig2mat_trial, nbins, FamilyStruct, W_1_trial);
            [allWs_trial, ~, ~] = FEMA_GLS(u_trial, X, W_1_trial, sig2mat_trial, ...
                                           RandomEffects, clusterinfo, nfamtypes, famtypevec, ...
                                           'GroupByFamType', GroupByFamType, 'precision', precision, ...
                                           'useLSQ', useLSQ, 'max_pages', max_pages, 'GLSflag', false, ...
                                           'binvec', binvec_trial, 'GLSCache', GLSCache, 'nbins', nbins);

            % calculate the probability according to the updated parameters
            [u_update_trial, prob_trial, deviance_trial] = compute_BLUP(u_trial, X, beta_trial, ...
                                                                        sig2mat_trial, allWs_trial, ...
                                                                        RandomEffects, RandomVar, ...
                                                                        ymat_current(:,ivec_pending));

            descent = 1e-3 * 2 .* step(ivec_pending) .* ...
                      sum((ymat_current(:,ivec_pending) - prob_old(:,ivec_pending)).^2 ./ ...
                      (prob_old(:,ivec_pending) .* (1 - prob_old(:,ivec_pending))), 1);
            [~, loc_previous] = ismember(ivec_pending, ivec_active);
            accept_local = deviance_trial < ...
                           previous_deviance(loc_previous) - descent;

            ivec_accept = ivec_pending(accept_local);
            loc_accept  = find(accept_local);
            for jj = 1:length(ivec_accept)
                yy = ivec_accept(jj);
                kk = loc_accept(jj);
                u(:,yy)                    = u_trial(:,kk);
                W_1(:,yy)                  = W_1_trial(:,kk);
                beta_hat_current(:,yy)     = beta_trial(:,kk);
                beta_cov_current(:,:,yy)   = beta_cov_trial(:,:,kk);
                sig2mat_current(:,yy)      = sig2mat_trial(:,kk);
                sig2tvec_current(yy)       = sig2tvec_trial(kk);
                allWsTerms{yy}             = allWs_trial{kk};
                u_update(:,yy)             = u_update_trial(:,kk);
                prob(:,yy)                 = prob_trial(:,kk);
                deviance(yy)               = deviance_trial(kk);
            end

            ivec_reject = ivec_pending(~accept_local);
            step(ivec_reject) = step(ivec_reject) * 0.5;

            % Match the original single-y fallback: when the step becomes too small, retain the complete previous state and stop trying.
            ivec_stop = ivec_reject(step(ivec_reject) < 1e-2);
            for yy = ivec_stop
                u(:,yy)                  = u_old(:,yy);
                u_update(:,yy)           = u_update_old(:,yy);
                prob(:,yy)               = prob_old(:,yy);
                W_1(:,yy)                = W_1_old(:,yy);
                beta_hat_current(:,yy)   = beta_old(:,yy);
                beta_cov_current(:,:,yy) = beta_cov_old(:,:,yy);
                sig2mat_current(:,yy)    = sig2mat_old(:,yy);
                sig2tvec_current(yy)     = sig2tvec_old(yy);
                allWsTerms{yy}           = allWsTerms_old{yy};
                deviance(yy)             = previous_deviance(ivec_active == yy);
            end

            ivec_pending = setdiff(ivec_reject, ivec_stop, 'stable');
        end

        beta_change = zeros(1, length(ivec_active));
        sig_change  = zeros(1, length(ivec_active));
        for jj = 1:length(ivec_active)
            yy = ivec_active(jj);
            beta_change(jj) = norm(beta_hat_current(:,yy) - beta_old(:,yy)) / ...
                              (norm(beta_old(:,yy)) + tol);
        end
        converged(ivec_active) = beta_change <= tol;

        beta_record(iter,:,:)   = reshape(beta_hat_current, [1 num_X num_y]);
        sigmat_record(iter,:,:) = reshape(sig2mat_current, [1 num_RFX num_y]);
        deviance_record(iter,:) = deviance;
        iter_last               = iter;
    end

    % Snap to the random-effects grid
    [binvec_current, nvec_bins_current] = ...
     FEMA_snap_to_sig2grid(sig2mat_current, nbins, FamilyStruct, W_1);

    if nbins == 0
         [allWsTerms, ~, beta_cov_current, tvec_bins_current] = ...
         FEMA_GLS(u, X, W_1, sig2mat_current, RandomEffects, ...
                  clusterinfo, nfamtypes, famtypevec, ...
                  'GroupByFamType', GroupByFamType, 'precision', precision, ...
                  'useLSQ', useLSQ, 'max_pages', max_pages, ...
                  'GLSflag', true, 'binvec', binvec_current, ...
                  'allWsTerms', allWsTerms, 'nbins', nbins);
    else     
        [allWsTerms, ~, beta_cov_current, tvec_bins_current] = ...
         FEMA_GLS(u, X, W_1, sig2mat_current, RandomEffects, ...
                  clusterinfo, nfamtypes, famtypevec, ...
                  'GroupByFamType', GroupByFamType, 'precision', precision, ...
                  'useLSQ', useLSQ, 'max_pages', max_pages, ...
                  'GLSflag', true, 'binvec', binvec_current, ...
                  'GLSCache', GLSCache, 'nbins', nbins);
    end

    beta_cov_current = beta_cov_current .* reshape(sig2tvec_current, [1 1 num_y]);

    residuals_GLS_current = u - X * beta_hat_current;


    % Laplace correction for all random effects
    % Currently resorting to loop because multiple things in
    % laplace_correct will need to be changed to make this compatible with
    % multiple y variables
    coeffCovar     = beta_cov_current;

    diag_idx = (1:num_X+1:num_X^2).' + (0:num_y-1) * num_X^2;

    beta_se = sqrt(coeffCovar(diag_idx));

    eta_fixed = X * beta_hat_current;

    [sig2mat_lap, intercept_lap] = laplace_correct(ymat_current, X, beta_hat_current, ...
                                                   eta_fixed, RandomEffects, RandomVar, ...
                                                   sig2tvec_current.* sig2mat_current, IC_fam, IC_subj);

    sig2mat_current(:,:)  = [sig2mat_lap; ones(1, num_y)];
    beta_hat_current(1,:) = intercept_lap;


    % Evaluate contrasts after all coefficient corrections.
    if isempty(contrasts)
        beta_hat_current_out = beta_hat_current;
        beta_se_current_out  = beta_se;
    else
        betacon_hat = contrasts * beta_hat_current;
        betacon_se  = zeros(size(contrasts,1), num_y, precision);
        for yy = 1:num_y
            for ci = 1:size(contrasts,1)
                betacon_se(ci,yy) = sqrt(contrasts(ci,:) * ...
                    coeffCovar(:,:,yy) * contrasts(ci,:)');
            end
        end
        beta_hat_current_out = [betacon_hat; beta_hat_current];
        beta_se_current_out  = [betacon_se; beta_se];
    end

    % z-statistics
    zmat_current    = double(beta_hat_current_out) ./ double(beta_se_current_out);
    logpmat_current = -log10(normcdf(-abs(zmat_current))*2);

    if logLikflag
        logLikvec_current = NaN(1, num_y, precision);
    end

    if permi == 0
        beta_hat = beta_hat_current_out;
        beta_se  = beta_se_current_out;
        zmat     = zmat_current;
        logpmat  = logpmat_current;
        sig2tvec = sig2tvec_current;
        sig2mat  = sig2mat_current;
        binvec_save = binvec_current;
        nvec_bins   = nvec_bins_current;
        tvec_bins   = tvec_bins_current;
        if returnResiduals
            residuals_GLS = residuals_GLS_current;
        end
        if logLikflag
            logLikvec = logLikvec_current;
        end
        info.converged           = converged;
        info.numIterations       = iter_last;
        info.deviance            = deviance;
        info.binning.binvec      = binvec_current;
        info.binning.nvec_bins   = nvec_bins_current;
        info.binning.tvec_bins   = tvec_bins_current;
    else
        beta_hat_perm(:,:,permi) = beta_hat_current_out;
        beta_se_perm(:,:,permi)  = beta_se_current_out;
        zmat_perm(:,:,permi)     = zmat_current;
        sig2tvec_perm(:,:,permi) = sig2tvec_current;
        sig2mat_perm(:,:,permi)  = sig2mat_current;
        if logLikflag
            logLikvec_perm(:,:,permi) = logLikvec_current;
        end
    end
end

info.timing.total = toc(tInit);

end