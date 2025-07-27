function [beta_hat,      beta_se,        zmat,        logpmat,              ...
          sig2tvec,      sig2mat,        Hessmat,     logLikvec,            ...
          beta_hat_perm, beta_se_perm,   zmat_perm,   sig2tvec_perm,        ...
          sig2mat_perm,  logLikvec_perm, binvec_save, nvec_bins,            ...
          tvec_bins,     FamilyStruct,   reusableVars] =                    ...
          FEMA_fit_binary(X, iid, eid, fid, agevec, ymat, niter, contrasts, nbins, ...
                   pihatmat, varargin)

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
% nbins           <num>            [1 x 1]    number of bins across Y for estimating random effects (default 20)
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
% SingleOrDouble  <char>           default 'double' --> other option: 'single' - for precision
% logLikflag      <boolean>        default true - compute log-likelihood
% PermType        <char>           permutation type:
%                                       * 'wildbootstrap':    residual boostrap --> creates null distribution by randomly flipping the sign of each observation
%                                       * 'wildbootstrap-nn': non-null boostrap --> estimates distribution around effect of interest using sign flipping (used for sobel test)
% returnReusable  <boolean>        default false - if true, additionally returns reusableVars as a structure with some variables that can be reused (primarily by FEMA-GWAS)
% MaxIter         <num>            default 200
% tol             <num>            default 1e-4


%% Outputs:
% beta_hat                         [c+p x v]  estimated beta coefficients
% beta_se                          [c+p x v]  estimated beta standard errors
% zmat                             [c+p x v]  z statistics
% logpmat                          [c+p x v]  log10 p-values
% sig2tvec                         [1   x v]  total residual error of model at each vertex/voxel
% sig2mat                          [r   x v]  normalized random effect variances
% binvec_save                      [1   x v]  bin number(s) for non-permuted ymat
% FamilyStruct                                structure type (can be passed as input to avoid re-parsing family structure etc.)
%

%%
starttime = now(); %#ok<*TNOW1>
logging('***Start***');

% Extremely quick sanity check on X and y variables
if logical(sum(any(isnan(X)))) || logical(sum(any(isinf(X))))       ||  ...
   logical(sum(any(isnan(ymat)))) || logical(sum(any(isinf(ymat))))
    error('X and/or ymat have NaN or Inf; please check your data');
end

p = inputParser;

if ~exist('niter', 'var') || isempty(niter)
    niter = 200;
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

if ~exist('nbins','var') || isempty(nbins)
    nbins = 20;
end

if ~exist('pihatmat','var')
    pihatmat = [];
end


% Should change to allow p to be passed in, so as to avoid having to
% duplicate input argument parsing in FEMA_wrapper and FEMA_fit
p = inputParser;
addParamValue(p,'CovType', 'analytic'); %#ok<*NVREPLA>
addParamValue(p,'RandomEstType', 'MoM');
addParamValue(p,'PermType', 'wildbootstrap');
addParamValue(p,'GroupByFamType', true);
addParamValue(p,'NonnegFlag', true); % Perform lsqnonneg on random effects estimation
addParamValue(p,'SingleOrDouble', 'double');
addParamValue(p,'RandomEffects', {'F' 'S' 'E'}); % Default to Family and Subject
addParamValue(p,'logLikflag', false);
addParamValue(p,'Hessflag', false);
addParamValue(p,'ciflag', false);
addParamValue(p,'nperms', 0);
addParamValue(p,'FatherID', {}); % Father ID, ordered same as pihatmat
addParamValue(p,'MotherID', {}); % Mother ID, ordered same as pihatmat
addParamValue(p,'PregID', {}); % Pregnancy effect (same ID means twins), ordered same as pihatmat
addParamValue(p,'HomeID', {}); % Home effect (defined as same address ID), ordered same as pihatmat
addParamValue(p,'FamilyStruct',{}); % Avoids recomputing family strucutre et al
addParamValue(p,'returnReusable',false); % Additionally returns a few useful variables
addParamValue(p,'synthstruct',''); % True / synthesized random effects
addParamValue(p,'MaxIter',200); % Maximum iteration numbers
addParamValue(p,'tol',1e-4); % tolerance set for both fixed and random effects - Line 448   

parse(p,varargin{:})

CovType              = p.Results.CovType; %#ok<*NASGU>
RandomEstType        = p.Results.RandomEstType;
GroupByFamType       = p.Results.GroupByFamType;
NonnegFlag           = p.Results.NonnegFlag;
SingleOrDouble       = p.Results.SingleOrDouble;
RandomEffects        = p.Results.RandomEffects;
MoMflag              = ismember(lower(RandomEstType), {'mom'});
MLflag               = ismember(lower(RandomEstType), {'ml'});
logLikflag           = p.Results.logLikflag;
Hessflag             = p.Results.Hessflag;
ciflag               = p.Results.ciflag;
nperms               = p.Results.nperms;
PermType             = p.Results.PermType;
FamilyStruct         = p.Results.FamilyStruct;
returnReusable       = p.Results.returnReusable;
synthstruct          = p.Results.synthstruct;
MaxIter              = p.Results.MaxIter;
tol                  = p.Results.tol;

% Grouping by family type is only supported for RandomEffects 'F' 'S' 'E'
if ~isempty(setdiff(RandomEffects,{'F' 'S' 'E'}))
    GroupByFamType = false;
end

fprintf(1,'ModelSingularityIndex = %g\n',cond(X'*X)/cond(diag(diag(X'*X))));

[logLikvec,     beta_hat_perm, beta_se_perm,   zmat_perm,   ...
 sig2tvec_perm, sig2mat_perm,  logLikvec_perm, binvec_save, ...
 nvec_bins,     tvec_bins] = deal([]);


if ~returnReusable
    reusableVars = [];
end


% Check if lsqminnorm can be used
if exist('lsqminnorm', 'file')
    useLSQ = true;
else
    useLSQ = false;
end

% Check if X is rank deficient
if rank(X) < size(X, 2)
    lowRank = true;
else
    lowRank = false;
end

% Save some variables for later
if returnReusable
    reusableVars.GroupByFamType = GroupByFamType;
    reusableVars.RandomEffects  = RandomEffects;
    reusableVars.SingleOrDouble = SingleOrDouble;
    reusableVars.useLSQ         = useLSQ;
    reusableVars.lowRank        = lowRank;
end


t0 = now;

RandomVar = cell(1,length(RandomEffects));

% Parse family structure, if necessary
if ~exist('FamilyStruct', 'var') || isempty(FamilyStruct)
    % tic
    [clusterinfo, Ss, iid, famtypevec, famtypelist, subj_famtypevec] =                  ...
     FEMA_parse_family(iid, eid, fid, agevec, pihatmat, 'RandomEffects', RandomEffects, ...
                       'FatherID', p.Results.FatherID,  'MotherID', p.Results.MotherID, ...
                       'PregID',   p.Results.PregID,    'HomeID',   p.Results.HomeID); %#ok<*ASGLU>
    
    [iid_list, IA, IC_subj] = unique(iid,'stable'); nsubj = length(iid_list); nobs = length(iid);
    [fid_list, IA, IC_fam]  = unique(fid,'stable'); nfam  = length(fid_list);
    nfamtypes = length(famtypelist);
    % toc
    
    RandomVar{1} = sparse(1:nobs, IC_fam, ones(nobs,1), nobs, nfam);
    RandomVar{2} = sparse(1:nobs, IC_subj, ones(nobs,1), nobs, nsubj);
    
    % Prepare generalized matrix version of MoM estimator
    % tic
    S_sum = Ss{1};
    for i = 2:length(Ss)
        S_sum = S_sum + Ss{i};
    end
    [subvec1, subvec2] = find(S_sum); % Use full matrix, to simplify IGLS -- should be possible to limit to tril
    %[subvec1 subvec2] = find(tril(S_sum)); % Should exclude diagonals: tril(S_sum,-1)
    indvec = sub2ind([nobs nobs],subvec1,subvec2);

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

    % Scale back to using tril on S_sum
    % [subvec1, subvec2] = find(tril(S_sum)); % Should exclude diagonals: tril(S_sum,-1)
    % indvec             = sub2ind([nobs nobs],subvec1,subvec2);

    M = zeros(length(indvec),length(Ss));
    for i = 1:length(Ss)
        M(:,i) = Ss{i}(indvec);
    end

    % Create grid of normalized random effects
    binvals_edges       = linspace(0,1,nbins+1); 
    binvals_edges(end)  = binvals_edges(end)+0.0001;

    % New ND version
    if length(RandomEffects) == 2
        sig2gridi = colvec(1:length(binvals_edges)-1);
        sig2gridl = colvec(binvals_edges(1:end-1));
        sig2gridu = colvec(binvals_edges(2:end));
    else
        sig2gridi = ndgrid_amd(repmat({1:length(binvals_edges)-1}, [1 length(RandomEffects)-1]));
        sig2gridl = ndgrid_amd(repmat({binvals_edges(1:end-1)},    [1 length(RandomEffects)-1]));
        sig2gridu = ndgrid_amd(repmat({binvals_edges(2:end)},      [1 length(RandomEffects)-1]));
    end
    sig2grid_ivec = find(sum(sig2gridl,2)<=1); % Get rid of "impossible" bins
    sig2gridl     = sig2gridl(sig2grid_ivec,:);
    sig2gridu     = sig2gridu(sig2grid_ivec,:);
    sig2gridi     = sig2gridi(sig2grid_ivec,:);
    sig2grid      = (sig2gridl+sig2gridu)/2;
    sig2gridind   = sub2ind_amd(nbins*ones(1,length(RandomEffects)-1),sig2gridi);
    nsig2bins     = size(sig2gridl,1); % Should handle case of no binning

    % Prepare FamilyStruct
    FamilyStruct = struct('clusterinfo', {clusterinfo}, 'M', {M},                     ...
                          'famtypevec',  {famtypevec},  'famtypelist', {famtypelist}, ...
                          'nfamtypes',   nfamtypes,     'iid', {iid},                 ...
                          'fid',         {fid},         'iid_list', {iid_list},       ...
                          'fid_list',    {fid_list},    'nfam', nfam,                 ...
                          'sig2grid',    sig2grid,      'sig2gridl', sig2gridl,       ...
                          'sig2gridu',   sig2gridu,     'sig2gridi', sig2gridi,       ...
                          'sig2gridind', sig2gridind,   'nsig2bins', nsig2bins,       ...
                          'subvec1',     subvec1,       'subvec2', subvec2,           ...
                          'Ss',          {Ss});
else
    clusterinfo = FamilyStruct.clusterinfo;
    M           = FamilyStruct.M;
    nsig2bins   = FamilyStruct.nsig2bins;
    nfam        = FamilyStruct.nfam;
    famtypevec  = FamilyStruct.famtypevec;
    nfamtypes   = FamilyStruct.nfamtypes;
    sig2grid    = FamilyStruct.sig2grid;
    sig2gridl   = FamilyStruct.sig2gridl;
    sig2gridu   = FamilyStruct.sig2gridu;
    subvec1     = FamilyStruct.subvec1;
    subvec2     = FamilyStruct.subvec2;
    Ss          = FamilyStruct.Ss;

    [iid_list, IA, IC_subj] = unique(iid,'stable'); nsubj = length(iid_list); nobs = length(iid);
    [fid_list, IA, IC_fam]  = unique(fid,'stable'); nfam  = length(fid_list);
    RandomVar{1} = sparse(1:nobs, IC_fam, ones(nobs,1), nobs, nfam);
    RandomVar{2} = sparse(1:nobs, IC_subj, ones(nobs,1), nobs, nsubj);
end

tshim = now-t0;

Mi = single(pinv(M));
Cov_MoM = Mi*Mi'; % Variance  / covariance of MoM estimates, per unit of residual error variance

% logging('size(M) = [%d %d]',size(M));
% logging('Cov_MoM:'); disp(Cov_MoM);
% logging('Mi*M:'); disp(Mi*M);

if ~isempty(synthstruct)
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

% Various initialization
beta_hat                                = zeros(size(X,2), size(ymat,2), class(ymat));
[beta_se, zmat, ymat_hat, ymat_res]     = deal(zeros(size(beta_hat), class(ymat)));
[betacon_hat, betacon_se]               = deal(zeros(size(contrasts,1), size(ymat,2), class(ymat)));
binvec                                  = NaN(1, size(ymat,2));
Randombeta                              = cell(1, size(ymat,2));

if Hessflag
    Hessmat = NaN([length(RandomEffects) length(RandomEffects) size(ymat,2)]);
else
    Hessmat = [];
end


%% Initialization
% beta_record = zeros(niter , size(X,2)); % used for convergence check
% sigmat_record = zeros(niter, 3);

iter = 1;

converged = false;

while ~converged && (iter < MaxIter)

    if iter == 1

        % initialization of working response and working weights
        u   = log((ymat + 0.5) ./ (1.5 - ymat));
        W   = (ymat + 0.5) .* (1.5 - ymat) ./4;
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
        beta_hat     = iXtX * (X' * u);
        u_hat        = X * beta_hat;
        u_res        = u - u_hat;
        sig2tvec     = sum(u_res.^2,1)/(size(u_res, 1) - size(X, 2)); % Adjust for the number of estimated parameters
    
        % Computate random variances (only fid and iid)
        [~,   sig2mat] =    ...
         FEMA_fit_simplified(X, iid, eid, fid, u_res, sig2tvec, pihatmat, W_1,  ...
         "MLflag", MLflag, "FamilyStruct", FamilyStruct, "NonnegFlag", NonnegFlag);
    
        % GLS updating fixed effects estimates
        [allWsTerms, beta_hat, ~, ~] = FEMA_GLS(u, X, W_1, sig2mat, RandomEffects, ...
                                                clusterinfo, nfamtypes, famtypevec, ...
                                                "GroupByFamType", GroupByFamType, "useLSQ", useLSQ);

        % recoring the beta and sig2mat in the previous estimation
        beta_old    = beta_hat;
        sig2mat_old = sig2mat;
   
        % beta_record(iter,:) = beta_hat(:,1);
        % sigmat_record(iter,:) = sig2mat;
        
    
    else % for iter>0
    
        u_res = u - (X * beta_hat);
    
        % BLUP for random coefficients gamma estimation
        for coli = 1:size(ymat,2)
    
            Wu                 = allWsTerms{coli} * u_res(:,coli);
            Randombeta{coli}   = zeros(nobs, length(RandomEffects)-1);
    
            for r = 1:(length(RandomEffects)-1) % excluding the error part
        
                Randombeta_tmp        = sig2mat(r,coli) * RandomVar{r}' * Wu; 
                Randombeta{coli}(:,r) = RandomVar{r} * Randombeta_tmp; % random coefficients for all observations
    
            end
    
        end
    
        u_update_r = cell2mat(cellfun(@(x) sum(x,2), Randombeta, 'UniformOutput', false));
        u_update = X * beta_hat + u_update_r;
    
        % logistic transformation
        p       = 1 ./ (1+exp(-u_update));
        p       = max(min(p, 1 - 0.005), 0.005);
        
    
        % update u and W
        W        = p .* (1-p);
        W_1      = 1 ./ W;
        u        = u_update + (ymat - p) .* W_1;

        % % GLS to update beta
        [allWsTerms, beta_hat, beta_se, beta_cov] = FEMA_GLS(u, X, W_1, sig2mat, RandomEffects, ...
                                                clusterinfo, nfamtypes, famtypevec, "allWsTerms", allWsTerms, ...
                                                "GroupByFamType", GroupByFamType, "useLSQ", useLSQ);



        % update sigma and V
        u_res = u - (X * beta_hat);
        sig2tvec     = sum(u_res.^2,1)/(size(u_res, 1) - size(X, 2)); 
        [~,   sig2mat] =    ...
         FEMA_fit_simplified(X, iid, eid, fid, u_res, sig2tvec, pihatmat, W_1,  ...
         "MLflag", MLflag, "FamilyStruct", FamilyStruct, "NonnegFlag", NonnegFlag);


        [allWsTerms, ~, ~, ~] = FEMA_GLS(u, X, W_1, sig2mat, RandomEffects, ...
                                                clusterinfo, nfamtypes, famtypevec, ...
                                                "GroupByFamType", GroupByFamType, "useLSQ", useLSQ, ...
                                                "GLSflag",false);


        if (norm(beta_hat - beta_old) <= tol)
            converged = true;
        end

        % recoring the beta and sig2mat in the previous estimation
        beta_old    = beta_hat;
        sig2mat_old = sig2mat;

        % beta_record(iter,:) = beta_hat(:,1);
        % sigmat_record(iter,:) = sig2mat;
    
    end
    
    iter = iter+1;

end

if converged
    fprintf('Converged at iteration %d.\n', iter);
else
    fprintf('Maximum iterations (%d) reached.\n', MaxIter);
    % display(beta_hat);
end

zmat = double(beta_hat) ./ double(beta_se);
logpmat = -sign(zmat) .* log10(normcdf(-abs(zmat))*2);

end


