function [beta_hat,          beta_se,          zmat,            logpmat,        sig2tvec,      sig2mat,      Hessmat,         ...
          logLikvec,         beta_hat_perm,    beta_se_perm,    zmat_perm,      sig2tvec_perm, sig2mat_perm, logLikvec_perm,  ...
          snp_beta_hat,      snp_beta_se,      snp_tStats,      snp_logpValue,                                                ...
          snp_beta_hat_perm, snp_beta_se_perm, snp_tStats_perm, snp_logpValue_perm] =                                         ...
          FEMA_fit_GWAS(X, iid, eid, fid, agevec, ymat, niter, contrasts, nbins, pihatmat, bFile, outDir, varargin)
% Function to fit fast and efficient linear mixed effects model and conduct 
% GWAS analyses across a series of phenotypes
%
% For notation below: 
% n = observations, 
% p = predictors (fixed effects), 
% v = imaging units (e.g. voxels/vertices)
% g = number of SNPs
%
% Fan et al., (2021) - FEMA: Fast and efficient mixed-effects algorithm for
%                      population-scale whole brain imaging data, BioRxiv 
%                      https://doi.org/10.1101/2021.10.27.466202
%
%% INPUTS:
%   X <num>                    :  design matrix (n x p), with intersept if needed
%   iid <char>                 :  subject IDs to match imaging data
%   eid <char>                 :  eventname (n x 1)
%   fid <num>                  :  family ID (members of the same family unit have same value)
%   agevec <num>               :  participants age (n x 1)
%   ymat <num>                 :  matrix of imaging data (n x v)
%   niter <num>                :  number of iterations (default 1)
%   contrasts <num> OR <path>  :  contrast matrix (c x p), where c is number of contrasts to compute,
%                                 OR path to file containing contrast matrix (readable by readtable)
%   nbins <num>                :  number of bins across Y for estimating random effects (default 20)
%   pihatmat <num>             :  matrix of genetic relatedness (n x n) --> already intersected to match X and Y sample
%   bFile <char>               :  full path to a PLINK file (without extension)
%   outDir <char>              :  full path to where the results should be saved; if the directory does not exist, it is created
%
% Optional input arguments:
%   RandomEffects <cell>       :  list of random effects to estimate (default {'F','S','E'}):
%                                   family relatedness (F)
%                                   subject - required for longitudinal analyses (S)
%                                   error (E) - always required
%                                   additive genetic relatedness (A) - must include file path to genetic relatedness data (pihat) for this option
%   nperms <num>               :  deault 0 --> if >0 will run and output permuted effects
%   CovType <char>             :  default 'analytic' --> no other options currently available
%   FixedEstType <char>        :  default 'GLS' --> other option: 'OLS'
%   RandomEstType <char>       :  default 'MoM' --> other option: 'ML' (much slower)
%   GroupByFamType <boolean>   :  default true
%   Parallelize <boolean>      :  default false
%   NonnegFlag <blooean>       :  default true - non-negativity constraint on random effects estimation
%   SingleOrDouble <char>      :  default 'double' --> other option: 'single' - for precision
%   logLikflag <boolean>       :  default true - compute log-likelihood
%   PermType <char>            :  input for FEMA_fit - options:
%                                   'wildbootstrap' - residual boostrap --> creates null distribution by randomly flipping the sign of each observation 
%                                   'wildbootstrap-nn' - non-null boostrap --> estimates distribution around effect of interest using sign flipping (used for sobel test)
%   chunkSize <num>            :  number of SNPs to be considered together when estimating the effects (default: 1000)
%   gStdType <char>            :  standardization for genotype data:
%                                   'none' (default)
%                                   'gcta'
%                                   'emperical'
%   pValType <char>            :  method for calculating p value for SNPs:
%                                   't'     (use t distribution with n-p degrees of freedom; default)
%                                   'z'     (use z distribution)
%                                   'chi'   (use chi square distribution with 1 degree of freedom)
%   GWASParallel <boolean>     :  default false - indicates whether GWAS should be performed by parallel computing
%   numParGWAS <num>           :  number of parallel processes to launch for GWAS; if 0 then GWASParallel is set to false (default: 10)
%   GWASParallelMode <char>    :  option indicating how GWAS parallel computing should be perfomed:
%                                   'chromosome' (SNPs are divided based on their chromosome number and numParGWAS chromosomes processed in parallel; default)
%                                   'chunk'      (SNPs are divided into chunkSize fragments and numParGWAS fragments are processed in parallel)
%   gImpute <boolean>          :  mean impute missing genotype data (default: true)
%
%% OUTPUTS:
%   beta_hat                   :  estimated beta coefficients (p x v)
%   beta_se                    :  estimated beta standard errors (p x v)
%   zmat                       :  z statistics (p x v)
%   logpmat                    :  log10 p-values (p x v)
%   sig2tvec                   :  total residual error of model at each vertex/voxel (1 x v)
%   sig2mat                    :  normalized random effect variances (length(random_effects) x v)
%   snp_beta_hat               :  estimated beta coefficients for SNPs (g x v)
%   snp_beta_se                :  estimated beta standard errors for SNPs (g x v)
%   snp_tStats:                :  one sample t-test statistics for SNPs (g x v)
%   snp_logpValue:             :  -log10 p value for the test statistics (g x v)

% This software is Copyright (c) 2021 The Regents of the University of California. All Rights Reserved.
% See LICENSE.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%% Check inputs and initialize
if ~exist('niter','var') || isempty(niter)
	niter = 0;
end

if ~exist('contrasts','var')
	contrasts = [];
end

if ~exist('nbins','var') || isempty(nbins)
    nbins = 20;
end

if ~exist('pihatmat','var')
    pihatmat = [];
end

% Ensure bFile and outDir are specified
if ~exist('bFile', 'var') || isempty(bFile)
    error('Please provide full path to a PLINK file');
end

if ~exist('outDir', 'var') || isempty(outDir)
    error('Please provide full path to where results should be saved');
else
    % Ensure directory exists
    if ~exist(outDir, 'dir')
        mkdir(outDir);
    end
end

% Should change to allow p to be passed in, so as to avoid having to 
% duplicate input argument parsing in FEMA_wrapper and FEMA_fit
p = inputParser;
addParamValue(p, 'CovType', 'analytic'); %#ok<*NVREPL>
addParamValue(p, 'FixedEstType', 'GLS');
addParamValue(p, 'RandomEstType', 'MoM');
addParamValue(p, 'PermType', 'wildbootstrap');
addParamValue(p, 'GroupByFamType', true);
addParamValue(p, 'Parallelize', false);
addParamValue(p, 'NonnegFlag', true);               % Perform lsqnonneg on random effects estimation
addParamValue(p, 'SingleOrDouble', 'double');
addParamValue(p, 'RandomEffects', {'F' 'S' 'E'});   % Default to Family, Subject, and eps
addParamValue(p, 'logLikflag', false);
addParamValue(p, 'Hessflag', false);
addParamValue(p, 'ciflag', false);
addParamValue(p, 'nperms', 0);
addParamValue(p, 'reverse_cols', 1);                % AMD in development
addParamValue(p, 'reverseinferenceflag', 0);        % AMD in development
addParamValue(p, 'FatherID', {});                   % Father ID, ordered same as pihatmat
addParamValue(p, 'MotherID', {});                   % Mother ID, ordered same as pihatmat
addParamValue(p, 'PregID', {});                     % Pregnancy effect (same ID means twins), ordered same as pihatmat
addParamValue(p, 'HomeID', {});                     % Home effect (defined as same address ID), ordered same as pihatmat
addParamValue(p, 'chunkSize', 1000);                % Default to 1000 as the chunk size for GWAS
addParamValue(p, 'gStdType', 'none');               % use raw genotype information - 0, 1, and 2
addParamValue(p, 'pValType', 't');                  % use t distribution for calculating p values
addParamValue(p, 'gImpute', true);                  % Mean impute missing genotyping data
addParamValue(p, 'GWASParallel', false);            % Whether to run GWAS in parallel
addParamValue(p, 'numParGWAS', 10);                 % Default is 10 GWAS in parallel
addParamValue(p, 'GWASParallelMode', 'chromosome'); % Default is to split by chromosomes

if ~isfinite(contrasts)
    fname_contrasts = p.Results.contrasts;
    logging('Reading contrast matrix from %s',fname_contrasts);
    contrasts = readtable(fname_contrasts);
end

if ~isempty(contrasts) && size(contrasts,2)<size(X,2) % Zeros-pad contrasts, if needed
    contrasts = cat(2,contrasts,zeros([size(contrasts,1) size(X,2)-size(contrasts,2)]));
end

parse(p,varargin{:})
CovType              = p.Results.CovType; %#ok<*NASGU>
FixedEstType         = p.Results.FixedEstType;
RandomEstType        = p.Results.RandomEstType;
GroupByFamType       = p.Results.GroupByFamType;
Parallelize          = p.Results.Parallelize;
NonnegFlag           = p.Results.NonnegFlag;
SingleOrDouble       = p.Results.SingleOrDouble;
RandomEffects        = p.Results.RandomEffects;
OLSflag              = ismember(lower(FixedEstType),{'ols'});
GLSflag              = ismember(lower(FixedEstType),{'gee' 'gls'});
MoMflag              = ismember(lower(RandomEstType),{'mom'});
MLflag               = ismember(lower(RandomEstType),{'ml'});
logLikflag           = p.Results.logLikflag;
Hessflag             = p.Results.Hessflag;
ciflag               = p.Results.ciflag;
nperms               = p.Results.nperms;
PermType             = p.Results.PermType;
chunkSize            = p.Results.chunkSize;
gStdType             = lower(p.Results.gStdType);
pValType             = lower(p.Results.pValType);
gImpute              = p.Results.gImpute;
GWASParallel         = p.Results.GWASParallel;
numParGWAS           = p.Results.numParGWAS;
GWASParallelMode     = lower(p.Results.GWASParallelMode);
reverse_cols         = p.Results.reverse_cols; % AMD
reverseinferenceflag = p.Results.reverseinferenceflag; % AMD

% Ensure gStdType is valid
if ~ismember(gStdType, {'none', 'gcta', 'emperical'})
    warning('Unknown parameter for gStdType; resorting to none');
    gStdType = 'none';
end

% Ensure pValType is valid
if ~ismember(pValType, {'t', 'z', 'chi'})
    warning('Unknown parameter for pValType; resorting to t');
    pValType = 't';
end

% Validate GWAS parallel settings
if GWASParallel
    % If numParGWAS is 0, do not run in parallel
    if numParGWAS == 0
        GWASParallel = false;
    end
    
    % Validate GWASParallelMode
    if ~ismember(GWASParallelMode, {'chromosome', 'chunk'})
        warning('Unknown parameter for GWASParallelMode; resorting to chromosome');
    end
end

if ~isempty(setdiff(RandomEffects,{'F' 'S' 'E'})) % Grouping by family type is only supported for RandomEffects 'F' 'S' 'E'
  GroupByFamType = false;
end

starttime = now();
logging('***Start***');
 
fprintf(1,'ModelSingularityIndex = %g\n',cond(X'*X)/cond(diag(diag(X'*X)))); % Should perhaps report a more standard measure of model singularity?

% Make sure all output params are defined
[logLikvec, beta_hat_perm, beta_se_perm, zmat_perm, sig2tvec_perm, sig2mat_perm, logLikvec_perm] = deal([]);

%% Parse family structure
[clusterinfo, Ss, iid, famtypevec, famtypelist, subj_famtypevec] = FEMA_parse_family(iid, eid, fid,  agevec, pihatmat,   ...
                                                                                    'RandomEffects', RandomEffects,      ...
                                                                                    'FatherID',      p.Results.FatherID, ...
                                                                                    'MotherID',      p.Results.MotherID, ...
                                                                                    'PregID',        p.Results.PregID,   ...
                                                                                    'HomeID',        p.Results.HomeID); %#ok<ASGLU>
[iid_list, ~, ~]  = unique(iid,'stable'); 
nsubj             = length(iid_list); 
nobs              = length(iid);
[fid_list, ~, ~]  = unique(fid,'stable'); 
nfam              = length(fid_list);
nfamtypes         = length(famtypelist);

%% Ensure that SNP data exists for all IIDs
[~, ~, ~, ~, check] = FEMA_parse_PLINK(bFile, iid, true);
if ~check
    error('One or more IIDs not found in PLINK files');
end

%% Prepare generalized matrix version of MoM estimator
S_sum = Ss{1};
for i = 2:length(Ss)
    S_sum = S_sum + Ss{i};
end
[subvec1, subvec2] = find(tril(S_sum));
indvec = sub2ind([nobs nobs], subvec1, subvec2);
M = zeros(length(indvec), length(Ss));
for i = 1:length(Ss)
    M(:,i) = Ss{i}(indvec);
end
Mi = single(pinv(M));
Cov_MoM = Mi*Mi'; % Variance  / covariance of MoM estimates, per unit of residual error variance

%% Higher order moments
if strcmpi('HMoM',RandomEstType)
  [Mrows, dummy, Mindvec] = unique(M,'rows','stable'); nMrows= size(Mrows,1); %#ok<ASGLU> % Should eliminate diagonal elements from M?; perfomr gridding of to handle continuous random effect weights (e.g., GRM)?
  n  = 1e5; % Takes about 12s, with FSE model
  tic
  mu = [0 0];
  morder = 4;
  mumat = NaN(nMrows,nsig2bins,morder); covmat = NaN(nMrows,nsig2bins,morder,morder); X_tmp = NaN(n,morder); %#ok<NODEF>
  for Mind = 1:nMrows
    for sig2bini = 1:nsig2bins
      rho = Mrows(Mind,end);
      Sigma = [1 rho; rho 1];
      ymat_tmp = sqrt(max(0,1-sum(sig2grid(sig2bini,:))))*mvnrnd(mu,Sigma,n); %#ok<NODEF>
      for j = 1:size(sig2grid,2)
        rho = Mrows(Mind,j);
        Sigma = [1 rho; rho 1];
        ymat_tmp = ymat_tmp + sqrt(sig2grid(sig2bini,j))*mvnrnd(mu,Sigma,n);
      end
      % Make sure ymat_tmp in fact has unit variance, and correct covariance
      LHS_tmp = ymat_tmp(:,1).*ymat_tmp(:,2);
      fprintf(1,'%d %d: %f %f mu=%f\n',Mind,sig2bini,sig2grid(sig2bini,:),mean(LHS_tmp));
      disp(cov(ymat_tmp))
      % Save distribution of LHS_tmp
      for oi = 1:morder
        X_tmp(:,oi) = LHS_tmp.^oi;
      end
      mumat(Mind,sig2bini,:) = mean(X_tmp);
      covmat(Mind,sig2bini,:,:) = cov(X_tmp);
    end
  end
  toc
end

logging('size(M) = [%d %d]',size(M));
logging('Cov_MoM:'); disp(Cov_MoM);
logging('Mi*M:'); disp(Mi*M);

%% Create grid of normalized random effects
% Currently supports only FSE models -- should generalize to arbitrary set of random effects
binvals_edges = linspace(0,1,nbins+1); binvals_edges(end) = binvals_edges(end)+0.0001;
% if 0 % Old 2D version
%   [tmpmat1l tmpmat2l] = ndgrid(binvals_edges(1:end-1),binvals_edges(1:end-1));
%   [tmpmat1u tmpmat2u] = ndgrid(binvals_edges(2:end),binvals_edges(2:end));
%   sig2gridl = cat(2,tmpmat1l(:),tmpmat2l(:));
%   sig2gridu = cat(2,tmpmat1u(:),tmpmat2u(:));
% else % New ND version
if length(RandomEffects)==2
    sig2gridl = colvec(binvals_edges(1:end-1));
    sig2gridu = colvec(binvals_edges(2:end));
else
    sig2gridl = ndgrid_amd(repmat({binvals_edges(1:end-1)},[1 length(RandomEffects)-1]));
    sig2gridu = ndgrid_amd(repmat({binvals_edges(2:end)},[1 length(RandomEffects)-1]));
end

ivec_tmp  = find(sum(sig2gridl,2)<=1); % Get rid of "impossible" bins
sig2gridl = sig2gridl(ivec_tmp,:);
sig2gridu = sig2gridu(ivec_tmp,:);
sig2grid  = (sig2gridl+sig2gridu)/2;
nsig2bins = size(sig2gridl,1); % Should handle case of no binning

%% Initialization
Vs_famtype  = cell(1, nfamtypes); 
Vis_famtype = cell(1, nfamtypes); 
Ws_famtype  = cell(1, nfamtypes);
Vs_fam      = cell(1, nfam); 
Vis_fam     = cell(1, nfam); 
Ws_fam      = cell(1, nfam);
numYvars    = size(ymat, 2);
beta_hat    = zeros(size(X,2), numYvars, class(ymat)); 
beta_se     = zeros(size(beta_hat),    class(ymat)); 
zmat        = zeros(size(beta_hat),    class(ymat)); 
ymat_hat    = zeros(size(ymat),        class(ymat)); 
ymat_res    = zeros(size(ymat),        class(ymat)); 
betacon_hat = zeros(size(contrasts,1), numYvars, class(ymat)); 
betacon_se  = zeros(size(contrasts,1), numYvars, class(ymat)); 
binvec      = NaN(1,                   numYvars, class(ymat));

if reverseinferenceflag
    ncols_ri          = numYvars;
    X_bak_ri          = X;
    ymat_bak_ri       = ymat;
    beta_hat_perm_ri  = NaN([size(beta_hat) nperms+1],class(ymat));
    beta_se_perm_ri   = NaN([size(beta_se) nperms+1],class(ymat));
    zmat_perm_ri      = NaN([size(zmat) nperms+1],class(ymat));
    sig2mat_perm_ri   = [];
    sig2tvec_perm_ri  = [];
    logLikvec_perm_ri = [];
else
    ncols_ri = 1;
end

if Hessflag
    Hessmat = NaN([length(RandomEffects) length(RandomEffects) numYvars]);
else
    Hessmat = [];
end

%% Main module
for coli_ri=1:ncols_ri

    if reverseinferenceflag
        ymat = X_bak_ri(:,reverse_cols);
        X = [ymat_bak_ri(:,coli_ri) X_bak_ri(:,setdiff(1:size(X_bak_ri,2),reverse_cols))];
    end

    % Control randseed here?
    digits_nperms = max(ceil(log10(nperms+1)),1);
    loop_timer_start = now();
    for permi = 0:nperms

        permstart = now();

        if permi == 1 % Initialize perm, based on initial fit
            sig2mat_bak  = sig2mat;
            sig2tvec_bak = sig2tvec;
            binvec_bak   = binvec;
            zmat_bak     = zmat;
            ymat_bak     = ymat;
            ymat_res_bak = ymat_res;

            if ismember(lower(PermType),{'wildbootstrap'}) % Residual bootstrap - DEFAULT
                ymat_hat_bak = zeros(size(ymat));
            elseif ismember(lower(PermType),{'wildbootstrap-nn'}) % Non-null wild boostrap
                ymat_hat_bak = ymat_hat;
            end
        end
            
        if permi>0 % Perform resampling
            if ismember(lower(PermType),{'wildbootstrap'}) || ismember(lower(PermType),{'wildbootstrap-nn'}) %DEFAULT
                for fi = 1:nfam
                    % ymat(clusterinfo{fi}.jvec_fam,:) = ymat_hat_bak(clusterinfo{fi}.jvec_fam,:) + (2*randi(2)-3)*ymat_res_bak(clusterinfo{fi}.jvec_fam,:); % Use Rademacher distribution (-1 or 1, with equal probability) for “wild weights”
                    ymat(clusterinfo{fi}.jvec_fam,:) = ymat_hat_bak(clusterinfo{fi}.jvec_fam,:) + randn*ymat_res_bak(clusterinfo{fi}.jvec_fam,:); % Use Normal distribution for “wild weights” -- gives really bad z-score estimates for zer-inflated covariates
                end
            elseif ~ismember(lower(PermType),{'wildbootstrap'}) || ~ismember(lower(PermType),{'wildbootstrap-nn'})
                error('Resampling scheme not available. PermType must equal wildbootstrap or wildbootstrap-nn')
            end
        end
            
        % Initial OLS estimate
        % At this stage, we only need ymat_res
        % Xi            = pinv(X);
        % beta_hat      = Xi*ymat;
        % ymat_hat      = X*beta_hat;
        % ymat_res      = ymat - ymat_hat;
        % ymat_hat_ols  = ymat_hat; ymat_res_ols = ymat_res;
        % sig2tvec      = mean(ymat_res.^2,1);
        % Cov_beta      = Xi*Xi';
        % beta_se       = sqrt(diag(Cov_beta)*sig2tvec);
        % for ci = 1:size(contrasts,1)
        %     betacon_hat(ci,:) = contrasts(ci,:)*beta_hat; 
        %     betacon_se(ci,:)  = sqrt(contrasts(ci,:)*Cov_beta*contrasts(ci,:)'*sig2tvec);
        % end
        ymat_res = ymat - (X * (X \ ymat));
        sig2tvec = sum(ymat_res.^2,1)/(size(ymat_res,1)-size(X,2)); % Adjust for the number of estimated parameters -- should use effective DOF instead?
  
        % Estimating random effects and updating fixed effects
        for iter = 1:max(1,niter)

            sig2tvec = sum(ymat_res.^2,1)/(size(ymat_res,1)-size(X,2));  % Should we use  ymat_res' * inv(V) * ymat_res instead, as ymat_res ~ N(0, sig_t * V), with inv(V)=Vis, as defined in FEMA_sig2binseg_parfeval?            
            LHS      = ymat_res(subvec1,:).*ymat_res(subvec2,:)./mean(ymat_res.^2,1); % use normalized residuals
            
            if ~NonnegFlag % Standard least squares and max(0,x)
                tmp     = Mi*LHS;
                sig2mat = max(0,tmp); % Variances must be non-negative
            else % Use new version of lsqnonneg_amd to enfoce non-negative variances
                sig2mat = lsqnonneg_amd2(M,LHS);
            end
            
            sig2mat     = sig2mat ./ max(eps,sum(sig2mat,1)); % Is this different from dividing by sig2tvec?
            sig2mat_bak = sig2mat;
            
            logLikvec = [];
            if MLflag
                logLikvec  = nan(1,size(ymat_res,2));
                sig2mat_ml = nan(size(sig2mat)); 
                sig2mat_ll = nan(size(sig2mat)); 
                sig2mat_ul = nan(size(sig2mat));
                
                for coli=1:size(ymat_res, 2)
                    f = @(x) (-1 * FEMA_logLik(exp(x),X, ymat_res(:, coli), clusterinfo, Ss));
                    g = @(x) (-1 * FEMA_logLik(x,     X, ymat_res(:, coli), clusterinfo, Ss));
                    sig2vec0 = double(sig2mat(:, coli) * sig2tvec(coli));
                    fprintf(1,'Optimizing using fmincon\n');
                    tic
                    [sig2vec_hat, cost] = fmincon(g, sig2vec0, [], [], [], [], 0*ones(size(sig2vec0)));
                    toc
                    
                    if Hessflag && permi==0
                        fprintf(1,'Computing Hessian\n');
                        tic
                        [sig2vec_hat, cost, exitflag, output, lambda, grad, H] = fmincon(g, sig2vec_hat, [], [], [], [], 0*ones(size(sig2vec0))); %#ok<ASGLU>
                        toc
                        Hessmat(:,:,coli) = H;
                        if any(~isfinite(H(:))), keyboard; end
                    end
                    
                    if ciflag % Compute confidence intervals on random effects?
                        fprintf(1,'Computing Confidence Interval\n');
                        tic
                        loglikthresh = chi2inv(1-0.05/2,1)/2;
                        sig2vec_ll   = nan(size(sig2vec0)); 
                        sig2vec_ul   = nan(size(sig2vec0));
                        
                        for ri = 1:length(sig2vec_hat)
                            tmpfun  = @(x) g(sig2vec_hat+(colvec([1:length(sig2vec_hat)])==ri)*x)-cost; %#ok<NBRAK>
                            dx0     = 0.01; 
                            y0      = tmpfun(dx0); 
                            dx1     = dx0*sqrt(2/y0); 
                            y1      = tmpfun(dx1); % Get scale
                            x       = [0 max([dx0 dx1])*[0.5 1]]; 
                            y       = [0 tmpfun(x(2)) tmpfun(x(3))]; 
                            p       = polyfit(x,y,2);
                            xvec    = linspace(0,max(x),101); 
                            yvec    = polyval(p,xvec);
                            % figure(coli*10); subplot(length(sig2vec_hat),2,(ri-1)*2+2); plot(xvec,yvec,x,y,'*','lineWidth',2); drawnow;
                            ul      = sig2vec_hat(ri) + (-p(2)+((p(2)^2-4*p(1)*(p(3)-loglikthresh)))^0.5)/(2*p(1));
                            ll      = 0;
                            if sig2vec_hat(ri)>0.01
                                if x(end)>sig2vec_hat(ri)
                                    x = x*sig2vec_hat(ri)/x(end); 
                                end
                                x    = -x; 
                                y    = [0 tmpfun(x(2)) tmpfun(x(3))]; 
                                p    = polyfit(x,y,2);
                                xvec = linspace(min(x),max(x),101); 
                                yvec = polyval(p,xvec);
                                % figure(coli*10);; subplot(length(sig2vec_hat),2,(ri-1)*2+1); plot(xvec,yvec,x,y,'*','lineWidth',2); drawnow;
                                ll   = max(0,sig2vec_hat(ri) + (-p(2)-((p(2)^2-4*p(1)*(p(3)-loglikthresh)))^0.5)/(2*p(1)));
                            end
                            sig2vec_ll(ri) = ll;
                            sig2vec_ul(ri) = ul;
                        end
                        toc
                    end
                    
                    disp(num2str(cost,'%0.6e') )
                    sig2mat_ml(:, coli) = sig2vec_hat;
                    if ciflag
                        sig2mat_ll(:, coli) = sig2vec_ll;
                        sig2mat_ul(:, coli) = sig2vec_ul;
                    end
                    disp(rowvec(sig2mat_ml(:, coli)/sum(sig2mat_ml(:, coli))))
                    logl_ml =  g(sig2mat_ml(:, coli)); % This takes ~0.13s per column
                    logl_mom = g(double(sig2mat(:, coli) * sig2tvec(coli)));
                    logging('pheno %i of %i, perm %i of %i: loglike(MoM)=%.2f, loglike(ML)=%.2f', coli, size(ymat_res, 2), permi, nperms, logl_mom, logl_ml);
                    logLikvec(coli) = -logl_ml;
                end
                sig2tvec_ml = sum(sig2mat_ml);
                sig2mat_ml  = sig2mat_ml ./ sig2tvec_ml;
                if ciflag
                    sig2mat_ci = cat(3,sig2mat_ll,sig2mat_ul) ./ sig2tvec_ml;
                end
                sig2mat  = sig2mat_ml;
                sig2tvec = sig2tvec_ml;
            end

            % Snap to random effects grid -- should copy to above
            nvec_bins = NaN(nsig2bins,1); tvec_bins = zeros(nsig2bins,1);
            for sig2bini = 1:nsig2bins
                tmpvec = true;
                for ri = 1:size(sig2mat,1)-1
                    tmpvec = tmpvec & sig2mat(ri,:)>=sig2gridl(sig2bini,ri) & sig2mat(ri,:)<sig2gridu(sig2bini,ri);
                end
                ivec_bin = find(tmpvec);
                nvec_bins(sig2bini) = length(ivec_bin);
                binvec(ivec_bin) = sig2bini;
            end

            if strcmpi('HMoM',RandomEstType)
                % Improved MoM,  using higher-order moments -- should iterate with binning
                nsig2bins2 = nsig2bins;
                [sv, si] = sort(nvec_bins,'descend'); si = si(sv>0); sv = sv(sv>0);
                for sig2bini = si
                    tic
                    costmat = zeros(nsig2bins2,nvec_bins(sig2bini));
                    for Mind = 1:nMrows
                        LHS_tmp = LHS(Mindvec==Mind,binvec==sig2bini); dims = size(LHS_tmp);  LHS_tmp = colvec(LHS_tmp); X_tmp = NaN(length(LHS_tmp),morder);
                        figure(Mind); clf; hist(colvec(LHS_tmp),100); drawnow; %#ok<HIST>
                        for oi = 1:morder
                            X_tmp(:,oi) = LHS_tmp.^oi;
                        end
                        for sig2bini2 = 1:nsig2bins2 % Should only compute around neighborhood of current grid point
                            costmat(sig2bini2,:) = costmat(sig2bini2,:) + sum(reshape(-mvnpdfln(X_tmp,rowvec(mumat(Mind,sig2bini2,:)),squeeze(covmat(Mind,sig2bini,:,:))),dims),1);
                        end
                    end
                    toc
                    [mv, mi]  = min(costmat,[],1); %#ok<ASGLU>
                    disp(sig2grid(mode(mi),:)) % Random effects estimates are off ([0.1250    0.6750] vs. [0.0750    0.7750]) -- should double-check distributions, moments in simularion code, above
                    disp(sig2grid(sig2bini,:))
                    fprintf(1,'%d %d: [%f %f] [%s] [%s] [%s]\n',Mind,sig2bini,sig2grid(sig2bini,:),num2str(mean(X_tmp,1),' %f'),num2str(rowvec(mumat(Mind,sig2bini,:)),' %f'),num2str(rowvec(mumat(Mind,mode(mi),:)),' %f'));
                end
                % Re-compute grid of random effects
                for sig2bini = 1:nsig2bins
                    tmpvec = true;
                    for ri = 1:size(sig2mat,1)-1
                        tmpvec = tmpvec & sig2mat(ri,:)>=sig2gridl(sig2bini,ri) & sig2mat(ri,:)<sig2gridu(sig2bini,ri);
                    end
                    ivec_bin = find(tmpvec);
                    nvec_bins(sig2bini) = length(ivec_bin);
                    binvec(ivec_bin) = sig2bini;
                end
            end

            % ToDo
            %   Test with synthesized data
            %   Perform inverse rank normalization for each moment, or use other distribution shape parameters

            if logLikflag && ~MLflag
                logLikvec = nan(1,size(ymat_res,2));
                for coli=1:size(ymat_res, 2)
                    logLikvec(coli) = FEMA_logLik(sig2tvec(coli) * sig2mat(:, coli), X, ymat_res(:, coli), clusterinfo, Ss); % Should modify FEMA_logLik to leverage gridding of random effects?
                end
            end

            if iter>niter, break; end

            if Parallelize
                pp   = gcp;
                nseg = pp.NumWorkers;
            else
                nseg = 1;
            end

            sig2mat_save  = sig2mat; % Ugly hack to save resampled random effects estimates
            sig2tvec_save = sig2tvec;
            binvec_save   = binvec;
            if permi>0
                sig2tvec  = sig2tvec_bak;
                sig2mat   = sig2mat_bak;
                binvec    = binvec_bak;
            end
            
            if nseg == 1
                [beta_hat, beta_se, betacon_hat, betacon_se] = FEMA_sig2binseg_parfeval(X, ymat, contrasts, clusterinfo,          ...
                                                                                        binvec, sig2mat, sig2tvec, RandomEffects, ...
                                                                                        GroupByFamType, nfamtypes, famtypevec, OLSflag, SingleOrDouble);
                                                                                            
            else % Should perhaps deprecate parfeval, since it doesn't seem to speed the computation?
                % Split into segments of similar total computation time
                b_bins          = [0.13 4.81e-5]; % Empirical computation time params
                tvec_bins_pred  = b_bins(1) + b_bins(2)*nvec_bins; tvec_bins_pred(nvec_bins==0) = NaN;
                [sv, si]        = sort(tvec_bins_pred);
                ivec_finite     = si(isfinite(sv));
                tvec_tmp        = [rowvec(cumsum(tvec_bins_pred(ivec_finite))/sum(tvec_bins_pred(ivec_finite)))]; %#ok<NBRAK>
                segnum          = 1+min(nseg,floor(nseg*tvec_tmp/(max(tvec_tmp)+eps)));
                segvec          = NaN(size(sig2tvec));
                
                for segi = 1:nseg
                    segvec(ismember(binvec,ivec_finite(segnum==segi))) = segi;
                end
                
                F = repmat([parallel.FevalFuture],[1 nseg]);
                ivecs = false(nseg,size(ymat,2));
                
                for segi = 1:nseg
                    ivec_tmp = segvec==segi;
                    F(segi) = parfeval(@FEMA_sig2binseg_parfeval,4,X,ymat(:,ivec_tmp),contrasts,clusterinfo,binvec(:,ivec_tmp),sig2mat(:,ivec_tmp),sig2tvec(:,ivec_tmp),RandomEffects,GroupByFamType,nfamtypes,famtypevec,OLSflag,SingleOrDouble);
                    ivecs(segi,:) = ivec_tmp;
                end
                
                for resi = 1:nseg
                    [completedIdx] = fetchNext(F);
                    ivec_tmp = ivecs(completedIdx,:);
                    beta_hat(:,ivec_tmp) = F(completedIdx).OutputArguments{1};
                    beta_se(:,ivec_tmp) = F(completedIdx).OutputArguments{2};
                    betacon_hat(:,ivec_tmp) = F(completedIdx).OutputArguments{3};
                    betacon_se(:,ivec_tmp) = F(completedIdx).OutputArguments{4};
                end
            end

            ymat_hat = X*beta_hat;
            ymat_res = ymat - ymat_hat;
        end
        
        if ~isempty(contrasts) % Handle non-empty betacon_hat
            beta_hat = cat(1,betacon_hat,beta_hat);
            beta_se = cat(1,betacon_se,beta_se);
        end
        
        zmat = beta_hat ./ beta_se;

        if nperms > 0
            
            if permi == 0
                beta_hat_perm   = NaN([size(beta_hat)   nperms+1], class(beta_hat));
                beta_se_perm    = NaN([size(beta_se)    nperms+1], class(beta_se));
                zmat_perm       = NaN([size(zmat)       nperms+1], class(zmat));
                sig2mat_perm    = NaN([size(sig2mat)    nperms+1], class(sig2mat));
                sig2tvec_perm   = NaN([size(sig2tvec)   nperms+1], class(sig2tvec));
                logLikvec_perm  = NaN([size(logLikvec)  nperms+1], class(logLikvec));
            end
            
            beta_hat_perm(:,:,permi+1)   = beta_hat;
            beta_se_perm(:,:,permi+1)    = beta_se;
            zmat_perm(:,:,permi+1)       = zmat;
            % sig2mat_perm(:,:,permi+1)  = sig2mat;
            % sig2tvec_perm(:,:,permi+1) = sig2tvec;
            sig2mat_perm(:,:,permi+1)    = sig2mat_save;
            sig2tvec_perm(:,:,permi+1)   = sig2tvec_save;
            if ~isempty(logLikvec)
                logLikvec_perm(:,:,permi+1) = logLikvec;
            end
      
            estimated_time_remaining = (now()-loop_timer_start)*3600*24/permi * (nperms - permi);
            logging('permi=%0*d/%d (%0.2fs - remaining %.0fs)',digits_nperms,permi,nperms,(now-permstart)*3600*24,estimated_time_remaining);
        end
            
        %% GWAS analyses
        if permi == 0
            % If it is the first time here, read genotype data, initialize,
            % and compute some important variables; additionally perform
            % GWAS for the non-permuted data
            
            % Update user
            logging('Finished estimating random effects and fixed effects; starting estimation of coefficients for SNPs');
            
            % Read genotype information
            gTic                            = tic;
            [genomat, Chr, SNPID, basePair] = FEMA_parse_PLINK(bFile, iid, false, gStdType, gImpute);
            logging(['Finished reading genotype file in ', num2str(toc(gTic), '%.2fs')]);
            
            % Initialize a few variables
            numSNPs = size(genomat, 2);
            [snp_beta_hat, snp_beta_se, snp_tStats, snp_logpValue] = deal(zeros(numSNPs, numYvars));
            if nperms > 0
                [snp_beta_hat_perm, snp_beta_se_perm, snp_tStats_perm, snp_logpValue_perm] = deal(zeros(numSNPs, numYvars, nperms, class(ymat)));
            else
                [snp_beta_hat_perm, snp_beta_se_perm, snp_tStats_perm, snp_logpValue_perm] = deal(NaN);
            end
            
            % Residualize phenotype - already computed as ymat_res
            % N.B: these residuals are calculated from fixed effects only
            % i.e., these are marginal residuals, not conditional residuals

            % Residualize genotype - overwrite genomat
            genomat = genomat - (X * (pinv(X) * genomat));
            % genomat = genomat - (X * (X \ genomat)); % seems to take more time
            
            % Degrees of freedom for GWAS: only adjusting for number of X
            % variables; may need additional adjustment for number of
            % random effects or accounting for the genomat being adjusted
            % for number of X variables too?
            df = size(ymat_res,1) - size(X,2);
            
            % Get Vs and Ws terms for GLS implementation
            [Vs_famtype, Ws_famtype, Vs_fam, Ws_fam] = FEMA_compileTerms(clusterinfo, binvec, nfamtypes, famtypevec, sig2mat, RandomEffects, GroupByFamType, SingleOrDouble, OLSflag);            
            
            % Divide SNPs into chunks based on GWASParallel settings
            if GWASParallel
                
                % Split by chromosomes
                if strcmpi(GWASParallelMode, 'chromosome')
                    
                    % Identify chromosomes
                    chr      = unique(Chr, 'stable');
                    numChr   = length(chr);
                
                    % Split genotype across chromosomes
                    splitInfo  = cell(numChr,1);
                    for chromo = 1:numChr
                        tmpLocs                      = strcmpi(Chr, chr{chromo});
                        splitInfo{chromo}.tmpLocs    = tmpLocs;
                        splitInfo{chromo}.Genomat    = genomat(:, tmpLocs);
                        splitInfo{chromo}.chrNumber  = Chr(tmpLocs);
                        splitInfo{chromo}.SNPList    = SNPID(tmpLocs);
                        splitInfo{chromo}.basePair   = basePair(tmpLocs);
                        splitInfo{chromo}.outName    = fullfile(outDir, ['Estimates_Chr', num2str(chr{chromo}, '%02d')]);
                    end
                else
                    % Split by chunks                    
                    allChunks = 1:chunkSize:numSNPs;
                    numChunks = length(allChunks);
                    splitInfo = cell(numChunks,1);
                    for chunk = 1:numChunks
                        if chunk == numChunks
                            tmpLocs = allChunks(chunk):numSNPs;
                        else
                            tmpLocs = allChunks(chunk):allChunks(chunk)+chunkSize-1;
                        end
                        splitInfo{chunk}.tmpLocs    = tmpLocs;
                        splitInfo{chunk}.Genomat    = genomat(:, tmpLocs);
                        splitInfo{chunk}.chrNumber  = Chr(tmpLocs);
                        splitInfo{chunk}.SNPList    = SNPID(tmpLocs);
                        splitInfo{chunk}.basePair   = basePair(tmpLocs);
                        splitInfo{chunk}.outName    = fullfile(outDir, ['Estimates_Chunk', num2str(chunk, '%05d')]);
                    end
                end
                
                % Create parallel pool, if necessary
                currPool = gcp('nocreate');
                if isempty(currPool)
                    parpool(numParGWAS);
                end
                
                % Execute in parallel
                numParts = size(splitInfo,1);
                parfor parts = 1:numParts
                    % Initialize timer
                    partTic = tic;
                
                    % Find estimates
                    FEMA_sig2binseg_parfeval_GWAS(splitInfo{parts}.Genomat, ymat_res,   ...
                                                  clusterinfo, binvec, sig2tvec,        ...
                                                  GroupByFamType, famtypevec,           ...
                                                  OLSflag, Vs_fam, Vs_famtype,          ...
                                                  Ws_fam, Ws_famtype, df,               ...
                                                  SingleOrDouble, pValType,             ...
                                                  splitInfo{parts}.outName,             ...
                                                  splitInfo{parts}.chrNumber,           ...
                                                  splitInfo{parts}.SNPList,             ...
                                                  splitInfo{parts}.basePair);

                    % Update user
                    logging(['Finished part: ', num2str(parts), ' of ', num2str(numParts), ' in ' num2str(toc(partTic), '%.2fs')]);
                end
                
                % Compile together as single variable
                temp_Chr   = cell(numSNPs, 1);
                temp_SNPID = cell(numSNPs, 1);
                for parts  = 1:numParts
                    temp   = load([splitInfo{parts}.outName, '.mat'], 'beta_hat', 'beta_se', 'tStats', 'logpValues', 'Chr', 'SNPID');
                    snp_beta_hat(splitInfo{parts}.tmpLocs, :) = temp.beta_hat;
                    snp_beta_se(splitInfo{parts}.tmpLocs,  :) = temp.beta_se;
                    snp_tStats(splitInfo{parts}.tmpLocs,   :) = temp.tStats;
                    snp_logpValue(splitInfo{parts}.tmpLocs,:) = temp.logpValues;
                    temp_Chr(splitInfo{parts}.tmpLocs,     :) = temp.Chr;
                    temp_SNPID(splitInfo{parts}.tmpLocs,   :) = temp.SNPID;
                    clear temp
                end
                
                % Ensure that SNPIDs are aligned with original order
                if sum(strcmpi(SNPID, temp_SNPID)) ~= numSNPs
                    [~, tmpOrd]   = ismember(SNPID, temp_SNPID);
                    snp_beta_hat  = snp_beta_hat(tmpOrd, :);
                    snp_beta_se   = snp_beta_se(tmpOrd,  :);
                    snp_tStats    = snp_tStats(tmpOrd,   :);
                    snp_logpValue = snp_logpValue(tmpOrd,:);
                end
            else
                % Split into chunks
                allChunks = 1:chunkSize:numSNPs;
                numChunks = length(allChunks);
                        
                % Loop over every chunk of SNPs and calculate slope, SE, t
                % statistics and p values
                for chunk = 1:numChunks
                
                    % Select SNPs within the chunk
                    if chunk == numChunks
                        tmpLocs = allChunks(chunk):numSNPs;
                    else
                        tmpLocs = allChunks(chunk):allChunks(chunk)+chunkSize-1;
                    end

                    % Initialize timer
                    chunkTic = tic;
                    
                    % Temporay save name
                    outName = fullfile(outDir, ['Estimates_Chunk', num2str(chunk, '%05d')]);

                    % Find estimates
                    [snp_beta_hat(tmpLocs,:), snp_beta_se(tmpLocs,:),               ...
                     snp_tStats(tmpLocs,  :), snp_logpValue(tmpLocs, :)] =          ...
                     FEMA_sig2binseg_parfeval_GWAS(genomat(:,tmpLocs), ymat_res,    ...
                                                   clusterinfo, binvec, sig2tvec,   ...
                                                   GroupByFamType, famtypevec,      ...
                                                   OLSflag, Vs_fam, Vs_famtype,     ...
                                                   Ws_fam, Ws_famtype, df,          ...
                                                   SingleOrDouble, pValType,        ...
                                                   outName, Chr(tmpLocs),           ...
                                                   SNPID(tmpLocs), basePair(tmpLocs));

                    % Update user
                    logging(['Finished part: ', num2str(chunk), ' of ', num2str(numChunks), ' in ' num2str(toc(chunkTic), '%.2fs')]);
                end
            end
            
            % Write out formatted GWAS table for every phenotype
            fmt = '%s \t %s \t %s \t %d \t %g \t %g \t %g \t %g \t %d \t %s \n';
            GWASvarNames = {'Chromosome', 'SNPID', 'BasePair', 'SampleSize', 'Beta', 'SE', 'TStatistics', 'log10pValue', 'DF', 'pValueType'};
            for phen = 1:numYvars
                % Prepare table
                GWASTable           = cell(numSNPs, length(GWASvarNames));
                GWASTable(:, 1)     = Chr;
                GWASTable(:, 2)     = SNPID;
                GWASTable(:, 3)     = basePair;
                GWASTable(:, 4:9)   = num2cell([repmat(nobs, numSNPs, 1), snp_beta_hat(:,phen), snp_beta_se(:,phen), snp_tStats(:,phen), snp_logpValue(:,phen), repmat(df, numSNPs, 1)]);
                GWASTable(:, 10)    = repmat(cellstr(pValType), numSNPs, 1);
                GWASTable           = GWASTable';
                
                % Write out the formatted GWAS table
                fid = fopen(fullfile(outDir, ['GWAS_Summary_Phenotype', num2str(phen, '%04d'), '.csv']), 'w');
                fprintf(fid, [repmat('%s \t ', 1, length(GWASvarNames)-1), ' %s \n'], GWASvarNames{:});
                fprintf(fid, fmt, GWASTable{:});
                fclose(fid);
            end
        else
        
            % Perform GWAS for permuted values
            if GWASParallel
                % Create parallel pool, if necessary
                currPool = gcp('nocreate');
                if isempty(currPool)
                    parpool(numParGWAS);
                end

                % Execute in parallel
                numParts = size(splitInfo,1);
                parfor parts = 1:numParts
                    % Initialize timer
                    partTic = tic;

                    % Find estimates
                    FEMA_sig2binseg_parfeval_GWAS(splitInfo{parts}.Genomat, ymat_res,    ...
                                                  clusterinfo, binvec, sig2tvec,   ...
                                                  GroupByFamType, famtypevec,      ...
                                                  OLSflag, Vs_fam, Vs_famtype,     ...
                                                  Ws_fam, Ws_famtype, df,          ...
                                                  SingleOrDouble, pValType,        ...
                                                  [splitInfo{parts}.outName, '_Permi_', num2str(permi, '%03d')], ...
                                                  splitInfo{parts}.chrNumber,      ...
                                                  splitInfo{parts}.SNPList,        ...
                                                  splitInfo{parts}.basePair);

                    % Update user
                    logging(['Finished part: ', num2str(parts), ' of ', num2str(numParts), ' in ' num2str(toc(partTic), '%.2fs')]);
                end
            else
                for chunk = 1:numChunks

                    % Select SNPs within the chunk
                    if chunk == numChunks
                        tmpLocs = allChunks(chunk):numSNPs;
                    else
                        tmpLocs = allChunks(chunk):allChunks(chunk)+chunkSize-1;
                    end

                    % Initialize timer
                    chunkTic = tic;
                    
                    % Temporay save name
                    outName = fullfile(outDir, ['Estimates_Chunk', num2str(chunk, '%05d'), '_Permi_', num2str(permi, '%03d')]);

                    % Find estimates
                    [snp_beta_hat_perm(tmpLocs, :, permi), snp_beta_se_perm(tmpLocs, :, permi),         ...
                     snp_tStats_perm(tmpLocs,   :, permi), snp_logpValue_perm(tmpLocs,  :, permi)] =    ...
                     FEMA_sig2binseg_parfeval_GWAS(genomat(:,tmpLocs), ymat_res,    ...
                                                   clusterinfo, binvec, sig2tvec,   ...
                                                   GroupByFamType, famtypevec,      ...
                                                   OLSflag, Vs_fam, Vs_famtype,     ...
                                                   Ws_fam, Ws_famtype, df,          ...
                                                   SingleOrDouble, pValType,        ...
                                                   outName, Chr(tmpLocs),           ...
                                                   SNPID(tmpLocs), basePair(tmpLocs));

                    % Update user
                    logging(['Finished (permutation ', num2str(permi)', ') chunk: ', num2str(chunk), ' of ', num2str(numChunks), ' in ' num2str(toc(chunkTic), '%.2fs')]);
                end
            end
        end
    end
      
    if reverseinferenceflag
        beta_hat_perm_ri(reverse_cols,coli_ri,:) = beta_hat_perm(reverse_cols, 1:length(reverse_cols),:);
        beta_se_perm_ri(reverse_cols,coli_ri,:)  = beta_se_perm(reverse_cols,  1:length(reverse_cols),:);
        zmat_perm_ri(reverse_cols,coli_ri,:)     = zmat_perm(reverse_cols,     1:length(reverse_cols),:);
    end
end

if reverseinferenceflag
    beta_hat_perm   = beta_hat_perm_ri;
    beta_se_perm    = beta_se_perm_ri;
    zmat_perm       = zmat_perm_ri;
    sig2mat_perm    = sig2mat_perm_ri;
    sig2tvec_perm   = sig2tvec_perm_ri;
    logLikvec_perm  = logLikvec_perm_ri;
end
    
% Note: should use same resampling scheme for all coli_ri --  use randseed trick?
    
if nperms > 0
    beta_hat        = double(beta_hat_perm(:,:,1));
    beta_se         = double(beta_se_perm(:,:,1));
    zmat            = double(zmat_perm(:,:,1));
    sig2mat         = double(sig2mat_perm(:,:,1));
    sig2tvec        = double(sig2tvec_perm(:,:,1));
    logLikvec       = double(logLikvec_perm(:,:,1));
elseif nperms == 0
    beta_hat_perm   = [];
    beta_se_perm    = [];
    zmat_perm       = [];
    sig2tvec_perm   = [];
    sig2mat_perm    = [];
    logLikvec_perm  = [];
    perms           = [];
end

zmat    = double(beta_hat)./double(beta_se); %CHECK WITH ANDERS
logpmat = -sign(zmat).*log10(normcdf(-abs(zmat))*2); % Should look for normcdfln function

if ciflag
  sig2mat = cat(3,sig2mat,sig2mat_ci);
end

% Save variables - use v7.3
save(fullfile(outDir, 'Variables_FEMA.mat'), 'beta_hat',          'beta_se',          'zmat',            'logpmat',         ...
                                             'sig2tvec',          'sig2mat',          'Hessmat',                            ...
                                             'logLikvec',         'beta_hat_perm',    'beta_se_perm',    'zmat_perm',       ...
                                             'sig2tvec_perm',     'sig2mat_perm',     'logLikvec_perm',                     ...
                                             'snp_beta_hat',      'snp_beta_se',      'snp_tStats',      'snp_logpValue',   ...
                                             'snp_beta_hat_perm', 'snp_beta_se_perm', 'snp_tStats_perm', 'snp_logpValue_perm', '-v7.3');

logging('***Done*** (%0.2f seconds)\n',(now-starttime)*3600*24);

%PrintMemoryUsage

return