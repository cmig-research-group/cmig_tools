function [beta_hat beta_se zmat logpmat sig2tvec sig2mat binvec logLikvec beta_hat_perm beta_se_perm zmat_perm sig2tvec_perm sig2mat_perm logLikvec_perm] = FEMA_fit(X,iid,eid,fid,agevec,ymat,niter,contrasts,nbins,pihatmat,varargin)
%
% Function to fit fast and efficient linear mixed effects model
%
% For notation below: n = observations, p = predictors (fixed effects), v = imaging units (e.g. voxels/vertices)
%
% Fan et al., (2021) - FEMA: Fast and efficient mixed-effects algorithm for population-scale whole brain imaging data, BioRxiv (https://doi.org/10.1101/2021.10.27.466202)
%
% INPUTS
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

%
% OUTPUTS
%   beta_hat                   :  estimated beta coefficients (p x v)
%   beta_se                    :  estimated beta standard errors (p x v)
%   zmat                       :  z statistics (p x v)
%   logpmat                    :  log10 p-values (p x v)
%   sig2tvec                   :  total residual error of model at each vertex/voxel (1 x v)
%   sig2mat                    :  normalized random effect variances (length(random_effects) x v)
%

%
% This software is Copyright (c) 2021 The Regents of the University of California. All Rights Reserved.
% See LICENSE.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

if ~exist('niter','var') || isempty(niter)
  niter = 0;
end

if ~exist('contrasts','var')
  contrasts = [];
end

if ~isfinite(contrasts)
      fname_contrasts = p.Results.contrasts;
      logging('Reading contrast matrix from %s',fname_contrasts);
      contrasts = readtable(fname_contrasts);
end

if ~isempty(contrasts) && size(contrasts,2)<size(X,2) % Zeros-pad contrasts, if needed
  contrasts = cat(2,contrasts,zeros([size(contrasts,1) size(X,2)-size(contrasts,2)]));
end

if ~exist('nbins','var') || isempty(nbins)
  nbins = 20;
end

if ~exist('pihatmat','var')
  pihatmat = []; 
end

p = inputParser;
addParamValue(p,'CovType','analytic');
addParamValue(p,'FixedEstType','GLS');
addParamValue(p,'RandomEstType','MoM');
addParamValue(p,'PermType','wildbootstrap');
addParamValue(p,'GroupByFamType',true);
addParamValue(p,'Parallelize',false);
addParamValue(p,'NonnegFlag',true); % Perform lsqnonneg on random effects estimation
addParamValue(p,'SingleOrDouble','double');
addParamValue(p,'RandomEffects',{'F' 'S' 'E'}); % Default to Family, Subject, and eps
addParamValue(p,'logLikflag',true);
addParamValue(p,'nperms',0);
addParamValue(p,'reverse_cols',1); % AMD in development
addParamValue(p,'reverseinferenceflag',0); % AMD in development

parse(p,varargin{:})
CovType = p.Results.CovType;
FixedEstType = p.Results.FixedEstType;
RandomEstType = p.Results.RandomEstType;
GroupByFamType = p.Results.GroupByFamType;
Parallelize = p.Results.Parallelize;
NonnegFlag = p.Results.NonnegFlag;
SingleOrDouble = p.Results.SingleOrDouble;
RandomEffects = p.Results.RandomEffects;
OLSflag = ismember(lower(FixedEstType),{'ols'});
GLSflag = ismember(lower(FixedEstType),{'gee' 'gls'});
MoMflag = ismember(lower(RandomEstType),{'mom'});
MLflag = ismember(lower(RandomEstType),{'ml'});
logLikflag = p.Results.logLikflag;
nperms = p.Results.nperms;
PermType = p.Results.PermType;

reverse_cols=p.Results.reverse_cols; % AMD
reverseinferenceflag=p.Results.reverseinferenceflag; % AMD

if length(setdiff(RandomEffects,{'F' 'S' 'E'}))>0 % Grouping by family type is only supported for RandomEffects 'F' 'S' 'E'
  GroupByFamType = false;
end


if 0 % ~isempty(setdiff(RandomEffects,{'F' 'S' 'E'})) % Thios is no longer needed, with ND random effects grid --  remove these three lines
  niter = 0;
end

starttime = now();
logging('***Start***');
 
fprintf(1,'ModelSingularityIndex = %g\n',cond(X'*X)/cond(diag(diag(X'*X)))); % Should perhaps report a more standard measure of model singularity?

binvec = []; logLikvec = []; beta_hat_perm = []; beta_se_perm = []; zmat_perm = []; sig2tvec_perm = []; sig2mat_perm = []; logLikvec_perm = []; perms = []; % Make sure all output params are defined

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%Parse family structure

[clusterinfo, Ss, iid, famtypevec, famtypelist, subj_famtypevec]=FEMA_parse_family(iid,eid,fid,agevec,pihatmat,'RandomEffects',RandomEffects);
[iid_list IA IC_subj] = unique(iid,'stable'); nsubj = length(iid_list); nobs = length(iid);
[fid_list IA IC_fam] = unique(fid,'stable'); nfam = length(fid_list);
nfamtypes = length(famtypelist);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Prepare generalized matrix version of MoM estimator
%tic
S_sum = Ss{1};
for i = 2:length(Ss)
  S_sum = S_sum + Ss{i};
end
[subvec1 subvec2] = find(tril(S_sum));
indvec = sub2ind([nobs nobs],subvec1,subvec2);
M = zeros(length(indvec),length(Ss));
for i = 1:length(Ss)
  M(:,i) = Ss{i}(indvec);
end
Mi = single(pinv(M));
%toc
Cov_MoM = Mi*Mi'; % Variance  / covariance of MoM estimates, per unit of residual error variance

logging('size(M) = [%d %d]',size(M));
logging('Cov_MoM:'); disp(Cov_MoM);
logging('Mi*M:'); disp(Mi*M);

% Create grid of normalized random effects -- currently supports only FSE models -- should generalize to arbitrary set of random effects
binvals_edges = linspace(0,1,nbins+1); binvals_edges(end) = binvals_edges(end)+0.0001;
if 0 % Old 2D version
  [tmpmat1l tmpmat2l] = ndgrid(binvals_edges(1:end-1),binvals_edges(1:end-1));
  [tmpmat1u tmpmat2u] = ndgrid(binvals_edges(2:end),binvals_edges(2:end));
  sig2gridl = cat(2,tmpmat1l(:),tmpmat2l(:));
  sig2gridu = cat(2,tmpmat1u(:),tmpmat2u(:));
else % New ND version
  if length(RandomEffects)==2
    sig2gridl = colvec(binvals_edges(1:end-1));
    sig2gridu = colvec(binvals_edges(2:end));
  else
    sig2gridl = ndgrid_amd(repmat({binvals_edges(1:end-1)},[1 length(RandomEffects)-1]));
    sig2gridu = ndgrid_amd(repmat({binvals_edges(2:end)},[1 length(RandomEffects)-1]));
  end
end
ivec_tmp = find(sum(sig2gridl,2)<=1); % Get rid of "impossible" bins
sig2gridl = sig2gridl(ivec_tmp,:);
sig2gridu = sig2gridu(ivec_tmp,:);

Vs_famtype = cell(1,nfamtypes); Vis_famtype = cell(1,nfamtypes); Ws_famtype = cell(1,nfamtypes);
Vs_fam = cell(1,nfam); Vis_fam = cell(1,nfam); Ws_fam = cell(1,nfam);

beta_hat = zeros(size(X,2),size(ymat,2),class(ymat)); beta_se = zeros(size(beta_hat),class(ymat)); zmat = zeros(size(beta_hat),class(ymat)); ymat_hat = zeros(size(ymat),class(ymat)); ymat_res = zeros(size(ymat),class(ymat)); 

if reverseinferenceflag
      ncols_ri = size(ymat,2);
      X_bak_ri = X;
      ymat_bak_ri = ymat;
      beta_hat_perm_ri = NaN([size(beta_hat) nperms+1],class(ymat));
      beta_se_perm_ri = NaN([size(beta_se) nperms+1],class(ymat));
      zmat_perm_ri = NaN([size(zmat) nperms+1],class(ymat));
      sig2mat_perm_ri = [];
      sig2tvec_perm_ri = [];
      logLikvec_perm_ri = [];
    else
      ncols_ri = 1;
    end

betacon_hat = zeros(size(contrasts,1),size(ymat,2),class(ymat)); betacon_se = zeros(size(contrasts,1),size(ymat,2),class(ymat)); binvec = NaN(1,size(ymat,2),class(ymat));

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
                  sig2mat_bak = sig2mat;
                  sig2tvec_bak = sig2tvec;
                  binvec_bak = binvec;
                  zmat_bak = zmat;
                  ymat_bak = ymat;
                  ymat_res_bak = ymat_res;

                  if ismember(lower(PermType),{'wildbootstrap'}) % Residual bootstrap - DEFAULT
                        ymat_hat_bak = zeros(size(ymat)); 
                  elseif ismember(lower(PermType),{'wildbootstrap-nn'}) % Non-null wild boostrap
                        ymat_hat_bak = ymat_hat;
                  end
            end
            
            if permi>0 % Perform resampling
                  if ismember(lower(PermType),{'wildbootstrap'}) | ismember(lower(PermType),{'wildbootstrap-nn'}) %DEFAULT
                        for fi = 1:nfam
                              % ymat(clusterinfo{fi}.jvec_fam,:) = ymat_hat_bak(clusterinfo{fi}.jvec_fam,:) + (2*randi(2)-3)*ymat_res_bak(clusterinfo{fi}.jvec_fam,:); % Use Rademacher distribution (-1 or 1, with equal probability) for “wild weights”
                              ymat(clusterinfo{fi}.jvec_fam,:) = ymat_hat_bak(clusterinfo{fi}.jvec_fam,:) + randn*ymat_res_bak(clusterinfo{fi}.jvec_fam,:); % Use Normal distribution for “wild weights” -- gives really bad z-score estimates for zer-inflated covariates
                        end
                  elseif ~ismember(lower(PermType),{'wildbootstrap'}) | ~ismember(lower(PermType),{'wildbootstrap-nn'})
                        error('Resampling scheme not available. PermType must equal wildbootstrap or wildbootstrap-nn')
                  end
            end
            
            Xi = pinv(X);
            beta_hat = Xi*ymat;  % Initially use OLS estimate
            ymat_hat = X*beta_hat; 
            ymat_res = ymat - ymat_hat; 
            ymat_hat_ols = ymat_hat; ymat_res_ols = ymat_res; 
            %  sig2tvec = mean(ymat_res.^2,1);
            sig2tvec = sum(ymat_res.^2,1)/(size(ymat_res,1)-size(X,2)); % Adjust for the number of estimated parameters -- should use effective DOF instead?
            Cov_beta = Xi*Xi';
            beta_se = sqrt(diag(Cov_beta)*sig2tvec);
            
            for ci = 1:size(contrasts,1)
                  betacon_hat(ci,:) = contrasts(ci,:)*beta_hat; betacon_se(ci,:) = sqrt(contrasts(ci,:)*Cov_beta*contrasts(ci,:)'*sig2tvec);
            end
            
            for iter = 1:max(1,niter)

                  %  sig2tvec = mean(ymat_res.^2,1);
                  sig2tvec = sum(ymat_res.^2,1)/(size(ymat_res,1)-size(X,2));  % Should we use  ymat_res' * inv(V) * ymat_res instead, as ymat_res ~ N(0, sig_t * V), with inv(V)=Vis, as defined in FEMA_sig2binseg_parfeval?

                  LHS = ymat_res(subvec1,:).*ymat_res(subvec2,:); % Look into using subasgn?
                  
                  if ~NonnegFlag % Standard least squares and max(0,x)
                        tmp = Mi*LHS;
                        sig2mat = max(0,tmp); % Variances must be non-negative
                  else % Use new version of lsqnonneg_amd to enfoce non-negative variances  
                        sig2mat = lsqnonneg_amd2(M,LHS);
                  end

                  sig2mat = sig2mat ./ max(eps,sum(sig2mat,1));

                  if MLflag
                    sig2mat_ml = nan(size(sig2mat));
                    for coli=1:size(ymat_res, 2)
                      f = @(x) (-1 * FEMA_logLik(exp(x),X,ymat_res(:, coli),clusterinfo,Ss));
                      %x0 = sig2mat(:, coli) * sig2tvec(coli);   % start from MoM estimate (the problem is that some of the estimates might be zero)
                      x0 = sig2tvec(coli) * ones(length(RandomEffects), 1) / length(RandomEffects);
                      sig2mat_ml(:, coli) = exp(fminsearch(f, log(x0), struct('Display', 'off', 'TolFun', 1e-2, 'TolX', 1e-2)));

                      logl_ml =  f(log(sig2mat_ml(:, coli)));
                      logl_mom = f(log(sig2mat(:, coli) * sig2tvec(coli)));
                      logging('pheno %i of %i, perm %i of %i: loglike(MoM)=%.2f, loglike(ML)=%.2f', coli, size(ymat_res, 2), permi, nperms, logl_mom, logl_ml);
                    end
                    sig2tvec_ml = sum(sig2mat_ml);
                    sig2mat_ml = sig2mat_ml ./ sig2tvec_ml;
                    %plot(sig2mat_ml(:), sig2mat(:), '*'); hold on;  plot([0, 1], [0, 1]); xlabel('sig2mat, ML'); ylabel('sig2mat, MoM');
                    sig2mat = sig2mat_ml;
                    sig2tvec = sig2tvec_ml;
                  end

                  if iter>niter, break; end

                  binvec = NaN(size(sig2tvec));
                  nsig2bins = size(sig2gridl,1); % Should handle case on no binning
                  nvec_bins = NaN(nsig2bins,1); tvec_bins = zeros(nsig2bins,1);
                        
                  for sig2bini = 1:nsig2bins
                        if 0 % Old 2D version
                          ivec_bin = find((sig2mat(1,:)>=sig2gridl(sig2bini,1) & sig2mat(1,:)<sig2gridu(sig2bini,1)) & (sig2mat(2,:)>=sig2gridl(sig2bini,2) & sig2mat(2,:)<sig2gridu(sig2bini,2)));
                        else % New ND version
                          tmpvec = true;
                          for ri = 1:size(sig2mat,1)-1
                            tmpvec = tmpvec & sig2mat(ri,:)>=sig2gridl(sig2bini,ri) & sig2mat(ri,:)<sig2gridu(sig2bini,ri);
                          end
                          ivec_bin = find(tmpvec);
                        end
                        nvec_bins(sig2bini) = length(ivec_bin);
                        binvec(ivec_bin) = sig2bini;
                  end

                  if Parallelize
                        pp = gcp;
                        nseg = pp.NumWorkers;
                  else
                        nseg = 1;
                  end

                  sig2mat_save = sig2mat; % Ugly hack to save resampled random effects estimates
                  sig2tvec_save = sig2tvec;
                  binvec_save = binvec;
                  if permi>0 
                        sig2tvec = sig2tvec_bak;
                        sig2mat = sig2mat_bak;
                        binvec = binvec_bak;
                  end

                  if nseg == 1
                        [beta_hat beta_se betacon_hat betacon_se] = FEMA_sig2binseg_parfeval(X,ymat,contrasts,clusterinfo,binvec,sig2mat,sig2tvec,RandomEffects,GroupByFamType,nfamtypes,famtypevec,OLSflag,SingleOrDouble);
                  else % Should perhaps deprecate parfeval, since it doesn't seem to speed the computation?
                        % Split into segments of similar total computation time
                        b_bins = [0.13 4.81e-5]; % Empirical computation time params
                        tvec_bins_pred = b_bins(1) + b_bins(2)*nvec_bins; tvec_bins_pred(nvec_bins==0) = NaN;
                        [sv si] = sort(tvec_bins_pred);
                        ivec_finite = si(isfinite(sv));
                        tvec_tmp = [rowvec(cumsum(tvec_bins_pred(ivec_finite))/sum(tvec_bins_pred(ivec_finite)))];
                        segnum = 1+min(nseg,floor(nseg*tvec_tmp/(max(tvec_tmp)+eps)));
                        segvec = NaN(size(sig2tvec));
                              
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
      
            if logLikflag
                logLikvec = nan(size(ymat_res,2), 1);
                for coli=1:size(ymat_res, 2)
                    logLikvec(coli) = FEMA_logLik(sig2tvec(coli) * sig2mat(:, coli), X, ymat_res, clusterinfo, Ss);
                end
            else
                  logLikvec=[];
            end

            if ~isempty(contrasts) % Handle non-empty betacon_hat
                  beta_hat = cat(1,betacon_hat,beta_hat);
                  beta_se = cat(1,betacon_se,beta_se);
            end

            zmat = beta_hat ./ beta_se;

            if nperms > 0
      
                  if permi == 0
                        beta_hat_perm = NaN([size(beta_hat) nperms+1],class(beta_hat));
                        beta_se_perm = NaN([size(beta_se) nperms+1],class(beta_se));
                        zmat_perm = NaN([size(zmat) nperms+1],class(zmat));
                        sig2mat_perm = NaN([size(sig2mat) nperms+1],class(sig2mat));
                        sig2tvec_perm = NaN([size(sig2tvec) nperms+1],class(sig2tvec));
                        logLikvec_perm = NaN([size(logLikvec) nperms+1],class(logLikvec));
                  end
            
                  beta_hat_perm(:,:,permi+1) = beta_hat;
                  beta_se_perm(:,:,permi+1) = beta_se;
                  zmat_perm(:,:,permi+1) = zmat;
%                  sig2mat_perm(:,:,permi+1) = sig2mat;
%                  sig2tvec_perm(:,:,permi+1) = sig2tvec;
                  sig2mat_perm(:,:,permi+1) = sig2mat_save;
                  sig2tvec_perm(:,:,permi+1) = sig2tvec_save;
                  if ~isempty(logLikvec)
                        logLikvec_perm(:,:,permi+1) = logLikvec;
                  end
      
                  estimated_time_remaining = (now()-loop_timer_start)*3600*24/permi * (nperms - permi);
                  logging('permi=%0*d/%d (%0.2fs - remaining %.0fs)',digits_nperms,permi,nperms,(now-permstart)*3600*24,estimated_time_remaining);
            end
      end
      if reverseinferenceflag
        beta_hat_perm_ri(reverse_cols,coli_ri,:) = beta_hat_perm(reverse_cols,1:length(reverse_cols),:);
        beta_se_perm_ri(reverse_cols,coli_ri,:) = beta_se_perm(reverse_cols,1:length(reverse_cols),:);
        zmat_perm_ri(reverse_cols,coli_ri,:) = zmat_perm(reverse_cols,1:length(reverse_cols),:);
      end
end

if reverseinferenceflag
  beta_hat_perm = beta_hat_perm_ri;
  beta_se_perm = beta_se_perm_ri;
  zmat_perm = zmat_perm_ri;
  sig2mat_perm = sig2mat_perm_ri;
  sig2tvec_perm = sig2tvec_perm_ri;
  logLikvec_perm = logLikvec_perm_ri;
end
    
% Note: should use same resampling scheme for all coli_ri --  use randseed trick?
    
if nperms>0
      beta_hat = double(beta_hat_perm(:,:,1));
      beta_se = double(beta_se_perm(:,:,1));
      zmat = double(zmat_perm(:,:,1));
      sig2mat = double(sig2mat_perm(:,:,1));
      sig2tvec = double(sig2tvec_perm(:,:,1));
      logLikvec = double(logLikvec_perm(:,:,1));
elseif nperms==0
      beta_hat_perm=[];
      beta_se_perm=[];
      zmat_perm=[];
      sig2tvec_perm=[];
      sig2mat_perm=[];
      logLikvec_perm=[];
      perms=[];
end

zmat = double(beta_hat)./double(beta_se); %CHECK WITH ANDERS
logpmat = -sign(zmat).*log10(normcdf(-abs(zmat))*2); % Should look for normcdfln function

logging('***Done*** (%0.2f seconds)\n',(now-starttime)*3600*24);

%PrintMemoryUsage

return

