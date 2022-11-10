function [beta_hat beta_se zmat logpmat sig2tvec sig2mat Hessmat logLikvec beta_hat_perm beta_se_perm zmat_perm sig2tvec_perm sig2mat_perm logLikvec_perm] = FEMA_fit(X,iid,eid,fid,agevec,ymat,niter,contrasts,nbins,pihatmat,varargin)
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

% Should change to allow p to be passed in, so as to avoid having to duplicate input argument parsing in FEMA_wrapper and FEMA_fit

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
addParamValue(p,'logLikflag',false);
addParamValue(p,'Hessflag',false);
addParamValue(p,'ciflag',false);
addParamValue(p,'nperms',0);
addParamValue(p,'reverse_cols',1); % AMD in development
addParamValue(p,'reverseinferenceflag',0); % AMD in development
addParamValue(p,'FatherID',{}); % Father ID, ordered same as pihatmat
addParamValue(p,'MotherID',{}); % Mother ID, ordered same as pihatmat
addParamValue(p,'PregID',{}); % Pregnancy effect (same ID means twins), ordered same as pihatmat
addParamValue(p,'HomeID',{}); % Home effect (defined as same address ID), ordered same as pihatmat

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
Hessflag = p.Results.Hessflag;
ciflag = p.Results.ciflag;
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

logLikvec = []; beta_hat_perm = []; beta_se_perm = []; zmat_perm = []; sig2tvec_perm = []; sig2mat_perm = []; logLikvec_perm = []; perms = []; % Make sure all output params are defined

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%Parse family structure

[clusterinfo, Ss, iid, famtypevec, famtypelist, subj_famtypevec]=FEMA_parse_family(iid,eid,fid,agevec,pihatmat,'RandomEffects',RandomEffects,...
  'FatherID',p.Results.FatherID,'MotherID',p.Results.MotherID,'PregID',p.Results.PregID,'HomeID',p.Results.HomeID);
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

if strcmpi('HMoM',RandomEstType)
  [Mrows dummy Mindvec] = unique(M,'rows','stable'); nMrows= size(Mrows,1); % Should eliminate diagonal elements from M?; perfomr gridding of to handle continuous random effect weights (e.g., GRM)?
  n  = 1e5; % Takes about 12s, with FSE model
  tic
  mu = [0 0];
  morder = 4;
  mumat = NaN(nMrows,nsig2bins,morder); covmat = NaN(nMrows,nsig2bins,morder,morder); X_tmp = NaN(n,morder);
  for Mind = 1:nMrows
    for sig2bini = 1:nsig2bins
      rho = Mrows(Mind,end);
      Sigma = [1 rho; rho 1];
      ymat_tmp = sqrt(max(0,1-sum(sig2grid(sig2bini,:))))*mvnrnd(mu,Sigma,n);
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
sig2grid = (sig2gridl+sig2gridu)/2;
nsig2bins = size(sig2gridl,1); % Should handle case of no binning

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

if Hessflag 
  Hessmat = NaN([length(RandomEffects) length(RandomEffects) size(ymat,2)]); 
else
  Hessmat = [];
end

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

%                  LHS = ymat_res(subvec1,:).*ymat_res(subvec2,:);
                  LHS = ymat_res(subvec1,:).*ymat_res(subvec2,:)./mean(ymat_res.^2,1); % use normalized residuals
                  
                  if ~NonnegFlag % Standard least squares and max(0,x)
                        tmp = Mi*LHS;
                        sig2mat = max(0,tmp); % Variances must be non-negative
                  else % Use new version of lsqnonneg_amd to enfoce non-negative variances  
                        sig2mat = lsqnonneg_amd2(M,LHS);
                  end

                  sig2mat = sig2mat ./ max(eps,sum(sig2mat,1)); % Is this different from dividing by sig2tvec?

                  sig2mat_bak = sig2mat;

                  logLikvec = [];
                  if MLflag
                    logLikvec = nan(1,size(ymat_res,2)); 
                    sig2mat_ml = nan(size(sig2mat)); sig2mat_ll = nan(size(sig2mat)); sig2mat_ul = nan(size(sig2mat));
                    for coli=1:size(ymat_res, 2)
                      f = @(x) (-1 * FEMA_logLik(exp(x),X,ymat_res(:, coli),clusterinfo,Ss));
                      g = @(x) (-1 * FEMA_logLik(x,X,ymat_res(:, coli),clusterinfo,Ss));
                      sig2vec0 = double(sig2mat(:, coli)*sig2tvec(coli));
                      fprintf(1,'Optimizing using fmincon\n');
                      tic
                      [sig2vec_hat cost] = fmincon(g,sig2vec0,[],[],[],[],0*ones(size(sig2vec0)));
                      toc
                      if Hessflag && permi==0
                        fprintf(1,'Computing Hessian\n');
                        tic 
                        [sig2vec_hat cost exitflag output lambda grad H] = fmincon(g,sig2vec_hat,[],[],[],[],0*ones(size(sig2vec0)));
                        toc
                        Hessmat(:,:,coli) = H;
                        if any(~isfinite(H(:))), keyboard; end
                      end
                      if ciflag % Compute confidence intervals on random effects?
                        fprintf(1,'Computing Confidence Interval\n');
                        tic
                        loglikthresh = chi2inv(1-0.05/2,1)/2;
                        sig2vec_ll = nan(size(sig2vec0)); sig2vec_ul = nan(size(sig2vec0));
                        for ri = 1:length(sig2vec_hat)
                          tmpfun = @(x) g(sig2vec_hat+(colvec([1:length(sig2vec_hat)])==ri)*x)-cost;
                          dx0 = 0.01; y0 = tmpfun(dx0); dx1 = dx0*sqrt(2/y0); y1 = tmpfun(dx1); % Get scale
                          x = [0 max([dx0 dx1])*[0.5 1]]; y = [0 tmpfun(x(2)) tmpfun(x(3))]; p = polyfit(x,y,2);
                          xvec = linspace(0,max(x),101); yvec = polyval(p,xvec);
%                            figure(coli*10); subplot(length(sig2vec_hat),2,(ri-1)*2+2); plot(xvec,yvec,x,y,'*','lineWidth',2); drawnow;
                          ul = sig2vec_hat(ri) + (-p(2)+((p(2)^2-4*p(1)*(p(3)-loglikthresh)))^0.5)/(2*p(1));
                          ll = 0;
                          if sig2vec_hat(ri)>0.01
                            if x(end)>sig2vec_hat(ri), x = x*sig2vec_hat(ri)/x(end); end 
                            x = -x; y = [0 tmpfun(x(2)) tmpfun(x(3))]; p = polyfit(x,y,2); 
                            xvec = linspace(min(x),max(x),101); yvec = polyval(p,xvec);
%                              figure(coli*10);; subplot(length(sig2vec_hat),2,(ri-1)*2+1); plot(xvec,yvec,x,y,'*','lineWidth',2); drawnow;
                            ll = max(0,sig2vec_hat(ri) + (-p(2)-((p(2)^2-4*p(1)*(p(3)-loglikthresh)))^0.5)/(2*p(1)));
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
                    sig2mat_ml = sig2mat_ml ./ sig2tvec_ml;
                    if ciflag
                      sig2mat_ci = cat(3,sig2mat_ll,sig2mat_ul) ./ sig2tvec_ml;
                    end
                    sig2mat = sig2mat_ml;
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
                    [sv si] = sort(nvec_bins,'descend'); si = si(sv>0); sv = sv(sv>0);
                    for sig2bini = si
                      tic
                      costmat = zeros(nsig2bins2,nvec_bins(sig2bini));
                      for Mind = 1:nMrows
                        LHS_tmp = LHS(Mindvec==Mind,binvec==sig2bini); dims = size(LHS_tmp);  LHS_tmp = colvec(LHS_tmp); X_tmp = NaN(length(LHS_tmp),morder);
                        figure(Mind); clf; hist(colvec(LHS_tmp),100); drawnow;
                        for oi = 1:morder
                          X_tmp(:,oi) = LHS_tmp.^oi;
                        end
                        for sig2bini2 = 1:nsig2bins2 % Should only compute around neighborhood of current grid point
                          costmat(sig2bini2,:) = costmat(sig2bini2,:) + sum(reshape(-mvnpdfln(X_tmp,rowvec(mumat(Mind,sig2bini2,:)),squeeze(covmat(Mind,sig2bini,:,:))),dims),1); 
                        end
                      end
                      toc
                      [mv mi]  = min(costmat,[],1);
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

                  if logLikflag & ~MLflag
                    logLikvec = nan(1,size(ymat_res,2));
                    for coli=1:size(ymat_res, 2)
                      logLikvec(coli) = FEMA_logLik(sig2tvec(coli) * sig2mat(:, coli), X, ymat_res(:, coli), clusterinfo, Ss); % Should modify FEMA_logLik to leverage gridding of random effects?
                    end
                  end

                  if iter>niter, break; end

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

if ciflag
  sig2mat = cat(3,sig2mat,sig2mat_ci);
end

logging('***Done*** (%0.2f seconds)\n',(now-starttime)*3600*24);

%PrintMemoryUsage

return

