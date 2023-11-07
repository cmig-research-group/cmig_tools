function [beta_hat,      beta_se,        zmat,        logpmat,              ...
          sig2tvec,      sig2mat,        Hessmat,     logLikvec,            ...
          beta_hat_perm, beta_se_perm,   zmat_perm,   sig2tvec_perm,        ...
          sig2mat_perm,  logLikvec_perm, binvec_save, nvec_bins,            ...
          tvec_bins,     FamilyStruct,   reusableVars] =                    ...
          FEMA_fit(X, iid, eid, fid, agevec, ymat, niter, contrasts, nbins, ...
                   pihatmat, varargin)
% Function to fit fast and efficient linear mixed effects model
%
% For notation below:
% n = observations,
% p = predictors (fixed effects),
% v = imaging units (e.g. voxels/vertices)
% c = number of contrasts to evaluate
% r = number of random effects
%
% Parekh et al., (2021) - Fast and efficient mixed-effects algorithm for 
%                         large sample whole-brain imaging data, bioRxiv 
%                         https://doi.org/10.1101/2021.10.27.466202
%
%% Inputs:
% X               <num>            [n x p]    design matrix, with intercept if needed
% iid             <cell>           [n x 1]    subject IDs to match imaging data
% eid             <cell>           [n x 1]    eventname
% fid             <num>            [n x 1]    family ID (members of the same family unit have same value)
% agevec          <num>            [n x 1]    participants age
% ymat            <num>            [n x v]    matrix of imaging data
% niter           <num>            [1 x 1]    number of iterations (default 1)
% contrasts       <num> OR <path>  [c x p]    contrast matrix, where c is number of contrasts to compute,
%                                             OR path to file containing contrast matrix (readable by readtable)
% nbins           <num>            [1 x 1]    number of bins across Y for estimating random effects (default 20)
% pihatmat        <num>            [n x n]    matrix of genetic relatedness --> already intersected to match X and Y sample
%
%% Optional input arguments:
% RandomEffects   <cell>           list of random effects to estimate (default {'F','S','E'}):
%                                       * F:  family relatedness
%                                       * S:  subject - required for longitudinal analyses
%                                       * E:  error - always required
%                                       * A:  additive genetic relatedness - must include file path to genetic relatedness data (pihat) for this option
%                                       * D:  dominant genetic relatedness - square of A
%                                       * M:  maternal effect - effect of having same mother
%                                       * P:  paternal effect  - effect of having same father
%                                       * H:  home effect - effect of living at the same address
%                                       * T:  twin effect - effect of having the same pregnancy ID
% nperms          <num>            deault 0 --> if >0 will run and output permuted effects
% CovType         <char>           default 'analytic' --> no other options currently available
% FixedEstType    <char>           default 'GLS' --> other option: 'OLS'
% RandomEstType   <char>           default 'MoM' --> other option: 'ML' (much slower)
% GroupByFamType  <boolean>        default true
% NonnegFlag      <blooean>        default true - non-negativity constraint on random effects estimation
% SingleOrDouble  <char>           default 'double' --> other option: 'single' - for precision
% logLikflag      <boolean>        default true - compute log-likelihood
% PermType        <char>           permutation type:
%                                       * 'wildbootstrap':    residual boostrap --> creates null distribution by randomly flipping the sign of each observation
%                                       * 'wildbootstrap-nn': non-null boostrap --> estimates distribution around effect of interest using sign flipping (used for sobel test)
% returnReusable  <boolean>        default false - if true, additionally returns reusableVars as a structure with some variables that can be reused (primarily by FEMA-GWAS)
%
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
% This software is Copyright (c) 2021 
% The Regents of the University of California. 
% All Rights Reserved. See LICENSE.

starttime = now(); %#ok<*TNOW1>
logging('***Start***');

% Extremely quick sanity check on X and y variables
if logical(sum(any(isnan(X)))) || logical(sum(any(isnan(ymat)))) || ...
   logical(sum(any(isinf(X)))) || logical(sum(any(isinf(ymat))))
    error('X and/or ymat have NaN or Inf; please check your data');
end

p = inputParser;

if ~exist('niter', 'var') || isempty(niter)
    niter = 0;
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
addParamValue(p,'FixedEstType', 'GLS');
addParamValue(p,'RandomEstType', 'MoM');
addParamValue(p,'PermType', 'wildbootstrap');
addParamValue(p,'GroupByFamType', true);
addParamValue(p,'NonnegFlag', true); % Perform lsqnonneg on random effects estimation
addParamValue(p,'SingleOrDouble', 'double');
addParamValue(p,'RandomEffects', {'F' 'S' 'E'}); % Default to Family, Subject, and eps
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

parse(p,varargin{:})
CovType              = p.Results.CovType; %#ok<*NASGU>
FixedEstType         = p.Results.FixedEstType;
RandomEstType        = p.Results.RandomEstType;
GroupByFamType       = p.Results.GroupByFamType;
NonnegFlag           = p.Results.NonnegFlag;
SingleOrDouble       = p.Results.SingleOrDouble;
RandomEffects        = p.Results.RandomEffects;
OLSflag              = ismember(lower(FixedEstType),  {'ols'});
GLSflag              = ismember(lower(FixedEstType),  {'gee' 'gls'});
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

% Grouping by family type is only supported for RandomEffects 'F' 'S' 'E'
if ~isempty(setdiff(RandomEffects,{'F' 'S' 'E'}))
    GroupByFamType = false;
end

% Should perhaps report a more standard measure of model singularity?
fprintf(1,'ModelSingularityIndex = %g\n',cond(X'*X)/cond(diag(diag(X'*X))));

% Make sure all output params are defined
[logLikvec,     beta_hat_perm, beta_se_perm,   zmat_perm, ...
 sig2tvec_perm, sig2mat_perm,  logLikvec_perm] = deal([]);

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
    reusableVars.OLSflag        = OLSflag;
    reusableVars.useLSQ         = useLSQ;
    reusableVars.lowRank        = lowRank;
end

t0 = now;

% Parse family structure, if necessary
if ~exist('FamilyStruct', 'var') || isempty(FamilyStruct)
    tic
    [clusterinfo, Ss, iid, famtypevec, famtypelist, subj_famtypevec] =                  ...
     FEMA_parse_family(iid, eid, fid, agevec, pihatmat, 'RandomEffects', RandomEffects, ...
                       'FatherID', p.Results.FatherID,  'MotherID', p.Results.MotherID, ...
                       'PregID',   p.Results.PregID,    'HomeID',   p.Results.HomeID); %#ok<*ASGLU>
    
    [iid_list, IA, IC_subj] = unique(iid,'stable'); nsubj = length(iid_list); nobs = length(iid);
    [fid_list, IA, IC_fam]  = unique(fid,'stable'); nfam  = length(fid_list);
    nfamtypes = length(famtypelist);
    toc
    
    % Prepare generalized matrix version of MoM estimator
    tic
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
    [subvec1, subvec2] = find(tril(S_sum)); % Should exclude diagonals: tril(S_sum,-1)
    indvec             = sub2ind([nobs nobs],subvec1,subvec2);

    M = zeros(length(indvec),length(Ss));
    for i = 1:length(Ss)
        M(:,i) = Ss{i}(indvec);
    end

    % Create grid of normalized random effects -- currently supports only FSE models -- should generalize to arbitrary set of random effects
    %binvals_edges = linspace(0,1,nbins+1); binvals_edges(end) = binvals_edges(end)+0.0001; % Should adjust max to include all values above sig2gridl
    binvals_edges       = linspace(0, 1, nbins+1);
    binvals_edges(end)  = Inf;
    if length(RandomEffects)==2 % Why is this needed? -- N-d version, below, should work for 2-d?
        sig2gridi = colvec(1:length(binvals_edges)-1);
        sig2gridl = colvec(binvals_edges(1:end-1));
        sig2gridu = colvec(binvals_edges(2:end));
    else
        sig2gridi = ndgrid_amd(repmat({1:length(binvals_edges)-1},[1 length(RandomEffects)-1]));
        sig2gridl = ndgrid_amd(repmat({binvals_edges(1:end-1)},   [1 length(RandomEffects)-1]));
        sig2gridu = ndgrid_amd(repmat({binvals_edges(2:end)},     [1 length(RandomEffects)-1]));
    end
    %sig2grid = (sig2gridl+sig2gridu)/2; % Should make sig2gridl+1/nbins
    sig2grid      = sig2gridl+(0.5/nbins);
    sig2grid_ivec = find(sum(sig2grid,2)<=1-0.5/nbins); % Get rid of "impossible" bins -- use middle of bin instead
    sig2gridl     = sig2gridl(sig2grid_ivec,:);
    sig2gridu     = sig2gridu(sig2grid_ivec,:);
    sig2grid      = sig2grid(sig2grid_ivec,:);
    sig2gridi     = sig2gridi(sig2grid_ivec,:);
    sig2gridind   = sub2ind_amd(nbins*ones(1,length(RandomEffects)-1),sig2gridi);
    nsig2bins     = size(sig2grid,1); % Should handle case of no binning

    % Prepare FamilyStruct
    FamilyStruct = struct('clusterinfo', {clusterinfo}, 'M', {M},                     ...
                          'famtypevec',  {famtypevec},  'famtypelist', {famtypelist}, ...
                          'nfamtypes',   nfamtypes,     'iid', {iid},                 ...
                          'fid',         {fid},         'iid_list', {iid_list},       ...
                          'fid_list',    {fid_list},    'nfam', nfam,                 ...
                          'sig2grid',    sig2grid,      'sig2gridl', sig2gridl,       ...
                          'sig2gridu',   sig2gridu,     'sig2gridi', sig2gridi,       ...
                          'sig2gridind', sig2gridind,   'nsig2bins', nsig2bins,       ...
                          'subvec1',     subvec1,       'subvec2', subvec2);
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
end

tshim = now-t0;

Mi = single(pinv(M));
Cov_MoM = Mi*Mi'; % Variance  / covariance of MoM estimates, per unit of residual error variance

logging('size(M) = [%d %d]',size(M));
logging('Cov_MoM:'); disp(Cov_MoM);
logging('Mi*M:'); disp(Mi*M);

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

if Hessflag
    Hessmat = NaN([length(RandomEffects) length(RandomEffects) size(ymat,2)]);
else
    Hessmat = [];
end

% Control randseed here?
digits_nperms = max(ceil(log10(nperms+1)),1);
loop_timer_start = now();
for permi = 0:nperms

    permstart = now();

    if permi == 1 % Initialize perm, based on initial fit
        sig2mat_bak   = sig2mat;
        sig2tvec_bak  = sig2tvec;
        binvec_bak    = binvec;
        zmat_bak      = zmat;
        ymat_bak      = ymat;
        ymat_res_bak  = ymat_res;

        if ismember(lower(PermType), {'wildbootstrap'}) % Residual bootstrap - DEFAULT
            ymat_hat_bak = zeros(size(ymat));
        elseif ismember(lower(PermType), {'wildbootstrap-nn'}) % Non-null wild boostrap
            ymat_hat_bak = ymat_hat;
        end
    end

    if permi > 0 % Perform resampling
        if ismember(lower(PermType), {'wildbootstrap'}) || ismember(lower(PermType), {'wildbootstrap-nn'}) %DEFAULT
            for fi = 1:nfam
                % Use Rademacher distribution (-1 or 1, with equal probability) for "wild weights"
                % ymat(clusterinfo{fi}.jvec_fam,:) = ymat_hat_bak(clusterinfo{fi}.jvec_fam,:) + ...
                %                                                (2*randi(2)-3) * ymat_res_bak(clusterinfo{fi}.jvec_fam,:);

                % Use Normal distribution for "wild weights" --
                % gives really bad z-score estimates for zer-inflated covariates
                ymat(clusterinfo{fi}.jvec_fam,:) = ymat_hat_bak(clusterinfo{fi}.jvec_fam,:) + ...
                                                                randn * ymat_res_bak(clusterinfo{fi}.jvec_fam,:);
            end
        elseif ~ismember(lower(PermType), {'wildbootstrap'}) || ~ismember(lower(PermType), {'wildbootstrap-nn'})
            error('Resampling scheme not available. PermType must equal wildbootstrap or wildbootstrap-nn')
        end
    end

    %% Initially use OLS estimate
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
    beta_hat     = iXtX * (X' * ymat);
    ymat_hat     = X * beta_hat;
    ymat_res     = ymat - ymat_hat;
    sig2tvec     = sum(ymat_res.^2,1)/(size(ymat_res, 1) - size(X, 2)); % Adjust for the number of estimated parameters -- should use effective DOF instead?
    beta_se      = sqrt(diag(iXtX) * sig2tvec);
    Cov_beta     = iXtX;

    for ci = 1:size(contrasts,1)
        betacon_hat(ci,:) = contrasts(ci,:)      * beta_hat;
        betacon_se(ci, :) = sqrt(contrasts(ci,:) * Cov_beta * contrasts(ci,:)' * sig2tvec);
    end

    % Save OLS residuals for future use
    % If ymat is huge, this will lead to a large use of RAM
    if returnReusable && permi == 0
        reusableVars.ymat_res_ols = ymat_res;
        reusableVars.MSE_OLS      = sum(ymat_res.^2,1);
    end

    for iter = 1:max(1,niter)
        %% Method of moments solution
        sig2tvec = sum(ymat_res.^2,1)/(size(ymat_res,1)-size(X,2));
        LHS      = ymat_res(subvec1,:) .* ymat_res(subvec2,:) ./ mean(ymat_res.^2,1); % use normalized residuals

        if ~NonnegFlag % Standard least squares and max(0,x)
            tmp     = Mi*LHS;
            sig2mat = max(0,tmp); % Variances must be non-negative
        else
            % Use new version of lsqnonneg_amd to enfoce non-negative variances
            sig2mat = lsqnonneg_amd3(M,LHS); % This doesn't actually ensure non-negative values! -- problem with complex ymat / LHS
            sig2mat = max(0,sig2mat); % This shouldn't be needed
        end

        sig2mat     = sig2mat ./ max(eps,sum(sig2mat,1)); % Is this different from dividing by sig2tvec?
        sig2mat_mom = sig2mat;
        logLikvec   = [];

        %% Using maximum likelihood solution
        if MLflag % Phenotypes should be pre-normalized! -- now, scale is all over the place
            options_fmincon = optimoptions('fmincon','Display','off');
            logLikvec       = nan(1,size(ymat_res,2));
            [sig2mat_ml, sig2mat_ll, sig2mat_ul] = deal(nan(size(sig2mat)));
            disp(var(ymat_res));

            for coli=1:size(ymat_res, 2)
                f = @(x) (-1 * FEMA_logLik(exp(x), X, ymat_res(:, coli), clusterinfo, Ss));
                g = @(x) (-1 * FEMA_logLik(x,      X, ymat_res(:, coli), clusterinfo, Ss));
                sig2vec0 = double(sig2mat(:, coli) * sig2tvec(coli));
                fprintf(1,'Optimizing using fmincon\n');
                tic
                [sig2vec_hat, cost, exitflag, output] = fmincon(g, sig2vec0, [], [], [], [], 0*ones(size(sig2vec0)), [], [], options_fmincon);
                toc
                if 1 % exitflag<0
                    fprintf(1,'fmincon exited with exitflag = %d:\n',exitflag);
                    disp(output)
                end
                if ciflag % Compute confidence intervals on random effects?
                    fprintf(1,'Computing Confidence Interval\n');
                    tic
                    loglikthresh = chi2inv(1-0.05/2,1)/2;
                    [sig2vec_ll, sig2vec_ul] = deal(nan(size(sig2vec0)));
                    for ri = 1:length(sig2vec_hat)
                        ivec   = double((1:length(sig2vec_hat))==ri);
                        tmpfun = @(x) g(fmincon(g, (1-ivec)' .* sig2vec_hat + ivec' * x, [], [], ivec, sig2vec_hat(ri)+x, 0*ones(size(sig2vec0)), [], [], options_fmincon))-cost; % allow other parameters to change
                        dx0    = 0.01 * sum(sig2vec_hat);
                        y0     = tmpfun(dx0); % Hack to scale initial step size by variance -- phenotypes should be pre-normalized (unity variance)
                        if y0 < 0.2 % Increase scale if change in cost is too small
                            dx0 = dx0 * 4;
                            y0  = tmpfun(dx0);
                        elseif y0 > 6 % Decrease scale if change in cost is too large
                            dx0 = dx0/4;
                            y0  = tmpfun(dx0);
                        end
                        dx1  = dx0 * sqrt(2/y0);
                        y1   = tmpfun(dx1); % Get scale -- should increase dx0, if y0 is too small (or negative)
                        x    = [0 max([dx0 dx1]) * [0.5 1]];
                        y    = [0 tmpfun(x(2)) tmpfun(x(3))];
                        p    = polyfit(x, y, 2);
                        xvec = linspace(0, max(x), 101);
                        yvec = polyval(p, xvec);

                        figure(coli*10);
                        subplot(length(sig2vec_hat), 2, (ri-1)*2+2);
                        plot(xvec, yvec, x, y, '*', 'lineWidth', 2);
                        drawnow;
                        ul = sig2vec_hat(ri) + (-p(2)+((p(2)^2-4*p(1)*(p(3)-loglikthresh)))^0.5)/(2*p(1));
                        ll = 0;

                        if sig2vec_hat(ri) > 0.01 * sum(sig2vec_hat)
                            if x(end) > sig2vec_hat(ri)
                                x = x*sig2vec_hat(ri)/x(end);
                            end
                            x     = -x;
                            y     = [0 tmpfun(x(2)) tmpfun(x(3))];
                            p     = polyfit(x, y, 2);
                            xvec  = linspace(min(x), max(x), 101);
                            yvec  = polyval(p, xvec);

                            figure(coli*10);
                            subplot(length(sig2vec_hat), 2, (ri-1)*2+1);
                            plot(xvec, yvec, x, y, '*', 'lineWidth', 2);
                            drawnow;

                            ll = max(0,sig2vec_hat(ri) + (-p(2)-((p(2)^2-4*p(1)*(p(3)-loglikthresh)))^0.5)/(2*p(1)));
                        end

                        sig2vec_ll(ri) = ll;
                        sig2vec_ul(ri) = ul;
                        fprintf(1,'ri=%d: ll=%f ul=%f (%s)\n',ri, ll, ul, char(datetime));
                        if ~isreal(ll+ul) || ~isfinite(ll+ul) % Stop if result is imaginary or not finite
                            fprintf(1,'Invalid confidence interval estimates\n');
                        end
                    end
                    toc
                end
                % disp(num2str(cost,'%0.6e') )
                sig2mat_ml(:, coli) = sig2vec_hat;
                if ciflag
                    sig2mat_ll(:, coli) = sig2vec_ll;
                    sig2mat_ul(:, coli) = sig2vec_ul;
                end
                % disp(rowvec(sig2mat_ml(:, coli)/sum(sig2mat_ml(:, coli))))
                logl_ml  = g(sig2mat_ml(:, coli)); % This takes ~0.13s per column
                logl_mom = g(double(sig2mat(:, coli) * sig2tvec(coli)));
                logging('pheno %i of %i, perm %i of %i: loglike(MoM)=%.2f, loglike(ML)=%.2f', coli, size(ymat_res, 2), permi, nperms, logl_mom, logl_ml);
                logLikvec(coli) = -logl_ml;
            end
            sig2tvec_ml = sum(sig2mat_ml);
            sig2mat_ml  = sig2mat_ml ./ sig2tvec_ml;
            if ciflag
                sig2mat_ci = cat(3, sig2mat_ll, sig2mat_ul) ./ sig2tvec_ml;
            end
            sig2mat  = sig2mat_ml;
            sig2tvec = sig2tvec_ml;
        end

        %% Snap to random effects grid -- should make this a script
        nvec_bins = NaN(nsig2bins, 1);
        tvec_bins = zeros(nsig2bins, 1);
        for sig2bini = 1:nsig2bins
            tmpvec = true;
            for ri = 1:size(sig2mat,1)-1
                tmpvec = tmpvec & sig2mat(ri,:) >= sig2gridl(sig2bini, ri) ...
                                & sig2mat(ri,:) <  sig2gridu(sig2bini, ri);
            end
            ivec_bin            = find(tmpvec);
            nvec_bins(sig2bini) = length(ivec_bin);
            binvec(ivec_bin)    = sig2bini;
        end

        % If by any chance the bin was not assigned, coerce to the next
        % nearest bin
        locNaN = find(isnan(binvec));
        if ~isempty(locNaN)
            % Set this flag to 1, to enable debugging
            if 0
                keyboard; %#ok<KEYBOARDFUN,UNRCH>
            end
            for tmpBin = 1:length(locNaN)
                % For this bin, find the closest neighboring bin on the
                % grid; this is the bin where the absolute difference
                % between the estimated effects and grid values are the 
                % smallest
                tmpVals = sig2mat(1:size(sig2mat,1)-1, locNaN(tmpBin))';
                absDiff = sum(abs(tmpVals - sig2grid), 2);
                binvec(locNaN(tmpBin)) = find(absDiff == min(absDiff), 1);
            end
            warning(['Bins for the following y variables were coerced to the next nearest: ', num2str(locNaN)]);
        end

        %% Using approach described by Goldstein and Lindquist
        if strcmpi(RandomEstType, 'IGLS')
            if 0 % ~isempty(sig2mat_true)
                % Initialize with true bin, to see if that improves
                % estimates -- doesn't seem to change things much
                binvec = binvec_true; %#ok<UNRCH>
            end

            lam = [0.1 0];
            % lam = [1.0 0]; % This should be equivalent to standard MoM -- seems to work
            sig2mat_igls  = NaN(size(sig2mat));
            binvec_unique = unique(binvec(isfinite(binvec)), 'stable');
            [msemat_igls, msemat_mom, muemat_igls, muemat_mom] = deal(NaN(size(sig2mat,1), nsig2bins));
            for sig2bini = binvec_unique
                fprintf(1, 'sig2bini=%d/%d (%s)\n', sig2bini, length(binvec_unique), char(datetime));
                t0 = now; % AMD save time at beginning of computation for given bin
                ivec_bin      = find(binvec      == sig2bini);
                ivec_bin_true = find(binvec_true == sig2bini);
                sig2vec       = mean(sig2mat(:, ivec_bin), 2);
                % sig2vec_true = mean(sig2mat_true(:,ivec_bin),2);

                if GroupByFamType % Compute Vs and Vis by family type
                    Vs_famtype  = cell(1,nfamtypes);
                    Ws_famtype  = cell(1,nfamtypes);
                    for fi = 1:nfamtypes
                        ivec = find(famtypevec==fi);
                        Vs_famtype{fi} = 0;
                        for ri = 1:length(RandomEffects)
                            Vs_famtype{fi} = Vs_famtype{fi} + sig2vec(ri) * clusterinfo{ivec(1)}.(['V_', RandomEffects{ri}]);
                        end
                        Ws_famtype{fi} = pinv(FEMA_reg(FEMA_kron(Vs_famtype{fi}), lam));
                    end
                else % Compute Vs and Vis for each family
                    Vs_fam  = cell(1,nfam);
                    Ws_fam  = cell(1,nfam);
                    for fi = 1:nfam
                        Vs_fam{fi} = 0;
                        for ri = 1:length(RandomEffects)
                            Vs_fam{fi} = Vs_fam{fi} + sig2vec(ri) * clusterinfo{fi}.(['V_', RandomEffects{ri}]);
                        end
                        Ws_fam{fi} = pinv(FEMA_reg(FEMA_kron(Vs_fam{fi}),lam));
                    end
                end

                if 0 % Compare theoretical vs. empirical random effect covariance matrices
                    for fi = 1:nfamtypes %#ok<UNRCH>
                        ivec = find(famtypevec==fi);
                        fprintf(1, 'fi=%d (%d)\n', fi, length(ivec));
                        tmp_cat = [];
                        for ii = 1:length(ivec)
                            if 1 % issorted(clusterinfo{ivec(ii)}.jvec_fam)
                                tmp = LHS(clusterinfo{ivec(ii)}.ivec_fam,ivec_bin);
                                % tmp = LHS(fnumvec==ivec(ii),ivec_bin);
                                tmp_cat = cat(2,tmp_cat,tmp);
                            end
                        end
                        V = Vs_famtype{fi};
                        disp(V)
                        disp(FEMA_kron(V))
                        disp(cov(tmp_cat'))
                        pause
                    end
                end

                XtWy = 0;
                XtWX = 0;
                if GroupByFamType
                    for fi = 1:length(clusterinfo)
                        ivec_tmp = clusterinfo{fi}.ivec_fam;
                        XtWy = XtWy + M(ivec_tmp,:)' * Ws_famtype{famtypevec(fi)} * LHS(ivec_tmp, ivec_bin);
                        XtWX = XtWX + M(ivec_tmp,:)' * Ws_famtype{famtypevec(fi)} * M(ivec_tmp, :);
                    end
                else
                    for fi = 1:length(clusterinfo)
                        ivec_tmp = clusterinfo{fi}.ivec_fam;
                        XtWy = XtWy + M(ivec_tmp,:)' * Ws_fam{fi} * LHS(ivec_tmp, ivec_bin);
                        XtWX = XtWX + M(ivec_tmp,:)' * Ws_fam{fi} * M(ivec_tmp, :);
                    end
                end

                sig2mat_igls(:,ivec_bin)   = max(0, pinv(XtWX) * XtWy);
                sig2mat_igls(end,ivec_bin) = 1 - sum(sig2mat_igls(1:end-1, ivec_bin));
                sig2mat_igls(:,ivec_bin)   = sig2mat_igls(:, ivec_bin) ./ sum(sig2mat_igls(:, ivec_bin));
                if ~isempty(sig2mat_true)
                    disp(sig2grid(sig2bini,:))
                    disp('IGLS')
                    disp(sqrt(mean((sig2mat_igls(:, ivec_bin) - sig2mat_true(:, ivec_bin)).^2, 2, 'omitnan')));
                    disp('MoM')
                    disp(sqrt(mean((sig2mat_mom(:, ivec_bin) - sig2mat_true(:, ivec_bin)).^2, 2, 'omitnan')));
                    muemat_igls(:, sig2bini) = mean((sig2mat_igls(:, ivec_bin) - sig2mat_true(:, ivec_bin)),    2, 'omitnan');
                    muemat_mom(:,  sig2bini) = mean((sig2mat_mom(:,  ivec_bin) - sig2mat_true(:, ivec_bin)),    2, 'omitnan');
                    msemat_igls(:, sig2bini) = mean((sig2mat_igls(:, ivec_bin) - sig2mat_true(:, ivec_bin)).^2, 2, 'omitnan');
                    msemat_mom(:,  sig2bini) = mean((sig2mat_mom(:,  ivec_bin) - sig2mat_true(:, ivec_bin)).^2, 2, 'omitnan');
                end

                t1 = now; % AMD save time at end of computation for given bin
                tvec_bins(sig2bini) = 24*3600*(t1-t0); % Save computation time in seconds
            end % sig2bini

            if ~isempty(sig2mat_true) % sig2mat_igls is much worse than sig2mat_mom for some bins E for bins with high E, F and S, for bins with high S, low F
                disp('IGLS')
                disp(sqrt(mean((sig2mat_igls-sig2mat_true).^2, 2, 'omitnan')))
                disp('MoM')
                disp(sqrt(mean((sig2mat_mom-sig2mat_true).^2, 2, 'omitnan')))
            end

            % ax = gca; ax.XTickLabel = num2str(binvals_edges(ax.XTick)) % Should replace voxel nuumbers with binvals

            % ToDo
            %   Check if greater error in igls estimates for high S & low F are associated with error in covariance estimate
            %   Does inverse normal transformation (preserving expectancy value of moment) improve estimates?
            %   Does  inclduing higher-order moments improve things?

            sig2mat = sig2mat_igls;

            % Re-compute grid of random effects
            for sig2bini = 1:nsig2bins
                tmpvec = true;
                for ri = 1:size(sig2mat,1)-1
                    tmpvec = tmpvec & sig2mat(ri, :) >= sig2gridl(sig2bini, ri) ...
                                    & sig2mat(ri, :) <  sig2gridu(sig2bini, ri);
                end
                ivec_bin            = find(tmpvec);
                nvec_bins(sig2bini) = length(ivec_bin);
                binvec(ivec_bin)    = sig2bini;
            end
        end

        if logLikflag && ~MLflag
            logLikvec = nan(1, size(ymat_res, 2));
            for coli = 1:size(ymat_res, 2)
                % Should modify FEMA_logLik to leverage gridding of random effects?
                logLikvec(coli) = FEMA_logLik(sig2tvec(coli) * sig2mat(:, coli), X, ...
                                              ymat_res(:, coli), clusterinfo, Ss); 
            end
        end

        if iter>niter, break; end

        % Ugly hack to save resampled random effects estimates
        sig2mat_save  = sig2mat;
        sig2tvec_save = sig2tvec;

        % Save bin info
        if permi == 0
            binvec_save = binvec;
            if returnReusable
                reusableVars.binvec = binvec_save;
            end
        end

        if permi>0
            sig2tvec = sig2tvec_bak;
            sig2mat  = sig2mat_bak;
            binvec   = binvec_bak;
        end

        %% Implement GLS solution 
        % Code imported from FEMA_sig2binseg_parfeval.m
        % Some initialization
        Ws_famtype  = cell(1, nfamtypes);
        Ws_fam      = cell(1, nfam);
        [betacon_hat, betacon_se] = deal(zeros(size(contrasts,1), size(ymat,2), class(ymat)));
        [beta_hat,    beta_se]    = deal(zeros(size(X,2), size(ymat,2), class(ymat)));

        % Get ordering of fields in clusterinfo - reasonable to assume that
        % fields are always ordered in the same way since clusterinfo is
        % created in the same way across all clusters
        ff           = fieldnames(clusterinfo{1});
        [a, b]       = ismember(ff, strcat('V_', RandomEffects));
        [~, RFX_ord] = sort(b(a));
        locJVec      = strcmpi(ff, 'jvec_fam');

        % Save this ordering info for future use
        if returnReusable && permi == 0
            reusableVars.RFX_ord = RFX_ord;
            reusableVars.locJVec = locJVec;
        end

        for sig2bini = unique(binvec(isfinite(binvec)), 'stable')
            t0         = now; %#ok<*TNOW1>
            ivec_bin   = find(binvec==sig2bini);
            nvec_bins(sig2bini) = length(ivec_bin);
            sig2vec    = mean(sig2mat(:, ivec_bin), 2);

            if ~isempty(ivec_bin)
                % Handle the case of OLS
                if OLSflag
                    XtX  = X' * X;
                    if lowRank
                        if useLSQ
                            iXtX = lsqminnorm(XtX, eye(size(XtX)));
                        else
                            iXtX = pinv(XtX);
                        end
                    else
                        iXtX = XtX \ eye(size(XtX));
                    end
                    beta_hat(:, ivec_bin) = iXtX * (X' * ymat(:, ivec_bin));
                    beta_se(:,  ivec_bin) = sqrt(diag(iXtX) * sig2tvec(ivec_bin));
                    Cov_beta              = iXtX;
                else
                    if GroupByFamType
                        % Compute Vs and Vis by family type
                        for fi = 1:nfamtypes
                            ivec       = find(famtypevec == fi);
                            currClus   = struct2cell(clusterinfo{ivec(1)});
                            tmpSize    = length(currClus{locJVec});
                            Vs_famtype = zeros(tmpSize);

                            % Compute V
                            for ri = 1:length(RandomEffects)
                                Vs_famtype = Vs_famtype + sig2vec(ri) * currClus{RFX_ord(ri)};
                            end

                            % Compute inverse of V
                            Vis_famtype = double(Vs_famtype) \ eye(tmpSize, SingleOrDouble);
                            if any(isnan(Vis_famtype))
                                Vis_famtype = cast(pinv(double(Vs_famtype)), SingleOrDouble);
                            end
                            Ws_famtype{fi} = Vis_famtype;
                        end
                    else
                        % Compute Vs and Vis for each family
                        for fi = 1:nfam
                            currClus   = struct2cell(clusterinfo{fi});
                            tmpSize    = length(currClus{locJVec});
                            Vs_fam     = zeros(tmpSize);

                            % Compute V
                            for ri = 1:length(RandomEffects)
                                Vs_fam = Vs_fam + sig2vec(ri) * currClus{RFX_ord(ri)};
                            end

                            % Compute inverse of V
                            Vis_fam = double(Vs_fam) \ eye(tmpSize, SingleOrDouble);
                            if any(isnan(Vis_fam))
                                Vis_fam = cast(pinv(double(Vs_fam)), SingleOrDouble);
                            end
                            Ws_fam{fi} = Vis_fam;
                        end
                    end

                    % Compute XtW
                    XtW   = zeros(fliplr(size(X)), class(X));
                    nClus = length(clusterinfo);

                    if GroupByFamType
                        for fi = 1:nClus
                            currClus = clusterinfo{fi};
                            XtW(:, currClus.jvec_fam) = X(currClus.jvec_fam,:)' * Ws_famtype{famtypevec(fi)};
                        end
                    else
                        for fi = 1:nClus
                            currClus = clusterinfo{fi};
                            XtW(:, currClus.jvec_fam) = X(currClus.jvec_fam,:)' * Ws_fam{fi};
                        end
                    end

                    % Compute XtWX
                    B  = XtW * X;

                    % Calculate inverse of XtWX
                    if lowRank
                        if useLSQ
                            Bi = lsqminnorm(B, eye(size(B)));
                        else
                            Bi = pinv(B);
                        end
                    else
                        Bi = B \ eye(size(B));
                    end

                    % Calculate beta coefficient
                    beta_hat_tmp          = Bi * (XtW * ymat(:, ivec_bin));
                    Cov_beta              = nearestSPD(Bi);
                    beta_hat(:, ivec_bin) = beta_hat_tmp; 
                    beta_se(:,  ivec_bin) = sqrt(diag(Cov_beta) * sig2tvec(ivec_bin));
                end

                % Evaluate contrasts
                for ci = 1:size(contrasts,1)
                    betacon_hat(ci, ivec_bin) = contrasts(ci,:) * beta_hat(:,ivec_bin);
                    betacon_se(ci,  ivec_bin) = sqrt(contrasts(ci,:) * Cov_beta * contrasts(ci,:)' * sig2tvec(ivec_bin));
                end
                tvec_bins(sig2bini) = (now-t0) * 3600 * 24; % Time in seconds
            end
        end

        ymat_hat = X * beta_hat;
        ymat_res = ymat - ymat_hat;

        % Save GLS residuals
        % If ymat is huge, this will take up quite a bit of RAM
        if returnReusable && permi == 0
            reusableVars.ymat_res_gls = ymat_res;
            reusableVars.MSE_GLS      = sum(ymat_res.^2,1);
        end
    end

    if ~isempty(contrasts) % Handle non-empty betacon_hat
        beta_hat = cat(1, betacon_hat, beta_hat);
        beta_se  = cat(1, betacon_se,  beta_se);
    end

    zmat = beta_hat ./ beta_se;

    if nperms > 0
        if permi == 0
            beta_hat_perm  = NaN([size(beta_hat)  nperms + 1], class(beta_hat));
            beta_se_perm   = NaN([size(beta_se)   nperms + 1], class(beta_se));
            zmat_perm      = NaN([size(zmat)      nperms + 1], class(zmat));
            sig2mat_perm   = NaN([size(sig2mat)   nperms + 1], class(sig2mat));
            sig2tvec_perm  = NaN([size(sig2tvec)  nperms + 1], class(sig2tvec));
            logLikvec_perm = NaN([size(logLikvec) nperms + 1], class(logLikvec));
        end

        beta_hat_perm(:,  :, permi+1) = beta_hat;
        beta_se_perm(:,   :, permi+1) = beta_se;
        zmat_perm(:,      :, permi+1) = zmat;
        sig2mat_perm(:,   :, permi+1) = sig2mat_save;
        sig2tvec_perm(:,  :, permi+1) = sig2tvec_save;
        if ~isempty(logLikvec)
            logLikvec_perm(:,:,permi+1) = logLikvec;
        end
        estimated_time_remaining = (now() - loop_timer_start) * 3600 * 24/permi * (nperms - permi);
        logging('permi=%0*d/%d (%0.2fs - remaining %.0fs)', digits_nperms, permi, nperms, (now-permstart) * 3600 * 24, estimated_time_remaining);
    end
end 

if nperms>0
    beta_hat  = double(beta_hat_perm(:,:,1));
    beta_se   = double(beta_se_perm(:,:,1));
    zmat      = double(zmat_perm(:,:,1));
    sig2mat   = double(sig2mat_perm(:,:,1));
    sig2tvec  = double(sig2tvec_perm(:,:,1));
    logLikvec = double(logLikvec_perm(:,:,1));
elseif nperms == 0
    beta_hat_perm  = [];
    beta_se_perm   = [];
    zmat_perm      = [];
    sig2tvec_perm  = [];
    sig2mat_perm   = [];
    logLikvec_perm = [];
    perms          = [];
end

zmat    = double(beta_hat) ./ double(beta_se); %CHECK WITH ANDERS
logpmat = -sign(zmat) .* log10(normcdf(-abs(zmat))*2); % Should look for normcdfln function

if ciflag
    sig2mat = cat(3, sig2mat, sig2mat_ci);
end

logging('***Done*** (%0.2f seconds)\n', (now-starttime-tshim) * 3600 * 24);

%PrintMemoryUsage

return