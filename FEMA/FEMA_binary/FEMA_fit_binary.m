function [beta_hat,      beta_se,        zmat,        logpmat,              ...
          sig2tvec,      sig2mat,        Hessmat,     logLikvec,            ...
          beta_hat_perm, beta_se_perm,   zmat_perm,   sig2tvec_perm,        ...
          sig2mat_perm,  logLikvec_perm, FamilyStruct,   reusableVars] =    ...
          FEMA_fit_binary(X, iid, eid, fid, agevec, ymat, niter, contrasts, ...
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
    reusableVars.useLSQ         = useLSQ;
    reusableVars.lowRank        = lowRank;
end

% permutation initialization
if nperms>0
    beta_hat_perm = zeros(size(X,2),size(ymat,2),nperms);
    beta_se_perm = zeros(size(X,2),size(ymat,2),nperms);
    zmat_perm = zeros(size(X,2),size(ymat,2),nperms);
    sig2tvec_perm = zeros(size(X,2),size(ymat,2),nperms);
    sig2mat_perm = zeros(length(RandomEffects), size(ymat,2), nperms);
    if logLikflag
        logLikvec_perm = zeros(size(X,2),size(ymat,2),nperms);
    else
        logLikvec_perm = [];
    end
else
    [beta_hat_perm, beta_se_perm, zmat_perm, sig2tvec_perm, sig2mat_perm, logLikvec_perm] = deal([]);
end


t0 = now;

[iid_list, IA, IC_subj] = unique(iid,'stable'); nsubj = length(iid_list); nobs = length(iid);
[fid_list, IA, IC_fam]  = unique(fid,'stable'); nfam  = length(fid_list);

% Parse family structure, if necessary
if ~exist('FamilyStruct', 'var') || isempty(FamilyStruct)
    % tic
    [clusterinfo, Ss, iid, famtypevec, famtypelist, subj_famtypevec] =                  ...
     FEMA_parse_family(iid, eid, fid, agevec, pihatmat, 'RandomEffects', RandomEffects, ...
                       'FatherID', p.Results.FatherID,  'MotherID', p.Results.MotherID, ...
                       'PregID',   p.Results.PregID,    'HomeID',   p.Results.HomeID); %#ok<*ASGLU>
    
    % remove before the if-else
    % [iid_list, IA, IC_subj] = unique(iid,'stable'); nsubj = length(iid_list); nobs = length(iid);
    % [fid_list, IA, IC_fam]  = unique(fid,'stable'); nfam  = length(fid_list);
    nfamtypes = length(famtypelist);
    % toc
    
    % remove after the if-else
    % RandomVar{1} = sparse(1:nobs, IC_fam, ones(nobs,1), nobs, nfam);
    % RandomVar{2} = sparse(1:nobs, IC_subj, ones(nobs,1), nobs, nsubj);
    
    % Prepare generalized matrix version of MoM estimator
    % tic
    S_sum = Ss{1};
    for i = 2:length(Ss)
        S_sum = S_sum + Ss{i};
    end
    [subvec1, subvec2] = find(S_sum); % Use full matrix, to simplify IGLS -- should be possible to limit to tril
    %[subvec1 subvec2] = find(tril(S_sum)); % Should exclude diagonals: tril(S_sum,-1)
    indvec = sub2ind([nobs nobs],subvec1,subvec2);

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
    % indvec             = sub2ind([nobs nobs],subvec1,subvec2);

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
                          'nfamtypes',   nfamtypes,     'iid', {iid},                 ...
                          'fid',         {fid},         'iid_list', {iid_list},       ...
                          'fid_list',    {fid_list},    'nfam', nfam,                 ...
                          'subvec1',     subvec1,       'subvec2', subvec2,           ...
                          'Ss',          {Ss});
else
    clusterinfo = FamilyStruct.clusterinfo;
    M           = FamilyStruct.M;
    % nsig2bins   = FamilyStruct.nsig2bins;
    % nfam        = FamilyStruct.nfam; % nfam defined by the fid? duplicate
    famtypevec  = FamilyStruct.famtypevec;
    nfamtypes   = FamilyStruct.nfamtypes;
    % sig2grid    = FamilyStruct.sig2grid;
    % sig2gridl   = FamilyStruct.sig2gridl;
    % sig2gridu   = FamilyStruct.sig2gridu;
    subvec1     = FamilyStruct.subvec1;
    subvec2     = FamilyStruct.subvec2;
    Ss          = FamilyStruct.Ss;

    % Remove outside the if-else - duplicated
    % [iid_list, IA, IC_subj] = unique(iid,'stable'); nsubj = length(iid_list); nobs = length(iid);
    % [fid_list, IA, IC_fam]  = unique(fid,'stable'); nfam  = length(fid_list);
    % RandomVar{1} = sparse(1:nobs, IC_fam, ones(nobs,1), nobs, nfam);
    % RandomVar{2} = sparse(1:nobs, IC_subj, ones(nobs,1), nobs, nsubj);
end

% Should be generalized ------- later
RandomVar = struct();
RandomVar.("V_F") = sparse(1:nobs, IC_fam, ones(nobs,1), nobs, nfam);
RandomVar.("V_S") = sparse(1:nobs, IC_subj, ones(nobs,1), nobs, nsubj);


tshim = now-t0;

Mi = single(pinv(M));
Cov_MoM = Mi*Mi'; % Variance  / covariance of MoM estimates, per unit of residual error variance

% logging('size(M) = [%d %d]',size(M));
% logging('Cov_MoM:'); disp(Cov_MoM);
% logging('Mi*M:'); disp(Mi*M);

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
% % binvec                                  = NaN(1, size(ymat,2));
% Randombeta                              = cell(1, size(ymat,2));


for permi = 0:nperms
    if permi == 0
        ymat_current = ymat;
    else
        beta_hat_null = zeros(size(X,2), 1);
        if ~exist("sig2mat_null","var")
            X_null = ones(nobs,1);
            [~, ~,   ~,   ~,              ...
             ~, sig2mat_null, ~, ~,       ...
             ~, ~,   ~,   ~,              ...
             ~, ~,   ~,   ~,] =           ...
                 FEMA_fit_binary(X_null, iid, eid, fid, agevec, ymat, niter, ones(1,size(X,2)),    ...
                                [], 'RandomEffects', {'F','S','E'}, 'returnReusable', true,           ...
                                'RandomEstType','MoM');
        end
        ymat_current = generate_null_sample(X, iid, fid, agevec, ...
            beta_hat_null, sig2mat_null, clusterinfo, RandomEffects);
    end

    if Hessflag
        Hessmat_current = NaN([length(RandomEffects) length(RandomEffects) size(ymat_current,2)]);
    else
        Hessmat = [];
    end
    
    beta_hat_current = zeros(size(X,2), size(ymat_current,2), class(ymat_current));
    [beta_se_current, zmat_current] = deal(zeros(size(beta_hat_current), class(ymat_current)));
    [betacon_hat_current, betacon_se_current] = deal(zeros(size(contrasts,1), size(ymat_current,2), class(ymat_current)));
    
    %% Initialization
    % record should be deleted after
    % beta_record = zeros(MaxIter , size(X,2)); % used for convergence check
    % sigmat_record = zeros(MaxIter, length(RandomEffects));
    deviance_record = zeros(MaxIter, 1);

    iter = 1;
    
    converged = false;

    while ~converged && (iter <= MaxIter)
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
            beta_hat_current      = iXtX * (X' * u);

            u_res        = u - X * beta_hat_current ;
            sig2tvec_current      = sum(u_res.^2,1)/(size(u_res, 1) - size(X, 2));
        
            % Computate random variances (only fid and iid)
            [~,   sig2mat_current] =    ...
             FEMA_fit_simplified(X, iid, eid, fid, u_res, sig2tvec_current, pihatmat, W_1, "MLflag", MLflag, ...
                                "FamilyStruct", FamilyStruct, "NonnegFlag", NonnegFlag);
        
            % GLS updating fixed effects estimates
            [allWsTerms, beta_hat_current, ~, ~] = FEMA_GLS(u, X, W_1, sig2mat_current, RandomEffects, ...
                                                    clusterinfo, nfamtypes, famtypevec, ...
                                                    "GroupByFamType", GroupByFamType, "useLSQ", useLSQ);
            
            % calculate the probability according to the updated parameters
            [u_update, p, deviance] = compute_BLUP(u, X, beta_hat_current, sig2mat_current, allWsTerms, ...
                                                   RandomEffects, RandomVar, ymat_current);
            deviance_record(iter,:) = deviance;
    
            % record the parameter estimation
            beta_record(iter,:) = beta_hat_current(:,1);
            sigmat_record(iter,:) = sig2mat_current;
        
        else % for iter>0
            step = 1;
            previous_deviance = deviance_record(iter-1,:); 
            u_old = u;
            u_update_old = u_update;
            p_old = p;
            beta_old = beta_hat_current;
            sig2mat_old = sig2mat_current;
            sig2tvec_old = sig2tvec_current;
            allWsTerms_old = allWsTerms;
            
    
            W = p .* (1 - p);

            sigma = sum(sig2mat_current(1:end-1));
            h = 2 * sigma.*W./(1+sigma.*W);
            vif = 1 - 2 * (1+h/2) .* h.*(2*p-1).^2 ./(1-(2*p-1).^2);

            W = W .* vif;
            W_1 = 1 ./ W;

            % Line searching
            while true

                u = u_update_old + step * (ymat_current - p_old) ./ (p_old .* (1-p_old));
    
                % reupdate the beta_hat and sig2mat
                [~, beta_hat_current, beta_se_current, beta_cov] = FEMA_GLS(u, X, W_1, sig2mat_old, RandomEffects, ...
                                               clusterinfo, nfamtypes, famtypevec, "allWsTerms", allWsTerms, ...
                                               "GroupByFamType", GroupByFamType, "useLSQ", useLSQ);
                u_res = u - (X * beta_hat_current);
                sig2tvec_current = sum(u_res.^2,1)/(size(u_res, 1) - size(X, 2)); 
                [~, sig2mat_current] = FEMA_fit_simplified(X, iid, eid, fid, u_res, sig2tvec_current, pihatmat, W_1,  ...
                                                   "MLflag", MLflag, "FamilyStruct", FamilyStruct, "NonnegFlag", NonnegFlag);
    
                [allWsTerms, ~, ~, ~] = FEMA_GLS(u, X, W_1, sig2mat_current, RandomEffects, ...
                                                clusterinfo, nfamtypes, famtypevec, ...
                                                "GroupByFamType", GroupByFamType, "useLSQ", useLSQ, ...
                                                'GLSflag',false);
    
                % calculate the probability according to the updated parameters
                [u_update, p, new_deviance] = compute_BLUP(u, X, beta_hat_current, ...
                                                           sig2mat_current, allWsTerms, ...
                                                           RandomEffects, RandomVar, ymat_current);
    
                if all(new_deviance < previous_deviance)
                    % accept current step
                    deviance_record(iter,:) = new_deviance;
                    break;
    
                else
    
                    step = step * 0.5;
    
                    if step < 1e-5
                        % take the (t-1) state
                        u = u_old;
                        u_update = u_update_old;
                        p = p_old;
                        beta_hat_current = beta_old;
                        sig2mat_current = sig2mat_old;
                        sig2tvec = sig2tvec_old;
                        allWsTerms = allWsTerms_old;
                        deviance_record(iter,:) = previous_deviance;
                        break;
                    end
                end
            end
    
            % Convergence Check
            beta_change = norm(beta_hat_current - beta_old) / (norm(beta_old) + tol);
            sig_change = norm(sig2mat_current - sig2mat_old) / (norm(sig2mat_old) + tol);
            if beta_change <= tol && sig_change <= tol
                converged = true;
            end
    
            % beta_record(iter,:) = beta_hat_current(:,1);
            % sigmat_record(iter,:) = sig2mat_current;
        
        end
        
        iter = iter + 1;
    
    end

    % if converged
    %     fprintf('Converged at iteration %d.\n', iter-1);
    % else
    %     fprintf('Maximum iterations (%d) reached.\n', MaxIter);
    % end

    % post-hoc adjustment
    % Pearson' chi-square statistics
    sigma2_total = sum(sig2mat_current(1:end-1));
    Phi = sum((ymat_current - p).^2 ./ (p .* (1-p))) / (nobs - size(X,2));

    % Relative information offered by random effects
    W = p .* (1-p);
    eta_marginal = X * beta_hat_current; 
    p_marginal = 1 ./ (1 + exp(-eta_marginal));
    W_marginal = p_marginal .* (1 - p_marginal);
    Info_Ratio = sum(W) / sum(W_marginal);

    Phi_corrected = Phi * Info_Ratio;
    beta_se_current = beta_se_current / sqrt(Phi_corrected);


    % z-statistics
    zmat_current = double(beta_hat_current) ./ double(beta_se_current);
    logpmat_current = -log10(normcdf(-abs(zmat_current))*2);


    if permi == 0
        beta_hat = beta_hat_current;
        beta_se = beta_se_current;
        zmat = zmat_current;
        logpmat = logpmat_current;
        sig2tvec = sig2tvec_current;
        sig2mat = sig2mat_current;
        if logLikflag
            logLikvec = logLikvec_current;
        end
    else
        beta_hat_perm(:,:,permi) = beta_hat_current;
        beta_se_perm(:,:,permi) = beta_se_current;
        zmat_perm(:,:,permi) = zmat_current;
        sig2tvec_perm(:,:,permi) = sig2tvec_current;
        sig2mat_perm(:,:,permi) = sig2mat_current;
        if logLikflag
            logLikvec_perm(:,:,permi) = logLikvec_current;
        end
    end
    
    % fprintf('Permutation at permi %d. Converged at iteration %d.\n', [permi, iter-1]);



end

end


