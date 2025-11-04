% This 'hello world' demo shows how to use FEMA_fit
%
% The script generates synthetic data with n = 30000 observations, up to 5
% observations per subject, and up to 5 individuals per family; then, it
% calls FEMA_fit; finally it visualize the estimates
%
%% Requirements
% To run this example, clone the following repos and add them to
% MATLAB path with 'addpath' command.
%
% - https://github.com/cmig-research-group/cmig_utils
% - https://github.com/andersonwinkler/PALM   
%   (PALM is optional - you only needed it if you play with restricted exchangeability blocks)

%% Set seed here (if required)
% seed = 100;
% rng(seed, 'twister');

%% Simulate data
% Setup parameters for FEMA_synthesize
n = 30000;
p = 5;
v = 10;
numVisits = 5;
maxNumInFamily = 5;

% X variables
X = randn(n, p);

% Generate individual IDs
iid_int = sort(repmat(1:n, 1, numVisits));

% Generate family IDs
fid_int = ceil(iid_int / maxNumInFamily);

% Using dec2base: zero prefixes shouldn't matter
iid = cellstr([repmat('I', length(iid_int),1), dec2base(iid_int,10)]);
fid = cellstr([repmat('F', length(fid_int),1), dec2base(fid_int,10)]);
clear('iid_int'); clear('fid_int');

% Generate event IDs
allEvents = (1:numVisits)';
eid       = repmat(allEvents, n, 1);

% Now, annihilate vists
toKeep   = randperm(length(iid), n);
toDelete = setdiff(1:length(iid), toKeep);
eid(toDelete) = [];
iid(toDelete) = [];
fid(toDelete) = [];

% Older solution
% % Create FID and IID
% iid_int = sort(randsample([1:n, 1:n, 1:n, 1:n, 1:n], n));  % up to 5 observations per individual
% fid_int = ceil(iid_int / 5);                               % up to 5 individuals per family
% iid     = cell(size(iid_int)); for i=1:length(iid_int), iid{i}=sprintf('I%i', iid_int(i)); end
% fid     = cell(size(fid_int)); for i=1:length(fid_int), fid{i}=sprintf('F%i', fid_int(i)); end

% %% Shuffle observations to validate FEMA_reorder_by_families
% jvec = randperm(n); 
% X    = X(jvec, :); 
% iid  = iid(jvec); 
% fid  = fid(jvec);
% iid  = iid(:, jvec); 
% fid  = fid(:, jvec);

% Additional variables
% eid    = ones(n, 1);  % not needed
agevec = zeros(n, 1); % not needed
ymat   = nan(n, v);   % only size of ymat is used by FEMA_synthesize

% pihatmat is a square matrix with size matching the number of individuals
% (not the number of observations) - should be in the order unique(iid, 'stable')
[~, IA]  = unique(iid, 'stable');
pihatmat = eye(length(IA));
[~, ~, IC_fam] = unique(fid(IA), 'stable'); nfam=max(IC_fam);
for fi = 1:nfam
  idxvec = rowvec(find(IC_fam==fi)); % Identify all subjects for a given family
  t=triu(0.5 + 0.05 * randn(length(idxvec)), 1); 
  pihatmat(idxvec, idxvec) = pihatmat(idxvec, idxvec) + t + t';
end
%histogram(pihatmat(pihatmat>0 & pihatmat <1))

% It's hard to distinguish 'F' vs 'A' effect, so with FASE model error bars can be fairly large
% (though 'F'+'A' effect is well determined). Alternatively, use {'F','S','E'} model.
RandomEffects     = {'F', 'A', 'S', 'E'};
[~, f_effect_idx] = ismember('F', RandomEffects);
[~, a_effect_idx] = ismember('A', RandomEffects);
[~, s_effect_idx] = ismember('S', RandomEffects);
[~, e_effect_idx] = ismember('E', RandomEffects);

% [X,iid,eid,fid,agevec,~,pihatmat]   = FEMA_reorder_by_families(X,iid,eid,fid,agevec,[],pihatmat);
[ymat, sig2tvec_true, sig2mat_true] = FEMA_synthesize(X,iid,eid,fid,agevec,ymat,pihatmat,'RandomEffects',RandomEffects);

% output of FEMA_FEMA_synthesize :
% sig2tvec_true      - total variance of random effects, one value per voxel
% sig2mat_true(1, :) - normalized variance for families
% sig2mat_true(2, :) - normalized variance for repeated observations
%                      (normalized means it's between 0 and 1, in units of sig2tvec_true)

% setup parameters for FEMA_fit
contrasts = ones(1, p);
% niter     = 1;
nbins     = 20;
nperms    = 20;

%% Fit the model
[beta_hat,      beta_se,        zmat,        logpmat,               ...
 sig2tvec,      sig2mat,        Hessmat,     logLikvec,             ...
 beta_hat_perm, beta_se_perm,   zmat_perm,   sig2tvec_perm,         ...
 sig2mat_perm,  logLikvec_perm, binvec_save, nvec_bins,             ...
 tvec_bins,     FamilyStruct,   coeffCovar,  unstructParams,        ...
 residuals_GLS, info] = FEMA_fit(X, iid, eid, fid, agevec, ymat,    ...
                                 contrasts, nbins, pihatmat,        ...
                                 'RandomEffects', RandomEffects, 'nperms', nperms, 'returnResiduals', true, 'CovType', 'unstructured', 'niter', 2);

%% Visualize estimates
figure(1); clf; hold on;
plot([sig2tvec_true; sig2tvec_true], [min(sig2tvec_perm, [], 3); max(sig2tvec_perm, [], 3)], 'k.-');
plot(sig2tvec_true, sig2tvec, '*');
plot(sig2tvec_true, sig2tvec_true, '-')
xlabel('sig2tvec(true)'); ylabel('sig2tvec');

effect_idx_vec = [f_effect_idx, a_effect_idx, s_effect_idx, e_effect_idx];
for figure_idx=1:4
    figure(figure_idx+1); clf; hold on;
    effect_idx = effect_idx_vec(figure_idx);
    if effect_idx
        x = sig2mat_true(effect_idx, :);
        y = sig2mat(effect_idx, :);
        y_perm = sig2mat_perm(effect_idx, :, :);
        plot([x; x], [min(y_perm, [], 3); max(y_perm, [], 3)], 'k.-');
        plot(x, y, '*');
        plot(x, x, '-');
        xlabel(sprintf('sig2mat(true), %s', RandomEffects{effect_idx}));
        ylabel(sprintf('sig2mat, %s', RandomEffects{effect_idx}));
    end
end