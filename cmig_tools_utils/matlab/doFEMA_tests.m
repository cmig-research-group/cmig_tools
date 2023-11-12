%% Compile mean squared error between FEMA and ground truth (similarly fitlmematrix and ground truth)
% Models: FE, SE, AE, FSE, FAE, SAE, FASE
%% Setup parameters 
rng(1000, 'twister');
nObservations    = 10000;
nXvars           = 5;
nyVars           = 50;
epsmin           = 0.2; 
sigmaLow         = 0.2;
numBins          = 20;
niter            = 1;
numPerms         = 100;
eid              = ones(nObservations, 1);
agevec           = zeros(nObservations, 1);
sig2tvec_true    = 1; % Total residual variance is 1
contrasts        = ones(1,nXvars);

% Beta is uniformally distributed in the range -0.02 and 0.02 => rand
betaLow     = -0.2;
betaHigh    = 0.2;
betaTrue    = betaLow + (betaHigh-betaLow).*rand(nXvars,nyVars);

% Start parallel pool
local            = parcluster('local');
local.NumThreads = 2;
pool             = local.parpool(20, 'IdleTimeout', 240);

%% FE model: no repeated observations, up to 5 members per family
% Settings
RandomEffects_FE    = {'F', 'E'};
nFamMembers_FE      = 5;
nRepObservations_FE = 1;

% Set seed
rng(202311101, 'twister');

% Genereate IIDs
iid_int_FE = sort(randsample(repmat(1:nObservations, 1, nRepObservations_FE), nObservations));

% Generate FIDs
fid_int_FE = ceil(iid_int_FE / nFamMembers_FE);

% Initialize
iid_FE = cell(size(iid_int_FE));
fid_FE = cell(size(fid_int_FE));

% Generate individual and family IDs
for i=1:length(iid_int_FE)
    iid_FE{i}=sprintf('I%i', iid_int_FE(i));
end
for i=1:length(fid_int_FE)
    fid_FE{i}=sprintf('F%i', fid_int_FE(i));
end

% Random effects are uniformally distributed in the range 0.2 and 0.8
sig2mat_true_FE = NaN(length(RandomEffects_FE)-1, nyVars); 
ivec            = 1:nyVars;
while ~isempty(ivec)
    sig2mat_true_FE(:,ivec) = rand([size(sig2mat_true_FE,1) length(ivec)]);
    cond                    = (sum(sig2mat_true_FE,1) >= sigmaLow) & (sum(sig2mat_true_FE,1) < 1-epsmin);
    ivec                    = find(not(cond));
    logging('length(ivec)=%d',length(ivec));
    if isempty(ivec)
        break;
    end
end
sig2mat_true_FE = cat(1, sig2mat_true_FE, 1-sum(sig2mat_true_FE,1));

% Initialize a few variables
y_RFX_FE = nan(nObservations, nyVars);

% Generate X variables: X is normally distributed => randn
% Make the first X variable as the intercept
X_FE     = [ones(nObservations,1), randn(nObservations, nXvars-1)];
y_FFX_FE = X_FE * betaTrue;

% Parse family structure
clusterinfo_FE = FEMA_parse_family(iid_FE, eid, fid_FE, agevec, [], 'RandomEffects', RandomEffects_FE);

% Generate y_RFX
nfam = length(clusterinfo_FE);
for fi = 1:nfam
    if mod(fi,100)==0
        logging('fi=%d/%d',fi,nfam);
    end
    tmp = 0;
    for ri = 1:length(RandomEffects_FE)
        tmp = tmp + sqrt(sig2mat_true_FE(ri,:)) .*                        ...
                    mvnrnd(zeros(length(clusterinfo_FE{fi}.jvec_fam), 1), ...
                    double(getfield(clusterinfo_FE{fi},                   ...
        sprintf('V_%s',RandomEffects_FE{ri}))),nyVars)'; %#ok<GFLD>
    end
    y_RFX_FE(clusterinfo_FE{fi}.jvec_fam,:) = sqrt(sig2tvec_true) .* tmp;
end

% Put together as a single y variable
ymat_FE = y_FFX_FE + y_RFX_FE;

% Initialize timer for FEMA
FEMAtimeInit_FE = tic;

% Eatimate model parameters using FEMA
[beta_hat_FEMA_FE,      beta_se_FEMA_FE,        zmat_FEMA_FE,        logpmat_FEMA_FE,       ...
 sig2tvec_FEMA_FE,      sig2mat_FEMA_FE,        Hessmat_FEMA_FE,     logLikvec_FEMA_FE,     ...
 beta_hat_perm_FEMA_FE, beta_se_perm_FEMA_FE,   zmat_perm_FEMA_FE,   sig2tvec_perm_FEMA_FE, ...
 sig2mat_perm_FEMA_FE,  logLikvec_perm_FEMA_FE, binvec_save_FEMA_FE, nvec_bins_FEMA_FE,     ...
 tvec_bins_FEMA_FE,     FamilyStruct_FEMA_FE,   reusableVars_FEMA_FE] =                     ...
 FEMA_fit(X_FE, iid_FE, eid, fid_FE, agevec, ymat_FE, niter, contrasts,           ...
          numBins, [], 'RandomEffects', RandomEffects_FE, 'nperms', numPerms, 'PermType', 'wildbootstrap-nn');

% End timer for FEMA
FEMAelapsed_FE = toc(FEMAtimeInit_FE);

% Estimate model using fitlmematrix
mdls_parFE      = cell(nyVars, 1);
lme_parFE_Init  = tic;
parfor yVar = 1:nyVars
    mdls_parFE{yVar} = fitlmematrix(X_FE, ymat_FE(:,yVar), {ones(length(fid_FE),1)}, {fid_FE'});
end
lme_parFE_Elapsed = toc(lme_parFE_Init);

% Compile results and get rid of mdls
[beta_hat_fitlme_FE, beta_se_fitlme_FE, sig2mat_fitlme_FE, sig2Low_fitlme_FE, sig2Upp_fitlme_FE] = compile_fitlmematrix(mdls_parFE, nXvars, nyVars, length(RandomEffects_FE));
clear mdls_parFE

% Compile FEMA
[sig2mat_FEMA_FE_lowCI, sig2mat_FEMA_FE_uppCI] = compile_FEMA(sig2mat_perm_FEMA_FE, sig2tvec_perm_FEMA_FE, length(RandomEffects_FE), nyVars);

%% SE model: no family structure (fid = iid), up to 5 repeated observations
% Settings
RandomEffects_SE    = {'S', 'E'};
nRepObservations_SE = 5;

% Set seed
rng(202311102, 'twister');

% Genereate IIDs
iid_int_SE = sort(randsample(repmat(1:nObservations, 1, nRepObservations_SE), nObservations));

% Generate FIDs
fid_int_SE = iid_int_SE;

% Initialize
iid_SE = cell(size(iid_int_SE));
fid_SE = cell(size(fid_int_SE));

% Generate individual and family IDs
for i=1:length(iid_int_SE)
    iid_SE{i}=sprintf('I%i', iid_int_SE(i));
end
for i=1:length(fid_int_SE)
    fid_SE{i}=sprintf('F%i', fid_int_SE(i));
end

% Random effects are uniformally distributed in the range 0.2 and 0.8
sig2mat_true_SE = NaN(length(RandomEffects_SE)-1, nyVars); 
ivec            = 1:nyVars;
while ~isempty(ivec)
    sig2mat_true_SE(:,ivec) = rand([size(sig2mat_true_SE,1) length(ivec)]);
    cond                    = (sum(sig2mat_true_SE,1) >= sigmaLow) & (sum(sig2mat_true_SE,1) < 1-epsmin);
    ivec                    = find(not(cond));
    logging('length(ivec)=%d',length(ivec));
    if isempty(ivec)
        break;
    end
end
sig2mat_true_SE = cat(1, sig2mat_true_SE, 1-sum(sig2mat_true_SE,1));

% Initialize a few variables
y_RFX_SE = nan(nObservations, nyVars);

% Generate X variables: X is normally distributed => randn
% Make the first X variable as the intercept
X_SE     = [ones(nObservations,1), randn(nObservations, nXvars-1)];
y_FFX_SE = X_SE * betaTrue;

% Parse family structure
clusterinfo_SE = FEMA_parse_family(iid_SE, eid, fid_SE, agevec, [], 'RandomEffects', RandomEffects_SE);

% Generate y_RFX
nfam = length(clusterinfo_SE);
for fi = 1:nfam
    if mod(fi,100)==0
        logging('fi=%d/%d',fi,nfam);
    end
    tmp = 0;
    for ri = 1:length(RandomEffects_SE)
        tmp = tmp + sqrt(sig2mat_true_SE(ri,:)) .*                        ...
                    mvnrnd(zeros(length(clusterinfo_SE{fi}.jvec_fam), 1), ...
                    double(getfield(clusterinfo_SE{fi},                   ...
        sprintf('V_%s',RandomEffects_SE{ri}))),nyVars)'; %#ok<GFLD>
    end
    y_RFX_SE(clusterinfo_SE{fi}.jvec_fam,:) = sqrt(sig2tvec_true) .* tmp;
end

% Put together as a single y variable
ymat_SE = y_FFX_SE + y_RFX_SE;

% Initialize timer for FEMA
FEMAtimeInit_SE = tic;

% Eatimate model parameters using FEMA
[beta_hat_FEMA_SE,      beta_se_SEMA_SE,        zmat_FEMA_SE,        logpmat_FEMA_SE,       ...
 sig2tvec_FEMA_SE,      sig2mat_FEMA_SE,        Hessmat_FEMA_SE,     logLikvec_FEMA_SE,     ...
 beta_hat_perm_FEMA_SE, beta_se_perm_FEMA_SE,   zmat_perm_FEMA_SE,   sig2tvec_perm_FEMA_SE, ...
 sig2mat_perm_FEMA_SE,  logLikvec_perm_FEMA_SE, binvec_save_FEMA_SE, nvec_bins_FEMA_SE,     ...
 tvec_bins_FEMA_SE,     FamilyStruct_FEMA_SE,   reusableVars_FEMA_SE] =                     ...
 FEMA_fit(X_SE, iid_SE, eid, fid_SE, agevec, ymat_SE, niter, contrasts,           ...
          numBins, [], 'RandomEffects', RandomEffects_SE, 'nperms', numPerms, 'PermType', 'wildbootstrap-nn');

% End timer for FEMA
FEMAelapsed_SE = toc(FEMAtimeInit_SE);

% Estimate model using fitlmematrix
mdls_parSE      = cell(nyVars, 1);
lme_parSE_Init  = tic;
parfor yVar = 1:nyVars
    mdls_parSE{yVar} = fitlmematrix(X_SE, ymat_SE(:,yVar), {ones(length(iid_SE),1)}, {iid_SE'});
end
lme_parSE_Elapsed = toc(lme_parSE_Init);

% Compile results and get rid of mdls
[beta_hat_fitlme_SE, beta_se_fitlme_SE, sig2mat_fitlme_SE, sig2Low_fitlme_SE, sig2Upp_fitlme_SE] = compile_fitlmematrix(mdls_parSE, nXvars, nyVars, length(RandomEffects_SE));
clear mdls_parSE

% Compile FEMA
[sig2mat_FEMA_SE_lowCI, sig2mat_FEMA_SE_uppCI] = compile_FEMA(sig2mat_perm_FEMA_SE, sig2tvec_perm_FEMA_SE, length(RandomEffects_SE), nyVars);

%% AE model: no repeated observations but A
% Settings
RandomEffects_AE    = {'A', 'E'};
nFamMembers_AE      = 5;
nRepObservations_AE = 1;

% Set seed
rng(202311103, 'twister');

% Genereate IIDs
iid_int_AE = sort(randsample(repmat(1:nObservations, 1, nRepObservations_AE), nObservations));

% Generate FIDs
fid_int_AE = ceil(iid_int_AE / nFamMembers_AE);

% Initialize
iid_AE = cell(size(iid_int_AE));
fid_AE = cell(size(fid_int_AE));

% Generate individual and family IDs
for i=1:length(iid_int_AE)
    iid_AE{i}=sprintf('I%i', iid_int_AE(i));
end
for i=1:length(fid_int_AE)
    fid_AE{i}=sprintf('F%i', fid_int_AE(i));
end

% Random effects are uniformally distributed in the range 0.2 and 0.8
sig2mat_true_AE = NaN(length(RandomEffects_AE)-1, nyVars); 
ivec            = 1:nyVars;
while ~isempty(ivec)
    sig2mat_true_AE(:,ivec) = rand([size(sig2mat_true_AE,1) length(ivec)]);
    cond                    = (sum(sig2mat_true_AE,1) >= sigmaLow) & (sum(sig2mat_true_AE,1) < 1-epsmin);
    ivec                    = find(not(cond));
    logging('length(ivec)=%d',length(ivec));
    if isempty(ivec)
        break;
    end
end
sig2mat_true_AE = cat(1, sig2mat_true_AE, 1-sum(sig2mat_true_AE,1));

% Initialize a few variables
y_RFX_AE = nan(nObservations, nyVars);

% Generate X variables: X is normally distributed => randn
% Make the first X variable as the intercept
X_AE     = [ones(nObservations,1), randn(nObservations, nXvars-1)];
y_FFX_AE = X_AE * betaTrue;

% pihatmat is a square matrix with size matching the number of
% individuals (not the number of observations)
[~, IA]         = unique(iid_AE, 'stable');
pihatmat_AE     = eye(length(IA));
[~, ~, IC_fam]  = unique(fid_AE(IA), 'stable');
nfam            = max(IC_fam);
for fi          = 1:nfam
    idxvec      = rowvec(find(IC_fam==fi)); % Identify all subjects for a given family
    t           = triu(0.5 + 0.05 * randn(length(idxvec)), 1);
    pihatmat_AE(idxvec, idxvec) = pihatmat_AE(idxvec, idxvec) + t + t';
end

% Parse family structure
clusterinfo_AE = FEMA_parse_family(iid_AE, eid, fid_AE, agevec, pihatmat_AE, 'RandomEffects', RandomEffects_AE);

% Generate y_RFX
nfam = length(clusterinfo_AE);
for fi = 1:nfam
    if mod(fi,100)==0
        logging('fi=%d/%d',fi,nfam);
    end
    tmp = 0;
    for ri = 1:length(RandomEffects_AE)
        tmp = tmp + sqrt(sig2mat_true_AE(ri,:)) .*                        ...
                    mvnrnd(zeros(length(clusterinfo_AE{fi}.jvec_fam), 1), ...
                    double(getfield(clusterinfo_AE{fi},                   ...
        sprintf('V_%s',RandomEffects_AE{ri}))),nyVars)'; %#ok<GFLD>
    end
    y_RFX_AE(clusterinfo_AE{fi}.jvec_fam,:) = sqrt(sig2tvec_true) .* tmp;
end

% Put together as a single y variable
ymat_AE = y_FFX_AE + y_RFX_AE;

% Initialize timer for FEMA
FEMAtimeInit_AE = tic;

% Eatimate model parameters using FEMA
[beta_hat_FEMA_AE,      beta_se_SEMA_AE,        zmat_FEMA_AE,        logpmat_FEMA_AE,       ...
 sig2tvec_FEMA_AE,      sig2mat_FEMA_AE,        Hessmat_FEMA_AE,     logLikvec_FEMA_AE,     ...
 beta_hat_perm_FEMA_AE, beta_se_perm_FEMA_AE,   zmat_perm_FEMA_AE,   sig2tvec_perm_FEMA_AE, ...
 sig2mat_perm_FEMA_AE,  logLikvec_perm_FEMA_AE, binvec_save_FEMA_AE, nvec_bins_FEMA_AE,     ...
 tvec_bins_FEMA_AE,     FamilyStruct_FEMA_AE,   reusableVars_FEMA_AE] =                     ...
 FEMA_fit(X_AE, iid_AE, eid, fid_AE, agevec, ymat_AE, niter, contrasts,                     ...
          numBins, pihatmat_AE, 'RandomEffects', RandomEffects_AE, 'nperms', numPerms, 'PermType', 'wildbootstrap-nn');

% End timer for FEMA
FEMAelapsed_AE = toc(FEMAtimeInit_AE);

% Compile FEMA
[sig2mat_FEMA_AE_lowCI, sig2mat_FEMA_AE_uppCI] = compile_FEMA(sig2mat_perm_FEMA_AE, sig2tvec_perm_FEMA_AE, length(RandomEffects_AE), nyVars);

%% FSE model: family design, repeated observations
% Settings
RandomEffects_FSE    = {'F', 'S', 'E'};
nFamMembers_FSE      = 5;
nRepObservations_FSE = 5;

% Set seed
rng(202311104, 'twister');

% Genereate IIDs
iid_int_FSE = sort(randsample(repmat(1:nObservations, 1, nRepObservations_FSE), nObservations));

% Generate FIDs
fid_int_FSE = ceil(iid_int_FSE / nFamMembers_FSE);

% Initialize
iid_FSE = cell(size(iid_int_FSE));
fid_FSE = cell(size(fid_int_FSE));

% Generate individual and family IDs
for i=1:length(iid_int_FSE)
    iid_FSE{i}=sprintf('I%i', iid_int_FSE(i));
end
for i=1:length(fid_int_FSE)
    fid_FSE{i}=sprintf('F%i', fid_int_FSE(i));
end

% Random effects are uniformally distributed in the range 0.2 and 0.8
sig2mat_true_FSE = NaN(length(RandomEffects_FSE)-1, nyVars); 
ivec             = 1:nyVars;
while ~isempty(ivec)
    sig2mat_true_FSE(:,ivec) = rand([size(sig2mat_true_FSE,1) length(ivec)]);
    cond                    = (sum(sig2mat_true_FSE,1) >= sigmaLow) & (sum(sig2mat_true_FSE,1) < 1-epsmin);
    ivec                    = find(not(cond));
    logging('length(ivec)=%d',length(ivec));
    if isempty(ivec)
        break;
    end
end
sig2mat_true_FSE = cat(1, sig2mat_true_FSE, 1-sum(sig2mat_true_FSE,1));

% Initialize a few variables
y_RFX_FSE = nan(nObservations, nyVars);

% Generate X variables: X is normally distributed => randn
% Make the first X variable as the intercept
X_FSE     = [ones(nObservations,1), randn(nObservations, nXvars-1)];
y_FFX_FSE = X_FSE * betaTrue;

% Parse family structure
clusterinfo_FSE = FEMA_parse_family(iid_FSE, eid, fid_FSE, agevec, [], 'RandomEffects', RandomEffects_FSE);

% Generate y_RFX
nfam = length(clusterinfo_FSE);
for fi = 1:nfam
    if mod(fi,100)==0
        logging('fi=%d/%d',fi,nfam);
    end
    tmp = 0;
    for ri = 1:length(RandomEffects_FSE)
        tmp = tmp + sqrt(sig2mat_true_FSE(ri,:)) .*                        ...
                    mvnrnd(zeros(length(clusterinfo_FSE{fi}.jvec_fam), 1), ...
                    double(getfield(clusterinfo_FSE{fi},                   ...
        sprintf('V_%s',RandomEffects_FSE{ri}))),nyVars)'; %#ok<GFLD>
    end
    y_RFX_FSE(clusterinfo_FSE{fi}.jvec_fam,:) = sqrt(sig2tvec_true) .* tmp;
end

% Put together as a single y variable
ymat_FSE = y_FFX_FSE + y_RFX_FSE;

% Initialize timer for FEMA
FEMAtimeInit_FSE = tic;

% Eatimate model parameters using FEMA
[beta_hat_FEMA_FSE,      beta_se_SEMA_FSE,        zmat_FEMA_FSE,        logpmat_FEMA_FSE,       ...
 sig2tvec_FEMA_FSE,      sig2mat_FEMA_FSE,        Hessmat_FEMA_FSE,     logLikvec_FEMA_FSE,     ...
 beta_hat_perm_FEMA_FSE, beta_se_perm_FEMA_FSE,   zmat_perm_FEMA_FSE,   sig2tvec_perm_FEMA_FSE, ...
 sig2mat_perm_FEMA_FSE,  logLikvec_perm_FEMA_FSE, binvec_save_FEMA_FSE, nvec_bins_FEMA_FSE,     ...
 tvec_bins_FEMA_FSE,     FamilyStruct_FEMA_FSE,   reusableVars_FEMA_FSE] =                     ...
 FEMA_fit(X_FSE, iid_FSE, eid, fid_FSE, agevec, ymat_FSE, niter, contrasts,           ...
          numBins, [], 'RandomEffects', RandomEffects_FSE, 'nperms', numPerms, 'PermType', 'wildbootstrap-nn');

% End timer for FEMA
FEMAelapsed_FSE = toc(FEMAtimeInit_FSE);

% Estimate model using fitlmematrix
mdls_parFSE      = cell(nyVars, 1);
lme_parFSE_Init  = tic;
parfor yVar = 1:nyVars
    mdls_parFSE{yVar} = fitlmematrix(X_FSE, ymat_FSE(:,yVar), [{ones(length(fid_FSE),1)}, {ones(length(iid_FSE),1)}], [{fid_FSE'}, {iid_FSE'}]);
end
lme_parFSE_Elapsed = toc(lme_parFSE_Init);

% Compile results and get rid of mdls
[beta_hat_fitlme_FSE, beta_se_fitlme_FSE, sig2mat_fitlme_FSE, sig2Low_fitlme_FSE, sig2Upp_fitlme_FSE] = compile_fitlmematrix(mdls_parFSE, nXvars, nyVars, length(RandomEffects_FSE));
clear mdls_parFSE

% Compile FEMA
[sig2mat_FEMA_FSE_lowCI, sig2mat_FEMA_FSE_uppCI] = compile_FEMA(sig2mat_perm_FEMA_FSE, sig2tvec_perm_FEMA_FSE, length(RandomEffects_FSE), nyVars);

%% ASE model: family design in the sense of A and repeated observations
% Settings
RandomEffects_SAE    = {'A', 'S', 'E'};
nRepObservations_SAE = 5;

% Set seed
rng(202311105, 'twister');

% Genereate IIDs
iid_int_SAE = sort(randsample(repmat(1:nObservations, 1, nRepObservations_SAE), nObservations));

% Generate FIDs = IIDs
fid_int_SAE = iid_int_SAE;

% Initialize
iid_SAE = cell(size(iid_int_SAE));
fid_SAE = cell(size(fid_int_SAE));

% Generate individual and family IDs
for i=1:length(iid_int_SAE)
    iid_SAE{i}=sprintf('I%i', iid_int_SAE(i));
end
for i=1:length(fid_int_SAE)
    fid_SAE{i}=sprintf('F%i', fid_int_SAE(i));
end

% Random effects are uniformally distributed in the range 0.2 and 0.8
sig2mat_true_SAE = NaN(length(RandomEffects_SAE)-1, nyVars); 
ivec             = 1:nyVars;
while ~isempty(ivec)
    sig2mat_true_SAE(:,ivec) = rand([size(sig2mat_true_SAE,1) length(ivec)]);
    cond                    = (sum(sig2mat_true_SAE,1) >= sigmaLow) & (sum(sig2mat_true_SAE,1) < 1-epsmin);
    ivec                    = find(not(cond));
    logging('length(ivec)=%d',length(ivec));
    if isempty(ivec)
        break;
    end
end
sig2mat_true_SAE = cat(1, sig2mat_true_SAE, 1-sum(sig2mat_true_SAE,1));

% Initialize a few variables
y_RFX_SAE = nan(nObservations, nyVars);

% Generate X variables: X is normally distributed => randn
% Make the first X variable as the intercept
X_SAE     = [ones(nObservations,1), randn(nObservations, nXvars-1)];
y_FFX_SAE = X_SAE * betaTrue;

% pihatmat is a square matrix with size matching the number of
% individuals (not the number of observations)
[~, IA]         = unique(iid_SAE, 'stable');
pihatmat_SAE    = eye(length(IA));
[~, ~, IC_fam]  = unique(fid_SAE(IA), 'stable');
nfam            = max(IC_fam);
for fi          = 1:nfam
    idxvec      = rowvec(find(IC_fam==fi)); % Identify all subjects for a given family
    t           = triu(0.5 + 0.05 * randn(length(idxvec)), 1);
    pihatmat_SAE(idxvec, idxvec) = pihatmat_SAE(idxvec, idxvec) + t + t';
end

% Parse family structure
clusterinfo_SAE = FEMA_parse_family(iid_SAE, eid, fid_SAE, agevec, pihatmat_SAE, 'RandomEffects', RandomEffects_SAE);

% Generate y_RFX
nfam = length(clusterinfo_SAE);
for fi = 1:nfam
    if mod(fi,100)==0
        logging('fi=%d/%d',fi,nfam);
    end
    tmp = 0;
    for ri = 1:length(RandomEffects_SAE)
        tmp = tmp + sqrt(sig2mat_true_SAE(ri,:)) .*                        ...
                    mvnrnd(zeros(length(clusterinfo_SAE{fi}.jvec_fam), 1), ...
                    double(getfield(clusterinfo_SAE{fi},                   ...
        sprintf('V_%s',RandomEffects_SAE{ri}))),nyVars)'; %#ok<GFLD>
    end
    y_RFX_SAE(clusterinfo_SAE{fi}.jvec_fam,:) = sqrt(sig2tvec_true) .* tmp;
end

% Put together as a single y variable
ymat_SAE = y_FFX_SAE + y_RFX_SAE;

% Initialize timer for FEMA
FEMAtimeInit_SAE = tic;

% Eatimate model parameters using FEMA
[beta_hat_FEMA_SAE,      beta_se_SEMA_SAE,        zmat_FEMA_SAE,        logpmat_FEMA_SAE,       ...
 sig2tvec_FEMA_SAE,      sig2mat_FEMA_SAE,        Hessmat_FEMA_SAE,     logLikvec_FEMA_SAE,     ...
 beta_hat_perm_FEMA_SAE, beta_se_perm_FEMA_SAE,   zmat_perm_FEMA_SAE,   sig2tvec_perm_FEMA_SAE, ...
 sig2mat_perm_FEMA_SAE,  logLikvec_perm_FEMA_SAE, binvec_save_FEMA_SAE, nvec_bins_FEMA_SAE,     ...
 tvec_bins_FEMA_SAE,     FamilyStruct_FEMA_SAE,   reusableVars_FEMA_SAE] =                     ...
 FEMA_fit(X_SAE, iid_SAE, eid, fid_SAE, agevec, ymat_SAE, niter, contrasts,           ...
          numBins, pihatmat_SAE, 'RandomEffects', RandomEffects_SAE, 'nperms', numPerms, 'PermType', 'wildbootstrap-nn');

% End timer for FEMA
FEMAelapsed_SAE = toc(FEMAtimeInit_SAE);

% Compile FEMA
[sig2mat_FEMA_SAE_lowCI, sig2mat_FEMA_SAE_uppCI] = compile_FEMA(sig2mat_perm_FEMA_SAE, sig2tvec_perm_FEMA_SAE, length(RandomEffects_SAE), nyVars);

%% FAE model: family design, A, and no repeated observations
% Settings
RandomEffects_FAE    = {'F', 'A', 'E'};
nFamMembers_FAE      = 5;
nRepObservations_FAE = 1;

% Set seed
rng(202311106, 'twister');

% Genereate IIDs
iid_int_FAE = sort(randsample(repmat(1:nObservations, 1, nRepObservations_FAE), nObservations));

% Generate FIDs
fid_int_FAE = ceil(iid_int_FAE / nFamMembers_FAE);

% Initialize
iid_FAE = cell(size(iid_int_FAE));
fid_FAE = cell(size(fid_int_FAE));

% Generate individual and family IDs
for i=1:length(iid_int_FAE)
    iid_FAE{i}=sprintf('I%i', iid_int_FAE(i));
end
for i=1:length(fid_int_FAE)
    fid_FAE{i}=sprintf('F%i', fid_int_FAE(i));
end

% Random effects are uniformally distributed in the range 0.2 and 0.8
sig2mat_true_FAE = NaN(length(RandomEffects_FAE)-1, nyVars); 
ivec             = 1:nyVars;
while ~isempty(ivec)
    sig2mat_true_FAE(:,ivec) = rand([size(sig2mat_true_FAE,1) length(ivec)]);
    cond                    = (sum(sig2mat_true_FAE,1) >= sigmaLow) & (sum(sig2mat_true_FAE,1) < 1-epsmin);
    ivec                    = find(not(cond));
    logging('length(ivec)=%d',length(ivec));
    if isempty(ivec)
        break;
    end
end
sig2mat_true_FAE = cat(1, sig2mat_true_FAE, 1-sum(sig2mat_true_FAE,1));

% Initialize a few variables
y_RFX_FAE = nan(nObservations, nyVars);

% Generate X variables: X is normally distributed => randn
% Make the first X variable as the intercept
X_FAE     = [ones(nObservations,1), randn(nObservations, nXvars-1)];
y_FFX_FAE = X_FAE * betaTrue;

% pihatmat is a square matrix with size matching the number of
% individuals (not the number of observations)
[~, IA]         = unique(iid_FAE, 'stable');
pihatmat_FAE    = eye(length(IA));
[~, ~, IC_fam]  = unique(fid_FAE(IA), 'stable');
nfam            = max(IC_fam);
for fi          = 1:nfam
    idxvec      = rowvec(find(IC_fam==fi)); % Identify all subjects for a given family
    t           = triu(0.5 + 0.05 * randn(length(idxvec)), 1);
    pihatmat_FAE(idxvec, idxvec) = pihatmat_FAE(idxvec, idxvec) + t + t';
end

% Parse family structure
clusterinfo_FAE = FEMA_parse_family(iid_FAE, eid, fid_FAE, agevec, pihatmat_FAE, 'RandomEffects', RandomEffects_FAE);

% Generate y_RFX
nfam = length(clusterinfo_FAE);
for fi = 1:nfam
    if mod(fi,100)==0
        logging('fi=%d/%d',fi,nfam);
    end
    tmp = 0;
    for ri = 1:length(RandomEffects_FAE)
        tmp = tmp + sqrt(sig2mat_true_FAE(ri,:)) .*                        ...
                    mvnrnd(zeros(length(clusterinfo_FAE{fi}.jvec_fam), 1), ...
                    double(getfield(clusterinfo_FAE{fi},                   ...
        sprintf('V_%s',RandomEffects_FAE{ri}))),nyVars)'; %#ok<GFLD>
    end
    y_RFX_FAE(clusterinfo_FAE{fi}.jvec_fam,:) = sqrt(sig2tvec_true) .* tmp;
end

% Put together as a single y variable
ymat_FAE = y_FFX_FAE + y_RFX_FAE;

% Initialize timer for FEMA
FEMAtimeInit_FAE = tic;

% Eatimate model parameters using FEMA
[beta_hat_FEMA_FAE,      beta_se_SEMA_FAE,        zmat_FEMA_FAE,        logpmat_FEMA_FAE,       ...
 sig2tvec_FEMA_FAE,      sig2mat_FEMA_FAE,        Hessmat_FEMA_FAE,     logLikvec_FEMA_FAE,     ...
 beta_hat_perm_FEMA_FAE, beta_se_perm_FEMA_FAE,   zmat_perm_FEMA_FAE,   sig2tvec_perm_FEMA_FAE, ...
 sig2mat_perm_FEMA_FAE,  logLikvec_perm_FEMA_FAE, binvec_save_FEMA_FAE, nvec_bins_FEMA_FAE,     ...
 tvec_bins_FEMA_FAE,     FamilyStruct_FEMA_FAE,   reusableVars_FEMA_FAE] =                     ...
 FEMA_fit(X_FAE, iid_FAE, eid, fid_FAE, agevec, ymat_FAE, niter, contrasts,           ...
          numBins, pihatmat_FAE, 'RandomEffects', RandomEffects_FAE, 'nperms', numPerms, 'PermType', 'wildbootstrap-nn');

% End timer for FEMA
FEMAelapsed_FAE = toc(FEMAtimeInit_FAE);

% Compile FEMA
[sig2mat_FEMA_FAE_lowCI, sig2mat_FEMA_FAE_uppCI] = compile_FEMA(sig2mat_perm_FEMA_FAE, sig2tvec_perm_FEMA_FAE, length(RandomEffects_FAE), nyVars);

%% FASE model: family design, GRM, and repeated observations
% Settings
RandomEffects_FASE    = {'F', 'A', 'S', 'E'};
nFamMembers_FASE      = 5;
nRepObservations_FASE = 5;

% Set seed
rng(202311107, 'twister');

% Genereate IIDs
iid_int_FASE = sort(randsample(repmat(1:nObservations, 1, nRepObservations_FASE), nObservations));

% Generate FIDs
fid_int_FASE = ceil(iid_int_FASE / nFamMembers_FASE);

% Initialize
iid_FASE = cell(size(iid_int_FASE));
fid_FASE = cell(size(fid_int_FASE));

% Generate individual and family IDs
for i=1:length(iid_int_FASE)
    iid_FASE{i}=sprintf('I%i', iid_int_FASE(i));
end
for i=1:length(fid_int_FASE)
    fid_FASE{i}=sprintf('F%i', fid_int_FASE(i));
end

% Random effects are uniformally distributed in the range 0.2 and 0.8
sig2mat_true_FASE = NaN(length(RandomEffects_FASE)-1, nyVars); 
ivec              = 1:nyVars;
while ~isempty(ivec)
    sig2mat_true_FASE(:,ivec) = rand([size(sig2mat_true_FASE,1) length(ivec)]);
    cond                    = (sum(sig2mat_true_FASE,1) >= sigmaLow) & (sum(sig2mat_true_FASE,1) < 1-epsmin);
    ivec                    = find(not(cond));
    logging('length(ivec)=%d',length(ivec));
    if isempty(ivec)
        break;
    end
end
sig2mat_true_FASE = cat(1, sig2mat_true_FASE, 1-sum(sig2mat_true_FASE,1));

% Initialize a few variables
y_RFX_FASE = nan(nObservations, nyVars);

% Generate X variables: X is normally distributed => randn
% Make the first X variable as the intercept
X_FASE     = [ones(nObservations,1), randn(nObservations, nXvars-1)];
y_FFX_FASE = X_FASE * betaTrue;

% pihatmat is a square matrix with size matching the number of
% individuals (not the number of observations)
[~, IA]         = unique(iid_FASE, 'stable');
pihatmat_FASE   = eye(length(IA));
[~, ~, IC_fam]  = unique(fid_FASE(IA), 'stable');
nfam            = max(IC_fam);
for fi          = 1:nfam
    idxvec      = rowvec(find(IC_fam==fi)); % Identify all subjects for a given family
    t           = triu(0.5 + 0.05 * randn(length(idxvec)), 1);
    pihatmat_FASE(idxvec, idxvec) = pihatmat_FASE(idxvec, idxvec) + t + t';
end

% Parse family structure
clusterinfo_FASE = FEMA_parse_family(iid_FASE, eid, fid_FASE, agevec, pihatmat_FASE, 'RandomEffects', RandomEffects_FASE);

% Generate y_RFX
nfam = length(clusterinfo_FASE);
for fi = 1:nfam
    if mod(fi,100)==0
        logging('fi=%d/%d',fi,nfam);
    end
    tmp = 0;
    for ri = 1:length(RandomEffects_FASE)
        tmp = tmp + sqrt(sig2mat_true_FASE(ri,:)) .*                        ...
                    mvnrnd(zeros(length(clusterinfo_FASE{fi}.jvec_fam), 1), ...
                    double(getfield(clusterinfo_FASE{fi},                   ...
        sprintf('V_%s',RandomEffects_FASE{ri}))),nyVars)'; %#ok<GFLD>
    end
    y_RFX_FASE(clusterinfo_FASE{fi}.jvec_fam,:) = sqrt(sig2tvec_true) .* tmp;
end

% Put together as a single y variable
ymat_FASE = y_FFX_FASE + y_RFX_FASE;

% Initialize timer for FEMA
FEMAtimeInit_FASE = tic;

% Eatimate model parameters using FEMA
[beta_hat_FEMA_FASE,      beta_se_SEMA_FASE,        zmat_FEMA_FASE,        logpmat_FEMA_FASE,       ...
 sig2tvec_FEMA_FASE,      sig2mat_FEMA_FASE,        Hessmat_FEMA_FASE,     logLikvec_FEMA_FASE,     ...
 beta_hat_perm_FEMA_FASE, beta_se_perm_FEMA_FASE,   zmat_perm_FEMA_FASE,   sig2tvec_perm_FEMA_FASE, ...
 sig2mat_perm_FEMA_FASE,  logLikvec_perm_FEMA_FASE, binvec_save_FEMA_FASE, nvec_bins_FEMA_FASE,     ...
 tvec_bins_FEMA_FASE,     FamilyStruct_FEMA_FASE,   reusableVars_FEMA_FASE] =                     ...
 FEMA_fit(X_FASE, iid_FASE, eid, fid_FASE, agevec, ymat_FASE, niter, contrasts,           ...
          numBins, pihatmat_FASE, 'RandomEffects', RandomEffects_FASE, 'nperms', numPerms, 'PermType', 'wildbootstrap-nn');

% End timer for FEMA
FEMAelapsed_FASE = toc(FEMAtimeInit_FASE);

% Compile FEMA
[sig2mat_FEMA_FASE_lowCI, sig2mat_FEMA_FASE_uppCI] = compile_FEMA(sig2mat_perm_FEMA_FASE, sig2tvec_perm_FEMA_FASE, length(RandomEffects_FASE), nyVars);

%% Close parallel pool
delete(pool);

%% Examine results: squared errors - Fixed effects
% FEMA family
SqErr_FEMA_FFX_FE   = (betaTrue - beta_hat_FEMA_FE(2:end,:)).^2;
SqErr_FEMA_FFX_SE   = (betaTrue - beta_hat_FEMA_SE(2:end,:)).^2;
SqErr_FEMA_FFX_AE   = (betaTrue - beta_hat_FEMA_AE(2:end,:)).^2;
SqErr_FEMA_FFX_FSE  = (betaTrue - beta_hat_FEMA_FSE(2:end,:)).^2;
SqErr_FEMA_FFX_SAE  = (betaTrue - beta_hat_FEMA_SAE(2:end,:)).^2;
SqErr_FEMA_FFX_FAE  = (betaTrue - beta_hat_FEMA_FAE(2:end,:)).^2;
SqErr_FEMA_FFX_FASE = (betaTrue - beta_hat_FEMA_FASE(2:end,:)).^2;

% LME errors
SqErr_LME_FFX_FE   = (betaTrue - beta_hat_fitlme_FE).^2;
SqErr_LME_FFX_SE   = (betaTrue - beta_hat_fitlme_SE).^2;
SqErr_LME_FFX_FSE  = (betaTrue - beta_hat_fitlme_FSE).^2;

% Difference between estimates from LME and FEMA
SqDiff_FEMA_LME_FFX_FE   = (beta_hat_FEMA_FE(2:end,:)  - beta_hat_fitlme_FE).^2;
SqDiff_FEMA_LME_FFX_SE   = (beta_hat_FEMA_SE(2:end,:)  - beta_hat_fitlme_SE).^2;
SqDiff_FEMA_LME_FFX_FSE  = (beta_hat_FEMA_FSE(2:end,:) - beta_hat_fitlme_FSE).^2;

%% Examine results: sum of squared errors (sum across y variables) - Random effects
% Do we really need to square in case of RFX? They are always going to be
% positive quantities...
% FEMA family (remember to scale using sig2tvec)
SqErr_FEMA_RFX_FE   = sum((sig2mat_true_FE   - sig2mat_FEMA_FE   .* sig2tvec_FEMA_FE).^2,   2);
SqErr_FEMA_RFX_SE   = sum((sig2mat_true_SE   - sig2mat_FEMA_SE   .* sig2tvec_FEMA_SE).^2,   2);
SqErr_FEMA_RFX_AE   = sum((sig2mat_true_AE   - sig2mat_FEMA_AE   .* sig2tvec_FEMA_AE).^2,   2);
SqErr_FEMA_RFX_FSE  = sum((sig2mat_true_FSE  - sig2mat_FEMA_FSE  .* sig2tvec_FEMA_FSE).^2,  2);
SqErr_FEMA_RFX_SAE  = sum((sig2mat_true_SAE  - sig2mat_FEMA_SAE  .* sig2tvec_FEMA_SAE).^2,  2);
SqErr_FEMA_RFX_FAE  = sum((sig2mat_true_FAE  - sig2mat_FEMA_FAE  .* sig2tvec_FEMA_FAE).^2,  2);
SqErr_FEMA_RFX_FASE = sum((sig2mat_true_FASE - sig2mat_FEMA_FASE .* sig2tvec_FEMA_FASE).^2, 2);

% LME errors (remember to square each RFX)
SqErr_LME_RFX_FE   = sum((sig2mat_true_FE  - sig2mat_fitlme_FE.^2).^2,  2);
SqErr_LME_RFX_SE   = sum((sig2mat_true_SE  - sig2mat_fitlme_SE.^2).^2,  2);
SqErr_LME_RFX_FSE  = sum((sig2mat_true_FSE - sig2mat_fitlme_FSE.^2).^2, 2);

% Difference between estimates from LME and FEMA
SqDiff_FEMA_LME_RFX_FE   = sum((sig2mat_FEMA_FE   .* sig2tvec_FEMA_FE  - sig2mat_fitlme_FE.^2).^2,  2);
SqDiff_FEMA_LME_RFX_SE   = sum((sig2mat_FEMA_SE   .* sig2tvec_FEMA_SE  - sig2mat_fitlme_SE.^2).^2,  2);
SqDiff_FEMA_LME_RFX_FSE  = sum((sig2mat_FEMA_FSE  .* sig2tvec_FEMA_FSE - sig2mat_fitlme_FSE.^2).^2, 2);

%% Visualize results - fixed effects
col_FEMA = [27,158,119]./255;
col_LME  = [217,95,2]./255;

x_FEMA = (1:nXvars) - 0.2;
x_LME  = (1:nXvars) + 0.2;
xlabs  = {'X0', 'X1', 'X2', 'X3', 'X4'};

% Plotting the sum (across y variables) of squared differences 
fig  = figure('Units', 'centimeters', 'Position', [10 10 16 16], 'Name', 'FixedEffects - Squared difference from ground truth');
allH = tight_subplot(3, 3, [0.08 0.04], [0.04 0.03], [0.1 0.02]);

% FE
hold(allH(1), 'on');
box(allH(1), 'off');
scatter(allH(1), x_FEMA, sum(SqErr_FEMA_FFX_FE,2), 30, 'MarkerFaceColor', col_FEMA, 'MarkerEdgeColor', col_FEMA, 'MarkerFaceAlpha', 1, 'Marker', 'o');
scatter(allH(1), x_LME,  sum(SqErr_LME_FFX_FE, 2), 30, 'MarkerFaceColor', col_LME,  'MarkerEdgeColor', col_LME,  'MarkerFaceAlpha', 1, 'Marker', 'o');
xticks(allH(1),  1:nXvars);
xticklabels(allH(1), xlabs);
doAxisStuff_SqErr(allH(1), true, nyVars);
title(allH(1), 'FE', 'FontName', 'Courier', 'FontSize', 10);

% SE
hold(allH(2), 'on');
box(allH(2), 'off');
scatter(allH(2), x_FEMA, sum(SqErr_FEMA_FFX_SE,2), 30, 'MarkerFaceColor', col_FEMA, 'MarkerEdgeColor', col_FEMA, 'MarkerFaceAlpha', 1, 'Marker', 'o');
scatter(allH(2), x_LME,  sum(SqErr_LME_FFX_SE, 2), 30, 'MarkerFaceColor', col_LME,  'MarkerEdgeColor', col_LME,  'MarkerFaceAlpha', 1, 'Marker', 'o');
xticks(allH(2),  1:nXvars);
xticklabels(allH(2), xlabs);
doAxisStuff_SqErr(allH(2), false, nyVars);
title(allH(2), 'SE', 'FontName', 'Courier', 'FontSize', 10);

% FSE
hold(allH(3), 'on');
box(allH(3), 'off');
scatter(allH(3), x_FEMA, sum(SqErr_FEMA_FFX_FSE,2), 30, 'MarkerFaceColor', col_FEMA, 'MarkerEdgeColor', col_FEMA, 'MarkerFaceAlpha', 1, 'Marker', 'o');
scatter(allH(3), x_LME,  sum(SqErr_LME_FFX_FSE,2),  30, 'MarkerFaceColor', col_LME,  'MarkerEdgeColor', col_LME,  'MarkerFaceAlpha', 1, 'Marker', 'o');
xticks(allH(3),  1:nXvars);
xticklabels(allH(3), xlabs);
doAxisStuff_SqErr(allH(3), false, nyVars);
title(allH(3), 'FSE', 'FontName', 'Courier', 'FontSize', 10);

% AE
hold(allH(4), 'on');
box(allH(4), 'off');
scatter(allH(4), 1:nXvars, sum(SqErr_FEMA_FFX_AE,2), 30, 'MarkerFaceColor', col_FEMA, 'MarkerEdgeColor', col_FEMA, 'MarkerFaceAlpha', 1, 'Marker', 'o');
xlim(allH(4), [0 6]);
xticks(allH(4),  1:nXvars);
xticklabels(allH(4), xlabs);
doAxisStuff_SqErr(allH(4), true, nyVars);
title(allH(4), 'AE', 'FontName', 'Courier', 'FontSize', 10);

% SAE
hold(allH(5), 'on');
box(allH(5), 'off');
scatter(allH(5), 1:nXvars, sum(SqErr_FEMA_FFX_SAE,2), 30, 'MarkerFaceColor', col_FEMA, 'MarkerEdgeColor', col_FEMA, 'MarkerFaceAlpha', 1, 'Marker', 'o');
xlim(allH(5), [0 6]);
xticks(allH(5),  1:nXvars);
xticklabels(allH(5), xlabs);
doAxisStuff_SqErr(allH(5), false, nyVars);
title(allH(5), 'SAE', 'FontName', 'Courier', 'FontSize', 10);

% FAE
hold(allH(6), 'on');
box(allH(6), 'off');
scatter(allH(6), 1:nXvars, sum(SqErr_FEMA_FFX_FAE,2), 30, 'MarkerFaceColor', col_FEMA, 'MarkerEdgeColor', col_FEMA, 'MarkerFaceAlpha', 1, 'Marker', 'o');
xlim(allH(6), [0 6]);
xticks(allH(6),  1:nXvars);
xticklabels(allH(6), xlabs);
doAxisStuff_SqErr(allH(6), false, nyVars);
title(allH(6), 'FAE', 'FontName', 'Courier', 'FontSize', 10);

% FASE
hold(allH(7), 'on');
box(allH(7), 'off');
scatter(allH(7), x_FEMA, sum(SqErr_FEMA_FFX_FASE,2), 30, 'MarkerFaceColor', col_FEMA, 'MarkerEdgeColor', col_FEMA, 'MarkerFaceAlpha', 1, 'Marker', 'o');
xlim(allH(7), [0 6]);
xticks(allH(7),  1:nXvars);
xticklabels(allH(7), xlabs);
doAxisStuff_SqErr(allH(7), true, nyVars);
title(allH(7), 'FASE', 'FontName', 'Courier', 'FontSize', 10);

% Delete axis 8 - used for legend
delete(allH(8));

% Make a common legend
ll = legend(allH(1), {'FEMA', 'fitlmematrix'}, 'Orientation', 'vertical', 'Box', 'off', 'FontSize', 10);
ll.Position(2) = 0.22;
ll.Position(1) = 0.4;

% Plot the time taken by different models + permutations
allTimes_FEMA = [FEMAelapsed_FE, FEMAelapsed_SE, FEMAelapsed_FSE, FEMAelapsed_AE, FEMAelapsed_SAE, FEMAelapsed_FAE, FEMAelapsed_FASE];
hold(allH(9), 'on');
box(allH(9), 'off');
plot(allH(9), 1:7, allTimes_FEMA, 'LineStyle', '-', 'LineWidth', 1, 'Color', [0 0 0], 'Marker', '.', 'MarkerSize', 10);
xlim(allH(9), [1 7]);
xticks(allH(9), 1:7);
xticklabels(allH(9), {'FE', 'SE', 'FSE', 'AE', 'SAE', 'FAE', 'FASE'});
allH(9).YAxis.TickLabelsMode = 'auto';
allH(9).XAxis.FontSize = 7;
allH(9).YAxis.FontSize = 8;
allH(9).XAxis.FontName = 'Courier';
allH(9).YAxis.FontName = 'Courier';
allH(9).YLabel.String  = 'Time taken (s)';  
title(allH(9), 'Fit + 100 permutations', 'FontName', 'Courier', 'FontSize', 10);

% Optionally save for future references
% outDir  = '/ess/p697/cluster/users/parekh/2023-02-02_FEMA-Experiments/2023-11-08_Redone';
% outName = fullfile(outDir, [char(datetime('now', 'Format', 'yyyy-MM-dd-hhmmss')), '-Results_FEMA_tests.png']);
% print(outName, '-dpng', '-r600');
% close(fig);

%% Supporting functions
function [beta_hat_fitlme, beta_se_fitlme, sig2mat_fitlme, sig2Low_fitlme, sig2Upp_fitlme] = compile_fitlmematrix(mdls, nX, nY, nRFX)
% Initialize
[beta_hat_fitlme, beta_se_fitlme]                   = deal(zeros(nX,   nY));
[sig2mat_fitlme,  sig2Low_fitlme, sig2Upp_fitlme]   = deal(zeros(nRFX, nY));

% Go over every cell and compile information
% Remember that these are in standard deviation units - to be squared
% before comparing with FEMA output
for mdl = 1:nY
    [~, ~, stats]               = covarianceParameters(mdls{mdl,1});
    beta_hat_fitlme(1:nX,  mdl) = mdls{mdl,1}.Coefficients.Estimate;
    beta_se_fitlme(1:nX,   mdl) = mdls{mdl,1}.Coefficients.SE;

    for rfx = 1:nRFX
        sig2mat_fitlme(rfx, mdl) = stats{rfx}.Estimate;
        sig2Low_fitlme(rfx, mdl) = stats{rfx}.Lower;
        sig2Upp_fitlme(rfx, mdl) = stats{rfx}.Upper;
    end
end
end

function [sig2mat_lowCI, sig2mat_uppCI] = compile_FEMA(sig2mat_perm, sig2tvec_perm, nRFX, nY)
[sig2mat_lowCI, sig2mat_uppCI] = deal(zeros(nRFX, nY));
for params  = 1:nRFX
    sig2mat_lowCI(params,:) = prctile(squeeze(sig2mat_perm(params,:,2:end)) .* squeeze(sig2tvec_perm(:,:,2:end)),  2.5, 2);
    sig2mat_uppCI(params,:) = prctile(squeeze(sig2mat_perm(params,:,2:end)) .* squeeze(sig2tvec_perm(:,:,2:end)), 97.5, 2);
end
end

function doAxisStuff_SqErr(h, doLabel, nY)
h.YAxis.TickLabelsMode = 'auto';
h.XAxis.FontSize = 8;
h.YAxis.FontSize = 8;
if doLabel
    ylabel(h, ['$$\bf\sum_{y=1}^{', num2str(nY), '} (\beta_{True} - \beta_{Estimate})^2$$'], 'FontSize', 10, 'FontName', 'Courier', 'Interpreter', 'latex');
    % ylabel(h, '\Sigma_{y=1} (\beta_{True} - \beta_{Estimate})^2', 'FontSize', 12, 'FontName', 'Courier');
end
h.XAxis.FontName = 'Courier';
h.YAxis.FontName = 'Courier';
end