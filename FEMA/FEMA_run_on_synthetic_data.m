% This 'hello world' demo example generates an artificial data with
% n=10000 observations, up to 3 measures per subject, up to 2  subjects
% per family. 
%
% To run this example clone the following repos and add them to MATLAB path
% https://github.com/andersonwinkler/PALM 
% https://github.com/cmig-research-group/cmig_utils

% setup parameters for FEMA_synthesize
n = 10000;
p = 5;
v = 10;

X = randn(n, p);
iid = sort(randsample([1:n, 1:n, 1:n], n));  % up to 3 observations per individual
fid = ceil(iid / 2);                         % up to 2 individuals per family
eid = ones(n, 1);                            % non needed
agevec = zeros(n, 1);                        % not needed
ymat = nan(n, v);   % only size of ymat is used by FEMA_synthesize
[ymat, sig2tvec_true, sig2mat_true] = FEMA_synthesize(X,iid,fid,agevec,ymat);

% output of FEMA_FEMA_synthesize :
% sig2tvec_true      - total variance of random effects, one value per voxel
% sig2mat_true(1, :) - normalized variance for families
% sig2mat_true(2, :) - normalized variance for repeated observations
%                      (normalized means it's between 0 and 1, in units of sig2tvec_true)

% setup parameters for FEMA_fit
contrasts = ones(1, p);
pihatmat  = eye(n);
RandomEffects = {'F', 'S', 'E'}; % 'A'
niter=1;
nbins=20;

[beta_hat, beta_se, zmat, logpmat, sig2tvec, sig2mat, binvec, logLikvec, ...
 beta_hat_perm, beta_se_perm, zmat_perm, sig2tvec_perm, sig2mat_perm, logLikvec_perm, perms] = FEMA_fit( ...
    X,iid,eid,fid,agevec,ymat,niter,contrasts,nbins,pihatmat, 'RandomEffects', RandomEffects);

% draw plots to illustrate that variance components are correctly recovered
figure(1); hold on;
plot(sig2tvec_true, sig2tvec, '*');
plot(sig2tvec_true, sig2tvec_true, '-')
xlabel('sig2tvec_true'); ylabel('sig2tvec');

figure(2); hold on;
plot(sig2mat_true(1, :), sig2mat(1, :), '*');
plot(sig2mat_true(1, :), sig2mat_true(1, :), '-')
xlabel('sig2mat_true, F'); ylabel('sig2mat, F');

figure(3); hold on;
plot(sig2mat_true(2, :), sig2mat(2, :), '*');
plot(sig2mat_true(2, :), sig2mat_true(2, :), '-')
xlabel('sig2mat_true, S'); ylabel('sig2mat, S');
