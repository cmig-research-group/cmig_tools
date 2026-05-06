function FEMA_DEAP_init(fstem_imaging, fname_design, dirname_out, dirname_imaging, datatype, dirname_cache, nfrac, varargin)

% Add optional required arguments
% addRequired(p, 'fstem_imaging');
% addRequired(p, 'fname_design');
% addRequired(p, 'dirname_out');
% addRequired(p, 'dirname_imaging');
% addRequired(p, 'datatype');
% addRequired(p, 'dirname_cache');
% addRequired(p, 'nfrac');
% Parse inputs
p                 = inputParser;
p.KeepUnmatched   = true;
addOptional(p, 'onlyReturn', true);
parse(p, varargin{:});

if ~exist('fstem_imaging', 'var') || isempty(fstem_imaging)
    error('fstem_imaging needs to be specified');
end

if ~exist('fname_design', 'var') || isempty(fname_design)
    error('fname_design needs to be specified');
else
    if ~exist(fname_design, 'file')
        % Will this handle multiple design matrices?
        error(['Unable to find file: ', fname_design]);
    end
end

if ~exist('dirname_out', 'var') || isempty(dirname_out)
    error('dirname_out needs to be specified');
else
    if ~exist(dirname_out, 'dir')
        mkdir(dirname_out);
    end
end

if ~exist('dirname_imaging', 'var') || isempty(dirname_imaging)
    error('dirname_imaging needs to be specified');
else
    if ~exist(dirname_imaging, 'dir')
        error(['Unable to find dirname_imaging: ', dirname_imaging]);
    end
end

if ~exist('datatype', 'var') || isempty(datatype)
    error('datatype needs to be specified');
end

if ~exist('dirname_cache', 'var') || isempty(dirname_cache)
    error('dirname_cache needs to be specified');
else
    if ~exist(dirname_cache, 'dir')
        fprintf('Creating cache directory %s\n', dirname_cache)
        mkdir(dirname_cache);
    end
end
fname_cache = sprintf('%s/%s_%s.mat',dirname_cache,datatype,fstem_imaging);

if ~exist('nfrac', 'var') || isempty(nfrac)
    error('nfrac needs to be specified');
else
    % Make sure nfrac is numeric
    if ~isnumeric(nfrac)
        nfrac = str2double(nfrac);
    end
end

% Call FEMA_wrapper
% [~, ~, ~, ~, ~, ~, ~, ~, ~, ~, ~, ~, ~, ~, ~, ~, ~, ~, ~, ~, ~, ~, nvec_bins, tvec_bins, FamilyStruct, ~, inputs_FEMA] = ...
 [beta_hat, beta_se, zmat, logpmat, sig2tvec, sig2mat, Hessmat, logLikvec,              ...
  beta_hat_perm, beta_se_perm, zmat_perm, sig2tvec_perm, sig2mat_perm, logLikvec_perm,  ...
  binvec_save, nvec_bins, tvec_bins, FamilyStruct, coeffCovar, unstructParams,          ...
  residuals_GLS, info_fit, inputs_FEMA, mask, tfce_perm] = ...
  FEMA_wrapper('fstem_imaging', fstem_imaging, 'fname_design', fname_design, 'dirname_out', dirname_out, ...
              'dirname_imaging', dirname_imaging, 'datatype', datatype, 'outputType', {'none', 'cache'});
% FEMA_wrapper(fstem_imaging, fname_design, dirname_out, dirname_imaging, datatype);

% Why do we need to call twice? Only to test?
% [fpaths_out, beta_hat, beta_se, zmat, logpmat, sig2tvec, sig2mat, ...
%  beta_hat_perm, beta_se_perm, zmat_perm, sig2tvec_perm, sig2mat_perm, ...
%  inputs, mask, tfce_perm, colnames_interest, save_params, logLikvec, Hessmat, coeffCovar] = ...
%  FEMA_wrapper(fstem_imaging,fname_design,dirname_out,dirname_imaging,datatype, 'save_params', []);

% Compute binvec based on binvec from random effects estimates
% Split ymat into "equal-time" subsets, save out to files
% Write code to add list of analyses  to perform (modify existing)
% Write code to pick up next analysis to do (modify existing)
% Write code for persistent workers to perform specific subset of total analysis segments (modify existing)
% Write code to check for completed analyses, gather results

% keyboard

% Get bin info

% PP: disabling pre-residualization option
% X_resid = X;
% ymat_resid = ymat;
% Pre-residualize both LHS and RHS
%  M = X(:,2:end); % Should specify "columns of interest", residualize LHS & RHS for others, only save stats for those
%  X_tmp = X(:,[1 end]);
%  X_resid = X_tmp - M*(M\X_tmp); X_resid(:,end) = 1;
%  ymat_resid = ymat - M*(M\ymat);

% Why do we need to run FEMA so many times?
% % Should update FEMA_fit to allow passing of pre-computed clusterinfo 
%   [beta_hat beta_se zmat logpmat sig2tvec sig2mat Hessmat logLikvec beta_hat_perm beta_se_perm zmat_perm sig2tvec_perm sig2mat_perm logLikvec_perm binvec nvec_bins tvec_bins FamilyStruct] = ...
%      FEMA_fit(X_resid,iid,eid,fid,agevec,ymat_resid,niter,contrasts,nbins, GRM,'RandomEffects',RandomEffects,'nperms',nperms,'CovType',CovType,'FixedEstType',FixedEstType,'RandomEstType',RandomEstType,...
%      'GroupByFamType',GroupByFamType,'NonnegFlag',NonnegFlag,'SingleOrDouble',SingleOrDouble,'logLikflag',logLikflag,'Hessflag',Hessflag,'ciflag',ciflag,...
%      'permtype',permtype,'PregID',PregID,'HomeID',HomeID,'synthstruct',synthstruct);
% 
%   [beta_hat beta_se zmat logpmat sig2tvec sig2mat Hessmat logLikvec beta_hat_perm beta_se_perm zmat_perm sig2tvec_perm sig2mat_perm logLikvec_perm binvec nvec_bins tvec_bins] = ...
%      FEMA_fit(X_resid,iid,eid,fid,agevec,ymat_resid,niter,contrasts,nbins, GRM,'RandomEffects',RandomEffects,'nperms',nperms,'CovType',CovType,'FixedEstType',FixedEstType,'RandomEstType',RandomEstType,...
%      'GroupByFamType',GroupByFamType,'NonnegFlag',NonnegFlag,'SingleOrDouble',SingleOrDouble,'logLikflag',logLikflag,'Hessflag',Hessflag,'ciflag',ciflag,...
%      'permtype',permtype,'PregID',PregID,'HomeID',HomeID,'synthstruct',synthstruct,'FamilyStruct',FamilyStruct);

% Should also keep track of time for rest of FEMA_fit?

defvec = isfinite(nvec_bins+tvec_bins);
P      = polyfit(nvec_bins(defvec),tvec_bins(defvec),1);
% P = [1.3839e-05 0.0624]; % Values for mmilcluster10
% P = [1.2792e-05 0.0360]; % Values for mmilcluster10, pre-resid
tvec_bins_pred = polyval(P,nvec_bins);
% figure; plot(nvec_bins,tvec_bins,'*',nvec_bins,tvec_bins_pred,'-'); ylim([0 max(ylim)]);

ivec_bins = find(isfinite(tvec_bins_pred));

% Split into N subsets of bins predicted  to take roughly equal time (intervals of tsum)
% nfrac = 20;
fraclist = linspace(0,1.001,nfrac+1);

% Save imaging data in segments here -- also check actual computation time (prediction doesn't seem accurate)
% dirname_cache = '~/FEMA_cache';

% Extract variables and save
% Should be a better way to do this
tic
fid             = inputs_FEMA.fid;
iid             = inputs_FEMA.iid;
eid             = inputs_FEMA.eid;
agevec          = inputs_FEMA.agevec;
GRM             = inputs_FEMA.GRM;
RandomEffects   = inputs_FEMA.RandomEffects;
PregID          = inputs_FEMA.PregID;
HomeID          = inputs_FEMA.HomeID;
ymat            = inputs_FEMA.ymat;
X               = inputs_FEMA.X;
niter           = inputs_FEMA.niter;
nbins           = inputs_FEMA.nbins;
contrasts       = inputs_FEMA.contrasts;
nperms          = inputs_FEMA.nperms;
CovType         = inputs_FEMA.CovType;
FixedEstType    = inputs_FEMA.FixedEstType;
RandomEstType   = inputs_FEMA.RandomEstType;
GroupByFamType  = inputs_FEMA.GroupByFamType;
NonnegFlag      = inputs_FEMA.NonnegFlag;
precision       = inputs_FEMA.precision;
logLikflag      = inputs_FEMA.logLikflag;
Hessflag        = inputs_FEMA.Hessflag;
ciflag          = inputs_FEMA.ciflag;
permtype        = inputs_FEMA.PermType;
synthstruct     = inputs_FEMA.synthstruct;

clear inputs_FEMA

save(fname_cache,'iid','eid','fid','agevec','GRM','RandomEffects','PregID','HomeID','FamilyStruct');
toc

%% Pass 1
nvec_frac  = NaN(nfrac,1); 
dtvec_frac = NaN(nfrac,1);

tsum = cumsum(tvec_bins_pred(ivec_bins))/sum(tvec_bins_pred(ivec_bins));
for fraci = 1:nfrac
    ivec_tmp  = find(tsum>=fraclist(fraci)&tsum<fraclist(fraci+1));
    ivec_frac = find(ismember(binvec_save,ivec_bins(ivec_tmp)));
    ymat_frac = ymat(:,ivec_frac);
    fprintf(1,'fraci=%d/%d: length(ivec_tmp)=%d length(ivec_frac)=%d sumt=%f\n',fraci,nfrac,length(ivec_tmp),length(ivec_frac),sum(tvec_bins_pred(binvec_save(ivec_frac))));

    % Run and time FEMA_fit; no output required
    t0 = now;
    FEMA_fit(X, iid, eid, fid, agevec, ymat_frac, contrasts, nbins, GRM,        ...
             'niter', niter, 'RandomEffects', RandomEffects, 'nperms', nperms,  ...
             'CovType', CovType, 'FixedEstType', FixedEstType,                  ...
             'RandomEstType', RandomEstType, 'GroupByFamType', GroupByFamType,  ...
             'NonnegFlag', NonnegFlag, 'precision', precision,                  ...
             'logLikflag', logLikflag, 'Hessflag', Hessflag, 'ciflag', ciflag,  ...
             'permtype', permtype, 'PregID', PregID, 'HomeID', HomeID,          ...
             'synthstruct', synthstruct, 'FamilyStruct', FamilyStruct);
    t1 = now;

    dt = (t1-t0)*3600*24; dt_pred = sum(tvec_bins_pred(ivec_bins(ivec_tmp)));
    fprintf(1,'fraci=%d/%d %0.1f seconds pred:  %0.1f seconds (diff:%0.1f seconds) length(ivec_frac)=%d\n',fraci,nfrac,dt,dt_pred,dt-dt_pred,length(ivec_frac));
    dtvec_frac(fraci) = dt-dt_pred; nvec_frac(fraci) = length(ivec_frac);
    M = cat(2,nvec_frac,ones(size(nvec_frac)));
    bvec = pinv(M)*dtvec_frac;
end

%% Pass 2
tvec_bins_pred = tvec_bins_pred + bvec(1)*nvec_bins;
tsum = cumsum(tvec_bins_pred(ivec_bins))/sum(tvec_bins_pred(ivec_bins));
for fraci = 1:nfrac
    ivec_tmp  = find(tsum>=fraclist(fraci)&tsum<fraclist(fraci+1));
    ivec_frac = find(ismember(binvec_save,ivec_bins(ivec_tmp)));
    ymat_frac = ymat(:,ivec_frac);
    fprintf(1,'fraci=%d/%d: length(ivec_tmp)=%d length(ivec_frac)=%d sumt=%f\n',fraci,nfrac,length(ivec_tmp),length(ivec_frac),sum(tvec_bins_pred(binvec_save(ivec_frac))));

    % Run and time FEMA_fit; no output required
    t0 = now;
    FEMA_fit(X, iid, eid, fid, agevec, ymat_frac, contrasts, nbins, GRM,        ...
             'niter', niter, 'RandomEffects', RandomEffects, 'nperms', nperms,  ...
             'CovType', CovType, 'FixedEstType', FixedEstType,                  ...
             'RandomEstType', RandomEstType, 'GroupByFamType', GroupByFamType,  ...
             'NonnegFlag', NonnegFlag, 'precision', precision,                  ...
             'logLikflag', logLikflag, 'Hessflag', Hessflag, 'ciflag', ciflag,  ...
             'permtype', permtype, 'PregID', PregID, 'HomeID', HomeID,          ...
             'synthstruct', synthstruct, 'FamilyStruct', FamilyStruct);
    t1 = now;

    dt = (t1-t0)*3600*24; dt_pred = sum(tvec_bins_pred(ivec_bins(ivec_tmp)));
    fprintf(1,'fraci=%d/%d %0.1f seconds pred:  %0.1f seconds (diff:%0.1f seconds) length(ivec_frac)=%d\n',fraci,nfrac,dt,dt_pred,dt-dt_pred,length(ivec_frac));
    fname_out = sprintf('%s/%s_%s_%02d_%02d.mat',dirname_cache,datatype,fstem_imaging,fraci,nfrac);

    % Should save random effects estimates from fit for ymat_frac / ivec_frac, for the current X
    save(fname_out,'ymat_frac','ivec_frac','-v7.3', '-nocompression'); 
end

% keyboard

% ToDos
%   Make code above iterative, predicting total time per set of bins, re-estimating predicted time as function of length(ivec_frac)
% Should print out sum(tvec_bins_pred(ivec_bins(ivec_tmp))) in loop above for cmparison with actual time

%% Save info for persistent concurrent workers: X_resid
jobname = 'job1';
taskdir = sprintf('%s/FEMA_pending/%s',dirname_cache,jobname);
jsystem(sprintf('mkdir -p %s',taskdir));
fname_task = sprintf('%s/taskinfo.mat',taskdir);
save(fname_task, 'X', 'iid', 'eid');

% Save pre-computed info needed for updated FEMA_fit:
%   Pre-resiudalized X_resid
%   iid, eid, fid, agevec, FamilyStruct
%   omit synthstruct=[], HomeID=[], PregID=[], permtype='wildbootstrap', ciflag=0, Hessflag=0, logLikflag=0, SingleOrDouble='double', NonnegFlag=1, Parallelize=0 (ignore), GroupByFamType=1,
%     RandomEstType='MoM', FixedEstType='GLS', CovType='analytic', nperms=0, RandomEffects={'F' 'S' 'E'}, nbins=20, contrasts=[], 