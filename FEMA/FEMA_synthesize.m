function [ymat sig2tvec_true sig2mat_true] = FEMA_synthesize(X,iid,eid,fid,agevec,ymat,pihatmat,varargin)

p = inputParser;
addParamValue(p,'sig2tval',NaN);
addParamValue(p,'sig2repval',NaN);
addParamValue(p,'sig2famval',NaN);
addParamValue(p,'RandomEffects',{'F' 'S' 'E'}); % Default to Family, Subject, and eps

parse(p,varargin{:})
sig2tval = p.Results.sig2tval;
sig2repval = p.Results.sig2repval;
sig2famval = p.Results.sig2famval;
RandomEffects = p.Results.RandomEffects;

%Parse family structure
[clusterinfo, Ss, iid, famtypevec, famtypelist, subj_famtypevec]=FEMA_parse_family(iid,eid,fid,agevec,pihatmat,'RandomEffects',RandomEffects);

nfam = length(clusterinfo); nvox = size(ymat,2);

single_or_double = 'double';
%single_or_double = 'single';

epsmin = 0.05; % This should be an input parameter

% Randomize random effects
sig2mat_true = NaN(length(RandomEffects)-1,size(ymat,2)); ivec = 1:size(ymat,2);
while length(ivec)>0
  sig2mat_true(:,ivec) = rand([size(sig2mat_true,1) length(ivec)]);
  ivec = find(sum(sig2mat_true,1)>1-epsmin);
  logging('length(ivec)=%d',length(ivec));
  if length(ivec)==0, break; end
end
sig2mat_true = cat(1,sig2mat_true,1-sum(sig2mat_true,1));

sig2tvec_true = logspace(-log10(10),log10(10),nvox); % Allow this to be passed in as an argument (also allow to vary spatially)
if isfinite(sig2tval)
  sig2tvec_true(:) = sig2tval;
end

for fi = 1:nfam
  if mod(fi,100)==0
    logging('fi=%d/%d',fi,nfam);
  end
  tmp = 0;
  for ri = 1:length(RandomEffects)
    tmp = tmp + sqrt(sig2mat_true(ri,:)).*mvnrnd(zeros(length(clusterinfo{fi}.jvec_fam),1),double(getfield(clusterinfo{fi},sprintf('V_%s',RandomEffects{ri}))),nvox)';
  end
  ymat(clusterinfo{fi}.jvec_fam,:) = sqrt(sig2tvec_true).*tmp;
end

% figure; plot(mean(ymat.^2,1),sig2tvec_true,'LineWidth',2)

return


% ToDo
%   Allow for simulation of non-nulli: add X*beta_true outside of function
%   Allow for simulation of non-Gaussian distribution: iapply transform to ymat outside of function


% Old slow version

for j = 1:nvox 
  if mod(j,100)==0
    logging('j=%d/%d',j,size(ymat,2));
  end
  for fi = 1:nfam
    tmp = 0;
    for ri = 1:length(RandomEffects)
      tmp = tmp + sqrt(sig2tvec_true(j)*sig2mat_true(ri,j))*mvnrnd(zeros(length(clusterinfo{fi}.jvec_fam),1),double(getfield(clusterinfo{fi},sprintf('V_%s',RandomEffects{ri}))));
    end
    ymat(clusterinfo{fi}.jvec_fam,j) = tmp;
  end
end

