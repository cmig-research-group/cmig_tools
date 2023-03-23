function [ymat sig2tvec_true sig2mat_true] = FEMA_synthesize(X,iid,eid,fid,agevec,ymat,pihatmat,varargin)

p = inputParser;
addParamValue(p,'sig2tval',1);
addParamValue(p,'RandomEffects',{'F' 'S' 'E'}); % Default to Family, Subject, and eps
addParamValue(p,'nbins',20);

parse(p,varargin{:})
sig2tval = p.Results.sig2tval;
RandomEffects = p.Results.RandomEffects;
nbins = str2num_amd(p.Results.nbins);

%Parse family structure
[clusterinfo, Ss, iid, famtypevec, famtypelist, subj_famtypevec]=FEMA_parse_family(iid,eid,fid,agevec,pihatmat,'RandomEffects',RandomEffects);

nfam = length(clusterinfo); nvox = size(ymat,2);

single_or_double = 'double';
%single_or_double = 'single';



%%%%%%%%%%%%%% Should separate the followng into a callable function

binvals_edges = linspace(0,1,nbins+1); binvals_edges(end) = inf;
if length(RandomEffects)==2 % Why is this needed? -- N-d version, below, should work for 2-d?
  sig2gridi = colvec(1:length(binvals_edges)-1);
  sig2gridl = colvec(binvals_edges(1:end-1));
  sig2gridu = colvec(binvals_edges(2:end));
else
  sig2gridi = ndgrid_amd(repmat({1:length(binvals_edges)-1},[1 length(RandomEffects)-1]));
  sig2gridl = ndgrid_amd(repmat({binvals_edges(1:end-1)},[1 length(RandomEffects)-1]));
  sig2gridu = ndgrid_amd(repmat({binvals_edges(2:end)},[1 length(RandomEffects)-1]));
end
%sig2grid = (sig2gridl+sig2gridu)/2; % Should make sig2gridl+1/nbins
sig2grid = sig2gridl+(0.5/nbins);
sig2grid_ivec = find(sum(sig2grid,2)<=1-0.5/nbins); % Get rid of "impossible" bins -- use middle of bin instead
sig2gridl = sig2gridl(sig2grid_ivec,:);
sig2gridu = sig2gridu(sig2grid_ivec,:);
sig2grid = sig2grid(sig2grid_ivec,:);
sig2gridi = sig2gridi(sig2grid_ivec,:);
nsig2bins = size(sig2grid,1); % Should handle case of no binning

%%%%%%%%%%%%%%


binvec_true = 1+floor(nsig2bins*[0:size(ymat,2)-1]/size(ymat,2));
sig2mat_true = NaN(length(RandomEffects),size(ymat,2));
sig2mat_true(1:end-1,:) = sig2grid(binvec_true,:)';
sig2mat_true(end,:) = 1-sum(sig2mat_true(1:end-1,:));

sig2tvec_true(:) = sig2tval;

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

return



% Old version -- should delete



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

