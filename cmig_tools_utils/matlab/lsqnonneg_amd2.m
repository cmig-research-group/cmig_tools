function [betahat_out resnorm residual] = lsqnonneg_amd(C,d,nnvec)

if ~exist('nnvec','var'), nnvec = true(1,size(C,2)); end

betahat_out = pinv(C)*d;

ivec_neg = (sum(betahat_out<0,1)>0);

d_bak = d;
d = d(:,ivec_neg);

nnlist = find(nnvec);
ncombos = 2^length(nnlist); nvox = size(d,2);
betamat = zeros(size(C,2),nvox,ncombos);
costmat = NaN(nvox,ncombos);  % costmat(:,1) = sum(d.^2,1);
for comboi = 1:ncombos % Should only re-do voxels that provide negative values for comboi = ncombos
%  fprintf(1,'comboi=%d of %d (now=%s)\n',comboi,ncombos,datestr(now));
  bitvec = colvec(~nnvec);
  bitvec(nnlist) = ~logical(bitget(comboi-1,[1:length(nnlist)]'));
  betahat = zeros(size(C,2),nvox);
  betahat(bitvec,:) = pinv(C(:,bitvec))*d;
  betahat(nnlist,:) = max(0,betahat(nnlist,:));
  costvec = sum((d-C*betahat).^2,1);
  betamat(:,:,comboi) = betahat;
  costmat(:,comboi) = costvec;
end

[mv mi] = min(costmat,[],2);
for j = 1:size(costmat,2)
  ivec = find(mi==j);
  betahat(:,ivec) = betamat(:,ivec,j);
end

betahat_out(:,ivec_neg) = betahat;

d_pred = C*betahat_out;
residual = d_pred-d_bak;
resnorm = sum(colvec(residual).^2)./sum(colvec(d).^2);
%resnorm = norm(residual); % Check why this takes so long

return

% ToDo
%   Make faster by re-compute only voxels that have any negative values for nnvec / nnlist params? -- doesn't seem to work
%     Try dropping columns with negative values, then recursively re-fitting
%   Perform regularization?

nnlist = find(nnvec);
ncombos = 2^length(nnlist); nvox = size(d,2);
betamat = zeros(size(C,2),nvox,ncombos);
costmat = NaN(nvox,ncombos);  % costmat(:,1) = sum(d.^2,1);
ivec_vox = true(1,nvox);
for comboi = 1:ncombos % Should only re-do voxels that provide negative values for comboi = ncombos
%  fprintf(1,'comboi=%d of %d (now=%s)\n',comboi,ncombos,datestr(now));
  tic
  bitvec = colvec(~nnvec);
  bitvec(nnlist) = ~logical(bitget(comboi-1,[1:length(nnlist)]'));
  if comboi==1
    betahat = zeros(size(C,2),nvox);
    costvec = NaN(1,nvox);
  else
    betahat = betamat(:,:,1);
  end
  betahat(bitvec,ivec_vox) = pinv(C(:,bitvec))*d(:,ivec_vox);
  if comboi==1
    ivec_neg = (sum(betahat<0,1)>0);
  end
  betahat(nnlist,ivec_vox) = max(0,betahat(nnlist,ivec_vox));
  costvec(ivec_vox) = sum((d(:,ivec_vox)-C*betahat(:,ivec_vox)).^2,1);
  betamat(:,:,comboi) = betahat;
  costmat(:,comboi) = costvec;
  ivec_vox = ivec_neg;
  toc
end

