function [sig2tvec,      sig2mat] =                                  ...
          FEMA_fit_simplified(X, iid, eid, fid, ymat_res, sig2tvec,  ...
                   pihatmat, W_1, varargin)

p = inputParser;
addParamValue(p,'ciflag', false);
addParamValue(p,'MLflag', false);
addParamValue(p,'FamilyStruct', {});
addParamValue(p,'NonnegFlag', true);

parse(p,varargin{:})

MLflag       = p.Results.MLflag;
FamilyStruct = p.Results.FamilyStruct;
NonnegFlag   = p.Results.NonnegFlag;
ciflag       = p.Results.ciflag;
clusterinfo = FamilyStruct.clusterinfo;
Ss          = FamilyStruct.Ss;
subvec1     = FamilyStruct.subvec1;
subvec2     = FamilyStruct.subvec2;

% LHS      = ymat_res(subvec1,:) .* ymat_res(subvec2,:) ./ mean(ymat_res.^2,1); % use normalized residuals
LHS      = ymat_res(subvec1,:) .* ymat_res(subvec2,:);
sig2mat  = NaN(size(Ss,2), size(ymat_res, 2));

% heterogeneity variance is considered, showed as E ~ N(0, diag(W)^{-1}).
M           = FamilyStruct.M;
subvec_e    = find(M(:,end)); % diagonal position of variance matrix for y 
LHS(subvec_e,:)  = LHS(subvec_e,:) - M(subvec_e,end) .*  W_1;
 
% no loop
for coli=1:size(ymat_res, 2)

    % % Use new version of lsqnonneg_amd to enfoce non-negative variances
    sig2mat_tmp     = lsqnonneg_amd3(M(:,1:end-1),LHS(:,coli)); % This doesn't actually ensure non-negative values! -- problem with complex ymat / LHS
    sig2mat(:,coli) = [sig2mat_tmp;1];


end

end