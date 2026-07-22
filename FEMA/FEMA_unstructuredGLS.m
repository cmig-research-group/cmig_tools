function [beta_hat, beta_se, coeffCovar, converge, betacon_hat, betacon_se] = ...
         FEMA_unstructuredGLS(X, currY, curr_sig2tvec, curr_sig2mat,          ...
                              clusterinfo, visitnum, allSz, allR, allC,       ...
                              nnz_max, num_obs, num_RFX, RFX_ord, locJVec,    ...
                              contrasts, SingleOrDouble, useLSQ)
% Initialize
allV  = zeros(nnz_max, 1);
count = 1;
nfam  = length(clusterinfo);

% Get warning statuses for singular and nearly singular cases;
% temporarily set their display off
statusSingular = warning('off', 'MATLAB:singularMatrix');
statusNearSing = warning('off', 'MATLAB:nearlySingularMatrix');

% Clear last warning
lastwarn('');

% Go over families and compute W
for fi = 1:nfam

    % Extract current cluster
    currClus   = struct2cell(clusterinfo{fi});

    tmpSize = allSz(fi);
    wchLocs = currClus{locJVec};

    % Compute V
    Vs_fam = zeros(tmpSize);
    b      = visitnum(wchLocs);
    for ri = 1:num_RFX-1
        Vs_fam = Vs_fam + curr_sig2mat(b, b, ri) .* currClus{RFX_ord(ri)};
    end

    % Compute inverse of V
    Vis_fam = cast(double(Vs_fam) \ eye(tmpSize), SingleOrDouble);
    msg     = lastwarn;
    if ~isempty(msg)
        Vis_fam = cast(pinv(double(Vs_fam)), SingleOrDouble);
        msg = ''; %#ok<NASGU>
        lastwarn('');
    end

    tmp                     = numel(Vis_fam);
    allV(count:count+tmp-1) = Vis_fam(:);
    count                   = count + tmp;
end
allWsFam = sparse(allR, allC, allV, num_obs, num_obs, nnz_max);

% Compute XtW
XtW = X' * allWsFam;

% Compute XtWX
B  = XtW * X;

% Calculate inverse of XtWX
if rank(B) < size(B,2)
    if useLSQ
        Bi = lsqminnorm(B, eye(size(B)));
    else
        Bi = pinv(B);
    end
else
    Bi = B \ eye(size(B));
end

% Calculate beta coefficient
[Cov_beta, converge]  = nearestSPD_timeout(Bi);
beta_hat              = Cov_beta * (XtW * currY);
beta_se               = sqrt(diag(Cov_beta) * curr_sig2tvec);
coeffCovar            = Cov_beta .* curr_sig2tvec;

% Evaluate contrasts
[betacon_hat, betacon_se] = deal(zeros(size(contrasts,1), 1, SingleOrDouble));
for ci = 1:size(contrasts,1)
    betacon_hat(ci, 1) = contrasts(ci,:) * beta_hat;
    betacon_se(ci,  1) = sqrt(contrasts(ci,:) * Cov_beta * contrasts(ci,:)' * curr_sig2tvec);
end

% Reset the status of warnings
warning(statusSingular);
warning(statusNearSing);
