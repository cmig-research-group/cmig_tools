function [beta_hat, beta_se, tStats, pValues] = FEMA_sig2binseg_parfeval_GWAS(genovec,    ymat,             clusterinfo,     ...
                                                                              binvec,     sig2tvec,         GroupByFamType,  ...
                                                                              famtypevec, OLSflag,          Vs_fam,          ...
                                                                              Vs_famtype, Ws_fam,           Ws_famtype,      ...
                                                                              df,         SingleOrDouble,   pValType,        ...
                                                                              saveName,   Chr,              SNPID)
% Function to estimate the effect of each (residualized) SNP on 
% (residualized) phenotype using OLS or GLS
%% Inputs:
% genovec:          [n x m]     matrix of n subjects and m SNPs
%
% ymat:             [n x v]     matrix of n subjects and v phenotypes
% 
% clusterinfo:      [1 x f]     cell type of f families/clusters, returned
%                               as an output from FEMA_parse_family
% 
% binvec:           [1 x b]     vector of bin values on which random
%                               effects are evaluated
% 
% sig2tvec:         [1 x v]     total residual variance of v phenotypes
% 
% GroupByFamType:   logical     true/false indicating if GroupByFamType
% 
% famtypevec:       [1 x f]     vector indicating family type
% 
% OLSflag:          logical     OLS or GLS solution
% 
% Vs_fam, 
% Vs_famtype,
% Ws_fam,
% Ws_famtype:                   various necessary variables returned as an
%                               output from FEMA_compileTerms
% 
% df:               [1 x 1]     degrees of freedom
% 
% SingleOrDouble:   character   single or double precision for genovec
% 
% pValType:         character   distribution for calculating p values:
%                                   * 't'
%                                   * 'z'
%                                   * 'chi'
% 
% saveName:         character   full path to where results will be saved
% 
% Chr:              [m x 1]     cell type of chromosome number for m SNPs
% 
% SNPID:            [m x 1]     cell type of SNP IDs for m SNPs
%
%% Outputs:
% beta_hat:         [m x v]     estimated beta coefficients for m SNPs and
%                               v phenotypes
%
% beta_se:          [m x v]     estimated standard error of the beta 
%                               coefficients for m SNPs and v phenotypes
%
% tStats:           [m x v]     ratio of beta_hat and beta_se
%
% pValues:          [m x v]     p values for m SNPs and v phenotypes

%% Initialize
beta_hat     = zeros(size(genovec,2), size(ymat,2), class(ymat)); 
beta_se      = zeros(size(beta_hat),  class(ymat));
tStats       = zeros(size(beta_hat),  class(ymat));
pValues      = zeros(size(beta_hat),  class(ymat));
binLoc       = 1;

% Cast genovec as double precision, if required
if strcmpi(SingleOrDouble, 'double')
    genovec = double(genovec);
end

for sig2bini = unique(binvec(isfinite(binvec)), 'stable')
    ivec_bin = find(binvec==sig2bini);
    XtW      = zeros(fliplr(size(genovec)), class(genovec)); 
    XtWVsWt  = zeros(fliplr(size(genovec)), class(genovec));

    if GroupByFamType
        for fi = 1:length(clusterinfo)
            XtW(:, clusterinfo{fi}.jvec_fam)    = genovec(clusterinfo{fi}.jvec_fam,:)' * Ws_famtype{binLoc, famtypevec(fi)};
            XtWVsWt(:,clusterinfo{fi}.jvec_fam) = XtW(:,clusterinfo{fi}.jvec_fam)      * Vs_famtype{binLoc, famtypevec(fi)} * Ws_famtype{binLoc, famtypevec(fi)}';
        end
    else
        for fi = 1:length(clusterinfo)
            XtW(:,clusterinfo{fi}.jvec_fam)     = genovec(clusterinfo{fi}.jvec_fam,:)' * Ws_fam{binLoc, fi};
            XtWVsWt(:,clusterinfo{fi}.jvec_fam) = XtW(:,clusterinfo{fi}.jvec_fam)      * Vs_fam{binLoc, fi} * Ws_fam{binLoc, fi}';
        end
    end
    
    % Calculating beta coefficient for every SNP
    % beta = inv(X' * W * X) * X' * W * y
    % where W is the inverse of the total variance V
    %
    % Calculating the inner term
    % X' * W has been already calculated ---> XtWVsWt
    % Therefore, the inner term is the inverse of XtWVsWt * X
    % However, since each SNP is being treated independently, we are only
    % interested in the diagonal term of the XtWVsWt * X multiplication
    % This is equivalent to a sum over element-wise multiplication:
    % sum(XtWVsWt' .* X)
    % 
    % XtWVsWt:              [chunk x subj] beta = inv(X' * W * X) * X' * W * yy
    % X or genovec:         [subj x chunk]
    % XtWVsWt' .* genovec:  [subj x chunk] 
    % Sum over the first dimension (subjects) will be [1 x chunk]
    %
    % If everything is double precision:
    % diag(XtWVsWt * genovec) approximately equal to sum(XtWVsWt' .* genovec)
    %
    % Since we are interested in the inverse of this inner term, and each
    % SNP is still being treated independently, therefore the inverse is
    % simply an inverse of a scaler, which is 1/value
    % For older MATLAB use: 1./(sum(bsxfun(@times, XtWVsWt', genovec)));
    Bi = 1./(sum(XtWVsWt' .* genovec));
    
    % XtWVsWtX = XtWVsWt * genovec; 
    % B        = XtWVsWtX; 
    % Bi       = pinv(B); % Note that Cov_beta = Bi = pinv(XtWVsWtX)*XtWVsWt*pinv(XtWVsWtX)

    % For the next set of multiplications:
    % For GLS:
    % beta      = Bi * X'  * W * y
    % or, beta  = Bi * XtW * y
    %
    % Bi:   [1     x chunk]
    % XtW:  [chunk x subj]
    % ymat: [subj  x bins]
    %
    % It is faster to implement XtW * ymat, resulting in [chunk x bins] 
    % which is then multiplied element-wise by Bi
    % 
    % A similar calculation is performed for OLS, excluding the XtW term;
    % in this case, Bi = 1./(sum(genovec .* genovec));
    
    if OLSflag
        beta_hat(:,ivec_bin) = (1./sum(genovec .* genovec))' .* (genovec' * ymat(:,ivec_bin));
        % beta_hat_tmp = pinv(genovec) * ymat(:,ivec_bin);
    else
        beta_hat(:,ivec_bin) = Bi' .* (XtW * ymat(:, ivec_bin));
        % beta_hat_tmp = bsxfun(@times, Bi', XtW) * ymat(:, ivec_bin);
        % beta_hat_tmp = Bi * XtW * ymat(:,ivec_bin); % Should generalize this to work with arbitrary W matrices
    end
    
    % Standard error of beta
    beta_se(:,ivec_bin) = sqrt(Bi' * sig2tvec(ivec_bin));
    
    % Calculate t statistics
    tStats(:,ivec_bin)  = beta_hat(:,ivec_bin) ./ beta_se(:,ivec_bin);
    
    % Calculate p value
    switch pValType
        case 't'
            pValues(:,ivec_bin) = 2 * tcdf(-abs(tStats(:,ivec_bin)), df);
    
        case 'z'
            pValues(:,ivec_bin) = 2 * normcdf(-abs(tStats(:,ivec_bin)), 0, 1);
            
        case 'chi'
            pValues(:,ivec_bin) = chi2cdf(tStats(:,ivec_bin).^2, 1, 'upper');
    end
    
    % Update binLoc
    binLoc = binLoc + 1;
end

% Save variables
% Get an estimate of size of the variables
tmpInfo = whos;
if sum([tmpInfo(ismember({tmpInfo(:).name}', {'beta_hat', 'beta_se', 'tStats', 'pValues', 'df', 'Chr', 'SNPID'})).bytes]) > 2^31
    save([saveName, '.mat'], 'beta_hat', 'beta_se', 'tStats', 'pValues', 'df', 'Chr', 'SNPID', '-v7.3');
else
    save([saveName, '.mat'], 'beta_hat', 'beta_se', 'tStats', 'pValues', 'df', 'Chr', 'SNPID');
end
end