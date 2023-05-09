function [beta_hat, beta_se, tStats, logpValues] = FEMA_sig2binseg_parfeval_GWAS_bf(genoCell,      ymat,           clusterinfo,     ...
                                                                                    binvec,        sig2tvec,       GroupByFamType,  ...
                                                                                    famtypevec,    Ws_famtype,     Ws_fam,          ...
                                                                                    basisFunction, pValType,       df,              ...
                                                                                    OLSflag,       SingleOrDouble, outDir,  outName)
% Function to estimate the effect of each (residualized) SNP on
% (residualized) phenotype using OLS or GLS
%% Inputs:
% genoCell:         [b x 1]     cell type of residualized genotype matrix 
%                               where each cell contains genotype data of
%                               n subjects and m SNPs; each cell
%                               corresponds to a bin value (see,
%                               FEMA_residualizeGenotype)
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
% Ws_famtype, Ws_fam:           inverses of the V term returned by 
%                               FEMA_compileTerms (for every bin)
% 
%% Optional inputs:
% basisFunction:    [n x q]     matrix of basis functions to use as an
%                               interaction term for (residualized) 
%                               genotyping matrix; defaults to an
%                               intercept term; the same set of basis
%                               functions are used for all bins
%
% pValType:         character   distribution for calculating p values:
%                                   * 't'
%                                   * 'z'
%                                   * 'chi'
%
% df:               [1 x 1]     degrees of freedom
%
% OLSflag:          logical     OLS or GLS solution
%
% SingleOrDouble:   character   single or double precision for genovec
%
% outDir:           character   full path to where the results should be
%                               saved
%
% outName:          character   name of the mat file to use for saving
%
%% Outputs:
% beta_hat:         [m x q x v] estimated beta coefficients for m SNPs and
%                               q basis functions for v phenotypes; if only
%                               one basis function is specified, then
%                               beta_hat is [m x v] matrix
%
% beta_se:          [m x q x v] estimated standard error of the beta 
%                               coefficients for m SNPs and q basis 
%                               functions for v phenotypes; if only
%                               one basis function is specified, then
%                               beta_se is [m x v] matrix
%
% tStats:           [m x q x v] ratio of beta_hat and beta_se
%
% logpValues:       [m x q x v] -log10 p values for the T statistics
%
%% Notes
% Results are only saved as a mat file if both outDir and outName are
% specified in the input
% 
% Notes for OLS:
% --------------
% Previously, using psuedo inverse to solve; however, using QR is a faster 
% solution. Therefore, replace:
%   beta_hat(snps, :, ivec_bin) = pinv(tmpX) * tmpY;
% with
%   beta_hat(snps, :, ivec_bin) = tmpX \ tmpY;
% and for standard error, replace:
%   beta_se(snps,  :, ivec_bin) = sqrt(diag(pinv(tmpX' * tmpX)) * sig2tvec(ivec_bin));
% with
%   beta_se(snps,  :, ivec_bin) = sqrt(diag((tmpX' * tmpX) \ eye(size(tmpX, 2)) * sig2tvec(ivec_bin)));
% 
% In the case of rank deficiency of tmpX, the use of mldivide or backslash 
% operator will result in zero estimate; in such cases, we could use pinv
%
% Since pinv(A) * B is the same as lsqminnorm(A,B), and that the behaviour 
% of lsqminnorm is the same as the mldivide or backslash operator, if no 
% rank deficiency. Therefore, using lsqminnorm for all solutions is easier.
% Note that lsqminnorm was introduced only in R2017b; therefore, for older 
% MATLAB check for rank deficiency every time and use appropriate solver
% (still faster than using pinv overall)
% -------------------------------------------------------------------------
% 
% An alternate solution using QR decomposition
% Replace 0 with 'econ' for newer MATLAB
% [Q, R, p]   = qr(Xvariables, 0);
% beta(p,:)   = R \ (Q \ yVariable);
% 
% Get variance-covariance matrix for coefficients or reuse R
% [~, R] = qr(Xvariables, 0);
% Rinv   = R \ eye(size(Xvariables, 2));
% RRt    = Rinv * Rinv';
% 
% Standard deviation of the residuals
% sigma2 = sum((yVariable - (Xvariables * coeff)).^2)./df;
% 
% Calculate standard error
% SE       = zeros(numX, numY);
% vcov     = zeros(numX, numX, numY);
% for yy   = 1:numY
%     vcov(:, :, yy) = sigma2(1,yy) * RRt;
%     SE(:,yy)       = sqrt(diag(vcov(:,:,yy)));
% end
% 
% Therefore, 
% [Q, R, p]                   = qr(tmpX, 0);
% beta_hat(snps, :, ivec_bin) = R \ (Q \ tmpY);
% Calculate variance covariance matrix of coefficients, or reuse R
% [~, R]  = qr(tmpX, 0);
% Rinv    = R \ eye(size(tmpX, 2));
% RRt     = Rinv * Rinv';
% vcov    = RRt * sig2tvec(ivec_bin); % or use sigma
% Calculate standard error
% beta_se(snps, :, ivec_bin)  = sqrt(diag(vcov));
% -------------------------------------------------------------------------
%
% Notes for GLS:
% --------------
% Previously, using pseudo inverse:
% XtWVsWtX                      = XtWVsWt(beginLocs(snps):endLocs(snps),:) * genovecNew(:, beginLocs(snps):endLocs(snps)); 
% Bi                            = pinv(XtWVsWtX);
% beta_hat(snps, :, ivec_bin)   = Bi * XtW(beginLocs(snps):endLocs(snps),:) * ymat(:,ivec_bin);
% beta_se(snps, :, ivec_bin)    = sqrt(diag(Bi) * sig2tvec(ivec_bin));
%
% Now, replaced the use of pinv for calculating Bi as:
% Bi = lsqminnorm(XtWVsWtX, eye(size(XtWVsWtX)));
%
%% Defaults:
% basisFunction:    intercept only (standard GWAS)
% OLSflag:          false
% SingleOrDouble:   double
% pValType:         'z'

%% Check inputs and assign defaults
% Check genoCell
if ~exist('genoCell', 'var') || isempty(genoCell)
    error('Please specify genoCell which contains genotyping data');
else
    if ~iscell(genoCell)
        genoCell = {genoCell};
    end
end

% Check y variables
if ~exist('ymat', 'var') || isempty(ymat)
    error('Please specify ymat which contains y variables to analyze');
end

% Check clusterinfo
if ~exist('clusterinfo', 'var') || isempty(clusterinfo)
    error('Please specify clusterinfo which is an output from FEMA_parse_family');
end

% Check binvec
if ~exist('binvec', 'var') || isempty(binvec)
    error('Please specify binvec');
end

% Check sig2tvec
if ~exist('sig2tvec', 'var') || isempty(sig2tvec)
    error('Please specify sig2tvec for every y variable');
else
    if size(sig2tvec, 2) ~= size(ymat, 2)
        error('Mismatch between number of y variables and number of sig2tvec entries');
    end
end

% Check GroupByFamType
if ~exist('GroupByFamType', 'var') || isempty(GroupByFamType)
    error('Please specify GroupByFamType');
end

% Check famtypevec
if ~exist('famtypevec', 'var') || isempty(famtypevec)
    if GroupByFamType
        error('Please specify famtypevec');
    end
end

% Check Ws_famtype
if ~exist('Ws_famtype', 'var') || isempty(Ws_famtype)
    if GroupByFamType
        error('Please specify Ws_famtype');
    end
end

% Check Ws_fam
if ~exist('Ws_fam', 'var') || isempty(Ws_fam)
    if ~GroupByFamType
        error('Please specify Ws_fam');
    end
end

% Check pValue type
if ~exist('pValType', 'var') || isempty(pValType)
    pValType = 'z';
else
    pValType = lower(pValType);
    if ~ismember(pValType, {'t', 'z', 'chi'})
        error('Unknown value for pValType specified');
    end
end

% Check degrees of freedom
if ~exist('df', 'var') || isempty(df)
    if strcmpi(pValType, 't')
        warning('Degrees of freedom not specified; using df = sample size - 1');
        df = size(genoCell{1}, 1) - 1;
    end
end

% Check basisFunction
if ~exist('basisFunction', 'var') || isempty(basisFunction)
    basisFunction = ones(size(genoCell{1}, 1), 1);
end

% Check OLSflag
if ~exist('OLSflag', 'var') || isempty(OLSflag)
    OLSflag = false;
end

% Check SingleOrDouble
if ~exist('SingleOrDouble', 'var') || isempty(SingleOrDouble)
    SingleOrDouble = 'double';
else
    SingleOrDouble =  lower(SingleOrDouble);
    if ~ismember(SingleOrDouble, {'single', 'double'})
        error('Unknown value for SingleOrDouble specified');
    end
end

% Check if output needs to be saved as mat files
if ~exist('outDir',  'var') || isempty(outDir) || ...
   ~exist('outName', 'var') || isempty(outName)
    writeResults = false;
else
    writeResults = true;
    if ~exist(outDir,  'dir')
        mkdir(outDir);
    end
end

% Determine if lsqminnorm can be used
if exist('lsqminnorm', 'file')
    useLSQ = true;
else
    useLSQ = false;
end

%% Initialize
numBasisFunc      = size(basisFunction, 2);
[numObs, numSNPs] = size(genoCell{1});
numYvars          = size(ymat, 2);

if numBasisFunc > 1
beta_hat          = zeros(numSNPs, numBasisFunc, numYvars, class(ymat));
else
    beta_hat      = zeros(numSNPs, numYvars, class(ymat));
end

beta_se           = zeros(size(beta_hat),  class(ymat));

%% Update degrees of freedom
% In the case of a single basis function (assuming intercept only), the
% actual GWAS df would be n - p where p is the number of covariates. For
% greater than 1 number of basis functions, then the df should be penalized
% by every increasing number of basis function beyond 1 - only relevant for
% when p values are calculated using t distribution
if numBasisFunc > 1
    df = df - numBasisFunc + 1;
end

%% Handle the case of only one basis function
if numBasisFunc == 1
    binLoc       = 1;
    for sig2bini = unique(binvec(isfinite(binvec)), 'stable')
        ivec_bin = find(binvec==sig2bini);

        % Genovec for this bin
        genovec  = genoCell{binLoc};

        % Cast genovec as double precision, if required
        if strcmpi(SingleOrDouble, 'double')
            genovec = double(genovec);
        end

        % Initialize
        XtW      = zeros(fliplr(size(genovec)), class(genovec)); 
        % XtWVsWt  = zeros(fliplr(size(genovec)), class(genovec));

        if GroupByFamType
            for fi = 1:length(clusterinfo)
                XtW(:, clusterinfo{fi}.jvec_fam)    = genovec(clusterinfo{fi}.jvec_fam,:)' * Ws_famtype{binLoc, famtypevec(fi)};
                % XtWVsWt(:,clusterinfo{fi}.jvec_fam) = XtW(:,clusterinfo{fi}.jvec_fam)      * Vs_famtype{binLoc, famtypevec(fi)} * Ws_famtype{binLoc, famtypevec(fi)}';
            end
        else
            for fi = 1:length(clusterinfo)
                XtW(:,clusterinfo{fi}.jvec_fam)     = genovec(clusterinfo{fi}.jvec_fam,:)' * Ws_fam{binLoc, fi};
                % XtWVsWt(:,clusterinfo{fi}.jvec_fam) = XtW(:,clusterinfo{fi}.jvec_fam)      * Vs_fam{binLoc, fi} * Ws_fam{binLoc, fi}';
            end
        end

        % Calculating beta coefficient for every SNP
        % beta = inv(X' * W * X) * X' * W * y
        % where W is the inverse of the total variance V
        %
        % Calculating the inner term
        % X' * W has been already calculated ---> XtWVsWt or XtW
        % Therefore, the inner term is the inverse of XtW * X
        % However, since each SNP is being treated independently, we are only
        % interested in the diagonal term of the XtW * X multiplication
        % This is equivalent to a sum over element-wise multiplication:
        % sum(XtW' .* X)
        % 
        % XtW:                  [chunk x subj] beta = inv(X' * W * X) * X' * W * yy
        % X or genovec:         [subj x chunk]
        % XtWVsWt' .* genovec:  [subj x chunk] 
        % Sum over the first dimension (subjects) will be [1 x chunk]
        %
        % If everything is double precision:
        % diag(XtW * genovec) approximately equal to sum(XtW' .* genovec)
        %
        % Since we are interested in the inverse of this inner term, and each
        % SNP is still being treated independently, therefore the inverse is
        % simply an inverse of a scaler, which is 1/value
        % For older MATLAB use: 1./(sum(bsxfun(@times, XtWVsWt', genovec)));
        % Bi = 1./(sum(XtWVsWt' .* genovec));
        Bi = 1./(sum(XtW' .* genovec));

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
        end

        % Standard error of beta
        beta_se(:,ivec_bin) = sqrt(Bi' * sig2tvec(ivec_bin));

        % Update binLoc
        binLoc = binLoc + 1;
    end
    
    % Calculate t statistics
    tStats = beta_hat ./ beta_se;
    
    % Calculate p value
    switch pValType
        case 't'
            logpValues = -log10(2 * tcdf(-abs(tStats), df));
            
        case 'z'
            logpValues = -log10(2 * normcdf(-abs(tStats), 0, 1));
            
        case 'chi'
            logpValues = -log10(chi2cdf(tStats.^2, 1, 'upper'));
    end
else
    binLoc       = 1;
    for sig2bini = unique(binvec(isfinite(binvec)), 'stable')
        ivec_bin = find(binvec==sig2bini);
        
        % Genovec for this bin
        genovec  = genoCell{binLoc};
        
        % Cast genovec as double precision, if required
        if strcmpi(SingleOrDouble, 'double')
            genovec = double(genovec);
        end
        
        % Use basis function to create new genotype matrix for this bin
        tmpGenotype = zeros(numObs, numBasisFunc, numSNPs);
        for basis = 1:numBasisFunc
            tmpGenotype(:, basis, :) = genovec .* basisFunction(:, basis);
        end
        
        % Reorder tmpGenotype such that for each SNP, the SNP*basisFucntion
        % columns are next to each other - possible that this might create
        % a very large matrix and potential RAM bottleneck?
        genovecNew = reshape(tmpGenotype, numObs, [], 1);
        
        % Get rid of tmpGenotype - might slightly increase execution time
        clear tmpGenotype
        
        % Initialize
        XtW      = zeros(fliplr(size(genovecNew)), class(genovec));
        % XtWVsWt  = zeros(fliplr(size(genovecNew)), class(genovec));
        
        if GroupByFamType
            for fi = 1:length(clusterinfo)
                XtW(:, clusterinfo{fi}.jvec_fam)    = genovecNew(clusterinfo{fi}.jvec_fam,:)' * Ws_famtype{binLoc, famtypevec(fi)};
                % XtWVsWt(:,clusterinfo{fi}.jvec_fam) = XtW(:,clusterinfo{fi}.jvec_fam)         * Vs_famtype{binLoc, famtypevec(fi)} * Ws_famtype{binLoc, famtypevec(fi)}';
            end
        else
            for fi = 1:length(clusterinfo)
                XtW(:,clusterinfo{fi}.jvec_fam)     = genovecNew(clusterinfo{fi}.jvec_fam,:)' * Ws_fam{binLoc, fi};
                % XtWVsWt(:,clusterinfo{fi}.jvec_fam) = XtW(:,clusterinfo{fi}.jvec_fam)         * Vs_fam{binLoc, fi} * Ws_fam{binLoc, fi}';
            end
        end
        
        % Figure out locations in genovecNew that need to be handled
        % together; essentially, these are X variables that should be
        % considered together
        beginLocs  = 1:numBasisFunc:size(genovecNew,2);
        endLocs    = 2:numBasisFunc:size(genovecNew,2);
        
        % Perform estimation - in this case, loop seems to be unavoidable
        % because parts of genovecNew need to be considered together - for
        % every SNP, there are as many columns as the number of basis 
        % functions. They cannot be considered together because that would
        % violate the assumption of each SNP acting independently (standard
        % GWAS assumption)
        
        % For this bin, the y value - easier than repeated indexing
        tmpY = ymat(:, ivec_bin);
        
        if OLSflag
            if useLSQ
                for snps = 1:numSNPs
                    % Find columns of genovecNew that correspond to this SNP
                    tmpX = genovecNew(:, beginLocs(snps):endLocs(snps));
                    
                    % For clarity, write out XtX
                    tmpXtX = tmpX' * tmpX;
                    
                    % Get estimates and standard error
                    beta_hat(snps, :, ivec_bin) = lsqminnorm(tmpX, tmpY);
                    beta_se(snps,  :, ivec_bin) = sqrt(diag(lsqminnorm(tmpXtX, eye(size(tmpX,2)))) * sig2tvec(ivec_bin));
                end
            else
                for snps = 1:numSNPs
                    % Find columns of genovecNew that correspond to this SNP
                    tmpX = genovecNew(:, beginLocs(snps):endLocs(snps));
                    
                    % For clarity, write out XtX
                    tmpXtX = tmpX' * tmpX;
                    
                    % Get estimates and standard error
                    if rank(tmpX) < size(tmpX, 2)
                        beta_hat(snps, :, ivec_bin) = pinv(tmpX) * tmpY;
                        beta_se(snps,  :, ivec_bin) = sqrt(diag(pinv(tmpXtX)) * sig2tvec(ivec_bin));
                     else
                         beta_hat(snps, :, ivec_bin) = tmpX \ tmpY;
                         beta_se(snps,  :, ivec_bin) = sqrt(diag(tmpXtX \ eye(size(tmpX, 2)) * sig2tvec(ivec_bin)));
                    end
                end
            end
        else
             if useLSQ
                 for snps = 1:numSNPs
                     % Find columns of genovecNew that correspond to this SNP
                     tmpX = genovecNew(:, beginLocs(snps):endLocs(snps));
                     
                     % Solve
                     term1                        = XtW(beginLocs(snps):endLocs(snps),:) * tmpX;
                     Bi                           = lsqminnorm(term1, eye(size(term1))); 
                     beta_hat(snps, :, ivec_bin)  = Bi * (XtW(beginLocs(snps):endLocs(snps), :) * ymat(:,ivec_bin));
                     beta_se(snps,  :, ivec_bin)  = sqrt(diag(Bi) * sig2tvec(ivec_bin));
                 end
             else
                 for snps = 1:numSNPs
                     % Find columns of genovecNew that correspond to this SNP
                     tmpX = genovecNew(:, beginLocs(snps):endLocs(snps));
                     
                     % Prepare XtW * X
                     term1 = XtW(beginLocs(snps):endLocs(snps),:) * tmpX;
                     
                     % Prepare inverse of term1
                     if rank(tmpX) < size(tmpX, 2)
                         Bi = pinv(term1);
                     else
                         Bi = term1 \ eye(size(term1));
                     end
                     
                     % Get coefficient and standard error
                     beta_hat(snps, :, ivec_bin) = Bi * (XtW(beginLocs(snps):endLocs(snps), :) * ymat(:,ivec_bin));
                     beta_se(snps,  :, ivec_bin) = sqrt(diag(Bi) * sig2tvec(ivec_bin));
                 end
             end
        end
        
        % Update binloc
        binLoc = binLoc + 1;
    end
    
    % Now compute t statistics
    tStats = beta_hat ./ beta_se;
    
    % Finally compute -log10 p values
    switch pValType
        case 't'
            logpValues = -log10(2 * tcdf(-abs(tStats), df));
            
        case 'z'
            logpValues = -log10(2 * normcdf(-abs(tStats), 0, 1));
            
        case 'chi'
            logpValues = -log10(chi2cdf(tStats.^2, 1, 'upper'));
    end
end

% Save variables, if necessary
if writeResults

    % Save name - make sure .mat is not already part of outName
    saveName = fullfile(outDir, [strrep(outName, '.mat', ''), '.mat']);

    % Get an estimate of size of the variables
    tmpInfo = whos;
    if sum([tmpInfo(ismember({tmpInfo(:).name}', {'beta_hat', 'beta_se', 'tStats', 'pValues', 'df'})).bytes]) > 2^31
        save(saveName, 'beta_hat', 'beta_se', 'tStats', 'logpValues', 'df', '-v7.3');
    else
        save(saveName, 'beta_hat', 'beta_se', 'tStats', 'logpValues', 'df');
    end
end