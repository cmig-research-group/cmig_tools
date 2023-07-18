function [beta_hat, beta_se, tStats, logpValues] = ...
          FEMA_sig2binseg_parfeval_GWAS_bf(genoMat, ymat, binvec, sig2tvec, varargin)
% Function to estimate the effect of each (residualized) SNP on
% (residualized) phenotype using OLS or GLS
%% Inputs:
% genoMat:          [n x m]     residualized genotype matrix containing
%                               genotype data of n subjects and m SNPs; 
%                               (see FEMA_residualizeGenotype and
%                               FEMA_OLSResiduals; assumes OLS residuals)
%
% ymat:             [n x v]     matrix of n subjects and v phenotypes
%                               (residual phenotype(s) from FEMA_fit)
%  
% binvec:           [1 x b]     vector of bin values on which random
%                               effects are evaluated
% 
% sig2tvec:         [1 x v]     total residual variance of v phenotypes
% 
%% Optional inputs:
% bfSNP:            [n x q]     matrix of basis functions to use as an
%                               interaction term for (residualized) 
%                               genotyping matrix; defaults to an
%                               intercept term; the same set of basis
%                               functions are used for all bins
%
% sex:              [n x 2]     dummy coded two-column variable indicating
%                               the sex for every observation (used if
%                               snp*sex interaction is required)
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
% SingleOrDouble:   character   single or double precision
%
% outDir:           character   full path to where the results should be
%                               saved
%
% outName:          character   name of the mat file to use for saving
%
% useShortcut:      logical     returned by FEMA_compileTerms indicating if
%                               a faster computation can be employed; if  
%                               true, allJVec and allWsTerms should be 
%                               present (default is false)
%
% Additional inputs if useShortcut is false; alternatively, pass
% reusableVars returned by FEMA_fit
%
% clusterinfo:      [1 x f]     cell type of f families/clusters, returned
%                               as an output from FEMA_parse_family
%
% GroupByFamType:   logical     true/false indicating if GroupByFamType
% 
% famtypevec:       [1 x f]     vector indicating family type
% 
% Ws_famtype, Ws_fam:           inverses of the V term returned by 
%                               FEMA_compileTerms (for every bin)
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
% sex               []
% OLSflag:          false
% SingleOrDouble:   double
% pValType:         'z'
% useShortcut       false

%% Checking mandatory inputs
% Check genoMat
if ~exist('genoMat', 'var') || isempty(genoMat)
    error('Please specify genoMat which contains genotyping data');
end

% Check y variables
if ~exist('ymat', 'var') || isempty(ymat)
    error('Please specify ymat which contains y variables to analyze');
end

% Check binvec
if ~exist('binvec', 'var') || isempty(binvec)
    error('Please specify binvec');
% else
    % Number of unique bins
    % numBins = length(unique(binvec, 'stable'));
end

% Check sig2tvec
if ~exist('sig2tvec', 'var') || isempty(sig2tvec)
    error('Please specify sig2tvec for every y variable');
else
    if size(sig2tvec, 2) ~= size(ymat, 2)
        error('Mismatch between number of y variables and number of sig2tvec entries');
    end
end

%% Checking optional inputs
p = inputParser;

% Add parameters to parsed inputs, if they are missing
p.addParameter('bfSNP',          ones(size(genoMat, 1), 1));
p.addParameter('sex',            []);
p.addParameter('pValType',       'z');
p.addParameter('df',             []);
p.addParameter('OLSflag',        false);
p.addParameter('SingleOrDouble', 'double');
p.addParameter('outDir',         '');
p.addParameter('outName',        '');
p.addParameter('useShortcut',    false);
p.addParameter('allJVec',        []);
p.addParameter('allWsTerms',     []);
p.addParameter('reusableVars',   []);
p.addParameter('clusterinfo',    []);
p.addParameter('GroupByFamType', []);
p.addParameter('famtypevec',     []);
p.addParameter('Ws_famtype',     []);
p.addParameter('Ws_fam',         []);

% Parse varargin
parse(p, varargin{:});

% Assign values to optional inputs
bfSNP           = p.Results.bfSNP;
sex             = p.Results.sex;
pValType        = p.Results.pValType;
df              = p.Results.df;
OLSflag         = p.Results.OLSflag;
SingleOrDouble  = p.Results.SingleOrDouble;
outDir          = p.Results.outDir;
outName         = p.Results.outName;
useShortcut     = p.Results.useShortcut;
allJVec         = p.Results.allJVec;
allWsTerms      = p.Results.allWsTerms;
reusableVars    = p.Results.reusableVars;
clusterinfo     = p.Results.clusterinfo;
GroupByFamType  = p.Results.GroupByFamType;
famtypevec      = p.Results.famtypevec;
Ws_famtype      = p.Results.Ws_famtype;
Ws_fam          = p.Results.Ws_fam;

% Deprecated slow solution
classicalSolution = false;

% Check whether SNP should interact with sex
if isempty(sex)
    interactSex = false;
else
    if size(sex, 2) ~= 2
        warning('Expected sex to have 0/1 coded two column matrix; skipping interaction with sex');
        interactSex = false;
    else
        interactSex = true;
    end
end

% Check pValue type
pValType = lower(pValType);
if ~ismember(pValType, {'t', 'z', 'chi'})
    error('Unknown value for pValType specified');
end

% Check degrees of freedom
if isempty(df) && strcmpi(pValType, 't')
    warning('Degrees of freedom not specified; using df = sample size - 1');
    df = size(genoMat, 1) - 1;
end

% Check SingleOrDouble
SingleOrDouble = lower(SingleOrDouble);
if ~ismember(SingleOrDouble, {'single', 'double'})
    error('Unknown value for SingleOrDouble specified');
end

% Determine if output needs to be saved as mat files
if isempty(outDir) || isempty(outName)
    writeResults = false;
else
    writeResults = true;
    if ~exist(outDir,  'dir')
        mkdir(outDir);
    end
end

% Check if required variables are present, if useShortcut
if useShortcut
    if isempty(allJVec) || isempty(allWsTerms)
        error('allJVec and allWsTerms are required inputs if useShortcut is true');
    end
else
    % If reusableVars is present, parse that
    if ~isempty(reusableVars)
        clusterinfo     = reusableVars.clusterinfo;
        GroupByFamType  = reusableVars.GroupByFamType;
        famtypevec      = reusableVars.famtypevec;
        Ws_famtype      = reusableVars.Ws_famtype;
        Ws_fam          = reusableVars.Ws_fam;
    else
        % Check individual required variables
        % Check clusterinfo
        if isempty(clusterinfo)
            error('Please specify clusterinfo');
        end

        % Check GroupByFamType
        if isempty(GroupByFamType)
            error('Please specify GroupByFamType');
        end

        % Check famtypevec
        if isempty(famtypevec) && GroupByFamType
            error('Please specify famtypevec');
        end

        % Check Ws_famtype
        if isempty(Ws_famtype) && GroupByFamType
            error('Please specify Ws_famtype');
        end

        % Check Ws_fam
        if isempty(Ws_fam) && ~GroupByFamType
            error('Please specify Ws_fam');
        end
    end
end

% Determine if lsqminnorm can be used
if exist('lsqminnorm', 'file')
    useLSQ = true;
else
    useLSQ = false;
end

%% Initialize
[numObs, numSNPs] = size(genoMat);
numYvars          = size(ymat, 2);

% Update basis function, if required
if interactSex
    bfSNP = [bfSNP .* sex(:,1), bfSNP .* sex(:,2)];
end
numBasisFunc      = size(bfSNP, 2);

if numBasisFunc > 1
    beta_hat      = zeros(numSNPs, numBasisFunc, numYvars, class(ymat));
else
    beta_hat      = zeros(numSNPs, numYvars, class(ymat));
end
beta_se           = zeros(size(beta_hat), class(ymat));

%% Update degrees of freedom
% In the case of the basis function only having an intercept term, the df
% is n - p - h, where p is the number of covariates and h is the number of
% columns in the timeMatrix. When additional basis functions are specified,
% then the df is n - p - q - h, where h are the number of columns in the
% basis function matrix. The user has already provided df = n - p.
% Therefore, the updated degrees of freedom are df - q - h + 1, where the
% additional +1 only accounts for the additional columns of the basis
% function (beyond the intercept). Since q + h = totalColumns, the updated
% df is simply df - totalColums + 1. Note that this is only relevant for
% calculating p values when using t distribution
% df = df - totalColumns + 1;
df = df - numBasisFunc + 1;

%% Handle the case of only one total column (standard GWAS)
if numBasisFunc == 1
    binLoc       = 1;
    for sig2bini = unique(binvec(isfinite(binvec)), 'stable')
        ivec_bin = find(binvec==sig2bini);

        % Cast genoMat as double precision, if required
        if strcmpi(SingleOrDouble, 'double')
            genoMat = double(genoMat);
        end

        % Account for OLS case
        if OLSflag
            % In this case, the solution is:
            % beta       = inv(X' * X) * X' * y
            % covariance = sigma2 * inv(X' * X)
            % SE         = sqrt(diag(covariance))
            % Since each SNP is being treated independently, we can use
            % element-wise multiplication (see notes in the GLS case below)

            % First compute X' * X
            term1 = sum(genoMat .* genoMat);

            % Inverse of term1
            iterm1 = 1./term1;

            % Calculate beta coefficient
            beta_hat(:, ivec_bin) = iterm1' .* (genoMat' * ymat(:,ivec_bin));

            % Calculate standard error
            beta_se(:, ivec_bin) = sqrt(iterm1' * sig2tvec(ivec_bin));
        else
            % Initialize
            XtW      = zeros(fliplr(size(genoMat)), class(genoMat));
            if ~useShortcut
                if GroupByFamType
                    for fi = 1:length(clusterinfo)
                        XtW(:, clusterinfo{fi}.jvec_fam) = genoMat(clusterinfo{fi}.jvec_fam,:)' * Ws_famtype{binLoc, famtypevec(fi)};
                    end
                else
                    for fi = 1:length(clusterinfo)
                        XtW(:,clusterinfo{fi}.jvec_fam)  = genoMat(clusterinfo{fi}.jvec_fam,:)' * Ws_fam{binLoc, fi};
                    end
                end
            else
                % Jvec and Ws terms for this bin
                currJVec    = allJVec{binLoc};
                currWsTerm  = allWsTerms{binLoc};
                XtW         = genoMat(currJVec, :)' * currWsTerm;
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
            Bi = 1./(sum(XtW' .* genoMat));

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
            % beta_hat = bsxfun(@times, Bi', XtW) * ymat(:, ivec_bin);
            beta_hat(:,ivec_bin) = Bi' .* (XtW * ymat(:, ivec_bin));

            % Standard error of beta
            beta_se(:,ivec_bin) = sqrt(Bi' * sig2tvec(ivec_bin));
        end

        % Update binLoc
        binLoc = binLoc + 1;
    end
else
    binLoc       = 1;
    for sig2bini = unique(binvec(isfinite(binvec)), 'stable')
        ivec_bin = find(binvec==sig2bini);
                
        % Cast genovec as double precision, if required
        if strcmpi(SingleOrDouble, 'double')
            genoMat = double(genoMat);
        end

        % Jvec and Ws terms for this bin
        if useShortcut
            currJVec    = allJVec{binLoc};
            currWsTerm  = allWsTerms{binLoc};
        else
            if ~classicalSolution
                % Use basis function to create new genotype matrix for this bin
                tmpGenotype = zeros(numObs, numBasisFunc, numSNPs);
                for basis   = 1:numBasisFunc
                    tmpGenotype(:, basis, :) = genoMat .* bfSNP(:, basis);
                end

                % Reorder tmpGenotype such that for each SNP, the SNP*basisFucntion
                % columns are next to each other - possible that this might create
                % a very large matrix and potential RAM bottleneck?
                genoMat = reshape(tmpGenotype, numObs, [], 1);

                % Get rid of tmpGenotype - might slightly increase execution time
                clear tmpGenotype

                % Compile XtW
                XtW = zeros(fliplr(size(genoMat)), class(genoMat));

                if GroupByFamType
                    for fi = 1:length(clusterinfo)
                        XtW(:, clusterinfo{fi}.jvec_fam) = genoMat(clusterinfo{fi}.jvec_fam,:)' * Ws_famtype{binLoc, famtypevec(fi)};
                    end
                else
                    for fi = 1:length(clusterinfo)
                        XtW(:,clusterinfo{fi}.jvec_fam) = genoMat(clusterinfo{fi}.jvec_fam,:)' * Ws_fam{binLoc, fi};
                    end
                end

                % Figure out locations in genoMat that need to be handled
                % together; essentially, these are X variables that should be
                % considered together
                beginLocs  = 1:numBasisFunc:size(genoMat,2);
                endLocs    = numBasisFunc:numBasisFunc:size(genoMat,2);
            end
        end
                
        % For this bin, the y value - faster than repeated indexing
        tmpY = ymat(:, ivec_bin);

        % Perform estimation - in this case, loop seems to be unavoidable
        % because parts of genovecNew need to be considered together - for
        % every SNP, there are as many columns as the number of basis
        % functions. They cannot be considered together because that would
        % violate the assumption of each SNP acting independently (standard
        % GWAS assumption)
        for snps = 1:numSNPs

            if useShortcut
                % Prepare matrix of X variables
                % tmpX = genoMat(:, snps) .* bfSNP;
                tmpX   = genoMat(currJVec, snps) .* bfSNP;
            else
                if ~classicalSolution
                    tmpX = genoMat(:, beginLocs(snps):endLocs(snps));
                end
            end

            % Handle case of OLS solution
            if OLSflag
                % beta       = inv(X' * X) * X' * y
                % covariance = sigma2 * inv(X' * X)
                % SE         = sqrt(diag(covariance))

                % First compute X' * X
                term1 = tmpX' * tmpX;

                % Inverse of term1
                % iterm1 = term1 \ I;
                I = eye(size(term1));
                if useLSQ
                    iterm1 = lsqminnorm(term1, I);
                else
                    if rank(tmpX) < size(tmpX, 2)
                        iterm1 = pinv(term1);
                    else
                        iterm1 = term1 \ I;
                    end
                end

                % Calculate beta coefficient
                beta_hat(snps, :, ivec_bin) = iterm1 * (tmpX' * ymat(:,ivec_bin));

                % Calculate standard error
                beta_se(snps, :, ivec_bin) = sqrt(diag(iterm1) * sig2tvec(ivec_bin));
            else
                % Prepare XtW
                if useShortcut
                    % XtW = tmpX(currJVec, :)' * currWsTerm;
                    XtW = tmpX' * currWsTerm;

                    % XtW is already transposed
                    term1 = XtW * tmpX;
                else
                    if classicalSolution
                        XtW = zeros(fliplr(size(tmpX)), class(genoMat)); %#ok<UNRCH>
                        if GroupByFamType 
                            for fi = 1:length(clusterinfo)
                                XtW(:, clusterinfo{fi}.jvec_fam) = tmpX(clusterinfo{fi}.jvec_fam,:)' * Ws_famtype{binLoc, famtypevec(fi)};
                            end
                        else
                            for fi = 1:length(clusterinfo)
                                XtW(:,clusterinfo{fi}.jvec_fam)  = tmpX(clusterinfo{fi}.jvec_fam,:)' * Ws_fam{binLoc, fi};
                            end
                        end

                        % XtW is already transposed
                        term1 = XtW * tmpX;
                    else
                        term1 = XtW(beginLocs(snps):endLocs(snps),:) * tmpX;
                    end
                end

                % Calculate inverse of term1
                % Since term1 dimensionality will typically be small
                % (numBasisFunc * numBasisFunc), it might be more 
                % efficient to use backslash instead of lsqminnorm
                I       = eye(size(term1));
                % iterm1  = term1 \ I;
                if useLSQ
                    iterm1 = lsqminnorm(term1, I);
                else
                    if rank(tmpX) < size(tmpX, 2)
                        iterm1 = pinv(term1);
                    else
                        iterm1 = term1 \ I;
                    end
                end

                % Calculate beta coefficient and standard error
                if ~useShortcut && ~classicalSolution
                    beta_hat(snps, :, ivec_bin) = iterm1 * (XtW(beginLocs(snps):endLocs(snps), :) * tmpY);
                else
                    beta_hat(snps, :, ivec_bin) = iterm1 * (XtW * tmpY);
                end
                beta_se(snps,  :, ivec_bin) = sqrt(diag(iterm1) * sig2tvec(ivec_bin));
            end
        end
        
        % Update binloc
        binLoc = binLoc + 1;
    end
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