function [beta_hat, beta_se, tStats, logpValues, Wald_F, Wald_p, coeffCovar, errFlag, timeStruct] = ...
          FEMA_fit_GWAS(genoMat, ymat_res, binvec, Xvars, allWsTerms, varargin)
% Function to perform GWAS using FEMA
%% Inputs:
% genoMat:          [n x m] OR  matrix containing genotype data of n
%                   structure   subjects and m SNPs; OR can be a structure
%                               which reads genotyping data for specified
%                               subjects and SNPs; should contain the
%                               following fields (see FEMA_parse_PLINK):
%                   * 'fname':      full path to a PLINK bed file
%                   * 'iid':        [n x 1] cell type having subject IDs
%                   * 'meanImpute': true or false
%                   * 'roundOff':   true or false
%                   * 'transform':  one of the following: 
%                                       * 'center'
%                                       * 'centre'
%                                       * 'std'
%                                       * 'none'
%                   * 'stdType':    if transform is 'std', then one of the
%                                   following: 'emperical' or 'gcta'
%                   * 'snpList':    list of SNPs to read
%
% ymat_res:         [n x v]     matrix of n subjects and v (GLS) residual
%                               phenotypes from FEMA_fit
%  
% binvec:           [1 x b]     vector of bin values on which random
%                               effects were evaluated in FEMA_fit
% 
% Xvars:            [n x p]     matrix containing X variables (or
%                               covariates) that were used for creating
%                               ymat_res when calling FEMA_fit
%
% allWsTerms:       cell        output from FEMA_compileTerms containing
%                               the inverse of the V term for each bin
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
% L:                [l x q]     contrast vector (l == 1) or matrix (l > 1)
%                               specifying the linear combination of
%                               weights that should be tested for null
%                               hypothesis; alternatively can be a cell
%                               type with each cell containing a vector or
%                               matrix of linear combination of weights
%                               (see FEMA_WaldTest)
%
% hypValue:         [1 x k]     numeric scalar or [1 x k] vector (where k
%                               is the number of cell in L) containing the
%                               hypothesised value against which the null
%                               hypothesis will be tested; defaults to zero
%                               (see FEMA_WaldTest)
%
% doF:              logical     true or false indicating if F test should
%                               be conducted instead of the default Wald
%                               test, and then use F distribution for
%                               calculating the p values (default: false) 
%                               (see FEMA_WaldTest)
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
% adjustMSE:        logical     specifies whether sig2tvec should be
%                               adjusted for the number of covariates and
%                               the number of basis functions or not
%                               (default is true)
%
% doCoeffCovar      logical     true or false indicating if coefficient
%                               covariance for every SNP needs to be
%                               computed (and optionally saved); this can
%                               be used, for example, in the future for
%                               Wald test but can take up a large amount of
%                               disk space; is not computed if there are no
%                               basis functions (i.e., for standard GWAS)
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
% Wald_F:           [m x v x l] Wald statistic (omnibus test across q basis
%                               functions) for m SNPs, v phenotypes, and l
%                               weights specified in L
%
% Wald_p:           [m x v x l] p value associated with omnibus Wald test
%
% coeffCovar:   [m x q x q x v] coefficient covariance for every SNP for q
%                               basis functions and v phenotypes; not
%                               computed if there are no basis functions
%
% timeStruct:       structure   contains various timing information
% 
%% Notes
% Results are only saved as a mat file if both outDir and outName are
% specified in the input
% 
% Coefficient covariance is only computed (and saved) if basis functions
% are specified
% 
%% Defaults:
% basisFunction:    intercept only (standard GWAS)
% sex               []
% OLSflag:          false
% SingleOrDouble:   'double'
% pValType:         'z'
% adjustMSE:        true
% L:                any basis functions eye(numBasisFunctions)
% hypValue:         0
% doF:              false
% doCoeffCovar:     false
% 
%% Additional defaults when genStruct is specified:
% iid:              [] (i.e., read all subjects)
% onlyCheck:        false
% lowMem:           false
% meanImpute:       true
% roundOff:         false
% transform:        'none'
% stdType:          'none'
% SNPList:          [] (i.e., read all SNPs)
%
%% Additional notes:
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

%% Checking mandatory inputs
tOver  = tic;
tCheck = tic;

% Check genoMat
if ~exist('genoMat', 'var') || isempty(genoMat)
    error('Please specify genoMat which contains genotyping data or specify a structure with required fields to read genotyping data');
else
    if isstruct(genoMat)
        genStruct = genoMat;
        toRead = true;
    else
        if ismatrix(genoMat)
            toRead = false;
        else
            error('Unknown genoMat data type');
        end
    end
end

% Check y variables
if ~exist('ymat_res', 'var') || isempty(ymat_res)
    error('Please specify ymat_res which contains residual y variables to analyze');
end

% Check binvec
if ~exist('binvec', 'var') || isempty(binvec)
    error('Please specify binvec');
end

% Check Xvars
if ~exist('Xvars', 'var') || isempty(Xvars)
    warning('Xvars not provided; we hope you know what you are doing!');
    Xvars = [];
    regressX = false;
else
    regressX = true;
end

% Check fields for genStruct
if toRead
    % Make sure file name is specified
    if ~isfield(genStruct, 'fname')
        error('genStruct should have fname field having full path to a PLINK bed file');
    end

    % If iid is missing, read full data
    if ~isfield(genStruct, 'iid')
        genStruct.iid = [];
    end

    % If meanImpute is not specified, set it to true
    if ~isfield(genStruct, 'meanImpute')
        genStruct.meanImpute = true;
    end

    % If roundOff is not specified, set it to false
    if ~isfield(genStruct, 'roundOff')
        genStruct.roundOff = false;
    end

    % If transform is not specified, set it to none
    if ~isfield(genStruct, 'transform')
        genStruct.transform = 'none';
    end

    % If stdType is not specified, set it to none
    if ~isfield(genStruct, 'stdType')
        genStruct.stdType = 'none';
    end
    
    % If snpList is not specified, read all SNPs
    if ~isfield(genStruct, 'snpList')
        genStruct.snpList = '';
    end

    % Set onlyCheck to false, unless user has overriden this
    if ~isfield(genStruct, 'onlyCheck')
        genStruct.onlyCheck = false;
    end

    % Set lowMem to false, unless user has overriden this
    if ~isfield(genStruct, 'lowMem')
        genStruct.lowMem = false;
    end

    % Check if genInfo structure is present
    if ~isfield(genStruct, 'genInfo')
        genStruct.genInfo = [];
    end
end

%% Checking optional inputs
p = inputParser;

% Add parameters to parsed inputs, if they are missing
p.addParameter('bfSNP',          []);
p.addParameter('sex',            []);
p.addParameter('pValType',       'z');
p.addParameter('df',             []);
p.addParameter('OLSflag',        false);
p.addParameter('SingleOrDouble', 'double');
p.addParameter('outDir',         '');
p.addParameter('outName',        '');
p.addParameter('useShortcut',    false);
p.addParameter('useLSQ',         []);
p.addParameter('lowRank',        []);
p.addParameter('L',              []);
p.addParameter('hypValue',       []);
p.addParameter('doF',            false);
p.addParameter('doCoeffCovar',   false);

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
useLSQ          = p.Results.useLSQ;
lowRank         = p.Results.lowRank;
L               = p.Results.L;
hypValue        = p.Results.hypValue;
doF             = p.Results.doF;
doCoeffCovar    = p.Results.doCoeffCovar;

% Check basis functions
if ~exist('bfSNP', 'var') || isempty(bfSNP)
    bfSNP = ones(size(ymat_res, 1), 1);
else
    % Make sure correct number of rows exist
    if size(bfSNP, 1) ~= size(ymat_res, 1)
        error('Mismatch between number of rows in ymat_res and in bfSNP');
    end
end

% Check whether SNPs (*age) should interact with sex
if ~isempty(sex)
    if size(sex, 2) ~= 2
        warning('Expected sex to have 0/1 coded two column matrix; skipping interaction with sex');
    else
        % Update basis function(s)
        bfSNP = [bfSNP .* sex(:,1), bfSNP .* sex(:,2)];
    end
end
numBasisFunc = size(bfSNP, 2);

% Check pValue type
pValType = lower(pValType);
if ~ismember(pValType, {'t', 'z', 'chi'})
    error('Unknown value for pValType specified');
end

% Check degrees of freedom
if isempty(df)
    df = size(ymat_res, 1) - (size(Xvars, 2) + numBasisFunc);
end

% Check SingleOrDouble
SingleOrDouble = lower(SingleOrDouble);
if ~ismember(SingleOrDouble, {'single', 'double'})
    error('Unknown value for SingleOrDouble specified');
end

% Check L
if ~exist('L', 'var') || isempty(L)
    L{1} = eye(numBasisFunc);
else
    if ~iscell(L)
        L = {L};
    end
end
numCon = length(L);

% Determine if output needs to be saved as mat files
if isempty(outDir) || isempty(outName)
    writeResults = false;
else
    writeResults = true;
    if ~exist(outDir,  'dir')
        mkdir(outDir);
    end
end

% Determine if lsqminnorm can be used
if isempty(useLSQ)
    if exist('lsqminnorm', 'file')
        useLSQ = true;
    else
        useLSQ = false;
    end
end

% Determine if Xvars matrix is low rank
if isempty(lowRank)
    if rank(Xvars) < size(Xvars, 2)
        lowRank = true;
    else
        lowRank = false;
    end
end

% Check allWsTerms
if ~exist('allWsTerms', 'var') || isempty(allWsTerms)
    if ~OLSflag
        error('allWsTerms should be provided when OLSflag is false');
    end
end

% Save timing information for checking inputs
timeStruct.tCheck = toc(tCheck);

%% Read genotyping data (if required)
if toRead
    [genoMat, ~, ~, ~, ~, ~, ~, timingRead] =                                          ...
     FEMA_parse_PLINK(genStruct.fname,     genStruct.iid,       genStruct.snpList,     ...
                      genStruct.onlyCheck, genStruct.lowMem,    genStruct.meanImpute,  ...
                      genStruct.roundOff,  genStruct.transform, genStruct.stdType,     ...
                      genStruct.genInfo);
    timeStruct.tRead = timingRead;
else
    genStruct = '';
end

%% Cast genoMat as double precision, if required
tDouble = tic;
if strcmpi(SingleOrDouble, 'double')
    genoMat = double(genoMat);
end
timeStruct.tDouble = toc(tDouble);

%% Initialize
t1                = tic;
[numObs, numSNPs] = size(genoMat);
numYvars          = size(ymat_res, 2);
binLoc            = 1;

% Initialize coefficients and standard error
if numBasisFunc > 1
    beta_hat = zeros(numSNPs, numBasisFunc, numYvars, class(ymat_res));
    beta_se  = zeros(numSNPs, numBasisFunc, numYvars, class(ymat_res));
else
    beta_hat = zeros(numSNPs, numYvars, class(ymat_res));
    beta_se  = zeros(numSNPs, numYvars, class(ymat_res));
end

%% Check if this is a case of standard GWAS
if numBasisFunc == 1
    if any(ones(numObs, 1) - bfSNP)
        standardGWAS = false;
    else
        standardGWAS = true;
    end
else
    standardGWAS = false;
end

%% Initialize omnibus test
if ~standardGWAS
    [Wald_F, Wald_p] = deal(zeros(numSNPs, numYvars, numCon));
else
    Wald_F = [];
    Wald_p = [];
end

%% Initialize coefficient covariance if required
if doCoeffCovar
    if standardGWAS
        coeffCovar = [];
    else
        coeffCovar = zeros(numSNPs, numBasisFunc, numBasisFunc, numYvars, class(ymat_res));
    end
else
    coeffCovar = [];
end

timeStruct.tinitilize  = toc(t1);

%% Do estimation
tEstimation = tic;

% Get warning statuses for singular and nearly singular cases;
% temporarily set their display off
statusSingular = warning('off', 'MATLAB:singularMatrix');
statusNearSing = warning('off', 'MATLAB:nearlySingularMatrix');

% Clear last warning
lastwarn('');

if standardGWAS
    if OLSflag
        % In case of OLS, there are no bins
        numY = size(ymat_res, 2);

        % Residualize genotype for X variables
        % genoMat can be overwritten as we don't need it later
        if regressX
            genoMat = FEMA_OLSResiduals(genoMat, Xvars);
        end

        % Calculate inverse of X'X: since each SNP is treated
        % independently, we can use element-wise multiplication for XtX
        % and then use 1/scalar for the inverse
        % term1 = sum(genoMat .* genoMat);
        iterm1 = 1./sum(genoMat .* genoMat);

        % Calculate beta coefficient
        tmpCoeff = iterm1' .* (genoMat' * ymat_res);
        beta_hat = tmpCoeff;

        % Calculate standard error for every phenotype
        for yy = 1:numY

            % Mean squared error - where mean is calculated using df
            tmp_MSE = sum((ymat_res(:,yy) - (genoMat .* tmpCoeff(:,yy)')).^2)/df;

            % Standard error
            beta_se(:, yy) = sqrt(iterm1 .* tmp_MSE);
        end
    else
        for sig2bini = unique(binvec(isfinite(binvec)), 'stable')
            % Find the bins
            ivec_bin = find(binvec == sig2bini);

            % Current value of y
            tmpY = ymat_res(:,ivec_bin);

            % First, get the inv(V) term for this bin
            currWsTerm = allWsTerms{binLoc};

            % Get GLS residuals for each SNP
            if regressX
                res_genoMat = FEMA_GLSResiduals(genoMat, Xvars, [], [], currWsTerm);
            else
                res_genoMat = genoMat; % Useless copy
            end

            % XtW term
            XtW = res_genoMat' * currWsTerm;

            for snps = 1:numSNPs
                Bi                       = 1/(XtW(snps,:) * res_genoMat(:,snps));
                tmpCoeff                 = Bi * XtW(snps,:) * tmpY;
                beta_hat(snps, ivec_bin) = tmpCoeff;
                tmp_MSE                  = sum((tmpY - (res_genoMat(:,snps) * tmpCoeff)).^2)/df;
                beta_se(snps, ivec_bin)  = sqrt(Bi .* tmp_MSE);
            end
            
            % % Calculating the inner term
            % % X' * W has been already calculated ---> XtWVsWt or XtW
            % % Therefore, the inner term is the inverse of XtW * X
            % % However, since each SNP is being treated independently, we are only
            % % interested in the diagonal term of the XtW * X multiplication
            % % This is equivalent to a sum over element-wise multiplication:
            % % sum(XtW' .* X)
            % %
            % % XtW:                  [chunk x subj] beta = inv(X' * W * X) * X' * W * yy
            % % X or genovec:         [subj x chunk]
            % % XtWVsWt' .* genovec:  [subj x chunk]
            % % Sum over the first dimension (subjects) will be [1 x chunk]
            % %
            % % If everything is double precision:
            % % diag(XtW * genovec) approximately equal to sum(XtW' .* genovec)
            % %
            % % Since we are interested in the inverse of this inner term, and each
            % % SNP is still being treated independently, therefore the inverse is
            % % simply an inverse of a scaler, which is 1/value
            % % For older MATLAB use: 1./(sum(bsxfun(@times, XtWVsWt', genovec)));
            % % Bi = 1./(sum(XtWVsWt' .* genovec));
            % Bi = 1./(sum(XtW' .* res_genoMat));
            % 
            % % For the next set of multiplications:
            % % beta      = Bi * X'  * W * y
            % % or, beta  = Bi * XtW * y
            % %
            % % Bi:   [1     x chunk]
            % % XtW:  [chunk x subj]
            % % ymat: [subj  x bins]
            % %
            % % It is faster to implement XtW * ymat, resulting in [chunk x bins]
            % % which is then multiplied element-wise by Bi
            % % beta_hat(:,ivec_bin) = Bi' .* (XtW * ymat_res(:, ivec_bin));
            % tmpCoeff             = Bi' .* (XtW * tmpY);
            % beta_hat(:,ivec_bin) = tmpCoeff;
            % 
            % % Calculate standard error - likely faster to loop over y 
            % % variables in ivec_bin than looping over columns of genoMat
            % for yy = 1:length(ivec_bin)
            %     % Mean squared error - where mean is calculated using df
            %     tmp_MSE = sum((tmpY(:,yy) - (res_genoMat .* tmpCoeff(:,yy)')).^2)/df;
            % 
            %     % Standard error
            %     beta_se(:, ivec_bin(yy)) = sqrt(Bi .* tmp_MSE);
            % end

            % Update bin counter
            binLoc = binLoc + 1;
        end
    end
else
    % Use basis function to create new genotype matrix
    tmpGenotype = zeros(numObs, numBasisFunc, numSNPs, SingleOrDouble);
    for basis   = 1:numBasisFunc
        tmpGenotype(:, basis, :) = genoMat .* bfSNP(:, basis);
    end

    % Reorder tmpGenotype such that for each SNP, the SNP*basisFucntion
    % columns are next to each other - possible that this might create
    % a very large matrix and potential RAM bottleneck?
    genoMat = reshape(tmpGenotype, numObs, [], 1);

    % Get rid of tmpGenotype - might slightly increase execution time
    clear tmpGenotype

    % Figure out locations in genoMat that need to be handled
    % together; essentially, these are X variables that should be
    % considered together
    beginLocs  = 1:numBasisFunc:size(genoMat,2);
    endLocs    = numBasisFunc:numBasisFunc:size(genoMat,2);

    if OLSflag
        % In case of OLS, there are no bins
        numY = size(ymat_res, 2);

        % Residualize genoMat
        % genoMat can be overwritten
        if regressX
            genoMat = FEMA_OLSResiduals(genoMat, Xvars);
        end

        % Loop over every SNP
        for snps = 1:numSNPs
            % Prepare matrix of X variables
            tmpX = genoMat(:, beginLocs(snps):endLocs(snps));

            % First compute X' * X
            term1 = tmpX' * tmpX;

            % Inverse of term1
            I      = eye(size(term1));
            iterm1 = term1 \ I;
            msg    = lastwarn;
            if ~isempty(msg)
                if useLSQ
                    iterm1 = lsqminnorm(term1, I);
                else
                    iterm1 = pinv(term1);
                end
                msg = ''; %#ok<NASGU>
                lastwarn('');
            end
            iterm1 = nearestSPD(iterm1);

            % Calculate beta coefficient
            tmpCoeff                          = iterm1 * (tmpX' * ymat_res);
            beta_hat(snps, 1:numBasisFunc, :) = tmpCoeff;

            % Mean squared error - where mean is calculated using df
            tmp_MSE = sum((ymat_res - (tmpX * tmpCoeff)).^2)/df;

            % Coefficient covariance
            if doCoeffCovar
                for ii = 1:size(ymat_res,2)
                    coeffCovar(snps, 1:numBasisFunc, 1:numBasisFunc, ii) = iterm1 * tmp_MSE(ii); %#ok<AGROW>
                end
            end

            % Calculate standard error
            beta_se(snps, 1:numBasisFunc, :) = sqrt(diag(iterm1) * tmp_MSE);

            % Do Wald test
            for ii = 1:numY
                coeffCovar = iterm1 * (tmp_MSE(ii) * eye(numBasisFunc));
                [Wald_F(snps, ii, :), Wald_p(snps, ii, :)] =  ...
                FEMA_WaldTest(L, tmpCoeff(:, ii), coeffCovar, hypValue, doF, numObs);
            end
        end
    else
        for sig2bini = unique(binvec(isfinite(binvec)), 'stable')
            ivec_bin = find(binvec == sig2bini);

            % Ws term for this bin
            currWsTerm = allWsTerms{binLoc};

            % For this bin, the y value - faster than repeated indexing
            tmpY = ymat_res(:, ivec_bin);

            % Residualize genoMat
            if regressX
                res_genoMat = FEMA_GLSResiduals(genoMat, Xvars, [], [], currWsTerm, [], useLSQ, lowRank);
            else
                res_genoMat = genoMat;
            end

            % Loop over every SNP
            for snps = 1:numSNPs

                % Prepare matrix of X variables
                tmpX = res_genoMat(:, beginLocs(snps):endLocs(snps));

                % Prepare XtW
                XtW = tmpX' * currWsTerm;

                % XtW is already transposed
                term1 = XtW * tmpX;

                % Calculate inverse of term1
                I       = eye(size(term1));
                iterm1  = term1 \ I;
                msg    = lastwarn;
                if ~isempty(msg)
                    if useLSQ
                        iterm1 = lsqminnorm(term1, I);
                    else
                        iterm1 = pinv(term1);
                    end
                    msg = ''; %#ok<NASGU>
                    lastwarn('');
                end
                iterm1 = nearestSPD(iterm1);

                % Calculate beta coefficient
                tmpCoeff                                 = iterm1 * (XtW * tmpY);
                beta_hat(snps, 1:numBasisFunc, ivec_bin) = tmpCoeff;

                % Mean squared error - where mean is calculated using df
                tmp_MSE = sum((tmpY - (tmpX * tmpCoeff)).^2)/df;

                % Calculate standard error
                beta_se(snps,  1:numBasisFunc, ivec_bin) = sqrt(diag(iterm1) * tmp_MSE);

                % Do Wald test
                for ii = 1:length(ivec_bin)
                    temp_coeffCovar = iterm1 * (tmp_MSE(ii) * eye(numBasisFunc));
                    [Wald_F(snps, ivec_bin(ii), :), Wald_p(snps, ivec_bin(ii), :)] =  ...
                    FEMA_WaldTest(L, tmpCoeff(:, ii), temp_coeffCovar, hypValue, doF, numObs);

                    if doCoeffCovar
                        coeffCovar(snps, 1:numBasisFunc, 1:numBasisFunc, ivec_bin(ii)) = temp_coeffCovar; %#ok<AGROW>
                    end
                end
            end

            % Update bin counter
            binLoc = binLoc + 1;
        end
    end
end

% Reset the status of warnings
warning(statusSingular);
warning(statusNearSing);

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

timeStruct.tEstimation = toc(tEstimation);
timeStruct.Overall     = toc(tOver);

%% Handle the situation where a complex valued standard error was generated
% Likely because of the random effects variance-covariance matrices not
% being semi-positive definite; this seemingly (and rarely) leads to a
% negative valued coefficient covariance inv(XtWX) term, for which the
% square root is a complex number
errFlag = imag(beta_se) ~= 0;
if numBasisFunc > 1
    errFlag = squeeze(any(errFlag, 2));
end

% Save variables, if necessary
if writeResults

    % Save name - make sure .mat is not already part of outName
    saveName = fullfile(outDir, [strrep(outName, '.mat', ''), '.mat']);

    % Variables to save
    toSave = {'beta_hat', 'beta_se', 'tStats', 'logpValues', 'df', 'Wald_F', 'Wald_p', 'coeffCovar', 'timeStruct', 'genStruct', 'errFlag'};

    % Get an estimate of size of the variables
    tmpInfo = whos;
    if sum([tmpInfo(ismember({tmpInfo(:).name}', toSave)).bytes]) > 2^31
        save(saveName, 'beta_hat', 'beta_se', 'tStats', 'logpValues', 'df', 'Wald_F', 'Wald_p', 'coeffCovar', 'timeStruct', 'genStruct', 'errFlag', '-v7.3');
    else
        save(saveName, 'beta_hat', 'beta_se', 'tStats', 'logpValues', 'df', 'Wald_F', 'Wald_p', 'coeffCovar', 'timeStruct', 'genStruct', 'errFlag');
    end
end