function [meta_beta_hat, meta_beta_SE, meta_tStats, meta_logpValues, meta_coeffCovar, ...
          Cochran_Q, Cochran_Q_logp, I2, tau2] = FEMA_metaAnalysis(estList, metaType, univariate)
% Function to perform meta analysis, given a list of FEMA estimate file
% names and the type of meta analysis to perform
%% Inputs:
% estList:             cell         full path to FEMA estimate files from
%                                   different sites
% 
% metaType:            char         should be one of the following:
%                                       * fixed
%                                       * random
% 
% univariate:          logical      true or false indicating if inverse
%                                   variance weighting ('univariate') or 
%                                   inverse coefficient covariance
%                                   weighting should be used (default:
%                                   true)
% 
%% Outputs:
% meta_beta_hat:        matrix      meta-analysis beta coefficients
% 
% meta_beta_SE:         matrix      meta-analysis standard errors
% 
% meta_tStats:          matrix      meta-analysis t/z statistics
% 
% meta_logpValues:      matrix      meta-analysis signed -log10(p) values
% 
% meta_coeffCovar:      matrix      meta-analysis coefficient covariance
%                                   matrix 
% 
% Cochran_Q:            matrix      Cochran's Q statistic for heterogeneity
% 
% Cochran_Q_logp:       matrix      unsigned -log10(p) values for Cochran's
%                                   Q statistic
% 
% I2:                   matrix      Alternate statistics for heterogeneity
% 
% tau2:                 matrix      random effects variance components
%                                   (only returned when metaType is
%                                   random)
% 
%% Notes:
% Assumes that the FEMA estimate files from each site contain two estimates:
%   * beta_hat:     p x v,      where p is the number of coefficients and v
%                               is the number of outcome variables
% 
%   * coeffCovar:   p x p x v,  where p is the number of coefficients and v
%                               is the number of outcome variables
% 
% Additionally, assumes that the coefficients are lined up in the same
% order across sites 
%
% In univariate case: all estimates are computed using the inverse of the
% variance of each coefficient (i.e., using diagonal entries from site-wise
% coefficient covariance).
% 
% In multivariate case: all estimate are computed using the inverse of the
% site-wise coefficient covariance matrix
% 
% I2 is clamped at zero
%
% All outputs are double precision regardless of how the site estimates
% were stored
%
% meta_logpValues is signed by the direction of the effect; Cochran_Q_logp
% is unsigned
% 
% In univariate case, the coefficient covariance matrix is a diagonal
% matrix
% 
% In the univariate case, for Cochran's Q p value calculation, we use chi2
% distribution with numFiles-1 degrees of freedom. In the multivariate
% case, because all coefficients (for each v) are tested simultaneously, we
% use chi2 distribution with numCoefficients*(numFiles-1) degrees of
% freedom
% 
% For multivariate case, we model the between-site heterogeneity as a
% single scalar per outcome, T = tau2 * I.
%
% Worthwhile to note that Cochran's Q, the associated p value, and the I2
% are computed from the fixed effects weights using inverse of the variance
% or covariance, and not re-computed even if metaType is random
% 
%% References:
% Mägi, R., Morris, A.P. 
% GWAMA: software for genome-wide association meta-analysis. 
% BMC Bioinformatics 11, 288 (2010). 
% https://doi.org/10.1186/1471-2105-11-288
% 
% Ioannidis J.P., Patsopoulos N.A., Evangelou E. 
% Heterogeneity in Meta-Analyses of Genome-Wide Association Investigations. 
% PLOS ONE 2(9): e841 (2007). 
% https://doi.org/10.1371/journal.pone.0000841
% 
% Cochran, W. G.
% The combination of estimates from different experiments.
% Biometrics 10 (1954): 101.
% https://www.jstor.org/stable/3001666
% 
% Huedo-Medina, T. B., Sánchez-Meca J., Marín-Martínez F., Botella J.
% Assessing heterogeneity in meta-analysis: Q statistic or I2 index?
% Psychological methods 11 2 (2006): 193-206.
% https://pubmed.ncbi.nlm.nih.gov/16784338/
% 
% Jackson, D., White, I.R. and Thompson, S.G. 
% Extending DerSimonian and Laird's methodology to perform multivariate random effects meta-analyses. 
% Statist. Med., 29 (2010): 1282-1297 
% https://doi.org/10.1002/sim.3602

%% Check inputs
% Check estList
if ~exist('estList', 'var') || isempty(estList)
    error('Please provide a list of estimate files that should be meta-analysed');
else
    % Check if all the files exist prior to starting meta-analysis
    if ischar(estList)
        estList = cellstr(estList);
    end
    estList = reshape(estList, numel(estList), 1);
    chkList = false(length(estList),1);
    for files = 1:length(estList)
        if exist(estList{files}, 'file')
            chkList(files) = true;
        end
    end
    if sum(chkList) ~= length(estList)
        disp(estList(~chkList));
        error('One or more files are missing');
    else
        numFiles = length(estList);

        if numFiles == 1
            error('Cannot perform meta-analysis with a single estimate file');
        end
    end
end

% Check metaType
if ~exist('metaType', 'var') || isempty(metaType)
    isRFX = false;
else
    metaType = lower(metaType);
    if ~ismember(metaType, {'fixed', 'random'})
        error(['Unknown meta-analysis type specified: ', metaType, ...
               '; should be either fixed or random']);
    else
        if strcmpi(metaType, 'random')
            isRFX = true;
        else
            isRFX = false;
        end
    end
end

% Check univariate flag
if ~exist('univariate', 'var') || isempty(univariate)
    univariate = true;
else
    if ~isscalar(univariate) || ~(islogical(univariate) || (isnumeric(univariate) && ismember(univariate, [0 1])))
        error('univariate flag should be either true or false');
    end
end

%% Check if page functions are available
if exist('pagemldivide', 'builtin')
    useMLDivide = true;
else
    useMLDivide = false;
end

if exist('pagemtimes', 'builtin')
    useMTimes = true;
else
    useMTimes = false;
end

%% Sanity checks
% First, get size of coefficients from all sites
all_sizes_beta_hat   = zeros(numFiles, 2);
all_sizes_coeffCovar = zeros(numFiles, 3);
for files = 1:numFiles
    info  = whos('-file', estList{files}, 'beta_hat', 'coeffCovar');
    if length(info) ~= 2
        error(['Site ', estList{files}, ' is missing beta_hat and/or coeffCovar']);
    else
        tmpNames                      = {info(:).name};
        idx_beta                      = strcmpi(tmpNames, 'beta_hat');
        idx_coeff                     = strcmpi(tmpNames, 'coeffCovar');
        tmp_sz_beta                   = info(idx_beta).size;
        tmp_sz_coeff                  = info(idx_coeff).size;
        all_sizes_beta_hat(files,:)   = [tmp_sz_beta,  ones(1, 2-numel(tmp_sz_beta))];
        all_sizes_coeffCovar(files,:) = [tmp_sz_coeff, ones(1, 3-numel(tmp_sz_coeff))];
    end
end

% Ensure the same number of coefficients from each site
if ~all(all_sizes_beta_hat == all_sizes_beta_hat(1,:), 'all')
    error('Mismatch in number of coefficients between sites');
end

% Ensure that the coefficient covariance sizes are the same (should be)
if ~all(all_sizes_coeffCovar == all_sizes_coeffCovar(1,:), 'all')
    error('Mismatch in size of coefficient covariance between sites');
end

% Ensure that the coefficient covariance sizes match beta coefficients
numCoefficients = all_sizes_beta_hat(1,1);
numYVars        = all_sizes_beta_hat(1,2);

if ~all(all_sizes_coeffCovar(1,1:2) == numCoefficients)
    error(['Expected the first two dimensions of coefficient covariance to be: ', num2str(numCoefficients)]);
end

if ~all(all_sizes_coeffCovar(1,3) == numYVars)
    error(['Expected the third dimension of coefficient covariance to be: ', num2str(numYVars)]);
end

%% Meta-analysis
% Set tau2 to be empty to start with
tau2 = [];

if univariate
    % Initialize
    accum_inv_variance     = zeros(numCoefficients, numYVars);
    accum_invw_BetaHatUniv = zeros(numCoefficients, numYVars);
    sitewise_betaHat       = zeros(numCoefficients, numYVars, numFiles);
    sitewise_invVariance   = zeros(numCoefficients, numYVars, numFiles);

    % Name-sake coefficient covariance
    meta_coeffCovar = zeros(numCoefficients, numCoefficients, numYVars);
    locDiag         = 1:(numCoefficients+1):numCoefficients^2;
    
    for files = 1:numFiles
        % Load parameters from the site
        tmpParameters = load(estList{files}, 'beta_hat', 'coeffCovar');
    
        % Cast site parameters as double precision
        tmpParameters.beta_hat   = double(tmpParameters.beta_hat);
        tmpParameters.coeffCovar = double(tmpParameters.coeffCovar);
    
        % Save a copy of these coefficients
        sitewise_betaHat(:,:,files) = tmpParameters.beta_hat;
    
        % inverse of the coefficient covariance and variance
        % Additionally calculate (inverse variance)^2 for random effects case
        if useMLDivide
            invVariance = 1./cell2mat(arrayfun(@(k) ...
                                      diag(tmpParameters.coeffCovar(:,:,k)), ...
                                           1:numYVars, 'UniformOutput', false));
        else
            invVariance = zeros(numCoefficients, numYVars);
    
            for yy = 1:numYVars
                invVariance(:,yy) = 1./diag(tmpParameters.coeffCovar(:,:,yy));
            end
        end
    
        % Save a copy of the inverse variance
        sitewise_invVariance(:,:,files) = invVariance;
    
        % Weighted coefficients: univariate
        weightedBeta_univariate = invVariance .* tmpParameters.beta_hat;
    
        % Add to accumulators
        accum_inv_variance     = accum_inv_variance     + invVariance;
        accum_invw_BetaHatUniv = accum_invw_BetaHatUniv + weightedBeta_univariate;
    end
    
    % Univariate beta coefficient using inverse of variances
    meta_beta_hat = accum_invw_BetaHatUniv./accum_inv_variance;
    
    % Standard error
    tmp_variance = 1./sum(sitewise_invVariance, 3);
    meta_beta_SE = sqrt(tmp_variance);
    
    % Meta-analysis t/z statistics
    meta_tStats = meta_beta_hat./meta_beta_SE;
    
    % Meta-analysis signed -log10(p) values
    meta_logpValues = -sign(meta_tStats) .* FEMA_get_logp_normcdf(meta_tStats, 'both', true);
    
    % Cochran's Q
    sqDiff    = (meta_beta_hat - sitewise_betaHat).^2;
    Cochran_Q = sum(sitewise_invVariance .* sqDiff, 3);
    
    if any(Cochran_Q < 0, 'all')
        error(['Cochran Q is negative for one or more coefficients; ', ...
                'check site estimates to ensure that diagonal entries of coeffCovar are non-negative']);
    end
    
    % Cochran's Q p value
    Cochran_Q_logp = FEMA_get_logp_chi2cdf(Cochran_Q, numFiles-1, 'upper', true);
    
    % I2
    I2 = (Cochran_Q - (numFiles - 1))./Cochran_Q;
    if any(I2 < 0, 'all')
        warning('I2 resulted in negative values; clamping them to zero');
        I2(I2 < 0) = 0;
    end
    if any(Cochran_Q <= 0, 'all')
        warning('Cochran Q was zero or negative for one or more coefficients; clamping those I2 to zero');
        I2(Cochran_Q <= 0) = 0;
    end

    % Update for random effects case
    if isRFX
        % Square of the inverse term
        accum_inv_variance2 = sum(sitewise_invVariance.^2, 3);

        % Calculate tau per coefficient
        numerator   = Cochran_Q - (numFiles - 1);
        denominator = accum_inv_variance - (accum_inv_variance2./accum_inv_variance);
        tau2        = max(0, numerator./denominator);

        % Inverse of the updated weights
        up_sitewise_invVariance = 1./(tau2 + (1./sitewise_invVariance));

        % Sum of the inverse of the udpated weights
        tmp_sum_up_sitewise_invVariance = sum(up_sitewise_invVariance, 3);

        % Updated beta coefficient (univariate)
        meta_beta_hat = sum(sitewise_betaHat .* up_sitewise_invVariance, 3) ./ ...
                             tmp_sum_up_sitewise_invVariance;

        % Updated beta SE
        tmp_variance = 1./tmp_sum_up_sitewise_invVariance;
        meta_beta_SE = sqrt(tmp_variance);

        % Updated meta-analysis t/z statistics
        meta_tStats = meta_beta_hat./meta_beta_SE;

        % Updated meta-analysis signed -log10(p) values
        meta_logpValues = -sign(meta_tStats) .* FEMA_get_logp_normcdf(meta_tStats, 'both', true);
    end
    
    % Name-sake coefficient covariance
    for yy = 1:numYVars
        tmp_coeffCovar          = zeros(numCoefficients, numCoefficients);
        tmp_coeffCovar(locDiag) = tmp_variance(:,yy);
        meta_coeffCovar(:,:,yy) = tmp_coeffCovar;
    end
    
else
    % Initialize
    accum_inv_coeffCovar   = zeros(numCoefficients, numCoefficients, numYVars);
    accum_invw_BetaHat     = zeros(numCoefficients, numYVars);
    % sitewise_betaHat       = zeros(numCoefficients, numYVars, numFiles);

    inverter = eye(numCoefficients);
    for files = 1:numFiles
        % Load parameters from the site
        tmpParameters = load(estList{files}, 'beta_hat', 'coeffCovar');

        % Cast site parameters as double precision
        tmpParameters.beta_hat   = double(tmpParameters.beta_hat);
        tmpParameters.coeffCovar = double(tmpParameters.coeffCovar);

        % Save a copy of these coefficients
        % sitewise_betaHat(:,:,files) = tmpParameters.beta_hat;

        % Inverse of the coefficient covariance
        if useMLDivide
            invCoeffCovar = pagemldivide(tmpParameters.coeffCovar, inverter);
        else
            invCoeffCovar = zeros(numCoefficients, numCoefficients, numYVars);
            for yy = 1:numYVars
                invCoeffCovar(:,:,yy) = tmpParameters.coeffCovar(:,:,yy) \ inverter;
            end
        end

        % Weighted coefficients: multivariate
        if useMTimes
            % Reshaping beta from p * v to p * 1 * v and then multiplying with
            % coefficient covariance p * p * v
            weightedBeta_multivariate = reshape(pagemtimes(invCoeffCovar, ...
                                        reshape(tmpParameters.beta_hat,   ...
                                                numCoefficients, 1, numYVars)), ...
                                                numCoefficients, numYVars);
        else
            weightedBeta_multivariate = zeros(numCoefficients, numYVars);
            for yy = 1:numYVars
                weightedBeta_multivariate(:,yy) = invCoeffCovar(:,:,yy) * tmpParameters.beta_hat(:,yy);
            end
        end

        % Add to accumulators
        accum_inv_coeffCovar   = accum_inv_coeffCovar   + invCoeffCovar;
        accum_invw_BetaHat     = accum_invw_BetaHat     + weightedBeta_multivariate;
    end

    % Compute statistics
    % Coefficient covariance
    if useMLDivide
        meta_coeffCovar = pagemldivide(accum_inv_coeffCovar, inverter);
    else
        % Initialize coefficientCovariance
        meta_coeffCovar = zeros(numCoefficients, numCoefficients, numYVars);
        for yy = 1:numYVars
            meta_coeffCovar(:,:,yy) = accum_inv_coeffCovar(:,:,yy) \ inverter;
        end
    end

    % Multivariate beta coefficient using coefficient covariance
    if useMTimes
        meta_beta_hat = reshape(pagemtimes(meta_coeffCovar, ...
                        reshape(accum_invw_BetaHat, numCoefficients, 1, numYVars)), ...
                                numCoefficients, numYVars);
    else
        meta_beta_hat = zeros(numCoefficients, numYVars);
        for yy = 1:numYVars
            meta_beta_hat(:,yy) = meta_coeffCovar(:,:,yy) * accum_invw_BetaHat(:,yy);
        end
    end

    % Standard error
    meta_beta_SE = sqrt(cell2mat(arrayfun(@(k) diag(meta_coeffCovar(:,:,k)), ...
                                               1:numYVars, 'UniformOutput', false)));

    % Meta-analysis t/z statistics
    meta_tStats = meta_beta_hat./meta_beta_SE;

    % Meta-analysis signed -log10(p) values
    meta_logpValues = -sign(meta_tStats) .* FEMA_get_logp_normcdf(meta_tStats, 'both', true);

    % Calculate Cochran's Q
    total_inv_coeffCovar = meta_coeffCovar;

    % Difference in the coefficients from the meta analysis estimate
    % coeffDiff = (meta_beta_hat - sitewise_betaHat);
    tmpQ = zeros(numFiles, numYVars);
    if isRFX
        accum_W_totalInv_W = zeros(numCoefficients, numCoefficients, numYVars);
    end

    % Now, go over each site, load coefficient covariance, and calculate Q
    for files = 1:numFiles
        tmpParameters            = load(estList{files}, 'beta_hat', 'coeffCovar');
        tmpParameters.beta_hat   = double(tmpParameters.beta_hat);
        tmpParameters.coeffCovar = double(tmpParameters.coeffCovar);

        % Difference from the pooled estimate
        coeffDiff = tmpParameters.beta_hat - meta_beta_hat;

        if useMLDivide && useMTimes
            invCoeffCovar = pagemldivide(tmpParameters.coeffCovar, inverter);
            tmp           = reshape(coeffDiff, numCoefficients, 1, numYVars);
            tmpQ(files,:) = reshape(pagemtimes(tmp, 'transpose',      ...
                                    pagemtimes(invCoeffCovar, tmp),   ...
                                               'none'), 1, numYVars);
            if isRFX
                accum_W_totalInv_W = accum_W_totalInv_W +                  ...
                                     pagemtimes(invCoeffCovar,             ...
                                     pagemtimes(total_inv_coeffCovar, invCoeffCovar));
            end
        else
            invCoeffCovar = zeros(numCoefficients, numCoefficients, numYVars);
            for yy = 1:numYVars
                invCoeffCovar(:,:,yy) = tmpParameters.coeffCovar(:,:,yy) \ inverter;
                tmpQ(files,yy)        = coeffDiff(:,yy)' * invCoeffCovar(:,:,yy) * coeffDiff(:,yy);
                if isRFX
                    accum_W_totalInv_W(:,:,yy) = accum_W_totalInv_W(:,:,yy) +     ...
                                                 invCoeffCovar(:,:,yy)      *     ...
                                                 total_inv_coeffCovar(:,:,yy) *   ...
                                                 invCoeffCovar(:,:,yy);
                end
            end
        end
    end

    % Generate an error if tmpQ was negative anywhere
    if any(tmpQ < 0, 'all')
        nonSPD_sites = find(any(tmpQ < 0, 2));
        if numel(nonSPD_sites) > 1
            disp(estList(nonSPD_sites));
            error('Multiple sites have non positive definite coefficient covariance matrices');
        else
            error(['Site: ', estList{nonSPD_sites}, ' has non positive definite coefficient covariance matrix']);
        end
    end
    Cochran_Q = sum(tmpQ,1);

    % Cochran's Q p value
    Cochran_Q_logp = FEMA_get_logp_chi2cdf(Cochran_Q, numCoefficients*(numFiles-1), 'upper', true);

    % I2
    I2 = (Cochran_Q - numCoefficients*(numFiles-1))./Cochran_Q;
    if any(I2 < 0, 'all')
        warning('I2 resulted in negative values; clamping them to zero');
        I2(I2 < 0) = 0;
    end
    if any(Cochran_Q <= 0, 'all')
        warning('Cochran Q was zero or negative for one or more outcomes; clamping those I2 to zero');
        I2(Cochran_Q <= 0) = 0;
    end

    % Update for random effects case
    if isRFX
        % Denominator term: C = sum(W_i) - sum(W_i * totalInv * W_i)
        C = accum_inv_coeffCovar - accum_W_totalInv_W;

        % Making the assumption that there is a common scalar term for
        % between site heterogeneity
        trC = zeros(1, numYVars);
        for yy = 1:numYVars
            trC(yy) = trace(C(:,:,yy));
        end

        % Tau squared (simplified)
        tau2 = max(0, (Cochran_Q - (numCoefficients * (numFiles - 1)))./trC);

        % Initialize
        accum_inv_coeffCovar_RE = zeros(numCoefficients, numCoefficients, numYVars);
        accum_invw_BetaHat_RE   = zeros(numCoefficients, numYVars);

        for files = 1:numFiles
            tmpParameters            = load(estList{files}, 'beta_hat', 'coeffCovar');
            tmpParameters.beta_hat   = double(tmpParameters.beta_hat);
            tmpParameters.coeffCovar = double(tmpParameters.coeffCovar);

            % Add tau2 to the diagonal of every outcome's covariance
            coeffCovar_RE = zeros(numCoefficients, numCoefficients, numYVars);
            for yy = 1:numYVars
                coeffCovar_RE(:,:,yy) = tmpParameters.coeffCovar(:,:,yy) + (tau2(yy) * inverter);
            end

            % Inverse of the updated covariance
            if useMLDivide
                invCoeffCovar_RE = pagemldivide(coeffCovar_RE, inverter);
            else
                invCoeffCovar_RE = zeros(numCoefficients, numCoefficients, numYVars);
                for yy = 1:numYVars
                    invCoeffCovar_RE(:,:,yy) = coeffCovar_RE(:,:,yy) \ inverter;
                end
            end

            % Weighted coefficients
            if useMTimes
                weightedBeta_RE = reshape(pagemtimes(invCoeffCovar_RE,        ...
                                  reshape(tmpParameters.beta_hat,             ...
                                          numCoefficients, 1, numYVars)),     ...
                                          numCoefficients, numYVars);
            else
                weightedBeta_RE = zeros(numCoefficients, numYVars);
                for yy = 1:numYVars
                    weightedBeta_RE(:,yy) = invCoeffCovar_RE(:,:,yy) * tmpParameters.beta_hat(:,yy);
                end
            end

            accum_inv_coeffCovar_RE = accum_inv_coeffCovar_RE + invCoeffCovar_RE;
            accum_invw_BetaHat_RE   = accum_invw_BetaHat_RE   + weightedBeta_RE;
        end

        % Updated coefficient covariance
        if useMLDivide
            meta_coeffCovar = pagemldivide(accum_inv_coeffCovar_RE, inverter);
        else
            meta_coeffCovar = zeros(numCoefficients, numCoefficients, numYVars);
            for yy = 1:numYVars
                meta_coeffCovar(:,:,yy) = accum_inv_coeffCovar_RE(:,:,yy) \ inverter;
            end
        end

        % Updated beta coefficients
        if useMTimes
            meta_beta_hat = reshape(pagemtimes(meta_coeffCovar,      ...
                            reshape(accum_invw_BetaHat_RE,           ...
                                    numCoefficients, 1, numYVars)),  ...
                                    numCoefficients, numYVars);
        else
            meta_beta_hat = zeros(numCoefficients, numYVars);
            for yy = 1:numYVars
                meta_beta_hat(:,yy) = meta_coeffCovar(:,:,yy) * accum_invw_BetaHat_RE(:,yy);
            end
        end

        % Updated standard error, t/z statistics, and signed -log10(p) values
        meta_beta_SE    = sqrt(cell2mat(arrayfun(@(k) diag(meta_coeffCovar(:,:,k)), ...
                               1:numYVars, 'UniformOutput', false)));
        meta_tStats     = meta_beta_hat./meta_beta_SE;
        meta_logpValues = -sign(meta_tStats) .* FEMA_get_logp_normcdf(meta_tStats, 'both', true);
    end
end