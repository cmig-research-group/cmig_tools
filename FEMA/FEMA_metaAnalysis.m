function [meta_beta_hat, meta_beta_SE, meta_tStats, meta_logpValues, meta_coeffCovar] = ...
          FEMA_metaAnalysis(estList, metaType)
% Function to perform meta analysis, given a list of FEMA estimate file
% names and the type of meta analysis to perform
%% Inputs:
% estList:          cell string     full path to FEMA estimate files from
%                                   different sites
% metaType:         character       should be one of the following:
%                                       * fixed
%                                       * random
% 
%% Outputs:
% meta_beta_hat:    matrix          meta-analysis beta coefficients
% meta_beta_SE:     matrix          meta-analysis standard errors
% meta_tStats:      matrix          meta-analysis t/z statistics
% meta_logpValues:  matrix          meta-analysis signed -log10(p) values
% meta_coeffCovar:  matrix          meta-analysis coefficient covariance matrix
%
%% Notes:
% Assumes that the FEMA estimate files from each site contain two estimates:
%   * beta_hat:     p x v,      where p is the number of coefficients and v
%                               is the number of outcome variables
%   * coeffCovar:   p x p x v,  where p is the number of coefficients and v
%                               is the number of outcome variables
% 
% Additionally, assumes that the coefficients are lined up in the same
% order across sites 

%% Check inputs
% Check estList
if ~exist('estList', 'var') || isempty(estList)
    error('Please provide a list of estimate files that should be meta-analysed');
else
    % Check if all the files exist prior to starting meta-analysis
    if ~iscellstr(estList) || ~isstring(estList)
        estList = cellstr(estList);
    end
    estList = reshape(estList, length(estList), 1);
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
    end
end

% Check metaType
if ~exist('metaType', 'var') || isempty(metaType)
    metaType = 'fixed';
else
    metaType = lower(metaType);
    if ~ismember(metaType, {'fixed', 'random'})
        error(['Unknown meta-analysis type specified: ', metaType, '; should be either fixed or random']);
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

%% Get the required parameters: beta hat and coeffCovar
siteEstimates = repmat(struct('beta_hat', [], 'coeffCovar', []), numFiles, 1);
for files = 1:numFiles
    siteEstimates(files) = load(estList{files}, 'beta_hat', 'coeffCovar');
end

%% Sanity check:
% Get size of estimates across all sites
all_sizes_beta_hat   = cell2mat(arrayfun(@(x) size(x.beta_hat),        siteEstimates, 'UniformOutput', false));
all_sizes_coeffCovar = cell2mat(arrayfun(@(x) size(x.coeffCovar, 1:3), siteEstimates, 'UniformOutput', false));

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

if numYVars > 1
    if ~all(all_sizes_coeffCovar(1,3) == numYVars)
        error(['Expected the third dimension of coefficient covariance to be: ', num2str(numYVars)]);
    end
end

%% Do meta-analysis
switch metaType
    case 'fixed'
        % Inverse of the coefficient covariance matrix from each site
        inverter = eye(numCoefficients);

        if useMLDivide
            for files = 1:numFiles
                siteEstimates(files).W = pagemldivide(siteEstimates(files).coeffCovar, inverter);
            end
        else
            for files = 1:numFiles
                % Initialize inverse of coefficientCovariance for this site
                siteEstimates(files).W = zeros(numCoefficients, numCoefficients, numYVars, ...
                                               class(siteEstimates(files).coeffCovar));
                for yy = 1:numYVars
                    siteEstimates(files).W(:,:,yy) = siteEstimates(files).coeffCovar(:,:,yy) \ inverter;
                end
            end
        end

        % Calculate the overall coefficient covariance matrix: inverse of
        % the sum of the site specific inverses
        meta_W = sum(cat(4,siteEstimates(:).W),4);

        if useMLDivide
            meta_coeffCovar = pagemldivide(meta_W, inverter);
        else
            % Initialize coefficientCovariance
            meta_coeffCovar = zeros(numCoefficients, numCoefficients, numYVars, class(meta_W));
            for yy = 1:numYVars
                meta_coeffCovar(:,:,yy) = meta_W(:,:,yy) \ inverter;
            end
        end

        % Compute the site-specific inverse coefficient covariance weighted
        % beta coefficient
        if useMTimes
            for files = 1:numFiles
                % Reshaping beta from p * v to p * 1 * v and then
                % multiplying with coefficient covariance p * p * v
                siteEstimates(files).weightedBeta = squeeze(pagemtimes(siteEstimates(files).W, ...
                                                    reshape(siteEstimates(files).beta_hat, numCoefficients, 1, numYVars)));
            end
        else
            for files = 1:numFiles
                % Initialize
                siteEstimates(files).weightedBeta = zeros(numCoefficients, numYVars, class(siteEstimates(files).beta_hat));
                for yy = 1:numYVars
                    siteEstimates(files).weightedBeta(:,yy) = siteEstimates(files).W(:,:,yy) * siteEstimates(files).beta_hat(:,yy);
                end
            end
        end

        % Pooled coefficients
        % Sum across the sites and then multiply with meta-analysis
        % coefficient covariance
        if useMTimes
            meta_beta_hat = squeeze(pagemtimes(meta_coeffCovar, ...
                                    reshape(sum(cat(3,siteEstimates(:).weightedBeta),3), numCoefficients, 1, numYVars)));
        else
            tmp_beta_hat  = sum(cat(3,siteEstimates(:).weightedBeta),3);
            meta_beta_hat = zeros(numCoefficients, numYVars, class(tmp_beta_hat));
            for yy = 1:numYVars
                meta_beta_hat(:,yy) = meta_coeffCovar(:,:,yy) * tmp_beta_hat(:,yy);
            end
        end
        
        % Pooled standard error
        meta_beta_SE = sqrt(cell2mat(arrayfun(@(k) diag(meta_coeffCovar(:,:,k)), ...
                                                   1:numYVars, 'UniformOutput', false)));

        % Meta-analysis t/z statistics
        meta_tStats = meta_beta_hat./meta_beta_SE;

        % Meta-analysis signed -log10(p) values
        meta_logpValues = -sign(z) .* FEMA_get_logp_normcdf(meta_tStats, 'both', true);

    case 'random'
        disp('Not yet supported; forthcoming');
        [meta_beta_hat, meta_beta_SE, meta_tStats, meta_logpValues, meta_coeffCovar] = deal([]);
end