function [phenotype,        effSize_SNPs, locCausal,        XVars,                      ...
          effSize_FFX,      effSize_RFX,  propTotVariance,  propTotNames,               ...
          propResVar_RFX,   propResNames] =                                             ...
          FEMA_simulateGWAS(numPhenotypes, genomat,     bimFileDat,     causalSNPs,     ...
                            effSize_SNPs,  V_A,         XVars,          effSize_FFX,    ...
                            V_FFX,        clusterinfo,  RandomEffects, effSize_RFX, numObs)
% Function to generate synthetic phenotype given a set of causal SNPs 
% Can additionally account for effects of covariates and random effects
%% Inputs:
% numPhenotypes:    [1 x 1]     integer indicating how many phenotypes
%                               should be simulated (v) (default: 1)
%
% genomat:          [n x m]     genotype matrix for n subjects and m SNPs;
%                               if empty, additive effect is not simulated
%
% bimFileDat:       [m x 1]     cell type having list of m rsIDs
%
% causalSNPs:       [1 x v]     cell type having v lists of causal SNPs;
%                               within each cell, there should be a list of 
%                               [k x 1] causal SNPs; each of these SNPs
%                               should be a subset of m; in case of only
%                               one phenotype, it can be a [k x 1] cell; in
%                               case of multiple phenotypes and only one
%                               cell/[k x 1] is provided, the same set of
%                               causal SNPs are used for all the phenotypes
%
% effSize_SNPs:     [1 x v]     cell type having v lists of vectors of 
%                               [k x 1] effect sizes for each k causal SNP;
%                               in case of only one phenotype, can be a 
%                               [k x 1] vector; in case of multiple
%                               phenotypes and only one cell/[k x 1], the
%                               same set of effect sizes are used for all
%                               the phenotypes; if not provided, generated
%                               from standard normal distribution
% 
% V_A:              [1 x v]     proportion of residual variance explained
%                               by additive genetic effect for each
%                               phenotype; if a scalar is provided, the
%                               same V_A is used for all phenotypes
%                               (default: 0.6)
% 
% XVars:            [n x p]     matrix of confounding covariates; can be
%                               left empty if no confounding covariates;
%                               alternatively, if it is a scaler p, then p
%                               X variables are created from standard normal
%                               distribution (the first X variable is the
%                               intercept term in this case)
%
% effSize_FFX:      [p x v]     vector of effect sizes for p X variables
%                               and v phenotypes; if only one column of p
%                               values are provided, the same are used for
%                               all the phenotypes; if not provided,
%                               generated from standard normal distribution
%
% V_FFX:            [1 x v]     variance of the fixed effects; note that
%                               this value is in raw units and not
%                               proportion of total or residual variance;
%                               if scaler is provided, the same V_FFX is
%                               used for all the phenotypes (default: 1)
%
% clusterinfo:      [1 x f]     cell type of f families/clusters, returned
%                               as an output from FEMA_parse_family
% 
% RandomEffects:    [1 x r]     cell type of names of the random effects
%                               default: F; NB: do not pass 'A' and 'E' as
%                               random effects or in the effSize_RFX -
%                               effect size for additive component should 
%                               be specified via the V_A parameter and V_E
%                               is automatically calculated
%
% effSize_RFX:      [r x v]     matrix of normalized effect sizes for r 
%                               random effects and v phenotypes; if only
%                               one column of r is provided, the same
%                               effect sizes are used for all the
%                               phenotypes
%
% numObs:           [1 x 1]     used in cases where number of observations
%                               to simulate cannot be determined
%                               (default: 1000)
%
%% Outputs:
% phenotype:        [n x v]     vector of simulated v phenotypes
% 
% effSize_SNPs:     [v x 1]     cell type with each row containing the
%                               SNP effect sizes for that phenotype; each
%                               cell contains [k x 1] vector of effect
%                               sizes for k causal SNP; k does not need to
%                               be consistent across phenotypes
% 
% locCausal:        [v x 1]     cell type with each row containing the SNP
%                               location within the bimFileDat for that
%                               phenotype; each cell contains [k x 1]
%                               vector of locations where the causal
%                               locations are found in bimFileDat
%
% XVars:            [n x p]     matrix of confounding variables
%
% effSize_FFX:      [p x v]     matrix of effect sizes for p X variables,
%                               for v phenotypes
%
% effSize_RFX:      [r x v]     matrix of normalized effect sizes for r
%                               random effects and v phenotypes
% 
% propTotVariance:  [6 x v]     matrix indicating the proportion of total
%                               variance corresponding to different parts
%                               of the v phenotype
%
% propTotNames:     [6 x 1]     vector indicating the names of the
%                               different components of propTotVariance
%
% propResVar_RFX:   [r+2 x v]   matrix indicating the porportion of
%                               residual variance that is being accounted
%                               by different random effects, with the first
%                               term being the additive effect and the last
%                               term being the error/noise effect for v
%                               phenotypes
%
% propResNames:     [r+2 x 1]   vector indicating the names of the
%                               different components of propResVar_RFX
%
%% Notes:
% We assume that the residual variance in the data is 1; the total variance
% in the data is 1 + variance of the fixed effects (V_FFX), which is input
% by the user
%
% Given, V_A, V_FFX, and V_RFX, V_E is calculated as:
% totVarFrac = V_A + sum(effSize_RFX);
% V_E        = 1 - totVarFrac;
%
% Case 1: no FFX, no A, no RFX: V_Total = V_E                       [V_E = 1]
% Case 2: only FFX:             V_Total = V_FFX + V_E               [V_E = 1]
% Case 3: only A:               V_Total = V_A + V_E                 [V_A + V_E = 1]
% Case 4: only RFX:             V_Total = V_RFX + V_E               [V_RFX + V_E = 1]
% Case 5: A + FFX:              V_Total = V_FFX + V_A + V_E         [V_A + V_E = 1]
% Case 6: A + RFX:              V_Total = V_A + V_RFX + V_E         [V_A + V_RFX + V_E = 1]
% Case 7: FFX + RFX:            V_Total = V_FFX + V_RFX + V_E       [V_RFX + V_E = 1]
% Case 8: A + FFX + RFX:        V_Total = V_FFX + V_A + V_RFX + V_E [V_A + V_RFX + V_E = 1]
%
% Calculation of modified effect sizes for additive and fixed effects:
% Since the y variable is being standardized, the original effect sizes
% will also change. The slopes get scaled by the standard deviation of y
% and the effect size for the intercept is (oriIntercept - mean(y))/std(y).
% In case of additive effects, only scaling by standard deviation is applied
%
%% TO DO
% Consider the case where the X variables are derived from the genotype -
% such as the first few principal components
%
% Add support for categorical X variables provided as dummy coded X
% variables - to exlucde when standardizing
%
% Add support for simulating effect of maternal and paternal genotype, over
% and above the offspring's genotype (additive effects on SNPs)
% 
% There is a small amount of covariance between noise and other RFX and
% therefore the overall reisdual variance can be slightly off - possible
% fix is to orthogonalize the random effects and then mean center +
% standardize the random effects
%
% Since the variables are being standardized at multiple points, the
% resulting data does not have the same effect sizes as effSize_FFX and
% effSize_SNP; need to figure out scaling of effect sizes which can be
% output for checking parameter recovery

%% Check inputs and assign defaults
% First figure out number of phenotypes
if ~exist('numPhenotypes', 'var') || isempty(numPhenotypes)
    numPhenotypes = 1;
    multiPhen = false;
else
    if numPhenotypes > 1
        multiPhen = true;
    else
        multiPhen = false;
    end
end

% Check genomat and other variables related to additive genetic effect
if ~exist('genomat', 'var') || isempty(genomat)
    doAddGen  = false;
    V_A       = 0;
    locCausal = [];
    disp('No genomat variables provided; skipping simulating genetic effects');
else
    doAddGen          = true;
    [numObs, numSNPs] = size(genomat);

    % Check bimFileDat
    if ~exist('bimFileDat', 'var') || isempty(bimFileDat)
        error('Please provide a list of rsIDs');
    else
        if length(bimFileDat) ~= numSNPs
            error(['Number of rsIDs in bim file (', num2str(numSNPs), ...
                   ') do not match number of SNPs in genotype matrix (', num2str(numSNPs), ')']);
        end
    end
   
    % Check list of causal variants
    if ~exist('causalSNPs', 'var') || isempty(causalSNPs)
        error('Please provide a list of causal variants');
    else
        % Check if input is [k x 1] list
        if ~iscell(causalSNPs{1})
            causalSNPs = {causalSNPs};
        end
        
        % Check if the same causal SNPs are to be used for all phenotypes
        if size(causalSNPs, 2) == 1 && multiPhen
            causalSNPs = repmat(causalSNPs, 1, numPhenotypes);
        end
        
        % Ensure v number of causal SNP sets are present
        if size(causalSNPs, 2) ~= numPhenotypes
            error(['Mismatch between list of causal SNPs (', num2str(size(causalSNPs, 2)), ...
                   ') and number of phenotypes (', num2str(numPhenotypes), ')']);
        end
    end
    
    % Make sure that all causal SNPs are present in bimFileDat
    [tmpCausal, locCausal] = cellfun(@(x) ismember(x, bimFileDat), causalSNPs, 'UniformOutput', false);
    numCausalSNPs          = cellfun(@length, locCausal);
    if sum(cellfun(@sum, tmpCausal) == numCausalSNPs) ~= numPhenotypes
        error('One or more causal SNPs not found in bimFileDat');
    end
    
    % Check effect size of causal variants
    if ~exist('effSize_SNPs', 'var') || isempty(effSize_SNPs)
        disp('Effect size for SNPs not provided; sampling from standard normal distribution');
        
        % effSize_SNPs is [1 x v] with contents being [k x 1]
        effSize_SNPs = cell(1, numPhenotypes);
        for phen     = 1:numPhenotypes
            effSize_SNPs{1,phen} = randn(numCausalSNPs(phen),1);
        end
    else
        % Check if input is a vector of effect sizes
        if ~iscell(effSize_SNPs)
            effSize_SNPs = {effSize_SNPs};
        end
        
        % Check if the same set of effect sizes are to be used
        if size(effSize_SNPs, 2) == 1 && multiPhen
            effSize_SNPs = repmat(effSize_SNPs, 1, numPhenotypes);
        end
        
        % Ensure v number of causal effect size sets are present
        if size(effSize_SNPs, 2) ~= numPhenotypes
            error(['Mismatch between set of causal SNPs effect sizes (', num2str(size(effSize_SNPs, 2)), ...
                   ') and number of phenotypes (', num2str(numPhenotypes), ')']);
        end
        
        % Ensure that number of causal effect sizes across sets are
        % consistent with the list of causal SNPs
        numCausalEffSizes = cellfun(@length, effSize_SNPs);
        if sum(numCausalEffSizes == numCausalSNPs) ~= numPhenotypes
            error('Mismatch between number of causal SNPs and number of causal SNP effect sizes');
        end
        
        % Ensure all sets are in k x 1 
        for phen = 1:numPhenotypes
            effSize_SNPs{1,phen} = reshape(effSize_SNPs{1,phen}, numCausalSNPs(phen), 1);
        end
    end
    
    % Check variance fraction of additive effect
    if ~exist('V_A', 'var') || isempty(V_A)
        warning('Fraction of residual variance explained by additive genetic effects not provided; setting to 0.6');
        V_A = repmat(0.6, 1, numPhenotypes);
    else
        if isscalar(V_A)
            V_A = repmat(V_A, 1, numPhenotypes);
        else
            if length(V_A) ~= numPhenotypes
                error(['Mismatch between number of elements in V_A (', num2str(length(V_A)), ...
                       ') and number of phenotypes (', num2str(numPhenotypes), ')']);
            end
        end
    end
end

% Check XVars and other variables related to fixed effects
if ~exist('XVars', 'var') || isempty(XVars)
    doFFX       = false;
    numFFX      = 0;
    XVars       = [];
    effSize_FFX = [];
    V_FFX       = 0;
    disp('No X variables provided; skipping simulating fixed effects');
else
    doFFX = true;
    if isscalar(XVars)
        numFFX = XVars;
    else
        if ~exist('numObs', 'var')
            numObs = size(XVars, 1);
        end
        [tmpLen, numFFX] = size(XVars);
        if tmpLen ~= numObs
            error(['Number of observations in X variables (', num2str(tmpLen), ...
                  ') do not match the number of observations in genotype matrix (', num2str(numObs), ')']);
        end
    end
    
    % Check for effect sizes for confounding fixed effects
    if ~exist('effSize_FFX', 'var') || isempty(effSize_FFX)
        disp('Effect size for X variables not provided; sampling from standard normal distribution');
        effSize_FFX = randn(numFFX, numPhenotypes);
    else
        % Check if the same effect sizes are to be used across phenotypes;
        % this part relies on the fact that user has correctly specified
        % p x v or p x 1 effect sizes for p X variables
        [tmpFFXnum, tmpFFXphen] = size(effSize_FFX);
        if tmpFFXphen == 1
            effSize_FFX = repmat(effSize_FFX, 1, numPhenotypes);
        end
        if tmpFFXnum ~= numFFX
            error('Please provide effect size for each X variable');
        end
    end
    
    % Check variance fraction for fixed effects
    if ~exist('V_FFX', 'var') || isempty(V_FFX)
        V_FFX = ones(1, numPhenotypes);
        warning('Variance of fixed effects not specified; setting to 1');
    else
        if isscalar(V_FFX)
            V_FFX = repmat(V_FFX, 1, numPhenotypes);
        else
            if length(V_FFX) ~= numPhenotypes
                error(['Mismatch between number of elements in V_FFX (', num2str(length(V_FFX)), ...
                   ') and number of phenotypes (', num2str(numPhenotypes), ')']);
            end
        end
    end
end

% Check clusterinfo and other variables related to random effects
if ~exist('clusterinfo', 'var') || isempty(clusterinfo)
    doRFX         = false;
    effSize_RFX   = zeros(1, numPhenotypes);
    RandomEffects = {''};
    numRFX        = 1; % Setting to 1 for initialization purposes
    disp('No clusterinfo variable provided; skipping simulating random effects');
else
    doRFX  = true;
    numFam = length(clusterinfo);
    temp   = [clusterinfo{:}];
    temp   = length(horzcat(temp(:).jvec_fam)');
    
    % Check if numObs are the same
    if exist('numObs', 'var')
        if temp ~= numObs
            error('Number of observations in random effects does not match number of observations in genomat');
        end
    else
        numObs = temp;
    end
    
    % Check remaining parameters relevant to random effects
    % First check RandomEffects
    if ~exist('RandomEffects', 'var') || isempty(RandomEffects)
        RandomEffects = {'F'};
    else
        % Make sure that the additive effect and error are not included
        loc_additive = ismember(RandomEffects, 'A');
        if sum(loc_additive) ~= 0
            warning(['Found additive effect specified as a random effect; ', ...
                    'ignoring this effect as well as associated effSize_RFX']);
            RandomEffects(loc_additive) = [];
        end
        loc_error = ismember(RandomEffects, 'E');
        if sum(loc_error) ~= 0
            warning(['Found error effect specified as a random effect; ', ...
                    'ignoring this effect as well as associated effSize_RFX']);
            RandomEffects(loc_error) = [];
        end
    end
    numRFX = length(RandomEffects);

    % Check effect sizes for random effects
    if ~exist('effSize_RFX', 'var') || isempty(effSize_RFX)
        error('Please provide effect sizes for random effects');
    else
        % Make sure that A and E are not part of the effect sizes
        effSize_RFX(loc_additive) = [];
        effSize_RFX(loc_error)    = [];
        
        % Check if same effect sizes are to be used for all phenotypes
        % this part relies on the fact that user has correctly specified
        % r x v or r x 1 effect sizes for r RFX variables
        if ~isempty(effSize_RFX)
            [tmpNumRFX, tmpNumPhenRFX] = size(effSize_RFX);
            if tmpNumRFX  ~= numRFX
                error(['Number of effect sizes for random effects (', num2str(tmpNumRFX), ...
                      ') does not match the number of random effects (', num2str(numRFX), ')']);
            end
            if tmpNumPhenRFX == 1
                effSize_RFX = repmat(effSize_RFX, 1, numPhenotypes);
            end
        else
            % Appears to be the case of noise only
            warning('Appears to be the case of RFX as noise only');
            effSize_RFX = zeros(1, numPhenotypes);
            numRFX      = 0;
        end
    end
end

% Check if number of observations was defined
if ~exist('numObs', 'var')
    disp('numObs not provided; simulating 1000 subjects');
    numObs = 1000;
end

% Figure out remaining variance in the data - attributed to error term
totVarFrac = V_A + sum(effSize_RFX,1);
V_E        = 1 - totVarFrac;
if any(V_E < 0)
    error('Ensure that V_A and V_RFX sum to 1 or less');
end

% Add noise/error as randomEffect, if random effects are being simulated
if doRFX
    % Handle the case of noise as RFX
    if isempty(RandomEffects)
        RandomEffects = {'E'};
        effSize_RFX   = V_E;
        numRFX        = 1;
    else
        RandomEffects = [RandomEffects, 'E'];
        effSize_RFX   = [effSize_RFX; V_E];
        numRFX        = numRFX + 1;
    end        
end

%% Main module
% Initialization
FFX_y    = zeros(numObs, numPhenotypes);
initPhen = zeros(numObs, numPhenotypes);

% Additive genetic effect, scaled to have unit variance
if doAddGen
    for phen = 1:numPhenotypes
        % This line can become a bottle neck if the number of phenotypes is
        % fairly large
        initPhen(:,phen) = double(genomat(:, locCausal{phen})) * effSize_SNPs{phen};
        
        % Scale the phenotype with the standard deviation
        stdPhen          = std(initPhen(:,phen));
        initPhen(:,phen) = initPhen(:,phen)./stdPhen;
    
        % Scaling the effect size by standard deviation - no intercepts
        effSize_SNPs{phen} = effSize_SNPs{phen}./stdPhen;
    end
end

% Fixed effects
if doFFX
    % Simulate XVariables, if necessary
    if isscalar(XVars)
        disp(['X variables not provided; drawing ', num2str(numFFX), ...
              ' X variables from standard normal distribution;',     ...
              ' the first X variable is the intercept']);
        XVars = [ones(numObs, 1), randn(numObs, numFFX-1)];
    end
    
    % Simulate data
    FFX_y = XVars * effSize_FFX;
    
    % Find out if there is an intercept
    locIntercept = std(XVars) == 0;
    
    % Leave effect size and FFX_y unchanged if only intercept model
    if ~(numFFX == 1 && sum(locIntercept) == 1)
        if sum(locIntercept) == 0
            % Only scale by standard deviation    
            std_FFX = std(FFX_y);
    
            % Standardize FFX_y: only scaling by standard deviation
            FFX_y = FFX_y ./ std_FFX;
    
            % Update effect sizes
            effSize_FFX = effSize_FFX ./ std_FFX;
        else
            % Mean center and scale by standard deviation    
            std_FFX  = std(FFX_y);
            mean_FFX = mean(FFX_y);
        
            % Standardize FFX_y: both mean centering and standardizing
            FFX_y = (FFX_y - mean_FFX)./ std_FFX;
    
            % Update effect sizes
            effSize_FFX(~locIntercept, :) = effSize_FFX(~locIntercept,:) ./ std_FFX;
            effSize_FFX(locIntercept,  :) = (effSize_FFX(locIntercept,:) - mean_FFX) ./ std_FFX;
        end
    end
end

% Loop over every family and add random effect term; if clusterinfo is
% present, then the noise ias added per cluster; otherwise noise is added
% to the overall data
% Adapted from FEMA_synthesize
if doRFX
    % At this stage, tmpRFX_y is 3D: Obs x RFX x numPhenotypes
    tmpRFX_y = zeros(numObs, numRFX, numPhenotypes);
    for ri   = 1:numRFX
        for fi  = 1:numFam
            tmp = 0;
            n   = length(clusterinfo{fi}.jvec_fam);
            tmp = tmp + sqrt(effSize_RFX(ri,:)) .* ...
                        mvnrnd(zeros(n, 1), double(clusterinfo{fi}.(['V_', RandomEffects{ri}])), numPhenotypes)';
        tmpRFX_y(clusterinfo{fi}.jvec_fam, ri, 1:numPhenotypes) = tmp;
        end
    end
    
    % Mean center and standardize RFX, including noise, for each phenotype
    tmpRFX_y = (tmpRFX_y - mean(tmpRFX_y,1))./std(tmpRFX_y,[],1);
    
    % The last column of random effects is noise - take it out as a
    % separate component; everything that remains are other random effects
    noise    = tmpRFX_y(:,end,:);
    tmpRFX_y = tmpRFX_y(:,1:end-1,:);
    
    % Drop the noise row from effSize_RFX and RandomEffects
    RandomEffects(end) = [];
    effSize_RFX(end,:) = [];
    numRFX             = numRFX - 1;
else
   % In this case, noise is added outside of random effects 
   noise = randn(numObs, numPhenotypes);
   noise = (noise - mean(noise))./std(noise);
end

% Put phenotype together
% First calculate piecewise variables by scaling to desired variance
y_initPhen  = sqrt(V_A)   .* initPhen;
y_noise     = sqrt(V_E)   .* squeeze(noise);
y_FFX       = sqrt(V_FFX) .* FFX_y;
if doRFX
    vec     = zeros(numObs, numRFX, numPhenotypes);
    for rfx = 1:numRFX
        vec(:,rfx,:) = sqrt(effSize_RFX(rfx,:)) .* squeeze(tmpRFX_y(:,rfx,:));
    end
    y_RFX = squeeze(sum(vec, 2));
else
    y_RFX = zeros(numObs, numPhenotypes);
    vec   = zeros(numObs, numPhenotypes);
end

% Put phenotype together
phenotype = double(y_FFX + y_initPhen + y_noise + y_RFX);

% Scale the effect sizes for fixed and additive effects
if doFFX
    effSize_FFX = sqrt(V_FFX) .* effSize_FFX;
end
if doAddGen
    for phen = 1:numPhenotypes
        effSize_SNPs{phen} = sqrt(V_A(phen)) .* effSize_SNPs{phen};
    end
end

% Recompute actual variances in components of the data
V_FFX = var(y_FFX);
V_A   = var(y_initPhen);
V_RFX = reshape(squeeze(var(vec)), numRFX, numPhenotypes);
V_E   = var(y_noise);

%% Summarize different variance components
% Total variance in the data
totalVariance = var(phenotype);

% Calculate residual variance after accounting for fixed effects
% This residual variance will be explained by additive, random, and noise
% effects; note that this includes the noise term even if it was modeled
% outside of random effects and should sum to 1
residualVariance = totalVariance - V_FFX;

% Calculate proportion of variances explained by different components, as a
% function of the total variance
propVar_fixed    = V_FFX./totalVariance;
propVar_additive = V_A./totalVariance;
propVar_random   = sum(V_RFX,1)./totalVariance;
propVar_noise    = V_E./totalVariance;

% Put these values together and return as a vector
propTotVariance = ([totalVariance;   residualVariance; propVar_fixed; ...
                   propVar_additive; propVar_random;   propVar_noise]);
propTotNames    = {'TotalVariance'; 'ResidualVariance'; 'Fixed'; ...
                   'Additive';      'RandomEffects';    'Noise'};

% Calculate proportion of residual variance that is accounted by components
% of additive, random, and noise effects: this will serve as the ground 
% truth for random effects estimates from FEMA_fit
propResVar_RFX = [V_A./residualVariance; V_RFX./residualVariance; V_E./residualVariance];
propResNames   = [{'AdditiveEffect'}; reshape(RandomEffects, numel(RandomEffects), 1); {'E'}];