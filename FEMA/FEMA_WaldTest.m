function [W, p, logp, df, L] = FEMA_WaldTest(L, beta_hat, coeffCovar, hypValue, doF, numObs)
% Function to compute multivariate Wald test given a set of parameters
%% Inputs:
% L:            numeric (vector or matrix) or cell type containing weights 
%               to be used for Wald test; if cell type, each cell can 
%               contain a vector or matrix of weights; should be a [1 x p] 
%               vector or a [l x p] matrix, where p is the number of 
%               covariates and l is the number of linear combinations to be 
%               tested); weights will be zero-padded on the right side
%
% beta_hat:     [p x v] matrix of estimated beta coefficients (such as from 
%               FEMA_fit) where p is the number of covariates and v is the 
%               number of voxels/vertices (y variable)
%
% coeffCovar:   [p x p x v] matrix of estimated beta coefficient covariance
%               (such as from FEMA_fit) where p is the number of covariates
%               and v is the number of voxels/vertices (y variable)
%
% hypValue:     numeric scalar or [1 x k] vector (where k is the number of 
%               cell in L) containing the hypothesised value against which 
%               the null hypothesis will be tested; defaults to zero (i.e., 
%               H0 = the linear combination of coefficients is zero); if 
%               a vector is specified, then the number of entries should be 
%               equal to number of cells in L (i.e., for every cell in L, 
%               a separate hypothesised value)
%
% doF:          true or false indicating if F test should be conducted
%               instead of the default Wald test for and then use F
%               distribution for calculating the p values (default: false)
%
% numObs:       [1 x 1] numeric value indicating how many observations were
%               present at the time of estimating beta_hat and coeffCovar;
%               used for calculating df2 (df2 = numObs - p); alternatively,
%               it is possible to provide a vector of [1 x k] values, where
%               k is the number of cells in L, indicating that for every
%               cell, a different number of observation / degree of freedom
%               should be used
%
%% Outputs:
% W:            [k x v] vector of Wald statistics, where k is the number of
%               cells in the input L; if L was numeric, F is scalar
%
% p:            [k x v] vector of p values, where k is the number of cells
%               in the input L; if L was numeric, p is scalar; this is the
%               p value under F or chi-squared distribution testing the
%               null hypothesis that the linear combination of coefficients
%               are equal to the hypothesised value (or zero)
%
% logp:         [k x v] vector of log10 p values, where k is the number of cells
%               in the input L; if L was numeric, logp is scalar; this is the
%               logp value under F or chi-squared distribution testing the
%               null hypothesis that the linear combination of coefficients
%               are equal to the hypothesised value (or zero)
%
% df:           [k x 1] vector of numerator degrees of freedom, where k is
%               the number of cells in the input L; if L was numeric, df1 
%               is a scalar; df1 = number of linearly independent rows in
%               L, determined as the rank of that weight matrix
%
% L:            [k x 1] cell type containing the validated and zero-padded
%               contrats / weights (useful sanity check)
%
%% Notes:
% MATLAB calculates the F statistics in the same way as the Wald statistics
% calculated in this function; however, MATLAB additionally penalizes the F
% by dividing it with the numerator degrees of freedom. Then, they perform
% a lookup using a F distribution.

%% Reference:
% Fitzmaurice, G. M., Laird, N. M., & Ware, J. H. (2011). Applied 
% longitudinal analysis (2nd ed.). Wiley Series in Probability and 
% Statistics. John Wiley & Sons. [Page 98]
%
%% See also:
% coefTest.m
% Lines 524-599 of StandardLinearLikeMixedModel.m 
% [23.2.0.2409890 (R2023b)Update 3]

%% Check inputs
% Check L
if ~exist('L', 'var') || isempty(L)
    error('Please provide a vector/matrix/cell of weights to test');
else
    % If input is not a cell, make it a 1x1 cell
    if ~iscell(L)
        L = {L};
    end
    numCells = length(L);
end

% Check beta_hat
if ~exist('beta_hat', 'var') || isempty(beta_hat)
    error('Please provide a vector/matrix of beta coefficients');
else
    % Determine number of X and y variables
    [numX, numY] = size(beta_hat);
end

% Check coeffCovar
if ~exist('coeffCovar', 'var') || isempty(coeffCovar)
    error('Please provide coefficient covariance matrix');
else
    % Should be numX by numX by numY matrix
    [a, b, c] = size(coeffCovar);

    if a ~= b | a ~= numX | c ~= numY
        error(['Expected coefficient covariance matrix to be: ', ...
              num2str(numX), ' x ', num2str(numX), ' x ', num2str(numY)]);
    end
end

% Check hypValue
if ~exist('hypValue', 'var') || isempty(hypValue)
    hypValue = zeros(1, numCells);
else
    if isscalar(hypValue)
        hypValue = repmat(hypValue, 1, numCells);
    else
        % Ensure that is it organized as row vector
        hypValue = reshape(hypValue, [1, numel(hypValue)]);

        % Check that the number of elements match number of cells
        if size(hypValue, 2) ~= numCells
            error(['Either provide a scalar as the hypothesised value or ', ...
                   'provide a vector of 1 x ', num2str(numCells), ' values']);
        end
    end
end

% Check doF flag
if ~exist('doF', 'var') || isempty(doF)
    doF = false;
else
    if ~islogical(doF)
        error('doF should be either true or false');
    end
end

% Check numObs
if ~exist('numObs', 'var') || isempty(numObs)
    if doF
        error('Please provide number of observations or degrees of freedom');
    else
        numObs = [];
    end
else
    if isscalar(numObs)
        numObs = repmat(numObs, numCells, 1);
    else
        % Ensure that it is organized as a column vector
        numObs = reshape(numObs, [numel(numObs), 1]);

        % Check that the number of elements match number of cells
        if size(numObs, 2) ~= numCells
            error(['numObs should either be a scalar or a vector of ', ...
                   num2str(numCells), ' x 1 values']);
        end
    end
end

%% Initialize
[W, p] = deal(zeros(numCells, numY));
df     = zeros(numCells, 1);
df2    = numObs - numX;

%% Calculate, for every cell
for cc = 1:numCells

    % Current weights - validated and padded
    L{cc}       = validateContrast(L{cc}, numX);
    currWeights = L{cc};

    % Numerator degree of freedom = number of linearly independent rows
    % df(cc,1) = size(currWeights, 1);
    df(cc,1) = rank(currWeights);

    for yy = 1:numY
        % Linear combination of beta
        LB = currWeights * beta_hat(:,yy);
    
        % Difference from hypothesised mean
        LB = LB - hypValue(cc);
    
        % Calculate W2 = (L * B)' * inv{L * Cov * L'} * (L * B)
        % W2 = LB' * ((L * coeffCovar(:,:,yy) * L') \ eye(numX)) * LB;
        innerTerm = (currWeights * coeffCovar(:,:,yy) * currWeights');
        if rank(innerTerm) < size(innerTerm, 2)
            W(cc, yy) = (LB' *  (pinv(innerTerm) * LB));
        else
            W(cc, yy) = (LB' *  (innerTerm \ LB));
        end
    end

    % Calculate p values
    % p(cc,:) = 1 - chi2cdf(F, df1);
    if doF
        W(cc,:) = W(cc,:)./df(cc,1);
        p(cc,:) = fcdf(W(cc,:), df(cc,1), df2(cc,1), 'upper');
    else
        p(cc,:) = chi2cdf(W(cc,:), df(cc,1), 'upper');
    end
end

% Calculate -log10 p-values
logp = -log10(max(1e-300,p));

function L = validateContrast(L, numX)
% Function that validates and pads a contrast

% Get size of contrast matrix/vector
[sz_r, sz_c] = size(L);

% If both dimensions are 1, zero pad
if sz_r == 1 && sz_c == 1
    L = [L, zeros(1, numX-1)];
else
    % If number of rows is 1, zero pad
    if sz_r == 1
        L = [L, zeros(1, numX - size(L,2))];
    else
        % If number of columns is 1, warn, make row vector, zero pad
        if sz_c == 1
            warning(['Contrast/Weights has 1 column and ', num2str(sz_r), ...
                     ' rows; expected 1 row, multiple columns; transposing']);
            L = L';
            L = [L, zeros(1, numX - size(L,2))];
        else
            % Matrix case, need to zero pad columns
            L = padarray(L, [0, numX - size(L,2)], 0, 'post');
        end
    end
end

%% Deprecated solution using cellfun
% Typically loop is more efficient than using cellfun
% LB = L * beta_hat
% cell2mat(cellfun(@(x,y) x' * y * x, num2cell(LB,1)', ...
%          cellfun(@(x) x \ eye(numX), ...
%          cellfun(@(x) L * x * L', squeeze(num2cell(coeffCovar, [1,2])), ...
%          'UniformOutput', false), 'UniformOutput', false), 'UniformOutput', false));