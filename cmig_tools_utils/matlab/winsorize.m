function [outputVar, lower_bound, upper_bound] = winsorize(inputVars, varargin)
% Function to winsorize a given vector or matrix with pre-specified lower
% and upper percentile range
% 
%% Input(s):
% inputVars:    [n x p]             vector or matrix where every column
%                                   needs to be winsorized
% 
%% Optional inputs:
% lower_bound:  [1 x 1] OR [1 x p]  scalar or vector determining the lower
%                                   percentile range (default: 0)
% 
% upper_bound:  [1 x 1] OR [1 x p]  scalar or vector determining the upper
%                                   percentile range (default: 100)
%
%% Output(s):
% outputVar:    [n x p]             vector or matrix with winsorized data

%% Check inputs
p = inputParser;
addRequired(p, 'inputVars', @(x) validateattributes(x, {'numeric'}, {'2d'}));
numCols = size(inputVars,2);

addOptional(p, 'lower_bound', 0,   @(x) isnumeric(x) && (isscalar(x) || (isvector(x) && size(x,2) == numCols)));
addOptional(p, 'upper_bound', 100, @(x) isnumeric(x) && (isscalar(x) || (isvector(x) && size(x,2) == numCols)));

parse(p, inputVars, varargin{:});
inputVars   = p.Results.inputVars;
lower_bound = p.Results.lower_bound;
upper_bound = p.Results.upper_bound;

% If lower and upper are scalars, winsorize in one fell swoop
if isscalar(lower_bound) && isscalar(upper_bound)

    % Compute the lower and upper bounds for winsorization
    lower_bound = prctile(inputVars, lower_bound, 1);
    upper_bound = prctile(inputVars, upper_bound, 1);
    
    % Apply winsorization
    outputVar = inputVars;
    outputVar(outputVar < lower_bound) = lower_bound;
    outputVar(outputVar > upper_bound) = upper_bound;
else
    % Use repmat if either is a scalar
    if isscalar(lower_bound)
        lower_bound = repmat(lower_bound, 1, numCols);
    end
    if isscalar(upper_bound)
        upper_bound = repmat(upper_bound, 1, numCols);
    end

    % Go over each column, compute percentile, replace values
    outputVar = inputVars;
    for cols = 1:numCols
        % Compute the lower and upper bounds for winsorization
        temp_lower = prctile(inputVars(:,cols), lower_bound(cols), 1);
        temp_upper = prctile(inputVars(:,cols), upper_bound(cols), 1);

        % Apply winsorization
        outputVar(outputVar(:,cols) < temp_lower, cols) = temp_lower;
        outputVar(outputVar(:,cols) > temp_upper, cols) = temp_upper;
    end
end