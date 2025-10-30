function [output, settings] = doTransformation(inputVars, transformType, varargin)
% Function that applies a specified transformation to input variable(s)
%% Inputs:
% inputVars:        [n x p]     vector or matrix of n observations and p
%                               (continuous) variables; transformations are
%                               applied independently on each of the p
%                               variables
% 
% transformType:    character   case-insensitive; specifies the type of
%                               transformation that should be applied;
%                               should be one of the following (see Notes):
%                                   * 'center' | 'centre' | 'demean'
%                                   * 'std' | 'standardize' | 'normalize'
%                                   * 'logn'
%                                   * 'log10'
%                                   * 'inverseranknorm' | 'ranknorm' | 'int'
%                                   * 'winsorize' | 'winsorise' | 'winsor'
%
%% Optional inputs:
% lower:    [1 x 1] OR [1 x p]  scalar or vector determining the lower
%                               percentile range (default: 0)
% 
% upper:    [1 x 1] OR [1 x p]  scalar or vector determining the upper
%                               percentile range (default: 100)
% 
%% Outputs:
% output:           [n x p]     vector or matrix of transformed p variables
%                               (returned in the same order as in inputVars)
% 
% settings:         structure   settings that were applied to create the
%                               transformed variables, including timing
%                               information
% 
%% Notes:
% transformType:    'center' | 'centre' | 'demean'
% column-wise mean subtraction (i.e., every column of output has zero
% mean); NaN values are ignored during mean calculation
% 
% transformType:    'std' | 'standardize' | 'normalize'
% column-wise mean subtraction and division by column-wise standard
% deviation (i.e., every column of output has zero mean and unit standard
% deviation); NaN values are ignored during the calculation of the mean and
% the standard deviation
% 
% transformType:    'logn'
% Natural log transformation of every column of the input; may result in
% complex valued output, in which case a warning is additionally displayed
% 
% 
% transformType:    'log10'
% log transformation (base 10) of every column of the input; may result in
% complex valued output, in which case a warning is additionally displayed
% 
% 
% transformType:    'inverseranknorm' | 'ranknorm' | 'int'
% inverse normal rank transformation of every column of the input; uses the
% function rank_based_INT
%
% transformType:    'winsorize' | 'winsorise' | 'winsor'
% winsorization of every column of the input; uses the function winsorize

p = inputParser;
addParameter(p, 'lower', 0);
addParameter(p, 'upper', 100);
parse(p, varargin{:});
lower_bound = p.Results.lower;
upper_bound = p.Results.upper;

%% Check inputs
tInit = tic;

if ~exist('inputVars', 'var') || isempty(inputVars)
    error('Please provide a vector or matrix of variables to work on');
end

if ~exist('transformType', 'var') || isempty(transformType)
    error('Please specify the transformation type');
else
    if ~ischar(transformType)
        error('transformType should be of character type');
    else
        transformType = lower(transformType);
        supportedTransforms = {'center', 'centre', 'demean',         ...
                               'std', 'standardize', 'normalize',    ...
                               'logn', 'log10',                      ...
                               'inverseranknorm', 'ranknorm', 'int', ...
                               'winsorize', 'winsorise', 'winsor'};

        if ~ismember(transformType, supportedTransforms)
            tmp = strcat(supportedTransforms, ',', {' '});
            tmp = horzcat(tmp{:});
            error(['Unknown transformType specified: ', transformType, ...
                   '; supported transformations are: ', tmp(1:end-2)]);
        end
    end
end

%% Do transformations
settings.transformType = transformType;

switch transformType

    case {'center', 'centre', 'demean'}
        meanVec = mean(inputVars, 1, 'omitnan');
        output  = inputVars - meanVec;
        settings.mean        = meanVec;
        settings.omitNaN     = true;
        settings.timeTaken   = toc(tInit);
        
    case {'std', 'standardize', 'normalize'}
        meanVec = mean(inputVars, 1, 'omitnan');
        stdVec  = std(inputVars, [], 1, 'omitnan');
        output  = (inputVars - meanVec)./stdVec;
        settings.mean        = meanVec;
        settings.std         = stdVec;
        settings.omitNaN     = true;
        settings.timeTaken   = toc(tInit);
        
    case 'logn'
        output = log(inputVars);
        if any(imag(output) ~= 0, 'all')
            warning('log transformed resulted in complex numbers');
            settings.complexWarn = true;
        else
            settings.complexWarn = false;
        end
        settings.timeTaken = toc(tInit);

    case 'log10'
        output = log10(inputVars);
        if any(imag(output) ~= 0, 'all')
            warning('log10 transformed resulted in complex numbers');
            settings.complexWarn = true;
        else
            settings.complexWarn = false;
        end
        settings.timeTaken = toc(tInit);
        
    case {'inverseranknorm', 'ranknorm', 'int'}
        output = rank_based_INT(inputVars);
        settings.timeTaken = toc(tInit);

    case {'winsorize', 'winsorise', 'winsor'}
        [output, settings.lower, settings.upper] = winsorize(inputVars, 'lower_bound', lower_bound, 'upper_bound', upper_bound);
        settings.timeTaken = toc(tInit);
end