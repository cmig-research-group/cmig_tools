function [basisFunction, bfRank, basisSubset] =                                ...
          FEMA_create_basisFunctions(age,       knots,  splineType, dfFlag,    ...
                                     intercept, method, toDrop,     ageSubset, ...
                                     addConst,  outDir, optCommand, optAppend, cleanUp)
% Function to create basis functions, given a set of age / time, knots,
% and optionally other parameters
%% Inputs:
% age:              vector of age/time values for which basis functions
%                   need to be created
%
% knots:            vector of values which serve as knots
%
% splineType:       character; one of the following:
%                       * 'ns'  (nautral cubic splines)
%                       * 'bs'  (B-splines)
%                       * 'nsk' (natural cubic with unit heights at knots)
%
% dfFlag:           logical; can be used to specify degrees of freedom
%                   instead of knots; enter a single number corresponding 
%                   to the degrees of freedom for "knots" and specify
%                   "dfFlag" as true
%
% intercept:        logical; indicates if intercept should be added (this
%                   parameter is specific to R - adding an intercept is not
%                   the same as adding a constant to the basis functions)
%
% method:           character; one of the following (see Notes):
%                       * 'default'
%                       * 'demean'
%                       * 'regress'
%
% toDrop:           numeric; indicate which basis function to drop
%                   (relevant if method is 'demean' or 'regress'); can be
%                   set to -1 in which case no columns will be dropped
%
% ageSubset:        a subset of age to be used for creating basis
%                   functions; the corresponding values for basis functions
%                   for other age values are then interpolated
% 
% addConst:         logical; if true, a vector of ones is added as the 
%                   first column of the output basisFunction
%
% outDir:           full path to an existing output directory where
%                   temporary files will be created; if empty, pwd is used
%
% optCommand:       optional command(s) to be issued prior to invoking
%                   R; useful, for example, on cluster environments where
%                   modules may need to be loaded before R is accessible
%
% optAppend:        optional command to be appended to the call to
%                   'Rscript' - for example, '/usr/local/bin/' can be
%                   appended prior to 'Rscript', if Rscript is not on PATH
%
% cleanUp:          logical; if true, deletes all temporary files that were
%                   created along the way
%
%% Output(s):
% basisFunction:    a matrix containing the spline basis functions
%
% bfRank:           rank of the basis functions (useful sanity check)
%
% basisSubset:      if ageSubset was specified, then basisSubset contains
%                   the results obtained from R prior to interpolation
%
%% Defaults:
% splineType:       'bs'
% dfFlag:           false
% intercept:        true
% method:           'default'
% toDrop:           last column (if method is not default)
% ageSubset:        []
% addConst:         false if method is default, otherwise true
% outDir:           pwd
% optCommand:       ''
% optAppend:        ''
% cleanUp:          true

%% Notes:
% There are three types of methods supported:
% default:  this returns the output as such from R
%
% demean:   in this case, first the basis functions are created; then, they
%           are column-wise mean centered; then, if addConst is true, 
%           an overall constant term is added; finally, one (or more) basis 
%           functions - by default the last column is dropped; the column 
%           which is being dropped can be overridden by passing a numeric
%           input to toDrop; if toDrop contains more than one number, that
%           many basis functions will be dropped
% 
% regress:  in this case, first the basis functions are created; then, the
%           effect of age is regressed from each of the basis functions;
%           the residuals are the new basis function (and by default the
%           last column is dropped; this can be overridden by toDrop)

%% Check inputs
% Check age
if ~exist('age', 'var') || isempty(age)
    error('Please provide a vector of age');
else
    % Ensure a n x 1 vector
    try
        age = reshape(age, length(age), 1);
    catch
        sz  = size(age);
        error(['Expected a vector for age but found: ', num2str(sz(1)), ' x ', num2str(sz(2))]);
    end
end

% Check knots
if ~exist('knots', 'var') || isempty(knots)
    error('Please provide either knots or df; if df is provided, set dfFlag to true');
else
    % Ensure a n x 1 vector
    try
        knots = reshape(knots, length(knots), 1);
    catch
        sz  = size(knots);
        error(['Expected a vector for knots but found: ', num2str(sz(1)), ' x ', num2str(sz(2))]);
    end
end

% Check splineType
if ~exist('splineType', 'var') || isempty(splineType)
    splineType = 'bs';
else
    splineType = lower(splineType);
    if ~ismember(splineType, {'ns', 'bs', 'nsk'})
        error(['Expected splineType to be one of: ns, bs, or nsk but found: ', splineType]);
    end
end

% Check dfFlag
if ~exist('dfFlag', 'var') || isempty(dfFlag)
    dfFlag = false;
else
    if ~islogical(dfFlag)
        error('dfFlag should be either true or false');
    end
end

% Check intercept
if ~exist('intercept', 'var') || isempty(intercept)
    intercept = true;
else
    if ~islogical(intercept)
        error('intercept should be either true or false');
    end
end

% Check method
if ~exist('method', 'var') || isempty(method)
    method = 'default';
else
    method = lower(method);
    if ~ismember(method, {'default', 'demean', 'regress'})
        error(['Expected method to be one of: default, demean, or regress but found: ', method]);
    end
end

% Check ageSubset
if ~exist('ageSubset', 'var') || isempty(ageSubset)
    ageSubset = [];
    useSubset = false;
else
    % Ensure a n x 1 vector
    try
        ageSubset = reshape(ageSubset, length(ageSubset), 1);
    catch
        sz  = size(ageSubset);
        error(['Expected a vector for ageSubset but found: ', num2str(sz(1)), ' x ', num2str(sz(2))]);
    end
    useSubset = true;
end

% Check addConst
if ~exist('addConst', 'var') || isempty(addConst)
    if strcmpi(method, 'default')
        addConst = false;
    else
        addConst = true;
    end
else
    if ~islogical(addConst)
        error('addConst should be either true or false');
    end
end

% Check outDir
if ~exist('outDir', 'var') || isempty(outDir)
    outDir = pwd;
else
    if ~exist(outDir, 'dir')
        mkdir(outDir);
    end
end

% Check optCommand
if ~exist('optCommand', 'var') || isempty(optCommand)
    optCommand = '';
else
    if ~ischar(optCommand)
        error('optCommand should be character type');
    else
        % Append a semicolon at the end of optCommand
        optCommand = [optCommand, '; '];
    end
end

% Check optAppend
if ~exist('optAppend', 'var') || isempty(optAppend)
    optAppend = '';
else
    if ~ischar(optAppend)
        error('optAppend should be character type');
    end
end

% Check cleanUp
if ~exist('cleanUp', 'var') || isempty(cleanUp)
    cleanUp = true;
else
    if ~islogical(cleanUp)
        error('cleanUp should be either true or false');
    end
end

%% Work out how many basis function columns should exist
if dfFlag
    nCols = knots;
else
    switch(splineType)
        case 'bs'
            nCols = length(knots) + 3 + double(intercept);
        case 'ns'
            nCols = length(knots) + 1 + double(intercept);
        case 'nsk'
            nCols = length(knots) + 1 + double(intercept);
    end
end

%% Check toDrop
if ~exist('toDrop', 'var') || isempty(toDrop)
    if ~strcmpi(method, 'default')
        % Default to last column
        toDrop = nCols;
    else
        toDrop = [];
    end
else
    toDrop = reshape(toDrop, 1, length(toDrop));

    % Check if nothing should be dropped
    if toDrop == -1
        toDrop = [];
    else
        % Check if the column number is within range
        if sum(ismember(toDrop, 1:nCols)) ~= length(toDrop)
            error(['toDrop values out of range; values specified were: ', num2str(toDrop), ...
                   ' and number of basis function columns are ', num2str(nCols)]);
        end
    end
end

%% Locate createBasis.R
% This function is assumed to be in the cmig_tools_utils/r folder
toCall = fullfile(fileparts(fileparts(which('colvec'))), 'r', 'caller_createBasis.R');
if ~exist(toCall, 'file')
    % Perform a search
    toCall = which('caller_createBasis.R');
    if ~exist(toCall, 'file')
        error(['Unable to locate createBasis.R; ', ...
               'make sure that createBasis.R is in cmig_tools_utils/r folder ', ...
               'or else on your MATLAB path']);
    end
end

%% Prepare data to be read in R
% This step can, of course, become a bit of a bottleneck if the data size
% is huge - unfortunately, MATLAB does not seem to have a supported engine
% that shares data with R (this solution exists for Python - maybe there is
% a spline equivalent in Python?)
txt_date = char(datetime('now', 'Format', 'uuuu-MM-dd_hhmmss'));
txt_age  = fullfile(outDir, [txt_date, '-forBF_age.txt']);
txt_knot = fullfile(outDir, [txt_date, '-forBF_knots.txt']);
txt_out  = fullfile(outDir, [txt_date, '-forBF_basis.txt']);

if useSubset
    nObs = length(ageSubset);
    fid  = fopen(txt_age, 'w');
    fprintf(fid, '%g\n', ageSubset);
    fclose(fid);
else
    nObs = length(age);
    fid  = fopen(txt_age, 'w');
    fprintf(fid, '%g\n', age);
    fclose(fid);
end

fid = fopen(txt_knot, 'w');
fprintf(fid, '%g\n', knots);
fclose(fid);

%% Create R command and execute
if intercept
    tmpIntercept = 'TRUE';
else
    tmpIntercept = 'FALSE';
end

if dfFlag
    tmpDF = 'TRUE';
else
    tmpDF = 'FALSE';
end

command       = [optCommand, optAppend, 'Rscript ', toCall, ' ', splineType, ' ', ...
                 txt_age, ' ', txt_knot, ' ', tmpIntercept, ' ', tmpDF, ' ',      ...
                 txt_out, ' ', fileparts(toCall)];
[status, txt] = system(command);
if logical(status)
    error(['Could not create basis function; please see output: ', txt]);
end

%% Read data back into MATLAB
fid           = fopen(txt_out, 'r');
basisFunction = fscanf(fid, [repmat('%f', 1, nCols-1), '%f\n'], [nCols, nObs]);
fclose(fid);
basisFunction = basisFunction';

%% Decide if interpolation is necessary
if useSubset
    basisSubset   = basisFunction;
    basisFunction = interp1(ageSubset, basisSubset, age);

    % Check if interpolation results in NaN or Inf
    if any(logical(sum(isinf(basisFunction) | isnan(basisFunction), 2)))
        warning('Interpolation led to NaN or Inf values; trying linear extrapolation; results may be incorrect');
        basisFunction = interp1(ageSubset, basisSubset, age, 'linear', 'extrap');
    end
else
    basisSubset = [];
end

%% Do we need to demean or regress?
if strcmpi(method, 'demean')
    basisFunction = basisFunction - mean(basisFunction);
else
    if strcmpi(method, 'regress')
        basisFunction = basisFunction - (age * (age \ basisFunction));
    end
end

%% Do we need to drop one or more columns?
basisFunction(:, toDrop) = [];

%% Do we need to add a constant?
if addConst
    basisFunction = [ones(length(age), 1), basisFunction];
end

%% Delete temporary files
if cleanUp
    delete(txt_age);
    delete(txt_knot);
    delete(txt_out);
end

%% Finally compute the rank of the basis functions
bfRank = rank(basisFunction);

% Additionally, give a warning if rank deficient
if bfRank < size(basisFunction, 2)
    warning('Rank deficient basis function');
end

% Deprecated solution:
% piecewisePoly  = spline(knots,        eye(length(knots)));
% basisFunction  = ppval(piecewisePoly, age);