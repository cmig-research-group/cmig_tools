function [basisFunction, bfRank, basisSubset, settings, timing] =          ...
          createBasisFunctions(age,       knots,  splineType, dfFlag,      ...
                               intercept, method, toDrop,     ageSubset,   ...
                               addConst,  outDir, optCommand, optAppend,   ...
                               cleanUp,   instance)
% Function to create basis functions, given a set of age / time and knots
% and other parameters - currently supports creation of natural cubic
% splines or B-splines using ns or bs functions from the splines package in
% R; or creating natural cubic splines with unit heights at knots using the 
% nsk function from the splines2 package in R
%
%% Inputs:
% age:              vector of age/time values for which basis functions
%                   need to be created
%
% knots:            vector of values which serve as knots
%
% splineType:       character; one of the following:
%                       * 'ns'  (nautral cubic splines)
%                       * 'bs'  (B-splines)
%                       * 'nsk' (natural cubic splines with unit heights at knots)
%
% dfFlag:           logical; can be used to specify degrees of freedom
%                   instead of knots; enter a single number corresponding 
%                   to the degrees of freedom for "knots" and specify
%                   "dfFlag" as true
%
% intercept:        logical; indicates if intercept should be used during
%                   the creation of the splines (this parameter is specific 
%                   to R - adding an intercept is not the same as adding a
%                   constant to the basis functions)
%
% method:           character; one of the following (see Notes):
%                       * 'default'
%                       * 'demean'
%                       * 'regress'
%                       * 'svd'
%
% toDrop:           numeric or character; indicates which basis functions
%                   to drop; can be set to -1 in which case no columns will
%                   be dropped; or more than one columns to be dropped can
%                   be specified as numeric vector; as character, the
%                   following can be specified:
%                       * 'first'
%                       * 'middle'
%                       * 'last'
%                       * 'none'
%
% ageSubset:        a subset of age to be used for creating basis
%                   functions; the corresponding values for basis functions
%                   for other age values are then interpolated (or
%                   extrapolated)
%
% addConst:         logical; if true, a vector of ones is added as the 
%                   first column of the output basisFunction
%
% outDir:           full path to where the output file should be 
%                   temporarily saved; if empty, pwd is used
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
% instance:         numeric; useful if function is being called in parallel
%                   to ensure independence of each call (and that files are
%                   not deleted, etc.)
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
% splineType:       'ns'
% dfFlag:           false
% intercept:        true
% method:           'default'
% toDrop:           'middle' (if method is not default or svd)
% ageSubset:        []
% addConst:         false if method is default, otherwise true
% outDir:           pwd
% optCommand:       ''
% optAppend:        ''
% cleanUp:          true
% instance:         1
% 
%% Notes:
% There are four types of methods supported:
% default:          this returns the output as such from R
%
% demean:           first, the basis functions are created; then, they are
%                   column-wise mean centered; next, one (or more) basis
%                   functions are dropped; finally a constant term is added
% 
% regress:          first, the basis functions are created; then, the
%                   effect of age is regressed from each of the basis
%                   functions; the residuals are the new basis function;
%                   next, one (or more) basis functions are dropped;
%                   finally, a constant term is added
%
% svd:              first, the basis functions are created; then, they are
%                   column-wise mean centered; next, their orthonormal
%                   basis is created (using singular value decomposition;
%                   this might drop one or more basis functions); finally,
%                   a constant term is added
%
% If toDrop is unspecified, depending on the method, the following happens:
% default:          no columns are dropped
%
% demean:           last column is dropped
%
% regress:          last column is dropped
%
% svd:              internally during the call to orth, one or more columns
%                   may be dropped; no columns are removed explicitly
% 
% Alternatively, if toDrop == -1, no columns will be dropped for 'demean'
% or 'regress' methods ('svd' may still drop column(s) implicitly)
%
% Alternatively, if toDrop > 0, one or more columns can be dropped; in this
% case, all methods will remove the specified columns including 'default'
% and 'svd' (over and above any implicit drops); same holds if toDrop is
% character type and is not equal to 'none'

%% Start global timer
tInit = tic;

%% Check inputs
tCheck = tic;

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
    splineType = 'ns';
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
    else
        if dfFlag
            if numel(knots) > 1
                error(['Expected one value for knots when dfFlag is true but found: ', num2str(numel(knots)), ' values']);
            end
        end
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
    if ~ismember(method, {'default', 'demean', 'regress', 'svd'})
        error(['Expected method to be one of: default, demean, regress, or svd but found: ', method]);
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

% Check instance
if ~exist('instance', 'var') || isempty(instance)
    instance = 1;
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
    % If method is default or SVD, do not drop
    if strcmpi(method, 'default') || strcmpi(method, 'svd')
        toDrop = [];
    else
        % Default to middle
        toDrop = round(median(1:nCols), 0);
    end
else
    if isnumeric(toDrop)
        toDrop = reshape(toDrop, 1, length(toDrop));

        % Check if nothing should be dropped
        if toDrop == -1
            toDrop = [];
        else
            % Check if the column number is within range
            if sum(ismember(toDrop, 1:nCols)) ~= length(toDrop)
                error(['toDrop values out of range; values specified were: ', ...
                       num2str(toDrop), ' and number of basis function columns are ', num2str(nCols)]);
            end
        end
    else
        if strcmpi(toDrop, 'first')
            toDrop = 1;
        else
            if strcmpi(toDrop, 'last')
                toDrop = nCols;
            else
                if strcmpi(toDrop, 'middle')
                    toDrop = round(median(1:nCols), 0);
                else
                    if strcmpi(toDrop, 'none')
                        toDrop = [];
                    else
                        error(['Unable to interpret , ', toDrop, ' as a value for for toDrop']);
                    end
                end
            end
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

% Record timing for end of checks
timing.tChecks = toc(tCheck);

%% Prepare data to be read in R
% This step can, of course, become a bit of a bottleneck if the data size
% is huge - unfortunately, MATLAB does not seem to have a supported engine
% that shares data with R (this solution exists for Python - maybe there is
% a spline equivalent in Python?)
tPrepR   = tic;
txt_date = [char(datetime('now', 'Format', 'uuuu-MM-dd_hhmmss')), '-', num2str(instance, '%03d')];
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

% Prepare command
command = [optCommand, optAppend, 'Rscript ', toCall,   ' ', ...
           splineType,   ' ', txt_age, ' ',   txt_knot, ' ', ...
           tmpIntercept, ' ', tmpDF,   ' ',   txt_out,  ' ', ...
           fileparts(toCall)];

% Record timing for end of preparing data for R
timing.tPrepR = toc(tPrepR);

% Call R
tCallR        = tic;
[status, txt] = system(command);
if logical(status)
    error(['Could not create basis function; please see output: ', txt]);
end

% Record timing for calling R
timing.tCallR = toc(tCallR);


%% Read data back into MATLAB
tRead         = tic;
fid           = fopen(txt_out, 'r');
basisFunction = fscanf(fid, [repmat('%f', 1, nCols-1), '%f\n'], [nCols, nObs]);
fclose(fid);
basisFunction = basisFunction';

% Make sure that the number of observations read back are correct
if nObs ~= size(basisFunction, 1)
    warning('Mismatch between number of observations for creating basis function and number of observations read back');
end

% Record timing for reading data back in R
timing.tRead = toc(tRead);

%% Decide if interpolation is necessary
tInterp = tic;
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

% Record timing for interpolation
timing.tInterpolation = toc(tInterp);

%% Do we need to modify the basis functions?
tModBF = tic;

if strcmpi(method, 'demean')
    basisFunction = basisFunction - mean(basisFunction);
else
    if strcmpi(method, 'regress')
        basisFunction = basisFunction - (age * (age \ basisFunction));
    else
        if strcmpi(method, 'svd')
            basisFunction = orth(basisFunction - mean(basisFunction));
        end
    end
end

%% Do we need to drop one or more columns?
basisFunction(:, toDrop) = [];

%% Do we need to add a constant?
if addConst
    basisFunction = [ones(length(age), 1), basisFunction];
end

% Record timing for modifying basis functions
timing.tModifyBF = toc(tModBF);

%% Delete temporary files
tCleanUp = tic;
if cleanUp
    delete(txt_age);
    delete(txt_knot);
    delete(txt_out);
end

% Record timing for deleting files
timing.tCleanUp = toc(tCleanUp);

%% Finally compute the rank of the basis functions and save settings
tRank = tic;
if nargout > 1
    bfRank = rank(basisFunction);
    
    % Additionally, give a warning if rank deficient
    if bfRank < size(basisFunction, 2)
        warning('Rank deficient basis function');
    end

    % Save settings
    settings.splineType = splineType;
    settings.knots      = knots;
    settings.dfFlag     = dfFlag;
    settings.intercept  = intercept;
    settings.method     = method;
    settings.toDrop     = toDrop;
    settings.addConst   = addConst;
    settings.instance   = instance;
    settings.toCall     = toCall;
    settings.optCommand = optCommand;
    settings.optAppend  = optAppend;
end

% Record timing for calculating rank and saving settings
timing.tRank_saveSettings = toc(tRank);

% Record timing for overall call
timing.tOverall = toc(tInit);

% Deprecated solution:
% piecewisePoly  = spline(knots,        eye(length(knots)));
% basisFunction  = ppval(piecewisePoly, age);