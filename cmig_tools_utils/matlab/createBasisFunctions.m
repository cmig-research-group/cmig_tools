function [basisFunction, basisSubset, bfRank, settings, timing] =             ...
          createBasisFunctions(valvec,     knots,      splineType, Xpowers,   ...
                               method,     outDir,     minMax,     intercept, ...
                               optCommand, optAppend,  cleanUp,    instance)
% Function to create basis functions, given a set of values, knot values 
% and other parameters - supports creation of natural cubic splines (ns), 
% B-splines (bs), and natural cubic splines with unit heights at knots
% (nsk): ns and bs functions are created using the splines package in R 
% while nsk are created either using an implementation in MATLAB or using
% the splines2 package in R
%
% The default approach is to create nsk splines using MATLAB
%
% The approach for creating basis functions:
% 1) Calculate min and max for the values in valvec
% 2) Create a linearly spaced vector of 100 values between the min and max
% 3) Create basis functions
% 4) Linearly regress out the powers of X variables (if required)
% 5) Perform SVD on the demeaned residuals / basis functions (if required)
% 6) Perform interpolation and scale data back to the values in valvec
% 7) Delete any temporary files created along the way
%
% The powers in step 4 refers to the powers of the vector created in step 2
% 
%% Inputs:
% valvec:           vector of values for which basis functions need to be
%                   created (for example, age / time)
%
% knots:            vector of values which serve as knots
%
% splineType:       character; one of the following:
%                       * 'nsk'   (natural cubic splines with unit heights at knots)
%                       * 'nsk-R' (natural cubic splines with unit heights at knots created using nsk in R)
%                       * 'ns'    (nautral cubic splines)
%                       * 'bs'    (B-splines)
%
% Xpowers:          vector indicating which powers of the linearly spaced
%                   variables should be regressed from the basis functions;
%                   for example, 0:1 would mean regressing out the
%                   intercept and the linear effect of the variable; 0:2
%                   would mean regressing out the intercept, the linear
%                   effect, and the quadratic effect of the variable
%
% method:           character; one of the following (see Notes):
%                       * 'svd'
%                       * 'raw'
% 
% outDir:           full path to where the output file should be 
%                   temporarily saved; if empty, pwd is used
%
% minMax:           numeric vector: if a single value is provided, this is
%                   the number of linearly spaced values that will be
%                   created; if two values are provided, these are the
%                   minimum and the maximum values which are used for
%                   creating a vector linearly spaced values; if three
%                   values are provided, the third value is used as the
%                   number of linearly spaced values; if left empty, the
%                   min and max of valvec is used
%
% intercept:        logical; indicates if intercept should be used during
%                   the creation of the splines (this parameter is specific 
%                   to R - adding an intercept is not the same as adding a
%                   constant to the basis functions)
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
% basisSubset:      the basis functions prior to interpolation to full data
% 
% bfRank:           rank of the basis functions (useful sanity check)
%
% settings:         structure containing some of the settings used for
%                   creating the basis functions (useful for posteriety)
% 
% timing:           structure containing time taken by various steps
%
%% Defaults:
% knots:            0, 25, 50, 75, and 100th percentiles (if using MATLAB);
%                   otherwise 25, 50, and 75th percentiles; note that
%                   percentiles are calculated of the valvec
% splineType:       'nsk'
% Xpowers:          []
% intercept:        true
% method:           'svd'
% outDir:           pwd
% minMax:           min of valvec, max of valvec
% optCommand:       ''
% optAppend:        ''
% cleanUp:          true
% instance:         1
% 
%% Notes:
% There are two supported methods:
%
% svd:              first, the basis functions are created; then, the
%                   powers of the uniformly spaced vector are regressed
%                   out; next, the residualized basis functions are
%                   column-wise mean centered, followed by creating their
%                   orthonormal basis (using singular value decomposition;
%                   this might drop one or more basis functions); these are
%                   then scaled to have values between 0 and 1; finally,
%                   the values are interpolated to generate the basis
%                   function for the full data
% 
% raw:              this returns the output as such; the basis functions
%                   are computed for the linearly spaced vector, powers are
%                   regressed out (if the user specifies so), and then the
%                   values are interpolated to generate the basis functions
%                   for the full data
%
%% Start global timer
tInit = tic;

%% Check inputs
tCheck = tic;

% Check valvec
if ~exist('valvec', 'var') || isempty(valvec)
    error('Please provide a vector of values to work with');
else
    % Ensure a n x 1 vector
    try
        valvec = reshape(valvec, length(valvec), 1);
    catch
        sz = size(valvec);
        error(['Expected a vector for valvec but found: ', num2str(sz(1)), ' x ', num2str(sz(2))]);
    end
end

% Check knots
if ~exist('knots', 'var') || isempty(knots)
    createKnots = true;
else
    % Ensure a n x 1 vector
    try
        knots = reshape(knots, length(knots), 1);
        createKnots = false;
    catch
        sz  = size(knots);
        error(['Expected a vector for knots but found: ', num2str(sz(1)), ' x ', num2str(sz(2))]);
    end
end

% Check splineType
if ~exist('splineType', 'var') || isempty(splineType)
    splineType = 'nsk';
else
    splineType = lower(splineType);
    if ~ismember(splineType, {'ns', 'bs', 'nsk', 'nsk-r', 'nskr'})
        error(['Expected splineType to be one of: ns, bs, nsk, nsk-R but found: ', splineType]);
    end
end

% Check Xpowers
if ~exist('Xpowers', 'var') || isempty(Xpowers)
    toRegress = false;
else
    % Ensure a vector of powers
    try
        Xpowers = reshape(Xpowers, 1, length(Xpowers));
        toRegress = true;
    catch
        sz  = size(Xpowers);
        error(['Expected a vector for Xpowers but found: ', num2str(sz(1)), ' x ', num2str(sz(2))]);
    end

    % Ensure powers are not negative
    if any(Xpowers < 0)
        error('One or more Xpowers indicarted are negative');
    end
end

% Check method
if ~exist('method', 'var') || isempty(method)
    method = 'svd';
else
    method = lower(method);
    if ~ismember(method, {'raw', 'svd'})
        error(['Expected method to be one of: raw or svd but found: ', method]);
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

% Check minMax
if ~exist('minMax', 'var') || isempty(minMax)
    calcMinMax = true;
    howMany    = 101;
else
    if isscalar(minMax)
        howMany = minMax;
        calcMinMax = true;
    else
        if numel(minMax) == 2
            howMany = 101;
            calcMinMax = false;
        else
            if numel(minMax) == 3
                howMany    = minMax(3);
                minMax     = minMax(1:2);
                calcMinMax = false;
            else
                error(['Expected a vector of either one, or two, or three values for minMax but found ', num2str(numel(minMax)), ' values']);
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

%% Determine if MATLAB or R should be used for making basis functions
if strcmpi(splineType, 'nsk')
    doR = false;
else
    doR = true;
    if ismember(splineType, {'nsk-r', 'nskr'})
        splineType = 'nsk';
    end
end

%% Generate a vector of linearly spaced values
if calcMinMax
    minX = min(valvec);
    maxX = max(valvec);
else
    minX = min(minMax);
    maxX = max(minMax);
end
Xvars = linspace(minX, maxX, howMany);

%% Determine knots
if createKnots
    if doR
        knots = prctile(valvec, [25 50 75]);
    else
        knots = prctile(valvec, [0 25 50 75 100]);
    end
else
    if ~doR
        % Add the min and max boundary knots if they are not already there
        % Does min knot exist?
        if ~any(knots == minX)
            knots = [minX; knots];
            disp(['Creating nsk splines using MATLAB; added knot at: ', num2str(minX)]);
        end

        % Does max knot exist?
        if ~any(knots == maxX)
            knots = [knots; maxX];
            disp(['Creating nsk splines using MATLAB; added knot at: ', num2str(maxX)]);
        end
    end
end

%% If using MATLAB, create basis functions now
if ~doR
    % Record timing for end of checks
    timing.tChecks = toc(tCheck);

    % Start timer for creating basis functions
    tCreate = tic;

    % Create basis function
    pp            = csape(knots, eye(length(knots)), 'variational'); 
    basisSubset = ppval(pp, Xvars)';

    % End timer for creating basis functions
    timing.tCreate = toc(tCreate);
else
    %% Work out how many basis function columns should exist
    if doR
        switch(splineType)
            case 'bs'
                nCols = length(knots) + 3 + double(intercept);
            case 'ns'
                nCols = length(knots) + 1 + double(intercept);
            case 'nsk'
                nCols = length(knots) + 1 + double(intercept);
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
    tPrepR   = tic;
    txt_date = [char(datetime('now', 'Format', 'uuuu-MM-dd_hhmmss')), '-', num2str(instance, '%03d')];
    txt_age  = fullfile(outDir, [txt_date, '-forBF_age.txt']);
    txt_knot = fullfile(outDir, [txt_date, '-forBF_knots.txt']);
    txt_out  = fullfile(outDir, [txt_date, '-forBF_basis.txt']);
    
    nObs = length(Xvars);
    fid  = fopen(txt_age, 'w');
    fprintf(fid, '%g\n', Xvars);
    fclose(fid);
    
    fid = fopen(txt_knot, 'w');
    fprintf(fid, '%g\n', knots);
    fclose(fid);
    
    %% Create R command and execute
    if intercept
        tmpIntercept = 'TRUE';
    else
        tmpIntercept = 'FALSE';
    end
    
    % Prepare command
    command = [optCommand, optAppend, 'Rscript ', toCall,     ' ',    ...
               txt_age,      ' ', txt_knot, ' ',  splineType, ' ',    ...
               tmpIntercept, ' ', txt_out,  ' ',  fileparts(toCall)];

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
    basisSubset   = fscanf(fid, [repmat('%f', 1, nCols-1), '%f\n'], [nCols, nObs]);
    fclose(fid);
    basisSubset   = basisSubset';
    
    % Make sure that the number of observations read back are correct
    if nObs ~= size(basisSubset, 1)
        warning('Mismatch between number of observations for creating basis function and number of observations read back');
    end
    
    % Record timing for reading data back in R
    timing.tRead = toc(tRead);
end

%% Do we need to modify the basis functions?
tModBF = tic;

% If Xpowers exist, regress them from basis functions
if toRegress
    % Ensure Xvars is a column vector
    Xvars = reshape(Xvars, [], 1);
    Xexpanded = Xvars .^ Xpowers;
    if rank(Xexpanded) == size(Xexpanded,2)
        basisSubset = basisSubset - (Xexpanded * (Xexpanded \ basisSubset));
    else
        try
            basisSubset = basisSubset - (Xexpanded * lsqminnorm(Xexpanded, eye(size(Xexpanded))) * basisSubset);
        catch
            basisSubset = basisSubset - (Xexpanded * pinv(Xexpanded) * basisSubset);
        end
    end
end

% Perform SVD on de-meaned basis functions; if regression was performed in 
% the previous step, the basis functions would already be mean centered 
% (assuming intercept was included as a power); additionally, scale the 
% orthonormal span to be between 0 and 1
if strcmpi(method, 'svd')

    % Demean
    basisSubset = basisSubset - mean(basisSubset);

    % Orthonormal span
    [U, S, V] = svd(basisSubset); %#ok<ASGLU>

    % Drop columns
    if toRegress
        basisSubset = U(:, 1:size(basisSubset, 2) - size(Xexpanded,2));
    else
        basisSubset = U(:, 1:size(basisSubset, 2) - size(Xvars,1));
    end
    % basisFunction = orth(basisFunction - mean(basisFunction));

    % Re-scale the values using column-wise maximum
    basisSubset = basisSubset ./ max(abs(basisSubset), [], 1);
end

% Record timing for modifying basis functions
timing.tModBF = tModBF;

%% Interpolate to full dataset
tInterp       = tic;
basisFunction = interp1(Xvars, basisSubset, valvec, 'linear');

% Check if interpolation results in NaN or Inf
if any(logical(sum(isinf(basisFunction) | isnan(basisFunction), 2)))
    warning('Interpolation led to NaN or Inf values; trying linear extrapolation; results may be incorrect');
    basisFunction = interp1(Xvars, basisSubset, valvec, 'linear', 'extrap');
end

% Record timing for interpolation
timing.tInterpolation = toc(tInterp);

%% Delete temporary files
tCleanUp = tic;
if cleanUp & doR
    delete(txt_age);
    delete(txt_knot);
    delete(txt_out);
end

% Record timing for deleting files
timing.tCleanUp = toc(tCleanUp);

%% Finally compute the rank of the basis functions and save settings
tRank = tic;
if nargout > 2
    bfRank = rank(basisFunction);
    
    % Additionally, give a warning if rank deficient
    if bfRank < size(basisFunction, 2)
        warning('Rank deficient basis function');
    end

    % Save settings
    settings.useR       = doR;
    settings.splineType = splineType;
    settings.knots      = knots;
    settings.dfFlag     = dfFlag;
    settings.intercept  = intercept;
    settings.method     = method;
    settings.instance   = instance;
    settings.toCall     = toCall;
    settings.optCommand = optCommand;
    settings.optAppend  = optAppend;
    settings.Xvars      = Xvars;
    settings.minX       = minX;
    settings.maxX       = maxX;

    if toRegress
        settings.Xexpanded  = Xexpanded;
    end
end

% Record timing for calculating rank and saving settings
timing.tRank_saveSettings = toc(tRank);

% Record timing for overall call
timing.tOverall = toc(tInit);