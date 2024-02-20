function basisFunction = createBasisFunctions(age, knots, intercept, addConst, ...
                                              outDir, optCommand, optAppend, cleanUp)
% Function to create basis functions, given a set of age / time and knots
% and other parameters - currently creates natural cubic splines using ns
% function from splines package in R
%% Inputs:
% age:              vector of age for which basis functions need to be
%                   created (should be sorted in ascending order)
%
% knots:            vector of values which serve as knots
%
% intercept:        logical; indicates if intercept should be added (this
%                   parameter is specific to R - adding an intercept is not
%                   the same as adding a constant to the basis functions)
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
%% Output(s):
% basisFunction:    a matrix of 1 + length(knots) + intercept + constant
%                   containing the natural cubic spline basis functions
%                   which can be used for non-linear expansion of age

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

    % Ensure age is sorted
    if ~issorted(age)
        error('age should be sorted in ascending order');
    end
end

% Check knots
if ~exist('knots', 'var') || isempty(knots)
    knots = [];
else
    % Ensure a n x 1 vector
    try
        knots = reshape(knots, length(knots), 1);
    catch
        sz  = size(knots);
        error(['Expected a vector for knots but found: ', num2str(sz(1)), ' x ', num2str(sz(2))]);
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

% Check addConst
if ~exist('addConst', 'var') || isempty(addConst)
    addConst = true;
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

%% Locate createBasisNS.R
% This function is assumed to be in the cmig_tools_utils/r folder
toCall = fullfile(fileparts(fileparts(which('colvec'))), 'r', 'createBasisNS.R');
if ~exist(toCall, 'file')
    % Perform a search
    toCall = which('createBasisNS.R');
    if ~exist(toCall, 'file')
        error(['Unable to locate createBasisNS.R; ', ...
               'make sure that createBasisNS.R is in cmig_tools_utils/r folder ', ...
               'or else on your MATLAB path']);
    end
end

%% Some useful variables 
% Number of output columns
nCols = 1 + length(knots) + double(intercept);

% Number of "observations"
nObs = length(age);

%% Prepare data to be read in R
% This step can, of course, become a bit of a bottleneck if the data size
% is huge - unfortunately, MATLAB does not seem to have a supported engine
% that shares data with R (this solution exists for Python - maybe there is
% a spline equivalent in Python?)
txt_date = char(datetime('now', 'Format', 'uuuu-MM-dd_hhmmss'));
txt_age  = fullfile(outDir, [txt_date, '-forBF_age.txt']);
txt_knot = fullfile(outDir, [txt_date, '-forBF_knots.txt']);
txt_out  = fullfile(outDir, [txt_date, '-forBF_basis.txt']);

fid = fopen(txt_age, 'w');
fprintf(fid, '%g\n', age);
fclose(fid);

fid = fopen(txt_knot, 'w');
fprintf(fid, '%g\n', knots);
fclose(fid);

%% Create R command and execute
command       = [optCommand, optAppend, 'Rscript ', toCall, ' ', txt_age,   ...
                 ' ', txt_knot, ' ', upper(char(string(intercept))), ' ', txt_out];
[status, txt] = system(command);
if logical(status)
    error(['Could not create basis function; please see output: ', txt]);
end

%% Read data back into MATLAB
fid           = fopen(txt_out, 'r');
basisFunction = fscanf(fid, [repmat('%f', 1, nCols-1), '%f\n'], [nCols, nObs]);
fclose(fid);
basisFunction = basisFunction';

%% Add a constant, if user wants
if addConst
    basisFunction = [ones(nObs, 1), basisFunction];
end

%% Delete temporary files
if cleanUp
    delete(txt_age);
    delete(txt_knot);
    delete(txt_out);
end

% Deprecated solution:
% piecewisePoly  = spline(knots,        eye(length(knots)));
% basisFunction  = ppval(piecewisePoly, age);