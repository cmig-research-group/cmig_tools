function caller_createBasisFunctions(file_valvec, dirOutput, outType, outPrefix, varargin)
% Caller function for createBasisFunctions (compiled version)
%% Mandatory inputs:
% ------------------
% file_valvec:      full path to a .csv file containing a vector of values 
%                   for which basis functions need to be created (for
%                   example, age / time)
%                       - option 1: values should be in a single column
%                                   without any column name
%                       - option 2: if the input is a full design matrix
%                                   (with variable names), then an
%                                   additional input argument 'varName'
%                                   must specify which column of the design
%                                   matrix needs to be used
% 
% dirOutput:        full path to where the results will be saved
% 
% outType:          which format should the results be saved; should be one
%                   of the following:
%                       * mat:  results are saved as .mat file containing
%                               the following variables: 
%                                   - basisFunction
%                                   - basisSubset
%                                   - Xvars
%                                   - bfRank
%                                   - settings
%                                   - timing
%                       * csv:  results are saved as separate .csv files;
%                               following variables are written out:
%                                   - basisFunction (with header)
%                                   - basisSubset   (with header)
%                                   - Xvars         (without header)
% 
%% Optional inputs (comma-separated name-value pairs):
% ----------------------------------------------------
% outPrefix:        prefix to the output file names (without extension;
%                   default: FEMA_createBasisFunctions-yyyyMMMdd-HHmmSS)
%
% varName:          only needed if file_valvec points to a full design
%                   matrix; in this case, varName must correspond to a
%                   column in the design matrix
% 
% file_knots:       full path to a .csv file containing a vector vector of
%                   values which serve as knots (single column without any
%                   variable names; default knot values are at the [0 25 50
%                   75 100] percentiles of valvec
%
% splineType:       should be one of the following (default: nsk):
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
%                   (default: [])
%
% method:           should be one of the following (default: svd):
%                       * 'svd'
%                       * 'raw'
% 
% minMax:           numeric vector: if a single value is provided, this is
%                   the number of linearly spaced values that will be
%                   created; if two values are provided, these are the
%                   minimum and the maximum values which are used for
%                   creating a vector linearly spaced values; if three
%                   values are provided, the third value is used as the
%                   number of linearly spaced values; if left empty, the
%                   min and max of valvec is used (default: span of the data)
%
% intercept:        logical; indicates if intercept should be used during
%                   the creation of the splines (this parameter is specific 
%                   to R - adding an intercept is not the same as adding a
%                   constant to the basis functions; default: true)
%
% optCommand:       optional command(s) to be issued prior to invoking
%                   R; useful, for example, on cluster environments where
%                   modules may need to be loaded before R is accessible
%                   (default: [])
%
% optAppend:        optional command to be appended to the call to
%                   'Rscript' - for example, '/usr/local/bin/' can be
%                   appended prior to 'Rscript', if Rscript is not on PATH
%                   (default: [])
%
% cleanUp:          logical; if true, deletes all temporary files that were
%                   created along the way (default: true)
%
% instance:         numeric; useful if function is being called in parallel
%                   to ensure independence of each call (and that files are
%                   not deleted, etc.; default: 1)

%% Start
updateString = [char(datetime('now')), ': caller_createBasisFunctions: job started'];
disp(updateString);

%% Check mandatory inputs
if ~exist('file_valvec', 'var') || isempty(file_valvec)
    error(['Please provide a full path to a file containing the list of values ', ...
           'for which basis functions need to be made OR ', ...
           'full path to a file containing the design matrix']);
else
    if ~exist(file_valvec, 'file')
        error(['Unable to find file: ', file_valvec]);
    end
end

if ~exist('dirOutput', 'var') || isempty(dirOutput)
    error('Please provide a full path to where results should be saved');
else
    if ~exist(dirOutput, 'dir')
        mkdir(dirOutput);
    end
end

if ~exist('outType', 'var') || isempty(outType)
    error('Please specify output type');
else
    if ~ismember(outType, {'mat', 'csv'})
        error(['Unknown outType specified: ', outType]);
    end
end

if ~exist('outPrefix', 'var') || isempty(outPrefix)
    outPrefix = ['FEMA_createBasisFunctions-', char(datetime('now', 'Format', 'yyyyMMMdd-HHmmSS'))];
end

%% Assign default optional inputs
p = inputParser;
addParameter(p, 'varName', []);
addParameter(p, 'file_knots', []);
addParameter(p, 'splineType', 'nsk');
addParameter(p, 'Xpowers', [])
addParameter(p, 'method', 'svd');
addParameter(p, 'minMax', []);
addParameter(p, 'intercept', true);
addParameter(p, 'optCommand', []);
addParameter(p, 'optAppend', []);
addParameter(p, 'cleanUp', true);
addParameter(p, 'instance', 1);

%% Parse optional inputs
parse(p, varargin{:})
varName     = p.Results.varName;
file_knots  = p.Results.file_knots;
splineType  = p.Results.splineType;
Xpowers     = p.Results.Xpowers;
method      = p.Results.method;
minMax      = p.Results.minMax;
intercept   = p.Results.intercept;
optCommand  = p.Results.optCommand;
optAppend   = p.Results.optAppend;
cleanUp     = p.Results.cleanUp;
instance    = p.Results.instance;

%% Load valvec
updateString = [char(datetime('now')), ': caller_createBasisFunctions: loading valvec'];
disp(updateString);

[~, ~, ext] = fileparts(file_valvec);
if isempty(ext)
    error('Please provide full path to file_valvec');
else
    if strcmpi(ext, '.csv')
        opts = detectImportOptions(file_valvec);

        % Is this a design matrix or a vector?
        if isscalar(opts.VariableNames)
            updateString = [char(datetime('now')), ': caller_createBasisFunctions: single column input detected'];
            disp(updateString);
            temp_X = readtable(file_valvec, 'ReadVariableNames', false);
            valvec = temp_X{:,1};
            clear temp_X
        else
            updateString = [char(datetime('now')), ': caller_createBasisFunctions: design matrix input detected'];
            disp(updateString);
            temp_X = readtable(file_valvec, 'Delimiter', '\t');

            % Check if varName was provided
            if ~exist('varName', 'var') || isempty(varName)
                error(['When specifying a design matrix, please additionally ', ...
                       'provide the name of the variable to be used to make basis functions']);
            else
                % Check that varName exists in the design matrix
                tmp = find(ismember(temp_X.Properties.VariableNames, varName));
                if isempty(tmp)
                    error(['Unable to find ', varName, ' in the column names of ', file_valvec]);
                else
                    valvec = temp_X{:,tmp};
                    clear temp_X
                end
            end
        end
    else
        error(['Unknown extension for file_X: ', ext]);
    end
end

%% Load knots, if required
if ~isempty(file_knots)
    updateString = [char(datetime('now')), ': caller_createBasisFunctions: loading knots'];
    disp(updateString);

    if ~exist(file_knots, 'file')
        error(['Unable to find: ', file_knots]);
    else
        [~, ~, ext] = fileparts(file_knots);
        if isempty(ext)
            error('Please provide full path to file_knots');
        else
            if strcmpi(ext, '.csv')
                temp_knots = readtable(file_knots, 'ReadVariableNames', false);

                % Make sure a single column exists
                if size(temp_knots, 2) > 1
                    error('More than one columns detected in the knots file');
                else
                    knots = temp_knots{:,1};
                end
            else
                error(['Unknown extension for file_knots: ', ext]);
            end
        end
    end
end

%% Pass everything to createBasisFunctions
% No need to check other optional arguments as they get checked as part of
% createBasisFunctions
updateString = [char(datetime('now')), ': caller_createBasisFunctions: calling createBasisFunctions'];
disp(updateString);

[basisFunction, basisSubset, Xvars, bfRank, settings, timing] =     ...
 createBasisFunctions(valvec,     knots,     splineType, Xpowers,   ...
                      method,     dirOutput, minMax,     intercept, ...
                      optCommand, optAppend, cleanUp,    instance);

%% Save results
updateString = [char(datetime('now')), ': caller_createBasisFunctions: saving results'];
disp(updateString);

% Make variable names for the output
varNames_basisFunction = strreplace(strcat({'BF_'}, num2str((1:size(basisFunction,2))')), ' ', '');
varNames_basisSubset   = strreplace(strcat({'BFsubset_'}, num2str((1:size(basisSubset,2))')), ' ', '');

if strcmpi(outType, 'mat')
    % Get a sense for all variables in the workspace
    tmpInfo  = whos;
    toSave   = {'basisFunction', 'basisSubset', 'Xvars', 'bfRank', 'settings', 'timing', ...
                'varNames_basisSubset', 'varNames_basisFunction'};
    saveName = fullfile(dirOutput, [outPrefix, '.mat']);

    if sum([tmpInfo(ismember({tmpInfo(:).name}', toSave)).bytes]) > 2^31
        save(saveName, 'basisFunction', 'basisSubset', 'Xvars', 'bfRank', 'settings', 'timing', ...
                       'varNames_basisSubset', 'varNames_basisFunction', '-v7.3');
    else
        save(saveName, 'basisFunction', 'basisSubset', 'Xvars', 'bfRank', 'settings', 'timing', ...
                       'varNames_basisSubset', 'varNames_basisFunction');
    end
else
    saveName = fullfile(dirOutput, [outPrefix, '-basisFunction.csv']);
    writetable(cell2table(num2cell(basisFunction), 'VariableNames', varNames_basisFunction), saveName);

    saveName = fullfile(dirOutput, [outPrefix, '-basisSubset.csv']);
    writetable(cell2table(num2cell(basisSubset), 'VariableNames', varNames_basisSubset), saveName);

    saveName = fullfile(dirOutput, [outPrefix, '-basis_Xvars.csv']);
    writetable(cell2table(num2cell(Xvars), 'VariableNames', {'X'}), saveName);
end