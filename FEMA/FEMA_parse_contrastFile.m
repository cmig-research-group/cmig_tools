function [contrasts, hypValues] = FEMA_parse_contrastFile(inName, colnames)
% Function to parse a minimally specified contrast csv file into different
% univariate and multivariate contrasts with appropriate zero padding
%% Inputs:
% --------
% inName:       full path to a csv file containing different contrasts
% 
% colnames:     cell type containing a list of column names of the design
%               matrix (critical that these are in the final order in which
%               the analysis will be performed; should not contain fid,
%               iid, eid, or agevec) 
%               OR 
%               full path to a .mat or a tab-delimited .csv file
%               containing the design matrix and other information
%                   - if .mat, it should contain 'colnames' as a variable
%                   - if .csv, the columns should be named as (in order):
%                       * fid:      the family ID for every observation
%                       * iid:      the individual ID for every observation
%                       * eid:      the event ID for every observation
%                       * agevec:   the age values for every observation
%                       * column 5 onwards should be design matrix
%                       variables (which will be used by this function)
%
%% Output:
% --------
% contrasts:    a cell type variable where each cell entry contains a
%               zero-padded contrast vector (univariate contrast) or matrix
%               (multivariate contrast)
% 
% hypValues:    a [1 x k] vector, where k are the number of cell entries in
%               contrasts, with each value corresponding to the
%               hypothesised value for the contrast in each cell
% 
%% Contrast file specification:
% -----------------------------
% The following section describes the contents of a contrast file:
% 
% Let there be six fixed effects: X1, X2, X3, X4, X5, and X6
% 
% Let the user want four contrasts:
% Contrast 1: X1 - X3 [univariate contrast]
% Contrast 2: X3 - X4 [univariate contrast]
% Contrast 3: omnibus across X3, X4, and X5
% Contrast 4: omnibus across X1, X3, and X5
% 
% Then, the contrast file can be defined as tab-separated file (extra space
% is added for visual clarity; in the file, these should be tab-separated
% values):
% First row:    Type    hypValue    X1  X3  X4  X5
% Second row:    u1        0        1   -1  0   0
% Third row:     u2        0        0   1   -1  0
% Fourth row:    m1        0        0   1   0   0
% Fifth row:     m1        0        0   0   1   0
% Sixth row:     m1        0        0   0   0   1
% Seventh row:   m2        0        1   0   0   0
% Eighth row:    m2        0        0   1   0   0
% Ninth row:     m2        0        0   0   0   1
% 
% In the above specification, the column 'Type' indicates whether the
% contrast type is 'univariate' or 'multivariate'; the 'univariate'
% contrast is one contrast per row while multivariate contrasts span
% multiple rows; the column 'hypValue' indicates the hypothesised value for
% null hypothesis testing (in this example, we are testing against zero)
% 
% Note that for multivariate contrasts, the code looks through all the
% corresponding rows for each multivariate contrast, and picks the unique
% nonzero value - if all values are zero, it does not matter; if any of
% these are non-zeros, we expect only one unique non-zero value; otherwise,
% an error is generated
% 
% For univariate contrast, the 'u' prefix indicates univariate; this is
% followed by a number which specifies the univariate contrast number; one
% contrast per line, therefore, there are two univariate contrasts in the
% example above
% 
% For multivariate contrast, the 'm' prefix indicates multivariate; this is
% followed by a number which specifics which rows should be read together.
% In the above example, m1 specifies a multivariate contrast where fourth,
% fifth, and sixth rows are read together to create a contrast matrix;
% similarly, seventh, eighth, and ninth rows are read together to create
% another multivariate contrast matrix
% 
% This function performs appropriate left and right zero padding and
% creates a full contrast vector/matrix, as required. The output is a cell
% type variable, where each entry in the cell contains a separate contrast
%                   X1  X2  X3  X4  X5
% First entry:      [1  0   -1  0   0]
% Second entry:     [0  0   1   -1  0]
% Third entry:      [0  0   1   0   0
%                    0  0   0   1   0
%                    0  0   0   0   1]
% Fourth entry:     [1  0   0   0   0
%                    0  0   1   0   0
%                    0  0   0   0   1]
% 
% The amount of zero padding that the user should provide is
% non-consequential as long as the variable names are correctly specified
% in the top row of the contrast file; the user may either entry the zero
% values or leave them empty

%% Check inputs
if ~exist('inName', 'var') || isempty(inName)
    error('Please provide full path to a file containing contrast information');
else
    if ~exist(inName, 'file')
        error(['Unable to find file: ', inName]);
    end
end

if ~exist('colnames', 'var') || isempty(colnames)
    error('Please provide a list of column names of the design matrix');
else
    if ischar(colnames) || isstring(colnames)
        if ~exist(colnames, 'file')
            error(['Unable to find: ', colnames]);
        else
            colnames = parse_designMatrix(colnames);
        end
    end
end

%% Read contrast file
contrasts_file = readtable(inName, 'ReadVariableNames', true, 'EmptyValue', 0);
conNames       = contrasts_file.Properties.VariableNames;

% Make sure we have the first column as 'Type'
if ~strcmpi(conNames{1}, 'type')
    error(['Expected the first column of contrast file to be Type but found ', conNames{1}]);
else
    % Make sure that the second column is hypValue
    if ~strcmpi(conNames{2}, 'hypValue')
        error(['Expected the second column of contrast file to be hypValue but found ', conNames{2}]);
    else
        conNames = conNames(3:end);
    end
end

% Ensure all variables in contrasts are part of colNames
chk_vars = ismember(conNames, colnames);
if sum(chk_vars) ~= length(conNames)
    missMsg = sprintf('%s, ', conNames{~chk_vars});
    error(['Following variables are present in the contrast file ', ...
           'but not in the design matrix: ', missMsg(1:end-2)]);
end

%% Get some basic information and initialize
wch_contrasts  = unique(contrasts_file.Type, 'stable');
num_contrasts  = length(wch_contrasts);
num_covariates = length(colnames);
contrasts      = cell(num_contrasts, 1);
hypValues      = zeros(1, num_contrasts);

%% How do conNames map to colNames
[~, ~, i_col] = intersect(conNames, colnames);

%% Loop over every entry in wch_contrasts and parse
for con = 1:num_contrasts
    if logical(regexpi(wch_contrasts{con}, '^u'))
        wch_row                  = strcmpi(wch_contrasts{con}, contrasts_file.Type);
        contrasts{con}           = zeros(1, num_covariates);
        contrasts{con}(1, i_col) = contrasts_file{wch_row, 3:end};
        hypValues(1, con)        = contrasts_file{wch_row, 2};
    else
        if logical(regexpi(wch_contrasts{con}, '^m'))
            wch_rows                 = strcmpi(wch_contrasts{con}, contrasts_file.Type);
            contrasts{con}           = zeros(sum(wch_rows), num_covariates);
            contrasts{con}(:, i_col) = contrasts_file{wch_rows, 3:end};
            temp_hypValues           = unique(nonzeros(contrasts_file{wch_rows, 2}));
            if isempty(temp_hypValues)
                hypValues(1, con)    = 0;
            else
                if isscalar(temp_hypValues)
                    hypValues(1, con) = temp_hypValues;
                else
                    error(['Found more than one hypValue for contrast ', wch_contrasts{con}]);
                end
            end
        else
            warning('Unknown contrast type detected; ignoring');
        end
    end
end
end

function colnames = parse_designMatrix(inFile)
[~, ~, ext] = fileparts(inFile);
if isempty(ext)
    error('Please provide full path to design matrix file');
else
    if strcmpi(ext, '.mat')
        temp_X = load(inFile);
        
        % Make sure colnames exists as a variable
        if ~ismember('colnames', fieldnames(temp_X))
            error(['Could not find variable colnames in: ', inFile]);
        else
            % Assign the right variables
            colnames = temp_X.colnames;
        end
    else
        if strcmpi(ext, '.csv')
            temp_X   = readtable(inFile, 'Delimiter', '\t');
            colnames = temp_X.Properties.VariableNames(5:end);
        else
            error(['Unknown extension for design matrix file: ', ext]);
        end
    end
end
end