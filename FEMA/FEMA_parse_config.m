function [FFX_names,        FFX_categorical,     FFX_vectorTransforms,   ...
          FFX_winsorize,    FFX_deltaTransforms, FFX_splineTransforms,   ...
          FFX_interactions, global_transform,    global_intercept,       ...
          loc_cont,         family_id] = FEMA_parse_config(configFile)
% Function that reads a user-created text configuration file and parses
% various information out of it which can then be used by FEMA_makeDesign
%% Input(s):
% configFile:               full path to a text configuration file
%
%% Output(s):
% FFX_names:                cell type having the names of all fixed effects
% 
% FFX_categorical:          structure type with three fields: 
%                               * name:             name of the fixed effect(s)
%                               * of_interest:      true/false
%                               * reference:        name of the level which
%                                                   serves as the reference 
% 
% FFX_vectorTransforms:     structure type with three fields:
%                               * name:             name of the fixed effect(s)
%                               * of_interest:      true/false
%                               * transformation:   local transform to apply
% 
% FFX_winsorize:            structure type with four fields:
%                               * name:             name of the fixed effect(s)
%                               * of_interest:      true/false
%                               * transformation:   'winsorize'
%                               * lower:            lower percentile
%                               * upper:            upper percentile
% 
% FFX_deltaTransforms:      structure type with three fields:
%                               * name:             name of the fixed effect(s)
%                               * of_interest:      true/false
%                               * transformation:   'delta'
% 
% FFX_splineTransforms:     structure with following fields:
%                               * name:             name of the fixed effect(s)
%                               * of_interest:      true/false
%                               * knots:            knot placement information
%                               * splineType:       type of splines to create
%                               * Xpowers:          number of powers to regress
%                               * method:           'svd' or 'raw'
%                               * minMax:           custom range for spline span
%                               * interecpt:        true/false
%                               * optCommand:       if additional command is to be run
%                               * optAppend:        if additional command is to be added
%                               * cleanUp:          true/false
%                               * instance:         instance number
% 
% FFX_global:               character type having the name of the global
%                           transformation (if any) that should be applied
%                           after all local transformations to continuous
%                           variables
% 
% global_intercept:         true or false indicating if an overall
%                           intercept should be added to the design matrix
% 
% loc_cont:                 logical vector of entries in FFX_names which
%                           are continuous variables
% 
% family_id:                character type having the name of the
%                           variable to use as the family_id

%% Decrypt text file
tmp_txt = fileread(configFile);

% Split into lines and separate into name-value pairs
cfg = cellfun(@(x) strsplit(x, ':'), strsplit(tmp_txt, '\n')', 'UniformOutput', false);

% Concatenate name-value pairs, get rid of leading or trailing spaces
cfg = strtrim(vertcat(cfg{:}));

%% Extract all fixed effects
loc_categorical = check_cfg_parameter('categorical', cfg(:,1), true);
loc_continuous  = check_cfg_parameter('continuous',  cfg(:,1), true);
FFX_names       = [];

if ~isempty(loc_continuous)
    names_cont = strtrim(strsplit(cfg{loc_continuous, 2}, ',')');
    FFX_names = [FFX_names; names_cont];
end

if ~isempty(loc_categorical)
    names_catg = strtrim(strsplit(cfg{loc_categorical, 2}, ',')');
    FFX_names = [FFX_names; names_catg];
end

if isempty(FFX_names)
    error('No continuous or categorical variables specified; nothing to do');
end

% Ensure names are unique
if length(unique(FFX_names)) ~= length(FFX_names)
    error('One or more fixed effects is duplicated in the config file');
end

%% Determine how many PCs do we need to retain?
% If PCs are to be added, add them to fixed effects
loc_geneticPC = check_cfg_parameter('n_gpcs', cfg(:,1), true);
if ~isempty(loc_geneticPC)
    numPCs = cfg{loc_geneticPC, 2};
    loc_basenamePC = check_cfg_parameter('gpc_basename', cfg(:,1), true);
    if isempty(loc_basenamePC)
        basename_PC = 'ab_g_stc__gen_pc__';
    else
        basename_PC = cfg{loc_basenamePC, 2};
    end

    % Make PC names and add them to fixed effects
    tmp_names = strcat(basename_PC, num2str((1:numPCs)', '%02d'));

    % Make sure that the PC names do not already exist
    tmp_exist = ismember(FFX_names, tmp_names);
    FFX_names = [FFX_names; tmp_names(~tmp_exist)];
end

%% Determine categorical and continuous
loc_cont = ismember(FFX_names, names_cont);
loc_catg = ~loc_cont;

%% Overall intercept
loc_intercept = check_cfg_parameter('intercept', cfg(:,1), true);
if ~isempty(loc_intercept)
    global_intercept = cfg{loc_intercept, 2};
    if ischar(global_intercept)
        global_intercept = logical(eval(global_intercept));
    end
else
    global_intercept = true;
end

%% Determine which fixed effects are variables of interest
loc_ofInterest = check_cfg_parameter('of_interest', cfg(:,1), true);
if loc_ofInterest
    vars_of_interest = strsplit(cfg{loc_ofInterest,2}, ' ');
else
    vars_of_interest = [];
end

%% Default reference level for categorical variables
loc_reference = check_cfg_parameter('reference', cfg(:,1), true);
if global_intercept
    global_refLevel = 'mode';
else
    global_refLevel = 'none';
end

%% Convert categorical into a standard structure format
% Split the reference line
% Format: <varName>=<string> <varName>=<string> ...
if ~isempty(loc_categorical)
    if ~isempty(loc_reference)
        cfg_ref_map = strsplit(cfg{loc_reference, 2}, ',')';
        cfg_ref_map = cellfun(@(x) strsplit(x, '='), cfg_ref_map, 'UniformOutput', false);
        try
            cfg_ref_map = vertcat(cfg_ref_map{:});
        catch
            error(['Please check the specification of reference level in the config file; ', ...
                   'should be: <varName>=<string> <varName>=<string> ...']);
        end
    else
        cfg_ref_map = [];
    end
    
    FFX_categorical = standardizeCategorical(FFX_names(loc_catg), cfg_ref_map, global_refLevel, vars_of_interest);
else
    FFX_categorical = [];
end

%% If global transformation is specified, extract that
loc_globalTransform = check_cfg_parameter({'global_transform', 'global', 'globalTransform'}, cfg(:,1), true);
if ~isempty(loc_globalTransform)
    global_transform = cfg{loc_globalTransform, 2};
else
    global_transform = 'none';
end

%% Family ID, if present
loc_FID = check_cfg_parameter('family_id', cfg(:,1), true);
if ~isempty(loc_FID)
    family_id = cfg{loc_FID,2};
else
    family_id = '';
end

if ~isempty(loc_continuous)
    %% Standardize vector transformations
    vectorTransforms = {'center', 'centre', 'demean', 'std', 'standardize', ...
                        'normalize', 'logn', 'log10', 'inverseranknorm',    ...
                        'ranknorm', 'int', 'none'};
    loc_vector = find(ismember(cfg(:,1), vectorTransforms));
    if ~isempty(loc_vector)
        FFX_vectorTransforms = standardizeVectorTransforms(cfg(loc_vector,:), vars_of_interest);
    else
        FFX_vectorTransforms = [];
    end

    %% Standardize winsorization
    winsorizeTransforms = {'winsorize', 'winsorise', 'winsor'};
    loc_winsorize = find(ismember(cfg(:,1), winsorizeTransforms));
    if ~isempty(loc_winsorize)
        FFX_winsorize = standardizeWinsorizeTransforms(cfg(loc_winsorize,2), vars_of_interest);
    else
        FFX_winsorize = [];
    end
    
    %% Standardize delta transformations
    loc_delta = check_cfg_parameter('delta', cfg, true);
    if ~isempty(loc_delta)
        FFX_deltaTransforms = standardizeDeltaTransforms(cfg{loc_delta,2}, vars_of_interest);
    else
        FFX_deltaTransforms = [];
    end
    
    %% Standardize spline transformations
    loc_splines = check_cfg_parameter('splines', cfg, false);
    if ~isempty(loc_splines)
        FFX_splineTransforms = standardizeSplineTransforms(cfg(loc_splines,2), vars_of_interest);
    else
        FFX_splineTransforms = [];
    end
    
    %% Update vector transforms with no transformation variables
    % Variables that are present but have no defined transformation, add 'none'
    % as the transformation type
    alreadyTransformed = [{FFX_categorical.name}, {FFX_vectorTransforms.name}, ...
                          {FFX_deltaTransforms.name}, {FFX_splineTransforms.name}]';
    remainingVariables = setdiff(FFX_names, alreadyTransformed);
    
    % Loop over these, and grow vectorTransforms
    count = length(FFX_vectorTransforms) + 1;
    for v = 1:length(remainingVariables)
        FFX_vectorTransforms(count).name        = remainingVariables{v};
        if strcmpi(remainingVariables{v}, vars_of_interest)
            FFX_vectorTransforms(count).of_interest = true;
        else
            FFX_vectorTransforms(count).of_interest = false;
        end
        FFX_vectorTransforms(count).transformation = 'none';
    
        count = count + 1;
    end
else
    FFX_vectorTransforms = [];
    FFX_deltaTransforms  = [];
    FFX_splineTransforms = [];
end

%% Handle interactions
loc_interact = check_cfg_parameter({'interaction', 'interact', 'interactions'}, cfg, true);
if ~isempty(loc_interact)
    FFX_interactions = standardizeInteractions(cfg{loc_interact,2});
end
end

function idx = check_cfg_parameter(paramName, cfg_cell, chkDup)
% Trivial function that performs a look up for a parameter in a param x 1
% cell configuration, checks to make sure that only one parameter exists,
% and then returns the location of this parameter in the cell
if ~exist('chkDup', 'var') || isempty(chkDup)
    chkDup = true;
end
idx = find(ismember(lower(cfg_cell), lower(paramName)));
if chkDup
    if length(idx) > 1
        error(['More than one parameter: ', paramName, ' found in configFile']);
    end
end
end


function output = standardizeCategorical(cfg, ref_map, defReference, of_interest)
% Output is a structure with fields: name, of_interest, reference level
output  = cell(length(cfg), 3);
outName = {'name', 'of_interest', 'reference'};

for lines = 1:length(cfg)

    % Name of this variable
    output{lines,1} = cfg{lines};
    
    % Is this variable of interest?
    if strcmpi(cfg{lines}, of_interest)
        output{lines,2} = true;
    else
        output{lines,2} = false;
    end

    % Reference level
    loc = find(strcmpi(cfg{lines}, ref_map(:,1)));
    if isempty(loc)
        output{lines,3} = defReference;
    else
        output{lines,3} = ref_map{loc,2};
    end
end

output = cell2struct(output', outName);
end


function output = standardizeVectorTransforms(cfg_vector, of_interest)
% Output is a structure with fields: name, of_interest, transformation

% First pass, work out how many variables to transform
numVars = sum(cellfun(@(x) length(strrep(strsplit(x, ','), ' ', '')), cfg_vector(:,2)));
output  = cell(numVars, 3);
outName = {'name', 'of_interest', 'transformation'};

% Now, go over every transform type, split variables, and update structure
count = 1;
for transform = 1:length(cfg_vector)

    % Name of this transformation
    transformName = cfg_vector{transform,1};

    % Which variables does this transformation apply to?
    tmp = strrep(strsplit(cfg_vector{transform,2}, ','), ' ', '')';

    for v = 1:length(tmp)
        % Name of this variable
        output{count,1} = tmp{v};

        % Is this variable of interest
        if strcmpi(tmp{v}, of_interest)
            output{count,2} = true;
        else
            output{count,2} = false;
        end

        % Name of the transform
        output{count,3} = transformName;

        % Update count
        count = count + 1;
    end
end
output = cell2struct(output', outName);
end


function output = standardizeWinsorizeTransforms(cfg_winsorize, of_interest)
% Output is a structure with fields: name, of_interest, transformation, lower, upper

% First pass, split into variables using comma
output  = cell(length(cfg_winsorize), 5);
outName = {'name', 'of_interest', 'transformation', 'lower', 'upper'};

% Assign defaults
output(:,3) = {'winsorize'};
output(:,4) = {0};
output(:,5) = {100};

% Now, go over every variable, split, and update structure
for lines = 1:length(cfg_winsorize)

    % Split on spaces
    temp = strsplit(cfg_winsorize{lines}, ' ');

    % Name of this variable
    output{lines,1} = temp{1};

    % Is this variable of interest?
    if strcmpi(temp{1}, of_interest)
        output{lines,2} = true;
    else
        output{lines,2} = false;
    end

    % Split using spaces
    if isscalar(temp)
        % User has not specified percentiles
        warning(['Winsorization selected for: ', output{lines,1}, ' but no percentiles specified; ignoring winsorization']);
    else
        % Go over every name-pair value, split on "=", and extract values
        temp = temp(2:end);
        for v = 1:length(temp)
            % Split on spaces to separate out name-pair values
            temp2 = strsplit(temp{v}, '=');

            if strcmpi(temp2{1}, 'lower')
                output{lines,4} = str2double(temp2{2});
            else
                if strcmpi(temp2{1}, 'upper')
                    output{lines,5} = str2double(temp2{2});
                else
                    error(['Unknown parameter: ', temp2{1}, ' specified for winsorization of ', output{lines,1}]);
                end
            end
        end
    end
end
output = cell2struct(output', outName);
end


function output = standardizeDeltaTransforms(cfg, of_interest)
% Output is a structure with fields: name, of_interest, transformation

% First pass, work out which variables to apply delta to
toApply = strrep(strsplit(cfg, ',')', ' ', '');
output  = cell(length(toApply), 3);
outName = {'name', 'of_interest', 'transformation'};

% Assign names of the variables
output(:,1) = toApply;
output(:,2) = num2cell(ismember(toApply, of_interest));
output(:,3) = {'delta'};

output = cell2struct(output', outName);
end


function output = standardizeSplineTransforms(cfg, of_interest)
% Output is a structure with fields: name, of_interest, knots, splineType,
%                                    Xpowers, method, minMax, interecpt,
%                                    optCommand, optAppend, cleanUp, instance

outName = {'name', 'of_interest', 'knots', 'splineType', 'Xpowers', 'method', ...
           'minMax', 'intercept', 'optCommand', 'optAppend', 'cleanUp', 'instance'};
output  = cell(length(cfg), length(outName));

% Default settings
def_knots       = {'percentiles'};
def_splineType  = {'nsk'};
def_Xpowers     = {[]};
def_method      = {'svd'};
def_minMax      = {[]};
def_intercept   = {true};
def_optCommand  = {''};
def_optAppend   = {''};
def_cleanUp     = {true};
def_instance    = {1};

% Assign defaults, then populate with what is found
output(:, 3)  = def_knots;
output(:, 4)  = def_splineType;
output(:, 5)  = def_Xpowers;
output(:, 6)  = def_method;
output(:, 7)  = def_minMax;
output(:, 8)  = def_intercept;
output(:, 9)  = def_optCommand;
output(:, 10) = def_optAppend;
output(:, 11) = def_cleanUp;
output(:, 12) = def_instance;

for lines = 1:length(cfg)

    % Parse this line to extract parameters
    tmp = cfg{lines};

    % Split by space to get name-pair values
    tmp = strsplit(tmp, ' ');

    % Name of this variable
    var_name = strtrim(tmp{1});
    output{lines,1} = var_name;

    % Is this variable of interest?
    if strcmpi(var_name, of_interest)
        output{lines,2} = true;
    else
        output{lines,2} = false;
    end

    % Now make a cell array of name-pair values only
    all_params = cellfun(@(x) strsplit(x, '='), tmp(2:end), 'UniformOutput', false);
    all_params = vertcat(all_params{:});

    % Go over every parameter, assign to relevant location
    for p = 1:size(all_params,1)

        % Name of the parameter
        tmp_name = all_params{p,1};

        % Values for this parameter
        tmp_value = all_params{p,2};

        % Which location to save the info to?
        loc = strcmpi(tmp_name, outName);

        % If loc is empty, some unspecified parameter was present
        if sum(loc) == 0
            tmp = strcat(outName(3:end), ',', {char(30)});
            tmp = horzcat(tmp{:});
            txt = sprintf('%s%s%s%s%s\n%s%s', 'Unknown parameter ', tmp_name, ' specified for splines of ', ...
                                              var_name, '; ', 'Valid parameters are: ', tmp(1:end-2));
            error(txt); %#ok<SPERR>
        else
            % Special handling
            if strcmpi(tmp_name, 'knots')
                if ismember(tmp_value, {'percentiles', 'percentile', 'quartiles', 'quartile'})
                    output{lines,loc} = tmp_value;
                else
                    output{lines,loc} = str2num(tmp_value); %#ok<ST2NM>
                end
            else
                if strcmpi(tmp_name, 'minMax')
                    output{lines,loc} = sort(str2num(tmp_value)); %#ok<ST2NM>
                else
                    if strcmpi(tmp_name, 'instance')
                        output{lines,loc} = eval(strtrim(tmp_value));
                    else
                        if strcmpi(tmp_name, 'Xpowers')
                            output{lines,loc} = str2double(tmp_value);
                        else
                            output{lines,loc} = strtrim(tmp_value);
                        end
                    end
                end
            end
        end
    end
end

output = cell2struct(output', outName);
end


function output = standardizeInteractions(cfg)
% Output is a structure with fields: interaction, of_interest
% Currently unable to handle interactions which may not be of interest
% Defaulting to of_interest being true for all interactions

% Initialize
% output  = cell(size(tmp,1), 2);
outName = {'interaction', 'of_interest'};

% Split based on comma
output(:,1) = strtrim(strsplit(cfg, ',')');
output(:,2) = {true};

output = cell2struct(output', outName);
end