function [FFX_names,            FFX_categorical,     FFX_vectorTransforms,    ...
          FFX_winsorize,        FFX_deltaTransforms, FFX_quadraticTransforms, ...
          FFX_splineTransforms, FFX_interactions,    global_transform,        ...
          global_intercept,     loc_cont,            names_of_interest,       ...
          RFX_names] = FEMA_parse_JSON(configFile)
% Function that reads a DEAP-created JSON specification files and parses
% various information out of it which can then be used by FEMA_makeDesign
%% Input(s):
% configFile:               full path to a JSON specification file
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
% FFX_global:               character/string type having the name of the 
%                           global transformation (if any) that should be 
%                           applied after all local transformations to 
%                           continuous variables
% 
% global_intercept:         true or false indicating if an overall
%                           intercept should be added to the design matrix
% 
% loc_cont:                 logical vector of entries in FFX_names which
%                           are continuous variables
% 

%% Decrypt JSON file
configFile = jsondecode(fileread(configFile));

if strcmpi(configFile.params.dependent.type_data, 'external')
    tmp_locs = find(cellfun(@(x) isfield(x, 'transform'), configFile.params.fixed.vars));
    for tmp = 1:length(tmp_locs)
        if strcmpi(configFile.params.fixed.vars{tmp}.transform, 'splines')
            configFile.params.fixed.vars{tmp}.of_interest = true;
        end
    end
end

%% Extract random effects
RFX_names = configFile.params.random;

%% Extract all fixed effects
cfg_fixed = configFile.params.fixed;
if isstruct(cfg_fixed.vars)
    cfg_fixed_back = cfg_fixed;
    tmp = num2cell(cfg_fixed.vars);
    cfg_fixed = rmfield(cfg_fixed, 'vars');
    cfg_fixed.vars = tmp;
end

% Handle a case where the fixed effects are not specified
skipSteps = false;
if isempty(cfg_fixed.vars)
    skipSteps = true;
end

if ~skipSteps
    % Are custom names specified for any fixed effects?
    loc_custom = find(cellfun(@(x) isfield(x, 'name_custom'), cfg_fixed.vars, 'UniformOutput', true));

    % If custom names are present, overwrite the "name" field so that it
    % has the name_custom value
    if ~isempty(loc_custom)
        % FFX_names(loc_custom) = cellfun(@(x) x.name_custom, cfg_fixed.vars(loc_custom), 'UniformOutput', false);
        for ii = 1:length(loc_custom)
            cfg_fixed.vars{loc_custom(ii)}.name = cfg_fixed.vars{loc_custom(ii)}.name_custom;
        end
    end    

    % Extract all fixed effects name
    FFX_names = cellfun(@(x) x.name, cfg_fixed.vars, 'UniformOutput', false);
    
    % Ensure names are unique
    if length(unique(FFX_names)) ~= length(FFX_names)
        error('One or more fixed effects is duplicated in the config file');
    end
else
    FFX_names = {};
end

%% Determine how many PCs do we need to retain?
if isfield(cfg_fixed, 'n_gpcs')
    numPCs      = cfg_fixed.n_gpcs;

    % Handle a case where no PCs are specified
    if numPCs > 0

        basename_PC = {'pc_'};
    
        % Make PC names and add them to fixed effects
        tmp_names = strcat(basename_PC, num2str((1:numPCs)', '%02d'));
    
        % Make sure that the PC names do not already exist
        tmp_exist = ismember(tmp_names, FFX_names);
        FFX_names = [FFX_names; tmp_names(~tmp_exist)];
    end
end

% % Expand the locs_of_interest to full size
% locs_of_interest = false(length(FFX_names),1);
% locs_of_interest(ismember(FFX_names, names_of_interest),1) = true;

%% Determine categorical and continuous
if ~skipSteps
    loc_cont = cellfun(@(x) strcmpi(x.type_var, 'continuous'), cfg_fixed.vars);
    loc_catg = ~loc_cont;
else
    loc_cont = [];
    loc_catg = [];
end

%% Overall intercept
if isfield(configFile.params.fixed, 'intercept')
    global_intercept = configFile.params.fixed.intercept;
else 
    global_intercept = true;
end
if ~islogical(global_intercept)
    global_intercept = logical(global_intercept);
end

%% Default reference level for categorical variables
if global_intercept
    global_refLevel = 'mode';
else
    global_refLevel = 'none';
end

%% Convert categorical into a standard structure format
if sum(loc_catg) == 0
    FFX_categorical = [];
else
    FFX_categorical = standardizeCategorical(cfg_fixed.vars(loc_catg), global_refLevel);
end

%% Subset continuous variables
if ~skipSteps
    cfg_continuous = cfg_fixed.vars(loc_cont);
end

%% If global transformation is specified, extract that
if isfield(cfg_fixed, 'transform_global')
    global_transform = cfg_fixed.transform_global;
else
    global_transform = 'none';
end

%% For every continuous variable, get a list of local transformations
if sum(loc_cont) == 0
    FFX_transforms = {};
else
    FFX_transforms = extractTransforms(FFX_names(loc_cont), cfg_continuous);
end

%% Standardize vector transformations
vectorTransforms = {'center', 'centre', 'demean', 'std', 'standardize', ...
                    'normalize', 'logn', 'log10', 'inverseranknorm',    ...
                    'ranknorm', 'int', 'none'};
loc_vector = find(ismember(FFX_transforms, vectorTransforms));
if ~isempty(loc_vector)
    FFX_vectorTransforms = standardizeVectorTransforms(cfg_continuous(loc_vector));
else
    FFX_vectorTransforms = {};
end

%% Standardize winsorization
winsorizeTransforms = {'winsorize', 'winsorise', 'winsor'};
loc_winsorize = find(ismember(FFX_transforms, winsorizeTransforms));
if ~isempty(loc_winsorize)
    FFX_winsorize = standardizeWinsorizeTransforms(cfg_continuous(loc_winsorize));
else
    FFX_winsorize = {};
end

%% Standardize delta transformations
loc_delta = find(strcmpi(FFX_transforms, 'delta'));
if ~isempty(loc_delta)
    FFX_deltaTransforms = standardizeDeltaTransforms(cfg_continuous(loc_delta));
else
    FFX_deltaTransforms = {};
end

%% Standardize quadratic transformations
loc_quadratic = find(strcmpi(FFX_transforms, 'quadratic'));
if ~isempty(loc_quadratic)
    FFX_quadraticTransforms = standardizeQuadraticTransforms(cfg_continuous(loc_quadratic));
else
    FFX_quadraticTransforms = {};
end


%% Standardize spline transformations
loc_splines = find(strcmpi(FFX_transforms, 'splines'));
if ~isempty(loc_splines)
    FFX_splineTransforms = standardizeSplineTransforms(cfg_continuous(loc_splines));
else
    FFX_splineTransforms = {};
end

%% Standardize interactions
if isfield(cfg_fixed, 'interaction') && ~isempty(cfg_fixed.interaction)
    FFX_interactions = standardizeInteractions(cfg_fixed.interaction);
else
    FFX_interactions = {};
end

%% Make a logical vector of all variables that are of interest
if skipSteps
    names_of_interest = {''};
else
    locs_of_interest  = cellfun(@(x) x.of_interest, cfg_fixed.vars);
    names_of_interest = FFX_names(locs_of_interest);
    
    % Additionally pad in the interaction name of interest
    if ~isempty(FFX_interactions)
        names_of_interest = [names_of_interest; colvec({cfg_fixed.interaction([cfg_fixed.interaction(:).of_interest]).vars})];
    end
end

%% Update vector transforms with no transformation variables
% Variables that are present but have no defined transformation, add 'none'
% as the transformation type
alreadyTransformed = {};
if ~isempty(FFX_categorical)
    alreadyTransformed = [alreadyTransformed, {FFX_categorical.name}];
end
if ~isempty(FFX_vectorTransforms)
    alreadyTransformed = [alreadyTransformed, {FFX_vectorTransforms.name}];
end
if ~isempty(FFX_deltaTransforms)
    alreadyTransformed = [alreadyTransformed, {FFX_deltaTransforms.name}];
end
if ~isempty(FFX_quadraticTransforms)
    alreadyTransformed = [alreadyTransformed, {FFX_quadraticTransforms.name}];
end
if ~isempty(FFX_winsorize)
    alreadyTransformed = [alreadyTransformed, {FFX_winsorize.name}];
end
if ~isempty(FFX_splineTransforms)
    alreadyTransformed = [alreadyTransformed, {FFX_splineTransforms.name}];
end
remainingVariables = setdiff(FFX_names, alreadyTransformed);

% Loop over these, and grow vectorTransforms
count = length(FFX_vectorTransforms) + 1;
for v = 1:length(remainingVariables)
    FFX_vectorTransforms(count).name = remainingVariables{v};
    % if strcmpi(remainingVariables{v}, vars_of_interest)
    %     FFX_vectorTransforms(count).of_interest = true;
    % else
        FFX_vectorTransforms(count).of_interest = false;
    % end
    FFX_vectorTransforms(count).transformation = 'none';

    count = count + 1;
end

end


function output = standardizeCategorical(cfg, defReference)
% Output is a structure with fields: name, of_interest, reference level
output  = cell(length(cfg), 3);
outName = {'name', 'of_interest', 'reference'};

for lines = 1:length(cfg)
    output{lines,1} = cfg{lines}.name;
    output{lines,2} = cfg{lines}.of_interest;
    if isfield(cfg{lines}, 'reference')
        output{lines,3} = cfg{lines}.reference;
    else
        output{lines,3} = defReference;
    end
end

output = cell2struct(output', outName);
end


function output = extractTransforms(allVarNames, inputCellStruct)
% inputCellStruct is a cell array with each entry containing a structure
% that may or may not have a field named "transformation"; not all values
% in allVarNames will be present in inputCellStruct - the output is all
% transformations corresponding to all entries in allVarNames (i.e., if an
% entry in allVarNames is not specified in inputCellStruct, the
% transformation is set to none)
%
% Throws an error if allVarNames has duplicates OR if more than one
% transform is specified for any variable

% Initialize
output = cell(length(allVarNames), 1);

% Go over every entry and extract transform
for vars = 1:length(inputCellStruct)
    tmp  = inputCellStruct{vars};
    loc  = find(strcmpi(allVarNames, tmp.name));
    if isempty(loc)
        % Transformation defined but variable unavailable in allVarNames
        warning(['Ignoring transformation for ', tmp.name, ' as variable not present']);
    else
        if isfield(tmp, 'transform')
            if isempty(output{loc})
                output{loc} = lower(inputCellStruct{vars}.transform);
            else
                % A transformation was already defined for this variable
                error(['More than one transformation defined for: ', inputCellStruct{vars}.name]);
            end
        else
            % Transformation not defined, check if duplicate, and assign
            % none as a transform
            if isempty(output{loc})
                output{loc} = 'none';
            else
                % A transformation was already defined for this variable
                error(['More than one transformation defined for: ', inputCellStruct{vars}.name]);
            end
        end
    end
end

% Any unassigned variable: variable does not exist in the spec file
output(cellfun(@isempty, output)) = {'none'};
end


function output = standardizeVectorTransforms(cfg)
% Output is a structure with fields: name, of_interest, transformation
output  = cell(length(cfg), 3);
outName = {'name', 'of_interest', 'transformation'};

for lines = 1:length(cfg)
    output{lines,1} = cfg{lines}.name;
    output{lines,2} = cfg{lines}.of_interest;
    if isfield(cfg{lines}, 'transform')
        output{lines,3} = cfg{lines}.transform;
    else
        output{lines,3} = 'none';
    end
end

output = cell2struct(output', outName);
end


function output = standardizeWinsorizeTransforms(cfg)
% Output is a structure with fields: name, of_interest, transformation, lower, upper

% First pass, split into variables using comma
output  = cell(length(cfg), 5);
outName = {'name', 'of_interest', 'transformation', 'lower', 'upper'};

for lines = 1:length(cfg)
    output{lines,1} = cfg{lines}.name;
    output{lines,2} = cfg{lines}.of_interest;
    output{lines,3} = 'winsorize';

    % Extract minmax field and sort
    tmp = cfg{lines}.winsorize.minmax;
    if isempty(tmp)
        warning('No winsorization bounds set; no winsorization applied');
        output{lines,4} = 0;
        output{lines,5} = 100;
    else
        tmp = sort(tmp);
        output{lines,4} = tmp(1);
        output{lines,5} = tmp(2);
    end
end
output = cell2struct(output', outName);
end


function output = standardizeDeltaTransforms(cfg)
% Output is a structure with fields: name, of_interest, transformation
output  = cell(length(cfg), 3);
outName = {'name', 'of_interest', 'transformation'};

for lines = 1:length(cfg)
    output{lines,1} = cfg{lines}.name;
    output{lines,2} = cfg{lines}.of_interest;
    output{lines,3} = 'delta';
end
output = cell2struct(output', outName);
end


function output = standardizeQuadraticTransforms(cfg)
% Output is a structure with fields: name, of_interest, transformation
output  = cell(length(cfg), 3);
outName = {'name', 'of_interest', 'transformation'};

for lines = 1:length(cfg)
    output{lines,1} = cfg{lines}.name;
    output{lines,2} = cfg{lines}.of_interest;
    output{lines,3} = 'quadratic';
end
output = cell2struct(output', outName);
end


function output = standardizeSplineTransforms(cfg)
% Output is a structure with fields: name, of_interest, knots, splineType,
%                                    Xpowers, method, minMax, interecpt,
%                                    optCommand, optAppend, cleanUp, instance
outName = {'name', 'of_interest', 'knots', 'splineType', 'Xpowers', 'method', ...
           'minMax', 'intercept', 'optCommand', 'optAppend', 'cleanUp', 'instance'};
output  = cell(length(cfg), length(outName));

% Default settings
def_knots       = {'quartiles'};
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
    output{lines,1} = cfg{lines}.name;

    ff = fieldnames(cfg{lines}.splines);
    for f = 1:length(ff)
        loc = find(strcmpi(outName, ff{f}));
        if ~isempty(loc)
            output{lines,loc} = cfg{lines}.splines.(ff{f});
        end
    end
end

output = cell2struct(output', outName);
end


function output = standardizeInteractions(cfg)
% Output is a structure with fields: interaction, of_interest

% Initialize
output  = cell(size(cfg,1), 2);
outName = {'interaction', 'of_interest'};

output(:,1) = {cfg(:).vars}';
output(:,2) = {cfg(:).of_interest}';

output = cell2struct(output', outName);
end