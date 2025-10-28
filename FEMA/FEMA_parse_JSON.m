function [FFX_names, FFX_categorical, FFX_vectorTransforms, ...
          FFX_deltaTransforms, FFX_splineTransforms] = FEMA_parse_JSON(configFile)
% Function that reads a DEAP-created JSON specification files and parses
% the fixed effects information out of it
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
%                               * transformation:   which transform to apply
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
%                               * instance:         1

% Decrpy JSON file
if ~isstruct(configFile)
    configFile = jsondecode(fileread(configFile));
end 

% Extract all fixed effects
cfg_fixed = configFile.params.fixed;
FFX_names = cellfun(@(x) x.name, cfg_fixed.vars, 'UniformOutput', false);

% Determine categorical and continuous
loc_cont = cellfun(@(x) strcmpi(x.type_var, 'continuous'), cfg_fixed.vars);
loc_catg = ~loc_cont;

% Convert categorical into a standard structure format
FFX_categorical = standardizeCategorical(cfg_fixed.vars(loc_catg), 'mode');

% Subset continuous variables
cfg_continuous = cfg_fixed.vars(loc_cont);
% if global transform specified apply to all continuous vars, error if 
% global and local transform specified
if isfield(cfg_fixed, 'transform_global')
    for vars = 1:length(cfg_continuous)
        if ~isfield(cfg_continuous{vars}, 'transform') || ...
            isempty(cfg_continuous{vars}.transform) || ...
            strcmpi(cfg_continuous{vars}.transform, 'none')
            cfg_continuous{vars}.transform = cfg_fixed.transform_global;
        else
            error('Error: cannot have two transformations for variable %s: %s and %s', ...
                  cfg_continuous{vars}.name, cfg_continuous{vars}.transform, cfg_fixed.transform_global);
        end
    end
end

% For every continuous variable, get a list of transformations
% need to include  global transform if it exists
FFX_transforms = extractTransforms(FFX_names(loc_cont), cfg_continuous);

% Standardize vector transformations
vectorTransforms = {'center', 'centre', 'demean', 'std', 'standardize', ...
                    'normalize', 'logn', 'log10', 'inverseranknorm',    ...
                    'ranknorm', 'int', 'none'};
loc_vector = find(ismember(FFX_transforms, vectorTransforms));
if ~isempty(loc_vector)
    FFX_vectorTransforms = standardizeVectorTransforms(cfg_continuous(loc_vector));
else
    FFX_vectorTransforms = [];
end

% Standardize delta transformations
loc_delta = find(strcmpi(FFX_transforms, 'delta'), 1);
if ~isempty(loc_delta)
    FFX_deltaTransforms = standardizeDeltaTransforms(cfg_continuous(loc_delta));
else
    FFX_deltaTransforms = [];
end

% Standardize spline transformations
loc_splines = find(strcmpi(FFX_transforms, 'splines'));
if ~isempty(loc_splines)
    FFX_splineTransforms = standardizeSplineTransforms(cfg_continuous(loc_splines));
else
    FFX_splineTransforms = [];
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
    % output(1:length(allVarNames),1) = deal({'none'});
    if isfield(inputCellStruct{1}, 'global_transform')
        output(1:length(allVarNames),1) = deal({lower(inputCellStruct.global_transform)});
    end 

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
        if isfield(cfg{lines}, 'transformation')
            output{lines,3} = cfg{lines}.transformation;
        else
            output{lines,3} = 'none';
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

function output = standardizeSplineTransforms(cfg)
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
        output{lines,1} = cfg{lines}.name;
        output{lines,2} = cfg{lines}.of_interest;

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