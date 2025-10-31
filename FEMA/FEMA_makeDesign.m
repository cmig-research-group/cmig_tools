function [designMatrix, vars_of_interest] = FEMA_makeDesign(configFile, varargin)
% Function that makes FEMA-compatible design matrix for ABCD data
%% Inputs:
% configFile:       character       full path to a configuration file
%                                   (see example configuration file in the
%                                   recipe folder), can either be a JSON or 
%                                   text file saved with the extension '.config'
%  
%% Optional inputs (name-pair values): 
% iid:              cell / char     individual IDs for all observations
%                                   that need to be retained OR full path
%                                   to a csv file with no header that has
%                                   individual IDs that need to be retained
% 
% eid:              cell / char     event IDs for all observations that
%                                   need to be retained OR full path to a
%                                   csv file with no header that has event
%                                   IDs that need to be retained
% 
% dataFile:         character       full path to a .csv/.tsv/.parquet/
%                                   file that has the variables specified
%                                   in configFile to read in / additional
%                                   variables to be included in the design matrix
%                                       - if .csv/.tsv/.parquet, the columns should be named as (in order):
%                                           * participant_id: the individual ID for every observation
%                                           * session_id:     the event ID for every observation
%                                           * column 3 onwards should be the variables to include
%                                           * if the `familyIDvar` specified in configFile exists, 
%                                             in the dataFile, this is respected; otherwise, an attempt
%                                             is made to read the family ID from the demographics table
%                                           * if a column named `agevec` exists, it is not 
%                                             considered to be a part of the design matrix
% 
%                                   N.B: if agevec is missing in dataFile,
%                                   it is assumed to be "0" throughout
%                                   (since it is not directly used for analysis)
% 
% dirTabulated:     character       full path to where the ABCD 6.0 onwards
%                                   tabulated data is saved
% 
% dropMissing:      logical         true or false indicating if missing
%                                   values in any of the variables should
%                                   lead to that observation being deleted
%                                   from the design matrix or not (default:
%                                   true)
% 
% isFE:             logical         true or false indicating if the
%                                   RandomEffects model type is F and E (in
%                                   which case iid = fid)
% 
% isSE:             logical         true or false indicating if the
%                                   RandomEffects model type is S and E (in
%                                   which case fid = iid)
% 
% outDir:           character       full path to where the design matrix
%                                   should be written out
% 
% outName:          character       name of the design matrix (if not
%                                   specified, defaults to:
%                                   DesignMatrix-yyyyMMMdd-HHmmSS
% 
% outType:          character/cell  one or more additional output formats
%                                   for the design matrix:
%                                       * 'csv'
%                                       * 'mat'
%                                       * 'compiled' | 'standalone'
%
%% To Dos
% Need additional checks for when this function will be called in compiled mode
% How to handle non-standard family IDs? For example, scanner ID?
% Need to read data if anything specified in configFile is not found in
% dataFile
% dataFile specification for mat file is really complicated - not currently
% handled (see end of file)
% Currently using session_id for sorting of events which is not a
% generalisable solution as session_id may be incorrectly sorted (if
% non-numeric); need to incorporate age-based sorting OR allow user to
% specify another flag for sorting of events for delta computation
% DP: 
% i think we always should save an output file so writeDesignMat = true and 
%     if no output directory is given either save to pwd or ~ 
%

%% Check mandatory inputs
if ~exist('configFile', 'var') || isempty(configFile)
    error('Please provide full path to a config file');
else
    if ~exist(configFile, 'file')
        error(['Unable to find: ', configFile]);
    end
end

%% Parse optional inputs
p = inputParser;
addParameter(p, 'iid',          '');
addParameter(p, 'eid',          '');
addParameter(p, 'dataFile',     '');
addParameter(p, 'dirTabulated', '');
addParameter(p, 'isFE',         false);
addParameter(p, 'isSE',         false);
addParameter(p, 'dropMissing',  true);
addParameter(p, 'outDir',       pwd);
addParameter(p, 'outName',      '');
addParameter(p, 'outType',      '');

parse(p, varargin{:})
iid          = p.Results.iid;
eid          = p.Results.eid;
dataFile     = p.Results.dataFile;
dirTabulated = p.Results.dirTabulated;
isFE         = p.Results.isFE;
isSE         = p.Results.isSE;
dropMissing  = p.Results.dropMissing;
outDir       = p.Results.outDir;
outName      = p.Results.outName;
outType      = p.Results.outType;

%% Make some decision based on inputs
if isempty(outName)
    writeDesignMat = false;
else
    if ~exist(outDir, 'dir')
        mkdir(outDir);
        writeDesignMat = true;
    end
end

if writeDesignMat
    if isempty(outName)
        outName = ['DesignMatrix-', char(datetime('now', 'Format', 'yyyyMMMdd-HHmmSS'))];
    end
end

%% Parse the configuration file
[~, ~, tmp_ext] = fileparts(configFile);
if strcmpi(tmp_ext, '.json')
    isJSON = true;
    [FFX_names, FFX_categorical, FFX_vectorTransforms, FFX_winsorize,  ...
     FFX_deltaTransforms, FFX_splineTransforms, FFX_interactions,      ...
     global_transform, global_intercept, loc_cont, names_of_interest] = ...
     FEMA_parse_JSON(configFile);
else
    isJSON = false;
    [FFX_names, FFX_categorical, FFX_vectorTransforms, FFX_winsorize, ...
     FFX_deltaTransforms, FFX_splineTransforms, FFX_interactions,     ...
     global_transform, global_intercept, loc_cont, family_id] = FEMA_parse_config(configFile);
end

% Make a copy of the fixed effects name
temp_FFXNames = FFX_names;

%% Determine family ID
if isJSON
    if isSE
        familyIDvar = 'participant_id';
    else
        familyIDvar = 'family_id';
    end
else
    if isempty(family_id)
        if isSE
            familyIDvar = 'participant_id';
        else
            familyIDvar = 'ab_g_stc__design_id__fam';
        end
        warning(['familyIDvar not specified in the configFile; using: ', familyIDvar, ' as family ID']);
    else
        familyIDvar = family_id;
    end
end

%% Determine if tables need to be read
% If the input configuration file was a JSON file, dataFile should exist
if isJSON
    readTabulated = false;
    if isempty(dataFile)
        error('dataFile is required if the configuration file is a JSON specification');
    else
        data = parquetread(dataFile);
        exists_agevec = false;
    end
else
    if isempty(dataFile)
        readTabulated = true;
    else
        readTabulated = false;
        if ~exist(dataFile, 'file')
            error(['Unable to find: ', dataFile]);
        else
            % User has provided a full path to a dataFile: could be
            % .csv/.tsv./.parquet/ file
            [~, ~, ext] = fileparts(dataFile);
            switch ext
                case {'.csv', '.tsv'}
                    data = readtable(dataFile, 'FileType', 'text');
    
                case '.parquet'
                    data = parquetread(dataFile);
    
                otherwise
                    error(['Unsupported dataFile format: ', ext]);
            end
    
            % Check if the read in data follows the specification
            toCheck    = {'participant_id', 'session_id'};
            temp_names = data.Properties.VariableNames;
            chk        = ismember(toCheck, temp_names);
            if sum(chk) ~= length(toCheck)
                missMsg = sprintf('%s, ', toCheck{~chk});
                error(['The following required variables are missing from ', file_X, ': ', missMsg(1:end-2)]);
            else
                % Prepare an ID column to be used for merging
                data.IDs_merge = strcat(data.participant_id, {'_'}, data.session_id);
            end
    
            % Check if familyIDvar exists in data
            if ~any(strcmpi(data.Properties.VariableNames, familyIDvar))
                % Check if a standard familyIDvar is requested
                if ~ismember(familyIDvar, {'ab_g_stc__design_id__fam', 'ab_g_stc__design_id__fam__gen', 'participant_id'})
                    tmp = familyIDvar;
                    if isSE
                        familyIDvar = 'participant_id';
                    else
                        familyIDvar = 'ab_g_stc__design_id__fam';
                    end
                    warning(['Unknown value specified for familyIDvar: ', tmp, '; defaulting to: ', familyIDvar]);
                    readTabulated = true;
                end
            end
    
            % Check if agevec exists in data
            if any(strcmpi(data.Properties.VariableNames, 'agevec'))
                exists_agevec = true;
            else
                exists_agevec = false;
            end
    
            % Check if all fixed effects already exist?
            vars_toRead = setdiff(temp_FFXNames, data.Properties.VariableNames);
            if isempty(vars_toRead)
                % Empty out temp_FFXNames
                temp_FFXNames = [];
            else
                readTabulated = true;
    
                % Overwrite temp_FFXNames to remaining variables to read
                temp_FFXNames = vars_toRead;
            end
        end
    end
end

if readTabulated
    if isempty(dirTabulated)
        error(['Please provide a full path to where the ABCD tabulated data is located ', ...
               'OR provide a full path to a saved dataFile']);
    else
        if ~exist(dirTabulated, 'dir')
            error(['Unable to find: ', dirTabulated]);
        end
    end
end

%% Read data if required
if readTabulated
    if ~isempty(vars_toRead)

        % Make a list of tables to read
        temp_tabNames = cellfun(@(x) x{1,1}, cellfun(@(x) strsplit(x, '__'), ...
                                temp_FFXNames, 'UniformOutput', false),      ...
                                'UniformOutput', false);
        tables_toRead = fullfile(dirTabulated, strcat(unique(temp_tabNames), 'stable'), '.parquet');

        % Make sure these tables exist
        chk_tables = not(cellfun(@(x) exist(x, 'file'), tables_toRead));
        if any(chk_tables)
            disp('Following tables could not be found: ');
            disp(tables_toRead(chk_tables));
            error('Cannot proceed with making design matrix as one or more tables could not be found');
        end
    end

    %% Read static and dynamic tables and put them together
    data_stc  = parquetread(fullfile(dirTabulated, 'ab_g_stc.parquet'));
    data_dyn  = parquetread(fullfile(dirTabulated, 'ab_g_dyn.parquet'));
    data_demo = outerjoin(data_dyn, data_stc, 'Keys', 'participant_id', 'MergeKeys', true);

    % Mandatory variables in data_demo to retain
    vars_always_retain = {'participant_id', 'session_id', 'ab_g_dyn__visit_age'};

    %% Put IDs together for merging
    % Only use iid_eid so that it is easier to join with other tables
    % Make sure to use participant_id for iid instead of subjectIDvar
    data_demo.IDs_merge = strcat(data_demo.participant_id, {'_'}, data_demo.session_id);

    %% Are there variables that we need from these two tables?
    tmp_demo_tab_locs = ismember(temp_tabNames, {'ab_g_stc', 'ab_g_dyn'});
    vars_demographics = temp_FFXNames(tmp_demo_tab_locs);

    % Subset demo
    data_demo = data_demo(:, [vars_always_retain, vars_demographics]);

    % What are the remaining tables?
    temp_tabNames = temp_tabNames(~tmp_demo_tab_locs);
    temp_FFXNames = temp_FFXNames(~tmp_demo_tab_locs);
    tables_toRead = tables_toRead(~tmp_demo_tab_locs);

    % Loop over remaining tables - innerjoin with data_demo; innerjoin will
    % keep dropping individuals who are not found across the board
    % Use participant_id and session_id for merging
    tables_toLoop = unique(temp_tabNames, 'stable');
    vars_toRead   = {'participant_id', 'session_id'};
    for tab = 1:length(tables_toLoop)
        % Which table to read?
        temp = parquetinfo(tables_toRead{tab});

        % Do we have all the variables in the table?
        tmp_chk = ismember(temp_FFXNames{tab}, temp.VariableNames);
        if sum(tmp_chk) ~= length(temp_FFXNames{tab})
            wch = strcat(temp_FFXNames{tab}(~tmp_chk), {', '});
            wch = horzcat(wch{:});
            error(['Could not find the following variables in table ', ...
                   tables_toRead{tab}, ': ', wch(1:end-2)]);
        else
            temp = parquetread(tables_toRead{tab}, 'SelectedVariableNames', ...
                               [vars_toRead, temp_FFXNames{tab}]);
            temp.IDs_merge = strcat(temp.participant_id,  {'_'}, temp.session_id);
            data_demo = innerjoin(data_demo, temp, 'Keys', 'IDs_merge', 'MergeKeys', true);
        end
    end
end

%% Merge data_demo with dataFile: make data_work table
if exist('data_demo', 'var') && exist('data', 'var')
    data_work = innerjoin(data, data_demo, 'Keys', 'IDs_merge', 'MergeKeys', true);
else
    if exist('data_demo', 'var')
        data_work = data_demo;
    else
        data_work = data;
    end
end

%% Make a copy of IID_EID to intersect with the left hand side
if strcmpi('IDs_merge', data_work.Properties.VariableNames)
    IDs_merge = data_work.IDs_merge;
    removevars(data_work, 'IDs_merge');
else
    IDs_merge = strcat(data_work.iid, {'_'}, data.eid);
end

%% If the user has provided iid_eid, filter data_work
iid_eid   = strcat(iid, {'_'}, eid);
[~, ia]   = intersect(IDs_merge, iid_eid);
data_work = data_work(ia, :);

%% Determine subject ID
if isFE
    subjectIDvar = familyIDvar;
else
    subjectIDvar = 'participant_id';
end

%% Assign agevec if required
if ~exists_agevec
    if readTabulated
        % If we read the data, we have agevec from demographics
        data_work.agevec = data_work.ab_g_dyn__visit_age;
    else
        % We do not have demographics table - assign zeros
        data_work.agevec = zeros(height(data_work), 1);
    end
end

%% Reorder so that the data is organised as fid, iid, eid, agevec
vars_to_move  = {familyIDvar, subjectIDvar, 'session_id', 'agevec'};
remainingVars = setdiff(data_work.Properties.VariableNames, vars_to_move);
data_work     = data_work(:, [vars_to_move, remainingVars]);

%% Renaming variables to fid, iid, and eid
data_work = renamevars(data_work, {familyIDvar, subjectIDvar, 'session_id'}, {'fid', 'iid', 'eid'});

%% At this stage, remove any missing (if required)
if dropMissing
    miss_indicators = {'n/a', 'NaN', 'NaT', '<undefined>', ' ', 'NA'};
    loc_missing     = logical(sum(ismissing(data_work, miss_indicators), 2));
    data_work(loc_missing, :) = [];
end

%% Initialize transformations
all_catg_variables = cell(length(FFX_categorical), 2);
all_cont_variables = cell(sum(loc_cont), 2);
count              = 1;

% Track mapping of original names and modified names
% Original name, new name
names_mapping = cell(length(FFX_names), 2);
count_mapping = 1;

%% Dummy code categorical variables and drop reference levels
if ~isempty(FFX_categorical)
    for cc = 1:length(FFX_categorical)
        % Which variable are we looking at?
        tmp_var_name = FFX_categorical(cc).name;

        % What is the specified reference level
        tmp_ref_level = FFX_categorical(cc).reference;

        % Create dummy coding
        [all_catg_variables{cc,1}, all_catg_variables{cc,2}] = ...
         expand_categorical(data_work.(tmp_var_name), tmp_ref_level, tmp_var_name);

        % Maintain mapping
        names_mapping{count_mapping, 1} = tmp_var_name;
        names_mapping{count_mapping, 2} = all_catg_variables{cc,1};
        count_mapping = count_mapping + 1;
    end

    % Put all categorical variables together
    X_name_categorical = horzcat(all_catg_variables{:,1});
    X_data_categorical = horzcat(all_catg_variables{:,2});
else
    X_data_categorical = [];
    X_name_categorical = '';
end

%% Perform vector transformations on continuous variables
if ~isempty(FFX_vectorTransforms)
    % First pass, handle transforms that generate a vector - these can be
    % handled together in a single iteration without a need for loop
    % Examples of these transforms are (all handled by doTransformation):
    %   * 'center' | 'centre' | 'demean'
    %   * 'std' | 'standardize' | 'normalize'
    %   * 'logn'
    %   * 'log10'
    %   * 'inverseranknorm' | 'ranknorm' | 'int'

    % Handle centering
    loc_center = find(ismember({FFX_vectorTransforms(:).transformation}, {'center', 'centre', 'demean'}));
    if ~isempty(loc_center)
        ori_names  = {FFX_vectorTransforms(loc_center).name};
        [all_cont_variables{count,1}, all_cont_variables{count,2}] = ...
            applyTransforms(data_work, ori_names, 'center');
    
        % Maintain mapping
        names_mapping(count_mapping:count_mapping+length(ori_names)-1, 1) = ori_names;
        names_mapping(count_mapping:count_mapping+length(ori_names)-1, 2) = all_cont_variables{count,1};
        count_mapping = count_mapping + 1;
        count = count + 1;
    end

    % Handle standardization
    loc_std = find(ismember({FFX_vectorTransforms(:).transformation}, {'std', 'standardize', 'normalize'}));
    if ~isempty(loc_std)
        ori_names  = {FFX_vectorTransforms(loc_std).name};
        [all_cont_variables{count,1}, all_cont_variables{count,2}] = ...
            applyTransforms(data_work, ori_names, 'std');
    
        % Maintain mapping
        names_mapping(count_mapping:count_mapping+length(ori_names)-1, 1) = ori_names;
        names_mapping(count_mapping:count_mapping+length(ori_names)-1, 2) = all_cont_variables{count,1};
        count_mapping = count_mapping + 1;
        count = count + 1;
    end

    % Handle logn
    loc_logn = find(ismember({FFX_vectorTransforms(:).transformation}, {'logn'}));
    if ~isempty(loc_logn)
        ori_names  = {FFX_vectorTransforms(loc_logn).name};
        [all_cont_variables{count,1}, all_cont_variables{count,2}] = ...
            applyTransforms(data_work, ori_names, 'logn');
    
        % Maintain mapping
        names_mapping(count_mapping:count_mapping+length(ori_names)-1, 1) = ori_names;
        names_mapping(count_mapping:count_mapping+length(ori_names)-1, 2) = all_cont_variables{count,1};
        count_mapping = count_mapping + 1;
        count = count + 1;
    end

    % Handle log10
    loc_log10 = find(ismember({FFX_vectorTransforms(:).transformation}, {'log10'}));
    if ~isempty(loc_log10)
        ori_names  = {FFX_vectorTransforms(loc_log10).name};
        [all_cont_variables{count,1}, all_cont_variables{count,2}] = ...
            applyTransforms(data_work, ori_names, 'log10');
    
        % Maintain mapping
        names_mapping(count_mapping:count_mapping+length(ori_names)-1, 1) = ori_names;
        names_mapping(count_mapping:count_mapping+length(ori_names)-1, 2) = all_cont_variables{count,1};
        count_mapping = count_mapping + 1;
        count = count + 1;
    end

    % Handle inverse rank norm
    loc_INT = find(ismember({FFX_vectorTransforms(:).transformation}, {'inverseranknorm', 'ranknorm', 'int'}));
    if ~isempty(loc_INT)
        ori_names  = {FFX_vectorTransforms(loc_INT).name};
        [all_cont_variables{count,1}, all_cont_variables{count,2}] = ...
            applyTransforms(data_work, ori_names, 'ranknorm');
    
        % Maintain mapping
        names_mapping(count_mapping:count_mapping+length(ori_names)-1, 1) = ori_names;
        names_mapping(count_mapping:count_mapping+length(ori_names)-1, 2) = all_cont_variables{count,1};
        count_mapping = count_mapping + 1;
        count = count + 1;
    end
end

%% Perform winsorization
if ~isempty(FFX_winsorize)
    ori_names  = {FFX_winsorize(:).name};
    [all_cont_variables{count,1}, all_cont_variables{count,2}] = ...
        applyWinsorization(data_work, ori_names, [FFX_winsorize(:).lower], [FFX_winsorize(:).upper]);

    % Maintain mapping
    names_mapping(count_mapping:count_mapping+length(ori_names)-1, 1) = ori_names;
    names_mapping(count_mapping:count_mapping+length(ori_names)-1, 2) = all_cont_variables{count,1};
    count_mapping = count_mapping + 1;
    count = count + 1;
end

%% Perform delta computation
if ~isempty(FFX_deltaTransforms)
    ori_names  = {FFX_deltaTransforms(:).name};
    [all_cont_variables{count,1}, all_cont_variables{count,2}] = computeDelta(data_work{:,ori_names}, ori_names, data_work.participant_id, data_work.session_id);

    % Maintain mapping
    names_mapping(count_mapping:count_mapping+length(ori_names)-1, 1) = ori_names;
    names_mapping(count_mapping:count_mapping+length(ori_names)-1, 2) = all_cont_variables{count,1};
    count_mapping = count_mapping + 1;
    count = count + 1;
end

%% Perform spline creation
if ~isempty(FFX_splineTransforms)
    % Initialize a variable for storing spline settings and other relevant
    % variables that are not part of design matrix but are necessary
    all_spline_settings = cell(size(FFX_splineTransforms,1), 1);

    % Go over every spline index
    for v = 1:size(FFX_splineTransforms,1)
        % Variable name
        tmp_var_name = FFX_splineTransforms(v).name;

        % Create basis functions
        [all_cont_variables{count,2}, all_spline_settings{v,1},   ...
         all_spline_settings{v,2},    all_spline_settings{v,3},   ...
         all_spline_settings{v,4},    all_spline_settings{v,5}] = ...
         createBasisFunctions(data_work.(tmp_var_name), FFX_splineTransforms(v).knots, ...
                              FFX_splineTransforms(v).splineType, ...
                              FFX_splineTransforms(v).Xpowers, ...
                              FFX_splineTransforms(v).method, outDir, ...
                              FFX_splineTransforms(v).minMax, ...
                              FFX_splineTransforms(v).intercept, ...
                              FFX_splineTransforms(v).optCommand, ...
                              FFX_splineTransforms(v).optAppend, ...
                              FFX_splineTransforms(v).cleanUp, FFX_splineTransforms(v).instance);

        % Append the name of the variable for which settings are saved
        all_spline_settings{v,6} = tmp_var_name;
    
        % Name of basis functions
        tmp_new_name = strcat({tmp_var_name}, {'_bf'}, num2str((1:size(all_cont_variables{count,2},2))'))';
        all_cont_variables{count,1} = tmp_new_name;

        % Maintain mapping
        names_mapping{count_mapping, 1} = tmp_var_name;
        names_mapping{count_mapping, 2} = tmp_new_name;
        count_mapping = count_mapping + 1;
        count = count + 1;
    end
end

%% Variables that were not transformed
vars_none_transform = {FFX_vectorTransforms(strcmpi({FFX_vectorTransforms(:).transformation}, 'none')').name};
if ~isempty(vars_none_transform)
    all_cont_variables{count,1} = vars_none_transform;
    all_cont_variables{count,2} = data_work{:,vars_none_transform};

    % Maintain mapping
    names_mapping(count_mapping:count_mapping+length(vars_none_transform)-1, 1) = vars_none_transform;
    names_mapping(count_mapping:count_mapping+length(vars_none_transform)-1, 2) = vars_none_transform;
    % count_mapping = count_mapping + 1;
    % count = count + 1;
end

%% Shrink all_cont_variables, if required
tmp = cellfun(@isempty, all_cont_variables);
all_cont_variables(tmp,:) = [];
    
%% Put all continuous variables together
if ~isempty(all_cont_variables)
    X_name_continuous = horzcat(all_cont_variables{:,1});
    X_data_continuous = horzcat(all_cont_variables{:,2});
else
    X_data_continuous = [];
    X_name_continuous = '';
end

%% Put categorical and continnuous variables together
if ~isempty(X_data_continuous) & ~isempty(X_data_categorical)
    designMatrix = [X_data_continuous, X_data_categorical];
    designNames  = [X_name_continuous, X_name_categorical];
else
    if isempty(X_data_continuous)
        designMatrix = X_data_categorical;
        designNames  = X_name_categorical;
    else
        if isempty(X_data_categorical)
            designMatrix = X_data_continuous;
            designNames  = X_name_continuous;
        else
            % Both cannot be empty; throw an error
            error('Something went wrong during design matrix creation; no variables remain');
        end
    end
end

%% Make a table
designMatrix = [data_work(:,1:4), cell2table(num2cell(designMatrix), 'VariableNames', designNames)];

%% Handle interactions
if ~isempty(FFX_interactions)

    % Initialize
    all_interact_variables = cell(size(FFX_interactions,1), 2);

    % For every interaction, split by "*", and create new variables
    for v = 1:size(FFX_interactions,1)
        % Variables to look for in this interaction
        toLook = strtrim(strsplit(FFX_interactions(v).interaction, '*'));

        % Using regular expressions, look for the entries in "toLook" in
        % the first column of the mapping variable - this accounts for the
        % user wanting to do interactions with some form of transformed or
        % derived variables
        locs_look = find(ismember(names_mapping(:,1), toLook));

        % Initialize
        col_count = 1;

        % Assign the first variable into tempData
        tempName = names_mapping{locs_look(1), 2};
        tempData = designMatrix{:, tempName};
        
        % Which variables remain?
        remVariables = setdiff(toLook, names_mapping{locs_look(1), 1});

        % Go over all remaining variables
        for vv = 1:length(remVariables)

            % Which transformed variable names in names_mapping do we need to look up
            wch_vars = names_mapping{ismember(names_mapping(:,1), remVariables{vv}),2};

            % Whatever variables we picked from names_mapping, interacts
            % with tempData and is saved out in temp
            temp = designMatrix{:, wch_vars};

            % Initialize where results will be temporarily saved
            tempSize    = size(tempData,2) * size(temp,2);
            tempResData = zeros(size(tempData,1), tempSize);
            tempResName = cell(tempSize, 1);

            % Every column in tempData needs to be multiplied with every
            % column in temp and assigned to tempResData
            for c1 = 1:size(tempData, 2)
                for c2 = 1:size(temp, 2)
                    % Name of source column
                    sName = tempName{c1};

                    % Name of target column
                    tName = names_mapping{locs_look(2),2}{c2};

                    % Interaction
                    tempResData(:,col_count) = tempData(:,c1) .* temp(:,c2);

                    % Interaction variable name
                    tempResName{col_count,1} =  [sName, '__', tName];
                    col_count = col_count + 1;
                end
            end

            % Now tempResData becomes the new tempData; tempResName becomes
            % the new tempName; and col_count is reset
            tempData = tempResData;
            tempName = tempResName;
            col_count = 1;
        end
        % Interaction results are in tempData - extract
        all_interact_variables{v,1} = tempData;
        all_interact_variables{v,2} = tempName;
    end

    % Put all interactions together
    X_data_interactions = horzcat(all_interact_variables{:,1});
    X_name_interactions = horzcat(all_interact_variables{:,2});
else
    X_data_interactions = [];
    X_name_interactions = '';
end

% Concatenate with design matrix
if ~isempty(X_data_interactions)
    X_data_interactions = cell2table(num2cell(X_data_interactions), 'VariableNames', matlab.lang.makeValidName(X_name_interactions));
    designMatrix = [designMatrix, X_data_interactions];
end

%% Perform any global transformation that may be required
if ~isempty(global_transform)
    designMatrix{:,5:end} = doTransformation(designMatrix{:,5:end}, global_transform);
end

%% Add global intercept if needed
if global_intercept
    designMatrix.Intercept = ones(size(designMatrix,1), 1);
end

%% Drop rank deficient columns
mdl    = fitlm(designMatrix{:,5:end}, randn(size(designMatrix,1),1), 'Intercept', not(global_intercept));
toDrop = find(isnan(mdl.Coefficients.tStat))+4;
designMatrix(:,toDrop) = [];
% designNames(:, toDrop) = [];

%% Prepare a list of variables of interest
locs             = ismember(names_mapping(:,1), names_of_interest);
vars_of_interest = names_mapping(locs,2);
vars_of_interest = horzcat(vars_of_interest{:})';

%% Convert to table for output
% designMatrix = cell2table(designMatrix, 'VariableNames', designNames);
end


function [colnames, expVariable] = expand_categorical(variable, refLevel, varName)
% This function takes a single categorical variable, a reference level, and
% the name of the categorical variable and returns the dummy coded expanded
% categorical variable and new column names
% 
% Possible reference levels:
% 'mode':   in this case, the category with the largest count is dropped (default)
% 'first':  in this case, the first category is dropped
% 'last':   in this case, the last category is dropped
% 'none':   in this case, no category is dropped
% specific reference level: in this case, the specified category is droped

if iscellstr(refLevel) || iscell(refLevel) %#ok<ISCLSTR>
    refLevel = refLevel{1};
end

% First, make into a categorical variable
if iscellstr(variable) | isstring(variable)
    tmp_variable = categorical(variable);
else
    tmp_variable = categorical(cellstr(variable));
end

% Second, make a list of categories that exist
tmp_categories = categories(tmp_variable);

% Third, expand into dummy variable
tmp_dummy = dummyvar(tmp_variable);

% Fourth, figure out which one is the reference level column and drop
% that column from tmp_dummy; additionally, make unique colnames
% colname generation scheme:
% a) make unique strings out of category names
% b) convert unique strings into valid names (NB prior to R2025a, this will
%    result in maximum 63 character long name)
% c) concatenate original variable name as a prefix

switch refLevel
    case 'mode'
        [~, ~, c]   = mode(tmp_variable);
        tmp_ref     = setdiff(1:length(tmp_categories), find(ismember(tmp_categories, c{1}), 1));
        expVariable = tmp_dummy(:, tmp_ref);
        colnames    = strcat(varName, {'_'}, matlab.lang.makeValidName(...
                                             matlab.lang.makeUniqueStrings(tmp_categories(tmp_ref))));

    case 'first'
        expVariable = tmp_dummy(:, 2:end);
        colnames    = strcat(varName, {'_'}, matlab.lang.makeValidName(...
                                             matlab.lang.makeUniqueStrings(tmp_categories(2:end))));

    case 'last'
        expVariable = tmp_dummy(:, 1:end-1);
        colnames    = strcat(varName, {'_'}, matlab.lang.makeValidName(...
                                             matlab.lang.makeUniqueStrings(tmp_categories(1:end-1))));

    case 'none'
        expVariable = tmp_dummy;
        colnames    = strcat(varName, {'_'}, matlab.lang.makeValidName(...
                                             matlab.lang.makeUniqueStrings(tmp_categories)));

    otherwise
        tmp_ref = find(ismember(tmp_categories, refLevel), 1);
        if isempty(tmp_ref)
            warning(['Could not find reference level: ', refLevel, '; defaulting to using mode']);
            [~, ~, c] = mode(tmp_variable);
            tmp_ref   = setdiff(1:length(tmp_categories), find(ismember(tmp_categories, c{1}), 1));
        else
            tmp_ref   = setdiff(1:length(tmp_categories), tmp_ref);
        end
        expVariable   = tmp_dummy(:, tmp_ref);
        colnames      = strcat(varName, {'_'}, matlab.lang.makeValidName(...
                                               matlab.lang.makeUniqueStrings(tmp_categories(tmp_ref))));
end

% If colnames is a string, convert to cellstr
if isstring(colnames)
    colnames = cellstr(colnames);
end

% Return as [1 x n]
colnames = colnames';
end

function [outName, outValue] = applyTransforms(data_work, var_names, transformName)
outName  = strcat(var_names, {'_'}, transformName);
outValue = doTransformation(data_work{:,var_names}, transformName);
end

function [outName, outValue] = applyWinsorization(data_work, var_names, lower, upper)
outName  = strcat(var_names, {'_winsorized'});
outValue = doTransformation(data_work{:,var_names}, 'winsorize', lower, upper);
end


function [outName, outValue] = computeDelta(var_delta, var_names, participant_id, session_id)
% Function to compute delta of a specified variable:
% This creates two columns from a single input: the baseline column is the
% value at the first event (whatever is the first event for this subject),
% and the delta column is the difference in the value from the baseline -
% for different events (after the baseline), the deltas would be different;
% they are all put together in a single "delta" column
% Alternatively, instead of using session_id, could use age for deciding
% the correct event order: this is more reliable than assuming that
% session_id will always be sorted BUT requires age to be present in the
% data_work which may not always be the case
% Currently uses session_id
% 
% var_delta can be multiple columns

% Simpler but (likely) slower solution: loop over every unique
% participant_id, find the first session, and compute delta
uqSubjects = unique(participant_id, 'stable');

% Initialize
[~, numCols] = size(var_delta);
outValue     = zeros(size(var_delta, 1), 2*numCols);

% Output names
% This solution only works for R2019 onwards
% tmp = append(var_names, {'_baseline', '_delta'}');
% outName = tmp(:);
outName = cell(1, size(outValue,2));
count = 1;
for v = 1:length(var_names)
    outName(1,count:count+1) = strcat(var_names{v}, {'_baseline', '_delta'});
    count = count + 2;
end

% Loop over every subject in uqSubjects, sort by event ID to find the first
% real event, then do the diff
for subj = 1:length(uqSubjects)
    locs = strcmpi(participant_id, uqSubjects{subj});
    [~, b] = sort(session_id(locs), 'ascend');

    % There are no missing values to handle, since those were dropped
    % before
    tmp_data = var_delta(locs,:);
    baseline = tmp_data(b,:);
    delta    = tmp_data - baseline;

    % Assign to outVar
    outValue(locs,:) = [baseline, delta];

    % % Temporary data sorted by events
    % tmp_data_ori  = var_delta(locs);
    % tmp_data_sort = tmp_data_ori(b);
    % 
    % % Find the first value ignoring NaN
    % baseline = tmp_data_ori(find(~isnan(tmp_data_sort), 1));
    % 
    % % Assign to outVar
    % outVar(locs,1) = baseline;
    % outVar(locs,2) = baseline - tmp_data_ori;
end
end


% function output = standardizeTransforms(allVarNames, inputCellStruct)
% % inputCellStruct is a cell array with each entry containing a structure
% % that may or may not have a field named "transformation"; not all values
% % in allVarNames will be present in inputCellStruct - the output is all
% % transformations corresponding to all entries in allVarNames (i.e., if an
% % entry in allVarNames is not specified in inputCellStruct, the
% % transformation is set to none)
% %
% % Throws an error if allVarNames has duplicates OR if more than one
% % transform is specified for any variable
% 
% % Initialize
% output = cell(length(allVarNames), 1);
% % output(1:length(allVarNames),1) = deal({'none'});
% 
% % Go over every entry and extract transform
% for vars = 1:length(inputCellStruct)
%     tmp  = inputCellStruct{vars};
%     loc  = find(strcmpi(allVarNames, tmp.name));
%     if isempty(loc)
%         % Transformation defined but variable unavailable in allVarNames
%         warning(['Ignoring transformation for ', tmp.name, ' as variable not present']);
%     else
%         if isfield(tmp, 'transformation')
%             if isempty(output{loc})
%                 output{loc} = lower(inputCellStruct{vars}.transformation);
%             else
%                 % A transformation was already defined for this variable
%                 error(['More than one transformation defined for: ', inputCellStruct{vars}.name]);
%             end
%         else
%             % Transformation not defined, check if duplicate, and assign
%             % none as a transform
%             if isempty(output{loc})
%                 output{loc} = 'none';
%             else
%                 % A transformation was already defined for this variable
%                 error(['More than one transformation defined for: ', inputCellStruct{vars}.name]);
%             end
%         end
%     end
% end
% 
% % Any unassigned variable: variable does not exist in the spec file
% output(cellfun(@isempty, output)) = {'none'};
% end

% Deprecated:
% Could use unstacking trick: re-organises data with participant_id as the
% first column, and then the variable split across sessions
% N.B: relies on MATLAB to "internally" figure out which way to unstack:
% should this be a concern?
% This solution does not work since we cannot restack the table to have the
% same missing session pattern

% % Unstack table
% temp = unstack(data_work, var_delta, session_id, 'GroupingVariables', participant_id);
% 
% % Initialize
% outVar = zeros(height(data_work), 2);
% 
% % Fill out missing entries in the first session
% % Unfortunately, no "next column" type specification exists, requiring
% % extraction of variables and taking the transpose; no way to specify that
% % only the first (transposed) row needs to be filled in, so takes a bit
% % more time than otherwise
% temp2     = fillmissing(temp{:,2:end}', 'next');
% temp{:,2} = temp2(1,:)';
% 
% % Now compute delta
% temp{:,3:end} = temp{:,2} - temp{:,3:end};

%% Currently unsupported dataFile .mat format
% - if .mat, there are two allowed specifications:
%     # Version 1:
%     Table type: contains a single table
%     type variable that has the
%     following columns (in this order):
%     * participant_id: the individual ID for every observation
%     * event_id:       the event ID for every observation
%     * column 3 onwards should be the variables to include
%     * if the `fid` specified in configFile exists as a column,
%         this is respected; otherwise, an attempt is made to read
%         the family ID from the demographics table
%         * if a column named `agevec` exists, it is not
%             considered to be a part of the design matrix
% 
%             # Version 2:
%             Free form: contains different vectors which are a single table
%             type variable that has the
%             following columns (in this order):
%             * participant_id: the individual ID for every observation
% 
%             * agevec:     the age at every event for every observation
%             * X_cont:     continuous variables to include (numeric matrix)
%             * X_catg:     categorical variables to include (cellstr)
%             * cont_names: the names of each column in X_cont
%             * catg_names: the names of each column in X_catg

% case '.mat'
%     temp = load(dataFile);
%     % Check to make sure mandatory variables are present
%     % agevec is optional
%     toCheck = {'fid', 'iid', 'eid'};
%     chk     = ismember(toCheck, fieldnames(temp));
%     if sum(chk) ~= length(toCheck)
%         missMsg = sprintf('%s, ', toCheck{~chk});
%         error(['The following required variables are missing from ', dataFile, ': ', missMsg(1:end-2)]);
%     else
%         % Now check for X_cont, X_catg, and corresponding column
%         % name fields in temp
%         chk_X_cont = isfield(temp, 'X_cont');
%         chk_X_catg = isfield(temp, 'X_catg');
%         if chk_X_cont
%             if ~isfield(temp, 'cont_names')
%                 error('X_cont specified in dataFile but missing cont_names field');
%             else
%                 num_X_cont = size(X_cont,2);
%             end
%         else
%             num_X_cont = 0;
%         end
%         if chk_X_catg
%             if ~isfield(temp, 'catg_names')
%                 error('X_catg specified in dataFile but missing catg_names field');
%             else
%                 num_X_catg = size(X_cont,2);
%             end
%         else
%             num_X_catg = 0;
%         end
%         num_X = num_X_cont + num_X_catg;
%         if num_X == 0
%             error('Missing X_cont and X_catg in the dataFile');
%         end
%
%         % Make a table out of it
%         table_dataFile      = cell(length(temp.fid), num_X + 4);
%         table_dataFile(:,1) = temp.fid;
%         table_dataFile(:,2) = temp.iid;
%         table_dataFile(:,3) = temp.eid;
%         if isfield(temp, 'agevec')
%             table_dataFile(:,4) = temp.agevec;
%         else
%             table_dataFile(:,4) = zeros(length(temp.fid), 1);
%         end
%
%         % Assign X variables
%         count = 5;
%         if num_X_cont > 0
%             table_dataFile(:,count:count+num_X_cont) = num2cell(temp.X_cont);
%             count = count + num_X_cont + 1;
%         end
%         if num_X_catg > 0
%             table_dataFile(:,count:count+num_X_catg) = temp.X_cont;
%         end
%
%         % Make table variable names
%         table_varNames = {'fid', 'iid', 'eid', 'agevec'};
%         if num_X_cont > 0
%             table_varNames = [table_varNames, rowvec(temp.cont_names)];
%         end
%         if num_X_catg > 0
%             table_varNames = [table_varNames, rowvec(temp.catg_names)];
%         end
%
%         % Convert to table
%         table_dataFile = cell2table(table_dataFile, 'VariableNames', table_varNames);
%
%         % IDs to use for merging with data_demo
%         table_dataFile.IDs_merge = strcat(table_dataFile.iid, {'_'}, table_dataFile.eid);
%     end

% Older deprecated code:
% if ~isempty(names_continuous)
%     all_cont_variables = cell(length(names_continuous), 2);
% 
% 
% 
%     % Use a counter for accessing all_cont_variables: shrink later
%     count = 1;
% 
%     % List of possible transforms (including spelling variations and
%     % shortened forms)
%     allowed_transforms = {'center', 'centre', 'demean', 'std', 'standardize', ...
%                           'normalize', 'logn', 'log10', 'inverseranknorm',    ...
%                           'ranknorm', 'int', 'splines', 'delta', 'none'};
% 
%     % Make a list of all transformations to be done
%     if isJSON
%         all_transforms     = lower({cfg.params.fixed.vars.transformation});
%         all_uq_transforms  = unique(all_transforms, 'stable');
%         tmp_diff_transform = setdiff(all_uq_transforms, allowed_transforms);
%         if ~isempty(tmp_diff_transform)
%             tmp = strcat(tmp_diff_transform, {', '});
%             tmp = horzcat(tmp{:});
%             warning(['The following transforms are not currently handled: ', tmp(1:end-2)]);
%         end
%     end
% 
%     % Handle centering
%     loc_center = check_cfg_parameter({'center', 'centre', 'demean'}, cfg, true);
%     if ~isempty(loc_center)
%         [all_cont_variables{count,1}, all_cont_variables{count,2}, tmp_var_name] = ...
%          applyTransforms(loc_center, names_continuous, cfg_FFXNames, cfg, isJSON, data_work, 'center');
% 
%         % Maintain mapping
%         names_mapping(count_mapping:count_mapping+length(tmp_var_name)-1, 1) = tmp_var_name;
%         names_mapping(count_mapping:count_mapping+length(tmp_var_name)-1, 2) = all_cont_variables{count,2};
%         count_mapping = count_mapping + 1;
%         count = count + 1;
%     end
% 
%     % Handle standardization
%     loc_std = check_cfg_parameter({'std', 'standardize', 'normalize'}, cfg, true);
%     if ~isempty(loc_std)
%         [all_cont_variables{count,1}, all_cont_variables{count,2}, tmp_var_name] = ...
%          applyTransforms(loc_center, names_continuous, cfg_FFXNames, cfg, isJSON, data_work, 'std');
% 
%         % Maintain mapping
%         names_mapping(count_mapping:count_mapping+length(tmp_var_name)-1, 1) = tmp_var_name;
%         names_mapping(count_mapping:count_mapping+length(tmp_var_name)-1, 2) = all_cont_variables{count,2};
%         count_mapping = count_mapping + 1;
%         count = count + 1;
%     end
% 
%     % Handle logn transformation
%     loc_logn = check_cfg_parameter('logn', cfg, true);
%     if ~isempty(loc_logn)
%         [all_cont_variables{count,1}, all_cont_variables{count,2}, tmp_var_name] = ...
%          applyTransforms(loc_center, names_continuous, cfg_FFXNames, cfg, isJSON, data_work, 'logn');
% 
%         % Maintain mapping
%         names_mapping(count_mapping:count_mapping+length(tmp_var_name)-1, 1) = tmp_var_name;
%         names_mapping(count_mapping:count_mapping+length(tmp_var_name)-1, 2) = all_cont_variables{count,2};
%         count_mapping = count_mapping + 1;
%         count = count + 1;
%     end
% 
%     % Handle log10 transformation
%     loc_log10 = check_cfg_parameter('log10', cfg, true);
%     if ~isempty(loc_log10)
%         [all_cont_variables{count,1}, all_cont_variables{count,2}, tmp_var_name] = ...
%          applyTransforms(loc_center, names_continuous, cfg_FFXNames, cfg, isJSON, data_work, 'log10');
% 
%         % Maintain mapping
%         names_mapping(count_mapping:count_mapping+length(tmp_var_name)-1, 1) = tmp_var_name;
%         names_mapping(count_mapping:count_mapping+length(tmp_var_name)-1, 2) = all_cont_variables{count,2};
%         count_mapping = count_mapping + 1;
%         count = count + 1;
%     end
% 
%     % Handle INT
%     loc_INT = check_cfg_parameter({'inverseranknorm', 'ranknorm', 'int'}, cfg, true);
%     if ~isempty(loc_INT)
%         [all_cont_variables{count,1}, all_cont_variables{count,2}, tmp_var_name] = ...
%          applyTransforms(loc_center, names_continuous, cfg_FFXNames, cfg, isJSON, data_work, 'ranknorm');
% 
%         % Maintain mapping
%         names_mapping(count_mapping:count_mapping+length(tmp_var_name)-1, 1) = tmp_var_name;
%         names_mapping(count_mapping:count_mapping+length(tmp_var_name)-1, 2) = all_cont_variables{count,2};
%         count_mapping = count_mapping + 1;
%         count = count + 1;
%     end
% 
%     % Second pass, now we need to handle delta and splines
%     % For each of these transforms, need to loop over - they also produce
%     % different number of outcome than just a single column
% 
%     % Handle delta transform
%     loc_delta = check_cfg_parameter('delta', cfg, true);
%     if ~isempty(loc_delta)
%         if isJSON
%             loc_delta = ismember(names_continuous, cfg_FFXNames(loc_delta));
%         else
%             % Replace space and split by comma
%             tmp_delta = strsplit(strrep(cfg{loc_delta, 2}, ' ', ''), ',');
%             loc_delta = ismember(names_continuous, tmp_delta);
%         end
% 
%         if ~isempty(loc_delta)
%             for v = 1:length(loc_delta)
%                 tmp_var_name = names_continuous{loc_delta(v)};
%                 tmp_new_name = [{[tmp_var_name, '_baseline']}, {[tmp_var_name, '_delta']}];
%                 all_cont_variables{count,1} = computeDelta(data_work{:,tmp_var_name}, ...
%                                                            data_work.participant_id, data_work.session_id);
%                 all_cont_variables{count,2} = tmp_new_name;
% 
%                 % Maintain mapping
%                 names_mapping{count_mapping, 1} = tmp_var_name;
%                 names_mapping{count_mapping, 2} = tmp_new_name;
%                 count_mapping = count_mapping + 1;
%                 count = count + 1;
%             end
%         end
%     end
% 
%     % Handle splines
%     loc_splines = check_cfg_parameter('splines', cfg, false);
% 
%     % Every value in loc_splines is a separate variable for which we need to
%     % make splines; loop over these and do what is necessary
%     if ~isempty(loc_splines)
%         % Initialize a variable for storing spline settings and other relevant
%         % variables that are not part of design matrix but are necessary
%         all_spline_settings = cell(length(loc_splines), 1);
% 
%         % Go over every spline index
%         for v = 1:length(loc_splines)
%             % Parse this line to extract parameters
%             tmp = cfg{loc_splines(v), 2};
% 
%             % Split by space to get name-pair values
%             tmp = strsplit(tmp, ' ');
% 
%             % First entry is the variable name
%             tmp_var_name = tmp{1};
% 
%             % For every remaining entry, need to split on "=" and determine what
%             % those parameters are; start by initializing every parameter with
%             % empty values so that parameters are never shared across variables
%             [tmp_knots, tmp_Xpowers, tmp_minMax, tmp_intercept, tmp_cleanUp, tmp_instance] = deal([]);
%             [tmp_splineType, tmp_method, tmp_outDir, tmp_optCommand, tmp_optAppend] = deal('');
% 
%             for p = 2:length(tmp)
%                 % Split
%                 tmp_param  = strsplit(tmp{p}, '=');
%                 tmp_config = tmp_value{1, 1};
%                 tmp_value  = tmp_param{1, 2};
% 
%                 % Depending on which parameter we find, we assign the values
%                 % Could probably do if-else
%                 if strcmpi(tmp_config, 'knots')
%                     % str2double will make comma separated values into one number
%                     tmp_knots = str2num(tmp_value); %#ok<*ST2NM>
%                 end
% 
%                 if strcmpi(tmp_config, 'splineType')
%                     tmp_splineType = strtrim(tmp_value);
%                 end
% 
%                 if strcmpi(tmp_config, 'Xpowers')
%                     tmp_Xpowers = str2num(tmp_value);
%                 end
% 
%                 if strcmpi(tmp_config, 'method')
%                     tmp_method = strtrim(tmp_value);
%                 end
% 
%                 if strcmpi(tmp_config, 'minMax')
%                     tmp_minMax = sort(str2num(tmp_value));
%                 end
% 
%                 if strcmpi(tmp_config, 'intercept')
%                     tmp_intercept = eval(strtrim(tmp_value));
%                 end
% 
%                 if strcmpi(tmp_config, 'optCommand')
%                     tmp_optCommand = strtrim(tmp_value);
%                 end
% 
%                 if strcmpi(tmp_config, 'optAppend')
%                     tmp_optAppend = strtrim(tmp_value);
%                 end
% 
%                 if strcmpi(tmp_config, 'cleanUp')
%                     tmp_cleanUp = eval(strtrim(tmp_value));
%                 end
% 
%                 if strcmpi(tmp_config, 'instance')
%                     tmp_instance = str2double(tmp_value);
%                 end
%             end
% 
%             % Create basis functions
%             [all_cont_variables{count,1}, all_spline_settings{v,1},   ...
%              all_spline_settings{v,2},    all_spline_settings{v,3},   ...
%              all_spline_settings{v,4},    all_spline_settings{v,5}] = ...
%              createBasisFunctions(data_work.(tmp_var_name), tmp_knots, tmp_splineType, ...
%                                   tmp_Xpowers, tmp_method, tmp_outDir, tmp_minMax,    ...
%                                   tmp_intercept, tmp_optCommand, tmp_optAppend,       ...
%                                   tmp_cleanUp, tmp_instance);
% 
%             % Append the name of the variable for which settings are saved
%             all_spline_settings{v,6} = tmp_var_name;
% 
%             % Name of basis functions
%             tmp_new_name = strcat({tmp_var_name}, {'_bf'}, num2str((1:size(all_cont_variables{count,1},2))'))';
%             all_cont_variables{count,2} = tmp_new_name;
% 
%             % Maintain mapping
%             names_mapping{count_mapping, 1} = tmp_var_name;
%             names_mapping{count_mapping, 2} = tmp_new_name;
%             count_mapping = count_mapping + 1;
%             count = count + 1;
%         end
%     end
% 
%     % Now track what variables were not transformed
%     vars_none_transform = setdiff(names_continuous, names_mapping(:,1));
% 
%     % Stack variables that were not transformed
%     for v = 1:length(vars_none_transform)
%         tmp_var_name = vars_none_transform{v};
%         all_cont_variables{count,1} = data_work{:,tmp_var_name};
%         all_cont_variables{count,2} = tmp_new_name;
%         count = count + 1;
%     end
% 
%     % Shrink all_cont_variables, if required
%     tmp = cellfun(@isempty, all_cont_variables);
%     all_cont_variables(tmp,:) = [];
% 
%     % Put all continuous variables together
%     X_data_continuous = horzcat(all_cont_variables{:,1});
%     X_name_continuous = horzcat(all_cont_variables{:,2});
% else
%     X_data_continuous = [];
%     X_name_continuous = '';
% end

% %% Split into categorical and continuous covariates
% names_continuous  = data_work.Properties.VariableNames(loc_cont);
% names_categorical = data_work.Properties.VariableNames(loc_catg);

% %% Dummy code categorical variables and drop reference levels
% if ~isempty(names_categorical)
%     all_catg_variables = cell(length(names_categorical), 2);
%     if ~isJSON
%         loc_reference = check_cfg_parameter('reference', cfg(:,1), true);
% 
%         if isempty(loc_reference)
%             % No reference levels specified
%             if intercept
%                 global_refLevel = 'mode';
%             else
%                 global_refLevel = 'none';
%             end
%         end
% 
%         % Split the reference line
%         % Format: <varName>=<string> <varName>=<string> ...
%         cfg_ref_map = strsplit(cfg{loc_reference, 2}, ' ')';
%         cfg_ref_map = cellfun(@(x) strsplit(x, '='), cfg_ref_map, 'UniformOutput', false);
%         try
%             cfg_ref_map = vertcat(cfg_ref_map{:});
%         catch
%             error(['Please check the specification of reference level in the config file; ', ...
%                    'should be: <varName>=<string> <varName>=<string> ...']);
%         end
%     end
% 
%     for cc = 1:length(names_categorical)
% 
%         % Which variable are we looking at?
%         tmp_var_name = names_categorical{cc};
% 
%         % Parse reference level from configuration file (if specified)
%         if isJSON
%             % Does this variable exist?
%             tmp_loc = find(strcmpi(tmp_var_name, cfg_FFXNames));
%             if ~isempty(tmp_loc)
%                 try
%                     tmp_ref_level = cfg.params.fixed.vars.dummy(tmp_loc);
%                 catch
%                     % Reference not specified; use defaults
%                     tmp_ref_level = global_refLevel;
%                 end
%             else
%                 % Variable is not part of the config file - use defaults
%                 tmp_ref_level = global_refLevel;
%             end
%         else
%             % Does this variable exist?
%             tmp_loc = find(strcmpi(tmp_var_name, cfg_ref_map));
%             if ~isempty(tmp_loc)
%                 tmp_ref_level = cfg_ref_map{tmp_loc,2};
%             else
%                 % Reference level not specified - use defaults
%                 tmp_ref_level = global_refLevel;
%             end
%         end
% 
%         % Create dummy coding
%         [all_catg_variables{cc,1}, all_catg_variables{cc,2}] = ...
%          expand_categorical(data_work.(tmp_var_name), tmp_var_name, tmp_ref_level);
% 
%         % Maintain mapping
%         names_mapping{count_mapping, 1} = tmp_var_name;
%         names_mapping{count_mapping, 2} = all_catg_variables{cc,2};
%         count_mapping = count_mapping + 1;
%     end
% 
%     % Put all categorical variables together
%     X_data_categorical = horzcat(all_catg_variables{:,1});
%     X_name_categorical = horzcat(all_catg_variables{:,2});
% else
%     X_data_categorical = [];
%     X_name_categorical = '';
% end

% %% Perform transformations on continuous variables
% if ~isempty(names_continuous)
%     all_cont_variables = cell(length(names_continuous), 2);
% 
%     % First pass, handle transforms that generate a vector - these can be
%     % handled together in a single iteration without a need for loop
%     % Examples of these transforms are (all handled by doTransformation):
%     %   * 'center' | 'centre' | 'demean'
%     %   * 'std' | 'standardize' | 'normalize'
%     %   * 'logn'
%     %   * 'log10'
%     %   * 'inverseranknorm' | 'ranknorm' | 'int'
% 
%     % Use a counter for accessing all_cont_variables: shrink later
%     count = 1;
% 
%     % List of possible transforms (including spelling variations and
%     % shortened forms)
%     allowed_transforms = {'center', 'centre', 'demean', 'std', 'standardize', ...
%                           'normalize', 'logn', 'log10', 'inverseranknorm',    ...
%                           'ranknorm', 'int', 'splines', 'delta', 'none'};
% 
%     % Make a list of all transformations to be done
%     if isJSON
%         all_transforms     = lower({cfg.params.fixed.vars.transformation});
%         all_uq_transforms  = unique(all_transforms, 'stable');
%         tmp_diff_transform = setdiff(all_uq_transforms, allowed_transforms);
%         if ~isempty(tmp_diff_transform)
%             tmp = strcat(tmp_diff_transform, {', '});
%             tmp = horzcat(tmp{:});
%             warning(['The following transforms are not currently handled: ', tmp(1:end-2)]);
%         end
%     end
% 
%     % Handle centering
%     loc_center = check_cfg_parameter({'center', 'centre', 'demean'}, cfg, true);
%     if ~isempty(loc_center)
%         [all_cont_variables{count,1}, all_cont_variables{count,2}, tmp_var_name] = ...
%          applyTransforms(loc_center, names_continuous, cfg_FFXNames, cfg, isJSON, data_work, 'center');
% 
%         % Maintain mapping
%         names_mapping(count_mapping:count_mapping+length(tmp_var_name)-1, 1) = tmp_var_name;
%         names_mapping(count_mapping:count_mapping+length(tmp_var_name)-1, 2) = all_cont_variables{count,2};
%         count_mapping = count_mapping + 1;
%         count = count + 1;
%     end
% 
%     % Handle standardization
%     loc_std = check_cfg_parameter({'std', 'standardize', 'normalize'}, cfg, true);
%     if ~isempty(loc_std)
%         [all_cont_variables{count,1}, all_cont_variables{count,2}, tmp_var_name] = ...
%          applyTransforms(loc_center, names_continuous, cfg_FFXNames, cfg, isJSON, data_work, 'std');
% 
%         % Maintain mapping
%         names_mapping(count_mapping:count_mapping+length(tmp_var_name)-1, 1) = tmp_var_name;
%         names_mapping(count_mapping:count_mapping+length(tmp_var_name)-1, 2) = all_cont_variables{count,2};
%         count_mapping = count_mapping + 1;
%         count = count + 1;
%     end
% 
%     % Handle logn transformation
%     loc_logn = check_cfg_parameter('logn', cfg, true);
%     if ~isempty(loc_logn)
%         [all_cont_variables{count,1}, all_cont_variables{count,2}, tmp_var_name] = ...
%          applyTransforms(loc_center, names_continuous, cfg_FFXNames, cfg, isJSON, data_work, 'logn');
% 
%         % Maintain mapping
%         names_mapping(count_mapping:count_mapping+length(tmp_var_name)-1, 1) = tmp_var_name;
%         names_mapping(count_mapping:count_mapping+length(tmp_var_name)-1, 2) = all_cont_variables{count,2};
%         count_mapping = count_mapping + 1;
%         count = count + 1;
%     end
% 
%     % Handle log10 transformation
%     loc_log10 = check_cfg_parameter('log10', cfg, true);
%     if ~isempty(loc_log10)
%         [all_cont_variables{count,1}, all_cont_variables{count,2}, tmp_var_name] = ...
%          applyTransforms(loc_center, names_continuous, cfg_FFXNames, cfg, isJSON, data_work, 'log10');
% 
%         % Maintain mapping
%         names_mapping(count_mapping:count_mapping+length(tmp_var_name)-1, 1) = tmp_var_name;
%         names_mapping(count_mapping:count_mapping+length(tmp_var_name)-1, 2) = all_cont_variables{count,2};
%         count_mapping = count_mapping + 1;
%         count = count + 1;
%     end
% 
%     % Handle INT
%     loc_INT = check_cfg_parameter({'inverseranknorm', 'ranknorm', 'int'}, cfg, true);
%     if ~isempty(loc_INT)
%         [all_cont_variables{count,1}, all_cont_variables{count,2}, tmp_var_name] = ...
%          applyTransforms(loc_center, names_continuous, cfg_FFXNames, cfg, isJSON, data_work, 'ranknorm');
% 
%         % Maintain mapping
%         names_mapping(count_mapping:count_mapping+length(tmp_var_name)-1, 1) = tmp_var_name;
%         names_mapping(count_mapping:count_mapping+length(tmp_var_name)-1, 2) = all_cont_variables{count,2};
%         count_mapping = count_mapping + 1;
%         count = count + 1;
%     end
% 
%     % Second pass, now we need to handle delta and splines
%     % For each of these transforms, need to loop over - they also produce
%     % different number of outcome than just a single column
% 
%     % Handle delta transform
%     loc_delta = check_cfg_parameter('delta', cfg, true);
%     if ~isempty(loc_delta)
%         if isJSON
%             loc_delta = ismember(names_continuous, cfg_FFXNames(loc_delta));
%         else
%             % Replace space and split by comma
%             tmp_delta = strsplit(strrep(cfg{loc_delta, 2}, ' ', ''), ',');
%             loc_delta = ismember(names_continuous, tmp_delta);
%         end
% 
%         if ~isempty(loc_delta)
%             for v = 1:length(loc_delta)
%                 tmp_var_name = names_continuous{loc_delta(v)};
%                 tmp_new_name = [{[tmp_var_name, '_baseline']}, {[tmp_var_name, '_delta']}];
%                 all_cont_variables{count,1} = computeDelta(data_work{:,tmp_var_name}, ...
%                                                            data_work.participant_id, data_work.session_id);
%                 all_cont_variables{count,2} = tmp_new_name;
% 
%                 % Maintain mapping
%                 names_mapping{count_mapping, 1} = tmp_var_name;
%                 names_mapping{count_mapping, 2} = tmp_new_name;
%                 count_mapping = count_mapping + 1;
%                 count = count + 1;
%             end
%         end
%     end
% 
%     % Handle splines
%     loc_splines = check_cfg_parameter('splines', cfg, false);
% 
%     % Every value in loc_splines is a separate variable for which we need to
%     % make splines; loop over these and do what is necessary
%     if ~isempty(loc_splines)
%         % Initialize a variable for storing spline settings and other relevant
%         % variables that are not part of design matrix but are necessary
%         all_spline_settings = cell(length(loc_splines), 1);
% 
%         % Go over every spline index
%         for v = 1:length(loc_splines)
%             % Parse this line to extract parameters
%             tmp = cfg{loc_splines(v), 2};
% 
%             % Split by space to get name-pair values
%             tmp = strsplit(tmp, ' ');
% 
%             % First entry is the variable name
%             tmp_var_name = tmp{1};
% 
%             % For every remaining entry, need to split on "=" and determine what
%             % those parameters are; start by initializing every parameter with
%             % empty values so that parameters are never shared across variables
%             [tmp_knots, tmp_Xpowers, tmp_minMax, tmp_intercept, tmp_cleanUp, tmp_instance] = deal([]);
%             [tmp_splineType, tmp_method, tmp_outDir, tmp_optCommand, tmp_optAppend] = deal('');
% 
%             for p = 2:length(tmp)
%                 % Split
%                 tmp_param  = strsplit(tmp{p}, '=');
%                 tmp_config = tmp_value{1, 1};
%                 tmp_value  = tmp_param{1, 2};
% 
%                 % Depending on which parameter we find, we assign the values
%                 % Could probably do if-else
%                 if strcmpi(tmp_config, 'knots')
%                     % str2double will make comma separated values into one number
%                     tmp_knots = str2num(tmp_value); %#ok<*ST2NM>
%                 end
% 
%                 if strcmpi(tmp_config, 'splineType')
%                     tmp_splineType = strtrim(tmp_value);
%                 end
% 
%                 if strcmpi(tmp_config, 'Xpowers')
%                     tmp_Xpowers = str2num(tmp_value);
%                 end
% 
%                 if strcmpi(tmp_config, 'method')
%                     tmp_method = strtrim(tmp_value);
%                 end
% 
%                 if strcmpi(tmp_config, 'minMax')
%                     tmp_minMax = sort(str2num(tmp_value));
%                 end
% 
%                 if strcmpi(tmp_config, 'intercept')
%                     tmp_intercept = eval(strtrim(tmp_value));
%                 end
% 
%                 if strcmpi(tmp_config, 'optCommand')
%                     tmp_optCommand = strtrim(tmp_value);
%                 end
% 
%                 if strcmpi(tmp_config, 'optAppend')
%                     tmp_optAppend = strtrim(tmp_value);
%                 end
% 
%                 if strcmpi(tmp_config, 'cleanUp')
%                     tmp_cleanUp = eval(strtrim(tmp_value));
%                 end
% 
%                 if strcmpi(tmp_config, 'instance')
%                     tmp_instance = str2double(tmp_value);
%                 end
%             end
% 
%             % Create basis functions
%             [all_cont_variables{count,1}, all_spline_settings{v,1},   ...
%              all_spline_settings{v,2},    all_spline_settings{v,3},   ...
%              all_spline_settings{v,4},    all_spline_settings{v,5}] = ...
%              createBasisFunctions(data_work.(tmp_var_name), tmp_knots, tmp_splineType, ...
%                                   tmp_Xpowers, tmp_method, tmp_outDir, tmp_minMax,    ...
%                                   tmp_intercept, tmp_optCommand, tmp_optAppend,       ...
%                                   tmp_cleanUp, tmp_instance);
% 
%             % Append the name of the variable for which settings are saved
%             all_spline_settings{v,6} = tmp_var_name;
% 
%             % Name of basis functions
%             tmp_new_name = strcat({tmp_var_name}, {'_bf'}, num2str((1:size(all_cont_variables{count,1},2))'))';
%             all_cont_variables{count,2} = tmp_new_name;
% 
%             % Maintain mapping
%             names_mapping{count_mapping, 1} = tmp_var_name;
%             names_mapping{count_mapping, 2} = tmp_new_name;
%             count_mapping = count_mapping + 1;
%             count = count + 1;
%         end
%     end
% 
%     % Now track what variables were not transformed
%     vars_none_transform = setdiff(names_continuous, names_mapping(:,1));
% 
%     % Stack variables that were not transformed
%     for v = 1:length(vars_none_transform)
%         tmp_var_name = vars_none_transform{v};
%         all_cont_variables{count,1} = data_work{:,tmp_var_name};
%         all_cont_variables{count,2} = tmp_new_name;
%         count = count + 1;
%     end
% 
%     % Shrink all_cont_variables, if required
%     tmp = cellfun(@isempty, all_cont_variables);
%     all_cont_variables(tmp,:) = [];
% 
%     % Put all continuous variables together
%     X_data_continuous = horzcat(all_cont_variables{:,1});
%     X_name_continuous = horzcat(all_cont_variables{:,2});
% else
%     X_data_continuous = [];
%     X_name_continuous = '';
% end
