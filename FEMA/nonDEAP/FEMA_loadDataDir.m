function [data, fname_data_out] = FEMA_loadDataDir(dirname_data, var_names, varargin)
% FEMA_LOADDATADIR  Load variables from a directory of tabular data files
% and merge them into a single table.
%
% Users store covariate / independent variable data across multiple
% tsv / csv / parquet files in a common directory. This function scans all
% such files, finds which file contains each requested variable (first match
% wins when a variable name appears in more than one file), reads only the
% needed columns from each file, merges everything by participant_id and
% session_id, and saves the result as a parquet file.
%
%% Usage:
%   [data, fname_data_out] = FEMA_loadDataDir(dirname_data, var_names)
%   [data, fname_data_out] = FEMA_loadDataDir(dirname_data, var_names, ...
%                               'id_vars', {'participant_id','session_id'}, ...
%                               'outDir',  tempdir, ...
%                               'outName', 'merged_data.parquet', ...
%                               'verbose', true)
%
%% Required inputs:
%   dirname_data  <char> or <cell of char>  One directory, or multiple roots
%                               scanned in order (e.g. dependent then data):
%                               all .tsv / .csv / .parquet files are listed and
%                               merged; same basename in later roots loses to an
%                               earlier file after format_priority dedup.
%   var_names     <cell of char> Variable names to collect (column names as
%                               they appear in the data files).
%
%% Optional name-value inputs:
%   id_vars         <cell of char> ID columns used as merge keys.
%                                 Default: {'participant_id', 'session_id'}
%   outDir          <char>        Directory to write the merged file.
%                                 Default: tempdir
%   outName         <char>        Filename for the merged parquet output.
%                                 Default: 'merged_data.parquet'
%   format_priority <cell of char> Format preference when the same file
%                                 exists in multiple formats. First entry
%                                 has highest priority.
%                                 Default: {'parquet', 'csv', 'tsv'}
%   name_map        containers.Map  Optional rename original -> short column name.
%   var_types       containers.Map  Optional; keys = final column names after
%                                 name_map, values = type_var from JSON
%                                 ('categorical'|'continuous'). Categorical
%                                 columns are written as MATLAB categorical.
%   verbose         <logical>     Print progress messages. Default: true
%
%% Outputs:
%   data          <table>       Merged table with id_vars + all found vars.
%   fname_data_out <char>       Full path to the saved merged parquet file.

    %% Parse optional inputs -----------------------------------------------
    p = inputParser();
    addParameter(p, 'id_vars',         {'participant_id', 'session_id'}, @iscell);
    addParameter(p, 'outDir',          pwd,                              @ischar);
    addParameter(p, 'outName',         'merged_data.parquet',            @ischar);
    addParameter(p, 'format_priority', {'parquet', 'csv', 'tsv'},        @iscell);
    addParameter(p, 'name_map',        containers.Map('KeyType', 'char', 'ValueType', 'char'), ...
                                                                         @(x) isa(x, 'containers.Map'));
    addParameter(p, 'var_types',     containers.Map('KeyType', 'char', 'ValueType', 'char'), ...
                                                                         @(x) isa(x, 'containers.Map'));
    addParameter(p, 'verbose',         true,                             @islogical);
    parse(p, varargin{:});

    id_vars         = p.Results.id_vars;
    outDir          = p.Results.outDir;
    outName         = p.Results.outName;
    format_priority = p.Results.format_priority;
    name_map        = p.Results.name_map;
    var_types       = p.Results.var_types;
    verbose         = p.Results.verbose;

    if ischar(var_names)
        var_names = {var_names};
    end

    % Always include the family ID variable
    fam_var = 'ab_g_stc__design_id__fam';
    if ~ismember(fam_var, var_names)
        var_names{end+1} = fam_var;
    end

    %% Normalize directory root(s) -------------------------------------------
    roots = normalizeDataRoots(dirname_data);
    roots_label = strjoin(roots, ', ');

    %% List data files across all roots (order preserved: earlier root first)
    all_files = [];
    for r = 1:numel(roots)
        d = roots{r};
        if ~isfolder(d)
            error('FEMA_loadDataDir: not a directory: %s', d);
        end
        parquet_files = dir(fullfile(d, '*.parquet'));
        csv_files     = dir(fullfile(d, '*.csv'));
        tsv_files     = dir(fullfile(d, '*.tsv'));
        all_files = [all_files; parquet_files; csv_files; tsv_files]; %#ok<AGROW>
    end

    if isempty(all_files)
        error('FEMA_loadDataDir: no .parquet, .csv, or .tsv files found under: %s', roots_label);
    end

    %% Deduplicate files with the same basename ----------------------------
    % If the same file exists in multiple formats (e.g. data.parquet and
    % data.csv), keep only the highest-priority format per basename.
    % Same basename in different roots: first occurrence in the combined list wins.
    all_files = deduplicateByFormat(all_files, format_priority, verbose);

    if verbose
        logging('FEMA_loadDataDir: found %d data files under %s\n', numel(all_files), roots_label);
    end

    %% Build variable -> file map (first match wins) -----------------------
    var_file_map = containers.Map('KeyType', 'char', 'ValueType', 'char');
    tBuildVarMap = tic;
    for f = 1:numel(all_files)
        fpath = fullfile(all_files(f).folder, all_files(f).name);
        [~, ~, fext] = fileparts(fpath);

        try
            colnames = getColumnNames(fpath, fext);
        catch ME
            warning('FEMA_loadDataDir: could not read column names from %s: %s', fpath, ME.message);
            continue
        end

        for v = 1:numel(var_names)
            vname = var_names{v};
            if any(strcmpi(colnames, vname))
                if isKey(var_file_map, vname)
                    existing = var_file_map(vname);
                    warning(['FEMA_loadDataDir: variable "%s" found in multiple files. ' ...
                             'Using first match: %s (ignoring %s).'], ...
                             vname, existing, fpath);
                else
                    var_file_map(vname) = fpath;
                end
            end
        end
    end
    toc(tBuildVarMap);
    
    %% Report missing variables --------------------------------------------
    found_vars   = keys(var_file_map);
    missing_vars = setdiff(var_names, found_vars);
    if ~isempty(missing_vars)
        error('FEMA_loadDataDir: the following requested variables were not found in any file:\n  %s', ...
                strjoin(missing_vars, ', '));
    end

    if isempty(found_vars)
        error('FEMA_loadDataDir: none of the requested variables were found in any file under %s', roots_label);
    end

    %% Group variables by source file --------------------------------------
    % For each unique file, collect the variables we need from it
    unique_files = unique(values(var_file_map));
    partial_tables = cell(numel(unique_files), 1);

    for f = 1:numel(unique_files)
        fpath = unique_files{f};
        [~, ~, fext] = fileparts(fpath);

        % Determine which variables come from this file
        file_vars = {};
        for v = 1:numel(found_vars)
            if strcmp(var_file_map(found_vars{v}), fpath)
                file_vars{end+1} = found_vars{v}; %#ok<AGROW>
            end
        end

        cols_to_read = [id_vars, file_vars];

        if verbose
            logging('FEMA_loadDataDir: reading %d variable(s) from %s\n', numel(file_vars), fpath);
        end

        try
            partial_tables{f} = readDataFile(fpath, fext, cols_to_read, id_vars);
        catch ME
            error('FEMA_loadDataDir: failed to read %s: %s', fpath, ME.message);
        end
    end

    %% Merge partial tables on id_vars -------------------------------------
    % Files without session_id (static variables) are joined on participant_id
    % only, broadcasting the value to all sessions for that participant.
    data = partial_tables{1};

    for f = 2:numel(partial_tables)
        join_keys = intersect(id_vars, partial_tables{f}.Properties.VariableNames, 'stable');
        try
            data = innerjoin(data, partial_tables{f}, 'Keys', join_keys);
        catch ME
            error('FEMA_loadDataDir: failed to merge table from file %s: %s', unique_files{f}, ME.message);
        end
    end

    if verbose
        logging('FEMA_loadDataDir: merged table has %d rows and %d columns.\n', height(data), width(data));
    end

    %% Reorder columns: id_vars first, then fam_var, then everything else --
    all_cols   = data.Properties.VariableNames;
    id_cols    = intersect(id_vars,   all_cols, 'stable');
    fam_col    = intersect({fam_var}, all_cols);
    other_cols = setdiff(all_cols, [id_cols, fam_col], 'stable');
    data = data(:, [id_cols, fam_col, other_cols]);

    %% Rename family ID column to family_id --------------------------------
    if ismember(fam_var, data.Properties.VariableNames)
        data.Properties.VariableNames{strcmp(data.Properties.VariableNames, fam_var)} = 'family_id';
    end

    %% Rename columns to short names ---------------------------------------
    if ~isempty(name_map)
        orig_names = data.Properties.VariableNames;
        for c = 1:numel(orig_names)
            if isKey(name_map, orig_names{c})
                data.Properties.VariableNames{c} = name_map(orig_names{c});
            end
        end
    end

    %% Apply variable types from job spec (type_var) before writing parquet
    if ~isempty(var_types) && var_types.Count > 0
        data = applyVarTypesForParquet(data, var_types);
    end

    %% Save merged table as parquet ----------------------------------------
    if ~exist(outDir, 'dir')
        mkdir(outDir);
    end

    fname_data_out = fullfile(outDir, outName);
    parquetwrite(fname_data_out, data);

    if verbose
        logging('FEMA_loadDataDir: merged data saved to %s\n', fname_data_out);
    end
end


%% =========================================================================
%  Helper functions
%% =========================================================================

function data = applyVarTypesForParquet(data, var_types)
% Set column types for merged output using JSON type_var (keys = final column names).
    vn = data.Properties.VariableNames;
    for k = 1:numel(vn)
        col = vn{k};
        if ~isKey(var_types, col)
            continue
        end
        typ = lower(strtrim(var_types(col)));
        if ~strcmpi(typ, 'categorical')
            continue
        end
        v = data.(col);
        if iscategorical(v)
            continue
        end
        if isnumeric(v) || islogical(v)
            data.(col) = categorical(v);
        elseif iscell(v) || isstring(v)
            data.(col) = categorical(v);
        elseif ischar(v) && size(v, 2) > 1
            data.(col) = categorical(cellstr(v));
        end
    end
end


function roots = normalizeDataRoots(dirname_data)
% Single char path, cell array of paths, or string array — unique paths, first wins.
    if ischar(dirname_data)
        roots = {dirname_data};
    elseif iscell(dirname_data)
        roots = dirname_data(:)';
    elseif isstring(dirname_data)
        roots = cellstr(dirname_data(:));
    else
        error('FEMA_loadDataDir: dirname_data must be char, cell, or string');
    end
    roots = roots(~cellfun(@isempty, roots));
    if isempty(roots)
        error('FEMA_loadDataDir: no directory paths provided');
    end
    seen = containers.Map('KeyType', 'char', 'ValueType', 'logical');
    keep = {};
    for i = 1:numel(roots)
        p = roots{i};
        if ~ischar(p)
            p = char(p);
        end
        if ~isKey(seen, p)
            seen(p) = true;
            keep{end+1} = p; %#ok<AGROW>
        end
    end
    roots = keep;
end


function files = deduplicateByFormat(files, format_priority, verbose)
% For files sharing the same basename, keep only the highest-priority format.
    basenames = cell(numel(files), 1);
    for i = 1:numel(files)
        [~, bn, ~] = fileparts(files(i).name);
        basenames{i} = bn;
    end

    unique_basenames = unique(basenames);
    keep = true(numel(files), 1);

    for i = 1:numel(unique_basenames)
        idx = find(strcmp(basenames, unique_basenames{i}));
        if numel(idx) == 1
            continue
        end

        % Multiple formats for the same basename — pick highest priority
        best_priority = numel(format_priority) + 1;
        best_idx      = idx(1);
        for j = 1:numel(idx)
            [~, ~, ext] = fileparts(files(idx(j)).name);
            ext = lower(ext(2:end));  % strip leading dot
            pri = find(strcmpi(format_priority, ext), 1);
            if isempty(pri)
                pri = numel(format_priority) + 1;
            end
            if pri < best_priority
                best_priority = pri;
                best_idx      = idx(j);
            end
        end

        discard_idx = idx(idx ~= best_idx);
        keep(discard_idx) = false;

        if verbose
            discard_names = strjoin({files(discard_idx).name}, ', ');
            logging(['FEMA_loadDataDir: "%s" exists in multiple formats. ' ...
                     'Using %s; ignoring %s.\n'], ...
                     unique_basenames{i}, files(best_idx).name, discard_names);
        end
    end

    files = files(keep);
end


function colnames = getColumnNames(fpath, fext)
% Return the column names of a file without loading all data.
    switch lower(fext)
        case '.parquet'
            info     = parquetinfo(fpath);
            colnames = info.VariableNames;
        case '.csv'
            opts = detectImportOptions(fpath);
            colnames = opts.VariableNames;
        case '.tsv'
            % Some MATLAB versions do not infer .tsv extensions automatically.
            opts = detectImportOptions(fpath, 'FileType', 'text', 'Delimiter', '\t');
            colnames = opts.VariableNames;
        otherwise
            error('Unsupported file extension: %s', fext);
    end
end


function tbl = readDataFile(fpath, fext, cols_to_read, id_vars)
% Read specified columns from a file.
% participant_id is always required. session_id is optional — files without
% it contain static variables that will be broadcast across sessions on merge.
    required_id = id_vars{1};   % participant_id

    switch lower(fext)
        case '.parquet'
            info      = parquetinfo(fpath);
            existing  = info.VariableNames;
            if ~ismember(required_id, existing)
                error('Required ID column "%s" not found in %s', required_id, fpath);
            end
            available = intersect(cols_to_read, existing, 'stable');
            tbl = parquetread(fpath, 'SelectedVariableNames', available);

        case '.csv'
            opts     = detectImportOptions(fpath);
            existing = opts.VariableNames;
            if ~ismember(required_id, existing)
                error('Required ID column "%s" not found in %s', required_id, fpath);
            end
            available = intersect(cols_to_read, existing, 'stable');
            opts.SelectedVariableNames = available;
            tbl = readtable(fpath, opts);

        case '.tsv'
            % Some MATLAB versions do not infer .tsv extensions automatically.
            opts     = detectImportOptions(fpath, 'FileType', 'text', 'Delimiter', '\t');
            existing = opts.VariableNames;
            if ~ismember(required_id, existing)
                error('Required ID column "%s" not found in %s', required_id, fpath);
            end
            available = intersect(cols_to_read, existing, 'stable');
            opts.SelectedVariableNames = available;
            tbl = readtable(fpath, opts);

        otherwise
            error('Unsupported file extension: %s', fext);
    end
end
