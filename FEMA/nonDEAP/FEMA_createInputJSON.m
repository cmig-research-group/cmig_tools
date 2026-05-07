function fname_json = FEMA_createInputJSON(fname_toml, varargin)
% FEMA_CREATEINPUTJSON  Convert a non-DEAP TOML config file into a FEMA
% job-spec JSON file compatible with FEMA_parseInputs and FEMA_parse_JSON.
%
% The generated JSON contains all fields read by the existing FEMA parsers
% plus two custom top-level fields (dirname_tabular, dirname_out) that are
% used by the non-DEAP FEMA_wrapper entry point and ignored by all other
% existing FEMA code.
%
%% Usage:
%   fname_json = FEMA_createInputJSON(fname_toml)
%   fname_json = FEMA_createInputJSON(fname_toml, 'outDir', '/my/dir')
%   fname_json = FEMA_createInputJSON(fname_toml, 'outName', 'my_spec.json')
%
%% Required input:
%   fname_toml  <char>   Full path to a TOML config file.
%
%% Optional name-value inputs:
%   outDir      <char>   Directory in which to write the JSON file.
%                        Default: [output].dirname from the TOML config
%   outName     <char>   Filename for the JSON file.
%                        Default: 'fema_spec_<timestamp>.json'
%   designOnly  <logical> If true, [dependent] may be omitted (design-matrix-only
%                TOML; FEMA_makeDesign passes this when converting .toml). Default false
%                for FEMA_wrapper / FEMA_readToml.
%
%% Output:
%   fname_json  <char>   Full path to the written JSON file.
%
%% TOML config format (all sections and keys):
%
%   [dependent]   % required unless designOnly is true (FEMA_makeDesign .toml path)
%   fstem    = "thickness-sm16"   % imaging phenotype / variable stem
%   datatype = "vertex"           % vertex | voxel | roi | corrmat | external
%   dirname  = "/path/to/data"    % path to imaging data dir or file
%
%   [output]
%   dirname  = "/path/to/results"
%
%   [data]
%   dirname  = "/path/to/tabular_data_dir"  % dir with tsv/csv/parquet files
%
%   [makeDesign]   % optional — same names as FEMA_makeDesign name-values (merged via FEMA_mergeArgs)
%   dropMissing = false
%   outName     = "MyDesign"
%
%   [model]
%   study      = "abcd"
%   release    = "6.0"
%   transformY = "inverseranknorm"
%   random     = ["F", "S", "E"]
%
%   [[model.vars]]
%   name        = "visit_age"
%   shortname   = "age"          % optional: display name in outputs
%   type        = "continuous"   % continuous | categorical
%   of_interest = true
%   transform   = "splines"      % optional inline transform for continuous vars
%   type_spline  = "nsk"          %   (splines only) nsk | bs | poly
%   method      = "svd"          %   (splines only) svd | raw
%   knots       = "quartiles"    %   (splines only) "quartiles" or [12, 14, 16]
%   xpowers     = [0, 1]         %   (splines only) optional
%   minmax      = [1, 99]        %   (winsorize only) percentile bounds
%
%   [[model.vars]]
%   name        = "site"
%   type        = "categorical"
%   reference   = "mode"         % mode | last | first | <specific level>
%   of_interest = false
%
%   [model.design]
%   transform_global = "demean"
%   intercept        = true
%   n_gpcs           = 0
%
%   % Alternative: transforms can also be specified in a separate section
%   % (fallback when no inline transform field is present on the var):
%   [[model.transformations]]
%   var        = "visit_age"     % must match a name in [[model.vars]]
%   type       = "splines"       % splines | delta | log10 | demean | center |
%                                %   std | inverseranknorm | ranknorm |
%                                %   quadratic | winsorize
%   type_spline = "nsk"
%   method     = "svd"
%   knots      = "quartiles"
%
%   [[model.interactions]]
%   vars        = "visit_age * sex"
%   of_interest = true
%
%   [advanced]
%   type_perm      = "none"      % none | wildbootstrap | wildbootstrap-nn
%   n_perm         = 0
%   type_cov       = "analytic"  % analytic | unstructured
%   type_fixed_est = "gls"       % gls | ols
%   n_bins         = 20
%   precision      = "double"    % double | single
%
%   [options]
%   % Any FEMA_wrapper parameter not covered by the sections above.
%   % Field names must match FEMA_wrapper parameter names (case-insensitive).
%   % Example:
%   outputType      = "mat"
%   corrvec_thresh  = 0.8
%   returnResiduals = false
%   ico             = 5
%
%   % Only needed if random includes "A" (additive genetic effect):
%   % [model.grm]
%   % dirname = "/path/to/grm.mat"

    %% Parse optional name-value arguments
    p = inputParser();
    addParameter(p, 'outDir',  '',  @ischar);
    addParameter(p, 'outName', '',  @ischar);
    addParameter(p, 'designOnly', false, @(x) islogical(x) || isnumeric(x));
    parse(p, varargin{:});

    outDir  = p.Results.outDir;
    outName = p.Results.outName;
    designOnly = logical(p.Results.designOnly);

    %% Parse TOML file
    cfg = toml_parse(fname_toml);

    % Default outDir to the output dirname in the TOML unless explicitly overridden
    if isempty(outDir)
        if isfield(cfg, 'output') && isfield(cfg.output, 'dirname')
            outDir = cfg.output.dirname;
        else
            outDir = pwd;
        end
    end

    if isempty(outName)
        outName = ['fema_spec.json'];
    end

    %% Validate required sections
    required = {'output', 'data', 'model'};
    for r = 1:numel(required)
        if ~isfield(cfg, required{r})
            error('FEMA_createInputJSON: required TOML section [%s] is missing.', required{r});
        end
    end

    if ~isfield(cfg, 'dependent')
        if designOnly
            cfg.dependent = struct('fstem', '', 'datatype', 'external');
        else
            error('FEMA_createInputJSON: required TOML section [dependent] is missing.');
        end
    end
    dep = cfg.dependent;

    if designOnly
        if ~isfield(dep, 'fstem') || isempty(dep.fstem)
            dep.fstem = '';
        end
        if ~isfield(dep, 'datatype') || isempty(dep.datatype)
            dep.datatype = 'external';
        end
    else
        if ~isfield(dep, 'fstem') || isempty(dep.fstem)
            error('FEMA_createInputJSON: [dependent].fstem is required.');
        end
        if ~isfield(dep, 'datatype') || isempty(dep.datatype)
            error('FEMA_createInputJSON: [dependent].datatype is required.');
        end
    end
    if ~isfield(cfg.output, 'dirname')
        error('FEMA_createInputJSON: [output].dirname is required.');
    end
    if ~isfield(cfg.data, 'dirname')
        error('FEMA_createInputJSON: [data].dirname is required.');
    end

    %% Build JSON struct ---------------------------------------------------

    % Top-level metadata
    ts = char(datetime('now', 'Format', 'yyyyMMdd_HHmmss'));
    s.study   = getField(cfg.model, 'study',   'abcd');
    s.release = getField(cfg.model, 'release', '6.0');
    s.type    = 'analysis';
    s.id      = ['nonDEAP_' ts];

    % Custom top-level fields (ignored by existing parsers)
    % External: scan both [dependent].dirname and [data].dirname (unique, order kept).
    % Other datatypes: single root from [data].dirname.
    if strcmpi(lower(dep.datatype), 'external')
        roots = {};
        if isfield(dep, 'dirname') && ~isempty(dep.dirname)
            ddn = dep.dirname;
            if ischar(ddn)
                roots{end+1} = ddn;
            elseif iscell(ddn)
                for ii = 1:numel(ddn)
                    if ~isempty(ddn{ii})
                        roots{end+1} = char(ddn{ii});
                    end
                end
            else
                roots{end+1} = char(ddn);
            end
        end
        if ~isempty(cfg.data.dirname)
            roots{end+1} = char(cfg.data.dirname);
        end
        roots = uniqueRootPaths(roots);
        if isempty(roots)
            if designOnly && ~isempty(cfg.data.dirname)
                roots = {char(cfg.data.dirname)};
            else
                error('FEMA_createInputJSON: for external, set [dependent].dirname and/or [data].dirname.');
            end
        end
        if numel(roots) == 1
            s.dirname_tabular = roots{1};
        else
            s.dirname_tabular = roots;
        end
    else
        s.dirname_tabular = cfg.data.dirname;
    end
    s.dirname_out     = cfg.output.dirname;
    if isfield(cfg.data, 'format_priority')
        fp = cfg.data.format_priority;
        s.data_format_priority = fp;
    end

    %% params.dependent ----------------------------------------------------
    s.params.dependent.name      = dep.fstem;
    s.params.dependent.type_data = lower(dep.datatype);

    if isfield(dep, 'dirname')
        s.params.dependent.dir_data = dep.dirname;
    end

    transformY = getField(cfg.model, 'transformY', 'none');
    if ~strcmpi(transformY, 'none')
        s.params.dependent.transform = transformY;
    end

    %% params.random -------------------------------------------------------
    if isfield(cfg.model, 'random')
        s.params.random = cfg.model.random;
        if ischar(s.params.random)
            s.params.random = {s.params.random};
        end
    else
        s.params.random = {'F', 'S', 'E'};
    end

    %% params.dir_grm (only if additive genetic effect requested) ----------
    if ismember('A', s.params.random)
        if isfield(cfg.model, 'grm') && isfield(cfg.model.grm, 'dirname')
            s.params.dir_grm = cfg.model.grm.dirname;
        else
            error(['FEMA_createInputJSON: [model.grm].dirname is required ' ...
                   'when random effects include "A" (additive genetic).']);
        end
    end

    %% params.fixed.vars ---------------------------------------------------
    vars_list  = {};   % cell array of var structs
    trans_map  = buildTransformMap(cfg);  % varname -> transform config struct

    if isfield(cfg.model, 'vars')
        raw_vars = cfg.model.vars;
        if isstruct(raw_vars)
            % toml_parse returns a struct array; convert to cell array
            raw_vars = num2cell(raw_vars);
        end
        for v = 1:numel(raw_vars)
            varStruct = struct();
            rv = raw_vars{v};
            varStruct.name = rv.name;

            if isfield(rv, 'shortname') && ~isempty(rv.shortname)
                varStruct.name_custom = rv.shortname;
            end

            % type_var is required by FEMA_parse_JSON (line 132)
            if ~isfield(rv, 'type')
                error('FEMA_createInputJSON: [[model.vars]] entry "%s" is missing the "type" field.', rv.name);
            end
            varStruct.type_var    = rv.type;
            varStruct.of_interest = getField(rv, 'of_interest', false);

            % Categorical: reference level
            if strcmpi(rv.type, 'categorical') && isfield(rv, 'reference')
                varStruct.reference = rv.reference;
            end

            % Continuous: attach transform if specified.
            % Inline fields in the [[model.vars]] block take precedence;
            % the separate [[model.transformations]] section is a fallback.
            if strcmpi(rv.type, 'continuous')
                if isfield(rv, 'transform')
                    % Inline transform — spline fields are siblings of transform
                    varStruct.transform = rv.transform;
                    if strcmpi(rv.transform, 'splines')
                        varStruct.splines = buildSplineStruct(rv);
                    elseif strcmpi(rv.transform, 'winsorize')
                        varStruct.winsorize = buildWinsorizeStruct(rv);
                    end
                elseif isKey(trans_map, rv.name)
                    % Fallback: separate [[model.transformations]] block
                    tc = trans_map(rv.name);
                    varStruct.transform = tc.type;
                    if strcmpi(tc.type, 'splines')
                        varStruct.splines = buildSplineStruct(tc);
                    elseif strcmpi(tc.type, 'winsorize')
                        varStruct.winsorize = buildWinsorizeStruct(tc);
                    end
                end
            end

            vars_list{end+1} = varStruct; %#ok<AGROW>
        end
    end

    s.params.fixed.vars = vars_list;

    % Build long→short name map for interaction substitution below
    name_sub = containers.Map('KeyType', 'char', 'ValueType', 'char');
    for vi = 1:numel(vars_list)
        vv = vars_list{vi};
        if isfield(vv, 'name_custom') && ~isempty(vv.name_custom)
            name_sub(vv.name) = vv.name_custom;
        end
    end

    %% params.fixed design options -----------------------------------------
    design = struct();
    if isfield(cfg.model, 'design')
        design = cfg.model.design;
    end

    % intercept — read by FEMA_parse_JSON from params.fixed.intercept
    s.params.fixed.intercept = logical(getField(design, 'intercept', true));

    % transform_global — read by FEMA_parse_JSON from params.fixed.transform_global
    s.params.fixed.transform_global = getField(design, 'transform_global', 'none');

    % n_gpcs — read by FEMA_parse_JSON from params.fixed.n_gpcs
    s.params.fixed.n_gpcs = getField(design, 'n_gpcs', 0);

    %% params.fixed.interaction --------------------------------------------
    if isfield(cfg.model, 'interactions')
        raw_int = cfg.model.interactions;
        if isstruct(raw_int)
            raw_int = num2cell(raw_int);
        end
        int_list = cell(numel(raw_int), 1);
        for k = 1:numel(raw_int)
            % Substitute long names with short names so the interaction
            % string matches the substituted variable names in vars_list
            raw_vars_str = raw_int{k}.vars;
            parts = strtrim(strsplit(raw_vars_str, '*'));
            for p = 1:numel(parts)
                if isKey(name_sub, parts{p})
                    parts{p} = name_sub(parts{p});
                end
            end
            int_list{k}.vars        = strjoin(parts, ' * ');
            int_list{k}.of_interest = getField(raw_int{k}, 'of_interest', false);
        end
        s.params.fixed.interaction = int_list;
    else
        s.params.fixed.interaction = {};
    end

    %% params.advanced -----------------------------------------------------
    adv = struct();
    if isfield(cfg, 'advanced')
        adv = cfg.advanced;
    end

    s.params.advanced.type_perm      = getField(adv, 'type_perm',      'none');
    s.params.advanced.n_perm         = getField(adv, 'n_perm',         0);
    s.params.advanced.type_cov       = getField(adv, 'type_cov',       'analytic');
    s.params.advanced.type_fixed_est = getField(adv, 'type_fixed_est', 'gls');
    s.params.advanced.n_bins         = getField(adv, 'n_bins',         20);
    s.params.advanced.precision      = getField(adv, 'precision',      'double');

    %% params.wrapper_args — extra FEMA_wrapper options from [options] -----
    % Any field in the TOML [options] section is forwarded verbatim as a
    % key-value pair.  FEMA_parseInputs reads this struct and appends each
    % entry to extraArgs so it reaches the FEMA_wrapper inputParser.
    % Empty-string values (used as placeholders in the TOML) are skipped.
    if isfield(cfg, 'options') && ~isempty(fieldnames(cfg.options))
        ff = fieldnames(cfg.options);
        wrapper_args = struct();
        for k = 1:numel(ff)
            val = cfg.options.(ff{k});
            % Skip empty-string placeholders (e.g. outPrefix = "")
            if ischar(val) && isempty(strtrim(val))
                continue
            end
            wrapper_args.(ff{k}) = val;
        end
        if ~isempty(fieldnames(wrapper_args))
            s.params.wrapper_args = wrapper_args;
        end
    end

    %% params.makeDesign — optional FEMA_makeDesign name-values (see FEMA_mergeArgs)
    if isfield(cfg, 'makeDesign') && isstruct(cfg.makeDesign) && ~isempty(fieldnames(cfg.makeDesign))
        s.params.makeDesign = cfg.makeDesign;
    end

    %% Boilerplate fields expected by the JSON format ----------------------
    s.params.model    = 'glmm';
    s.params.contrast = {};
    s.params.query    = 'non-DEAP: variables loaded from dirname_tabular';

    %% Write JSON to file --------------------------------------------------
    if ~exist(outDir, 'dir')
        mkdir(outDir);
    end

    fname_json = fullfile(outDir, outName);
    json_txt   = save_jsonencode(s, 'PrettyPrint', true);

    fid = fopen(fname_json, 'w');
    if fid == -1
        error('FEMA_createInputJSON: cannot open file for writing: %s', fname_json);
    end
    fwrite(fid, json_txt);
    fclose(fid);

    logging('FEMA_createInputJSON: JSON written to %s\n', fname_json);
end


%% =========================================================================
%  Helper functions
%% =========================================================================

function roots = uniqueRootPaths(rootsIn)
% Unique directory paths; first occurrence wins (dependent before data when built that way).
    roots = {};
    if isempty(rootsIn)
        return
    end
    if ischar(rootsIn)
        rootsIn = {rootsIn};
    elseif ~iscell(rootsIn)
        rootsIn = cellstr(rootsIn);
    end
    seen = containers.Map('KeyType', 'char', 'ValueType', 'logical');
    for i = 1:numel(rootsIn)
        if isempty(rootsIn{i}), continue; end
        p = char(rootsIn{i});
        if ~isKey(seen, p)
            seen(p) = true;
            roots{end+1} = p; %#ok<AGROW>
        end
    end
end


function v = getField(s, fname, default)
% Return s.(fname) if it exists, otherwise return default.
    if isstruct(s) && isfield(s, fname)
        v = s.(fname);
    else
        v = default;
    end
end


function trans_map = buildTransformMap(cfg)
% Build a containers.Map from variable name -> transformation config struct.
    trans_map = containers.Map('KeyType', 'char', 'ValueType', 'any');

    if ~isfield(cfg.model, 'transformations')
        return
    end

    raw = cfg.model.transformations;
    if isstruct(raw)
        raw = num2cell(raw);
    end

    for k = 1:numel(raw)
        tc = raw{k};
        if ~isfield(tc, 'var') || ~isfield(tc, 'type')
            warning('FEMA_createInputJSON: skipping transformation entry %d — missing "var" or "type" field.', k);
            continue
        end
        varname = tc.var;
        if isKey(trans_map, varname)
            warning('FEMA_createInputJSON: duplicate transformation for variable "%s"; using first definition.', varname);
        else
            trans_map(varname) = tc;
        end
    end
end


function sp = buildSplineStruct(tc)
% Build the splines sub-object for a variable JSON entry.
%
% Maps TOML field names to the field names that FEMA_parse_JSON's
% standardizeSplineTransforms expects (case-insensitive via strcmpi):
%   knots, type_spline, Xpowers, method, minMax, intercept
%
% The TOML field type_spline -> JSON field type_spline (matched by strcmpi)

    sp = struct();

    % type_spline (TOML: type_spline) — maps to outName{4} = 'type_spline'
    if isfield(tc, 'type_spline')
        sp.type_spline = tc.type_spline;
    end

    % knots
    if isfield(tc, 'knots')
        sp.knots = tc.knots;
    end

    % method
    if isfield(tc, 'method')
        sp.method = tc.method;
    end

    % Xpowers (TOML: xpowers)
    if isfield(tc, 'xpowers')
        sp.Xpowers = tc.xpowers;
    end

    % minMax (TOML: minmax)
    if isfield(tc, 'minmax')
        sp.minMax = tc.minmax;
    end

    % intercept
    if isfield(tc, 'intercept')
        sp.intercept = logical(tc.intercept);
    end
end


function wz = buildWinsorizeStruct(tc)
% Build the winsorize sub-object expected by FEMA_parse_JSON:
%   winsorize.minmax
    wz = struct();
    if isfield(tc, 'minmax')
        wz.minmax = tc.minmax;
    elseif isfield(tc, 'lower') && isfield(tc, 'upper')
        % Backward-compatible fallback for older configs.
        wz.minmax = [tc.lower, tc.upper];
    else
        wz.minmax = [];
    end
end
