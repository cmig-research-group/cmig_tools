function [fname_json, fname_data, dirname_out_val, study_val, timing_rt] = FEMA_readToml(fname_toml, varargin)
% FEMA_READTOML  Convert non-DEAP TOML to job JSON and merged tabular parquet.
%
% Shared by FEMA_wrapper('toml', ...) and FEMA_makeDesign when given a .toml config.
% Optional name-value arguments are forwarded to FEMA_createInputJSON (e.g. 'outDir',
% 'outName', 'designOnly', true). FEMA_wrapper calls this with no extras; FEMA_makeDesign
% passes designOnly=true so [dependent] may be omitted in the TOML.
%
% See FEMA_createInputJSON for TOML format.
%
% Outputs:
%   fname_json       Path to written JSON specification
%   fname_data       Path to merged parquet from FEMA_loadDataDir
%   dirname_out_val  [output].dirname from the TOML (also used as default JSON outDir)
%   study_val        Top-level study string from the generated JSON (e.g. 'abcd')
%   timing_rt        struct with .tCreateInputJSON, .tLoadDataDir, .tReadTomlTotal (seconds)

    timing_rt = struct('tCreateInputJSON', [], 'tLoadDataDir', [], 'tReadTomlTotal', []);

    tReadTomlAll = tic;

    if ~exist(fname_toml, 'file')
        error('FEMA_readToml: TOML file does not exist: %s', fname_toml);
    end

    tJSON = tic;
    fname_json = FEMA_createInputJSON(fname_toml, varargin{:});
    timing_rt.tCreateInputJSON = toc(tJSON);

    tmp_cfg = jsondecode(fileread(fname_json));
    dirname_tabular = tmp_cfg.dirname_tabular;
    dirname_out_val = tmp_cfg.dirname_out;

    vars_json = tmp_cfg.params.fixed.vars;
    if isstruct(vars_json)
        vars_json = num2cell(vars_json);
    end
    FFX_names_tmp = cellfun(@(v) v.name, vars_json, 'UniformOutput', false);
    FFX_names_tmp = reshape(FFX_names_tmp, 1, []);

    dep_type = strrep(tmp_cfg.params.dependent.type_data, 'wise', '');
    if strcmpi(dep_type, 'external')
        dn = tmp_cfg.params.dependent.name;
        if isempty(dn)
            dep_cells = {};
        elseif ischar(dn)
            dep_cells = {dn};
        elseif iscell(dn)
            dep_cells = dn(:)';
        elseif isstring(dn)
            dep_cells = cellstr(dn(:)');
        else
            dep_cells = {char(dn)};
        end
        FFX_names_tmp = [dep_cells, FFX_names_tmp];
        FFX_names_tmp = unique(FFX_names_tmp, 'stable');
    end

    name_map = containers.Map('KeyType', 'char', 'ValueType', 'char');
    var_types = containers.Map('KeyType', 'char', 'ValueType', 'char');
    for v = 1:numel(vars_json)
        vv = vars_json{v};
        outCol = vv.name;
        if isfield(vv, 'name_custom') && ~isempty(vv.name_custom)
            name_map(vv.name) = vv.name_custom;
            outCol = vv.name_custom;
        end
        if isfield(vv, 'type_var')
            var_types(outCol) = lower(strtrim(vv.type_var));
        end
    end

    n_gpcs_tmp = 0;
    if isfield(tmp_cfg.params.fixed, 'n_gpcs')
        n_gpcs_tmp = tmp_cfg.params.fixed.n_gpcs;
    end
    if n_gpcs_tmp > 0
        pc_names = arrayfun(@(k) sprintf('ab_g_stc__gen_pc__%02d', k), ...
                            1:n_gpcs_tmp, 'UniformOutput', false);
        FFX_names_tmp = [FFX_names_tmp, pc_names];
        for k = 1:n_gpcs_tmp
            name_map(sprintf('ab_g_stc__gen_pc__%02d', k)) = sprintf('pc_%02d', k);
        end
    end

    loadDataArgs = {};
    if isfield(tmp_cfg, 'data_format_priority')
        fp = tmp_cfg.data_format_priority;
        if ischar(fp)
            fp = {fp};
        end
        loadDataArgs = {'format_priority', fp};
    end
    tLoad = tic;
    [~, fname_data] = FEMA_loadDataDir(dirname_tabular, FFX_names_tmp, 'outDir', dirname_out_val, ...
        'name_map', name_map, 'var_types', var_types, loadDataArgs{:});
    timing_rt.tLoadDataDir = toc(tLoad);

    if isfield(tmp_cfg, 'study')
        study_val = tmp_cfg.study;
        if isstring(study_val)
            study_val = char(study_val);
        end
    else
        study_val = '';
    end

    timing_rt.tReadTomlTotal = toc(tReadTomlAll);
end
