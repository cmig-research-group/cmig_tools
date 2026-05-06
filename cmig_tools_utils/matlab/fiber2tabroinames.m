function [roinames_out, roicodes] = fiber2tabroinames(roinames_in, direction)
    % Maps fiber tract ROI names between atlas and tabulated naming
    %
    % Inputs:
    %   roinames_in - cell array of ROI names
    %   direction   - mapping direction (optional):
    %                 'atlas2tab' (default) or 'tab2atlas'
    %
    % Outputs:
    %   roinames_out - cell array of mapped ROI names
    %   roicodes     - numeric ROI codes from support_files/fiber_lut.csv
    %
    % Long column-style names (e.g. mr_y_rsi__hnd__at__atr__lh_wmean) are
    % resolved via the same embedded-token logic as aparc/aseg mappers; tab
    % tokens may be followed by _wmean / _mean-style suffixes.

    if nargin < 2 || isempty(direction)
        direction = 'atlas2tab';
    end
    validatestring(direction, {'atlas2tab', 'tab2atlas'}, mfilename, 'direction', 2);

    [fiber_map, map_tab2atlas] = local_loadFiberMaps();
    map_use = fiber_map;
    if strcmp(direction, 'tab2atlas')
        map_use = map_tab2atlas;
    end

    roinames_out = cell(size(roinames_in));
    roicodes = nan(size(roinames_in));
    code_map = local_loadCodeMap();

    for i = 1:length(roinames_in)
        name = roinames_in{i};
        if isempty(name)
            roinames_out{i} = '';
            roicodes(i) = NaN;
            continue
        end

        if isKey(map_use, name)
            roinames_out{i} = map_use(name);
        else
            [mapped_name, tf] = local_mapEmbeddedToken(name, map_use, direction);
            if tf
                roinames_out{i} = mapped_name;
            else
                roinames_out{i} = local_mapFiberName(name, direction);
            end
        end

        if isempty(roinames_out{i})
            warning('Unknown fiber ROI name for direction %s: %s', direction, name);
            roinames_out{i} = '';
        end

        atlas_name = local_atlasNameForCode(name, roinames_out{i}, direction, fiber_map, map_tab2atlas);
        roicodes(i) = local_codeFromName(atlas_name, code_map);
    end
end

function tab = local_atlasToTab(atlas_name)
    tab = '';
    switch atlas_name
        case 'Fmaj'
            tab = 'fmaj';
            return
        case 'Fmin'
            tab = 'fmin';
            return
        case 'CC'
            tab = 'cc';
            return
        case 'AllFibers'
            tab = 'at_wmean';
            return
    end
    if startsWith(atlas_name, 'L_')
        fiber_name = atlas_name(3:end);
        if strcmp(fiber_name, 'AllFib')
            tab = 'lh';
        elseif strcmp(fiber_name, 'AllFibnoCC')
            tab = 'nocc__lh';
        else
            tab = [lower(fiber_name) '__lh'];
        end
        return
    end
    if startsWith(atlas_name, 'R_')
        fiber_name = atlas_name(3:end);
        if strcmp(fiber_name, 'AllFib')
            tab = 'rh';
        elseif strcmp(fiber_name, 'AllFibnoCC')
            tab = 'nocc__rh';
        else
            tab = [lower(fiber_name) '__rh'];
        end
        return
    end
end

function name_out = local_mapFiberName(name_in, direction)
    % Fallback for short forms not produced by LUT build; no warning here.
    name_out = '';
    switch direction
        case 'atlas2tab'
            name_out = local_atlasToTab(name_in);
        case 'tab2atlas'
            if ~isempty(regexp(lower(name_in), '(^|__)at_wmean(?=_|$)', 'once'))
                name_out = 'AllFibers';
                return
            end
            switch lower(name_in)
                case 'fmaj'
                    name_out = 'Fmaj';
                    return
                case 'fmin'
                    name_out = 'Fmin';
                    return
                case 'cc'
                    name_out = 'CC';
                    return
                case 'all'
                    name_out = 'AllFibers';
                    return
                case 'at_wmean'
                    name_out = 'AllFibers';
                    return
                case 'lh'
                    name_out = 'L_AllFib';
                    return
                case 'rh'
                    name_out = 'R_AllFib';
                    return
                case 'nocc__lh'
                    name_out = 'L_AllFibnoCC';
                    return
                case 'nocc__rh'
                    name_out = 'R_AllFibnoCC';
                    return
            end
            if endsWith(name_in, '__lh')
                fiber_name = name_in(1:end-4);
                name_out = ['L_' fiber_name];
                return
            end
            if endsWith(name_in, '__rh')
                fiber_name = name_in(1:end-4);
                name_out = ['R_' fiber_name];
                return
            end
    end
end

function atlas_name = local_atlasNameForCode(name_in, name_out, direction, map_atlas2tab, map_tab2atlas)
    atlas_name = '';
    if strcmp(direction, 'tab2atlas')
        atlas_name = name_out;
        return
    end
    if isKey(map_atlas2tab, name_in)
        atlas_name = name_in;
        return
    end
    if isKey(map_tab2atlas, name_out)
        atlas_name = map_tab2atlas(name_out);
        return
    end
    [tmp_name, tf] = local_mapEmbeddedToken(name_out, map_tab2atlas, 'tab2atlas');
    if tf
        atlas_name = tmp_name;
    else
        atlas_name = name_in;
    end
end

function [fiber_map, map_tab2atlas] = local_loadFiberMaps()
    persistent fm m2a
    if isempty(fm)
        this_dir = fileparts(mfilename('fullpath'));
        repo_root = fileparts(fileparts(this_dir));
        lut_path = fullfile(repo_root, 'support_files', 'fiber_lut.csv');
        lut = readtable(lut_path, 'TextType', 'string');
        atlas_list = cellstr(lut.roinames);
        fiber_map = containers.Map('KeyType', 'char', 'ValueType', 'char');
        for ii = 1:numel(atlas_list)
            a = atlas_list{ii};
            t = local_atlasToTab(a);
            if isempty(t)
                error('fiber2tabroinames:LUT', ...
                    'No atlas2tab rule for fiber ROI ''%s'' in %s.', a, lut_path);
            end
            fiber_map(a) = t;
        end
        ks = keys(fiber_map);
        vs = values(fiber_map);
        map_tab2atlas = containers.Map('KeyType', 'char', 'ValueType', 'char');
        for ii = 1:numel(vs)
            map_tab2atlas(vs{ii}) = ks{ii};
        end
        % Whole-brain / aggregate aliases used in some Y column naming (vs at_wmean from atlas2tab)
        map_tab2atlas('at_wmean') = 'AllFibers';
        map_tab2atlas('all') = 'AllFibers';
        fm = fiber_map;
        m2a = map_tab2atlas;
    end
    fiber_map = fm;
    map_tab2atlas = m2a;
end

function code_map = local_loadCodeMap()
persistent map_cache
if isempty(map_cache)
    this_dir = fileparts(mfilename('fullpath'));
    repo_root = fileparts(fileparts(this_dir));
    lut_path = fullfile(repo_root, 'support_files', 'fiber_lut.csv');
    lut = readtable(lut_path, 'TextType', 'string');
    map_cache = containers.Map(lower(cellstr(lut.roinames)), num2cell(double(lut.roicodes)));
end
code_map = map_cache;
end

function code = local_codeFromName(name, code_map)
code = NaN;
if isempty(name)
    return
end
key = char(strtrim(lower(string(name))));
if isKey(code_map, key)
    code = code_map(key);
end
end
