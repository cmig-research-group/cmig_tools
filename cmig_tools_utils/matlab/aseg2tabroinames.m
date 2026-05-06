function [roinames_out, roicodes] = aseg2tabroinames(roinames_in, direction)
    % Maps Freesurfer aseg ROI names between atlas and tabulated naming
    %
    % Inputs:
    %   roinames_in - cell array of ROI names 
    %   direction   - mapping direction (optional):
    %                 'atlas2tab' (default) or 'tab2atlas'
    %
    % Outputs:
    %   roinames_out - cell array of mapped ROI names 
    %                  (e.g., 'th__lh', 'hc__rh', '3rdv')
    %   roicodes     - numeric ROI codes from aseg_lut.csv

    if nargin < 2 || isempty(direction)
        direction = 'tab2atlas';
    end
    validatestring(direction, {'atlas2tab', 'tab2atlas'}, mfilename, 'direction', 2);
    
    % Define mapping from aseg names to tabulated abbreviations
    aseg_map = containers.Map();
    
    % Left hemisphere structures
    aseg_map('Left-Cerebral-White-Matter') = 'cwm__lh';
    aseg_map('Left-Cerebral-Cortex') = 'cc__lh';
    aseg_map('Left-Lateral-Ventricle') = 'lv__lh';
    aseg_map('Left-Inf-Lat-Vent') = 'ilv__lh';
    aseg_map('Left-Cerebellum-White-Matter') = 'cbwm__lh';
    aseg_map('Left-Cerebellum-Cortex') = 'cbc__lh';
    aseg_map('Left-Thalamus') = 'th__lh';
    aseg_map('Left-Caudate') = 'cd__lh';
    aseg_map('Left-Putamen') = 'pt__lh';
    aseg_map('Left-Pallidum') = 'pl__lh';
    aseg_map('Left-Hippocampus') = 'hc__lh';
    aseg_map('Left-Amygdala') = 'ag__lh';
    aseg_map('Left-Accumbens-area') = 'ab__lh';
    aseg_map('Left-VentralDC') = 'vdc__lh';
    
    % Right hemisphere structures
    aseg_map('Right-Cerebral-White-Matter') = 'cwm__rh';
    aseg_map('Right-Cerebral-Cortex') = 'cc__rh';
    aseg_map('Right-Lateral-Ventricle') = 'lv__rh';
    aseg_map('Right-Inf-Lat-Vent') = 'ilv__rh';
    aseg_map('Right-Cerebellum-White-Matter') = 'cbwm__rh';
    aseg_map('Right-Cerebellum-Cortex') = 'cbc__rh';
    aseg_map('Right-Thalamus') = 'th__rh';
    aseg_map('Right-Caudate') = 'cd__rh';
    aseg_map('Right-Putamen') = 'pt__rh';
    aseg_map('Right-Pallidum') = 'pl__rh';
    aseg_map('Right-Hippocampus') = 'hc__rh';
    aseg_map('Right-Amygdala') = 'ag__rh';
    aseg_map('Right-Accumbens-area') = 'ab__rh';
    aseg_map('Right-VentralDC') = 'vdc__rh';
    
    % Midline/bilateral structures
    aseg_map('3rd-Ventricle') = '3rdv';
    aseg_map('4th-Ventricle') = '4thv';
    aseg_map('Brain-Stem') = 'bs';
    aseg_map('CSF') = 'csf';
    
    % Corpus Callosum
    aseg_map('CC_Posterior') = 'ccp';
    aseg_map('CC_Mid_Posterior') = 'ccpm';
    aseg_map('CC_Central') = 'ccc';
    aseg_map('CC_Mid_Anterior') = 'ccam';
    aseg_map('CC_Anterior') = 'cca';
    
    % Whole-cortex / aggregate metric token (e.g. ...gm__dsk_mean) — not __dsk__<parcel>__...
    aseg_map('all') = 'dsk_mean';
    
    % Structures to exclude (not in tabulated data)
    exclude_set = {'Left-vessel', 'Right-vessel', 'Left-choroid-plexus', ...
                   'Right-choroid-plexus', '5th-Ventricle', 'WM-hypointensities', ...
                   'non-WM-hypointensities', 'Optic-Chiasm'};
    
    map_use = aseg_map;
    map_tab2atlas = containers.Map(values(aseg_map), keys(aseg_map));
    exclude_use = exclude_set;
    if strcmp(direction, 'tab2atlas')
        map_use = map_tab2atlas;
        exclude_use = {};
    end

    roinames_out = cell(size(roinames_in));
    roicodes = nan(size(roinames_in));
    code_map = local_loadCodeMap();
    
    for i = 1:length(roinames_in)
        name = roinames_in{i};
        
        % Handle empty cells {0x0 char}
        if isempty(name)
            roinames_out{i} = '';
            roicodes(i) = NaN;
            continue
        end
        
        % Check if name is in exclusion list
        if ismember(name, exclude_use)
            roinames_out{i} = '';
            roicodes(i) = NaN;
            continue
        end
        
        % Look up in mapping
        if isKey(map_use, name)
            roinames_out{i} = map_use(name);
        else
            [mapped_name, tf] = local_mapEmbeddedToken(name, map_use, direction);
            if tf
                roinames_out{i} = mapped_name;
            else
                % Unknown name - return empty
                warning('Unknown aseg ROI name for direction %s: %s', direction, name);
                roinames_out{i} = '';
            end
        end

        atlas_name = local_atlasNameForCode(name, roinames_out{i}, direction, aseg_map, map_tab2atlas);
        roicodes(i) = local_codeFromName(atlas_name, code_map);
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

function code_map = local_loadCodeMap()
persistent map_cache
if isempty(map_cache)
    this_dir = fileparts(mfilename('fullpath'));
    repo_root = fileparts(fileparts(this_dir));
    lut_path = fullfile(repo_root, 'support_files', 'aseg_lut.csv');
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
