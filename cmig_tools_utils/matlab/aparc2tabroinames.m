function [roinames_out, roicodes] = aparc2tabroinames(roinames_in, direction)
    % Maps Freesurfer aparc (Desikan-Killiany) ROI names between atlas and tabulated naming
    %
    % Inputs:
    %   roinames_in - cell array of ROI names
    %   direction   - mapping direction (optional):
    %                 'atlas2tab' (default) or 'tab2atlas'
    %
    % Outputs:
    %   roinames_out - cell array of mapped ROI names 
    %                  (e.g., 'bstmps__lh', 'cn__rh')
    %   roicodes     - numeric ROI codes from support_files/aparc_lut.csv

    if nargin < 2 || isempty(direction)
        direction = 'tab2atlas';
    end
    validatestring(direction, {'atlas2tab', 'tab2atlas'}, mfilename, 'direction', 2);
    
    % Define mapping from aparc names to tabulated abbreviations
    % Ordered alphabetically by output abbreviation names
    aparc_map = containers.Map();
    
    % Left hemisphere structures (sorted alphabetically by output abbreviation)
    aparc_map('bankssts_lh') = 'bstmps__lh';
    aparc_map('caudalanteriorcingulate_lh') = 'cac__lh';
    aparc_map('corpuscallosum_lh') = 'cc__lh';
    aparc_map('caudalmiddlefrontal_lh') = 'cmfrt__lh';
    aparc_map('cuneus_lh') = 'cn__lh';
    aparc_map('entorhinal_lh') = 'er__lh';
    aparc_map('fusiform_lh') = 'ff__lh';
    aparc_map('isthmuscingulate_lh') = 'ic__lh';
    aparc_map('insula_lh') = 'ins__lh';
    aparc_map('inferiorparietal_lh') = 'iprt__lh';
    aparc_map('inferiortemporal_lh') = 'itmp__lh';
    aparc_map('lingual_lh') = 'lg__lh';
    aparc_map('lateralorbitofrontal_lh') = 'lobfrt__lh';
    aparc_map('lateraloccipital_lh') = 'locc__lh';
    aparc_map('medialorbitofrontal_lh') = 'mobfrt__lh';
    aparc_map('middletemporal_lh') = 'mtmp__lh';
    aparc_map('paracentral_lh') = 'pactr__lh';
    aparc_map('pericalcarine_lh') = 'pcc__lh';
    aparc_map('posteriorcingulate_lh') = 'pcg__lh';
    aparc_map('frontalpole_lh') = 'pfrt__lh';
    aparc_map('parahippocampal_lh') = 'ph__lh';
    aparc_map('parsorbitalis_lh') = 'pob__lh';
    aparc_map('postcentral_lh') = 'poctr__lh';
    aparc_map('parsopercularis_lh') = 'pop__lh';
    aparc_map('precuneus_lh') = 'prcn__lh';
    aparc_map('precentral_lh') = 'prctr__lh';
    aparc_map('parstriangularis_lh') = 'ptg__lh';
    aparc_map('temporalpole_lh') = 'ptmp__lh';
    aparc_map('rostralanteriorcingulate_lh') = 'rac__lh';
    aparc_map('rostralmiddlefrontal_lh') = 'rmfrt__lh';
    aparc_map('superiorfrontal_lh') = 'sfrt__lh';
    aparc_map('supramarginal_lh') = 'sm__lh';
    aparc_map('superiorparietal_lh') = 'sprt__lh';
    aparc_map('superiortemporal_lh') = 'stmp__lh';
    aparc_map('transversetemporal_lh') = 'ttmp__lh';
    
    % Right hemisphere structures (matching lh, sorted alphabetically by output abbreviation)
    aparc_map('bankssts_rh') = 'bstmps__rh';
    aparc_map('caudalanteriorcingulate_rh') = 'cac__rh';
    aparc_map('corpuscallosum_rh') = 'cc__rh';
    aparc_map('caudalmiddlefrontal_rh') = 'cmfrt__rh';
    aparc_map('cuneus_rh') = 'cn__rh';
    aparc_map('entorhinal_rh') = 'er__rh';
    aparc_map('fusiform_rh') = 'ff__rh';
    aparc_map('isthmuscingulate_rh') = 'ic__rh';
    aparc_map('insula_rh') = 'ins__rh';
    aparc_map('inferiorparietal_rh') = 'iprt__rh';
    aparc_map('inferiortemporal_rh') = 'itmp__rh';
    aparc_map('lingual_rh') = 'lg__rh';
    aparc_map('lateralorbitofrontal_rh') = 'lobfrt__rh';
    aparc_map('lateraloccipital_rh') = 'locc__rh';
    aparc_map('medialorbitofrontal_rh') = 'mobfrt__rh';
    aparc_map('middletemporal_rh') = 'mtmp__rh';
    aparc_map('paracentral_rh') = 'pactr__rh';
    aparc_map('pericalcarine_rh') = 'pcc__rh';
    aparc_map('posteriorcingulate_rh') = 'pcg__rh';
    aparc_map('frontalpole_rh') = 'pfrt__rh';
    aparc_map('parahippocampal_rh') = 'ph__rh';
    aparc_map('parsorbitalis_rh') = 'pob__rh';
    aparc_map('postcentral_rh') = 'poctr__rh';
    aparc_map('parsopercularis_rh') = 'pop__rh';
    aparc_map('precuneus_rh') = 'prcn__rh';
    aparc_map('precentral_rh') = 'prctr__rh';
    aparc_map('parstriangularis_rh') = 'ptg__rh';
    aparc_map('temporalpole_rh') = 'ptmp__rh';
    aparc_map('rostralanteriorcingulate_rh') = 'rac__rh';
    aparc_map('rostralmiddlefrontal_rh') = 'rmfrt__rh';
    aparc_map('superiorfrontal_rh') = 'sfrt__rh';
    aparc_map('supramarginal_rh') = 'sm__rh';
    aparc_map('superiorparietal_rh') = 'sprt__rh';
    aparc_map('superiortemporal_rh') = 'stmp__rh';
    aparc_map('transversetemporal_rh') = 'ttmp__rh';
    
    % Special cases
    aparc_map('unknown') = 'unknown';
    % Whole-cortex DK metric (e.g. ...gm__dsk_mean) — not __dsk__<parcel>__...
    aparc_map('all') = 'dsk_mean';
    aparc_map('lh') = 'lh';
    aparc_map('rh') = 'rh';
    
    map_use = aparc_map;
    map_tab2atlas = containers.Map(values(aparc_map), keys(aparc_map));
    if strcmp(direction, 'tab2atlas')
        map_use = map_tab2atlas;
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
        
        % Look up in mapping
        if isKey(map_use, name)
            roinames_out{i} = map_use(name);
        else
            [mapped_name, tf] = local_mapEmbeddedToken(name, map_use, direction);
            if tf
                roinames_out{i} = mapped_name;
            else
                % Unknown name - return empty or 'unknown'
                warning('Unknown aparc ROI name for direction %s: %s', direction, name);
                roinames_out{i} = '';
            end
        end

        atlas_name = local_atlasNameForCode(name, roinames_out{i}, direction, aparc_map, map_tab2atlas);
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
    lut_path = fullfile(repo_root, 'support_files', 'aparc_lut.csv');
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
