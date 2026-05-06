function [roinames_out, roicodes] = aparc2009s2tabroinames(roinames_in, direction)
    % Maps Freesurfer aparc.a2009s ROI names between atlas and tabulated naming
    %
    % Inputs:
    %   roinames_in - cell array of ROI names
    %   direction   - mapping direction (optional):
    %                 'atlas2tab' (default) or 'tab2atlas'
    %
    % Outputs:
    %   roinames_out - cell array of mapped ROI names 
    %   roicodes     - numeric ROI codes from support_files/aparc_a2009s_lut.csv

    if nargin < 2 || isempty(direction)
        direction = 'tab2atlas';
    end
    validatestring(direction, {'atlas2tab', 'tab2atlas'}, mfilename, 'direction', 2);
    
    % Define mapping from aparc.a2009s names to tabulated abbreviations
    % Sorted alphabetically by output abbreviation
    aparc2009s_map = containers.Map();
    
    % Left hemisphere structures (alphabetically by output abbreviation)
    aparc2009s_map('G_cuneus_lh') = 'gcn__lh';
    aparc2009s_map('G_cingul-Post-dorsal_lh') = 'gcpd__lh';
    aparc2009s_map('G_cingul-Post-ventral_lh') = 'gcpv__lh';
    aparc2009s_map('G_front_inf-Opercular_lh') = 'gfio__lh';
    aparc2009s_map('G_front_inf-Orbital_lh') = 'gfiob__lh';
    aparc2009s_map('G_front_inf-Triangul_lh') = 'gfit__lh';
    aparc2009s_map('G_front_middle_lh') = 'gfm__lh';
    aparc2009s_map('G_front_sup_lh') = 'gfs__lh';
    aparc2009s_map('G_Ins_lg_and_S_cent_ins_lh') = 'gilsci__lh';
    aparc2009s_map('G_insular_short_lh') = 'gis__lh';
    aparc2009s_map('G_orbital_lh') = 'go__lh';
    aparc2009s_map('G_occipital_middle_lh') = 'gom__lh';
    aparc2009s_map('G_occipital_sup_lh') = 'gos__lh';
    aparc2009s_map('G_oc-temp_lat-fusifor_lh') = 'gotlf__lh';
    aparc2009s_map('G_oc-temp_med-Lingual_lh') = 'gotml__lh';
    aparc2009s_map('G_oc-temp_med-Parahip_lh') = 'gotmp__lh';
    aparc2009s_map('G_precentral_lh') = 'gpc__lh';
    aparc2009s_map('G_pariet_inf-Angular_lh') = 'gpia__lh';
    aparc2009s_map('G_pariet_inf-Supramar_lh') = 'gpis__lh';
    aparc2009s_map('G_precuneus_lh') = 'gprcn__lh';
    aparc2009s_map('G_postcentral_lh') = 'gprct__lh';
    aparc2009s_map('G_parietal_sup_lh') = 'gps__lh';
    aparc2009s_map('G_rectus_lh') = 'gr__lh';
    aparc2009s_map('G_subcallosal_lh') = 'gs__lh';
    aparc2009s_map('G_and_S_cingul-Ant_lh') = 'gsca__lh';
    aparc2009s_map('G_and_S_cingul-Mid-Ant_lh') = 'gscma__lh';
    aparc2009s_map('G_and_S_cingul-Mid-Post_lh') = 'gscmp__lh';
    aparc2009s_map('G_and_S_frontomargin_lh') = 'gsfm__lh';
    aparc2009s_map('G_and_S_occipital_inf_lh') = 'gsoi__lh';
    aparc2009s_map('G_and_S_paracentral_lh') = 'gspc__lh';
    aparc2009s_map('G_and_S_subcentral_lh') = 'gssc__lh';
    aparc2009s_map('G_and_S_transv_frontopol_lh') = 'gstf__lh';
    aparc2009s_map('G_temporal_inf_lh') = 'gti__lh';
    aparc2009s_map('G_temporal_middle_lh') = 'gtm__lh';
    aparc2009s_map('G_temp_sup-G_T_transv_lh') = 'gtsgtt__lh';
    aparc2009s_map('G_temp_sup-Lateral_lh') = 'gtsl__lh';
    aparc2009s_map('G_temp_sup-Plan_polar_lh') = 'gtspp__lh';
    aparc2009s_map('G_temp_sup-Plan_tempo_lh') = 'gtspt__lh';
    aparc2009s_map('Lat_Fis-ant-Horizont_lh') = 'lfah__lh';
    aparc2009s_map('Lat_Fis-ant-Vertical_lh') = 'lfav__lh';
    aparc2009s_map('Lat_Fis-post_lh') = 'lfp__lh';
    aparc2009s_map('Medial_wall_lh') = 'mw__lh';
    aparc2009s_map('Pole_occipital_lh') = 'pocc__lh';
    aparc2009s_map('Pole_temporal_lh') = 'ptmp__lh';
    aparc2009s_map('S_calcarine_lh') = 'scc__lh';
    aparc2009s_map('S_central_lh') = 'sc__lh';
    aparc2009s_map('S_circular_insula_ant_lh') = 'scia__lh';
    aparc2009s_map('S_circular_insula_inf_lh') = 'scii__lh';
    aparc2009s_map('S_circular_insula_sup_lh') = 'scis__lh';
    aparc2009s_map('S_cingul-Marginalis_lh') = 'scm__lh';
    aparc2009s_map('S_collat_transv_ant_lh') = 'scta__lh';
    aparc2009s_map('S_collat_transv_post_lh') = 'sctp__lh';
    aparc2009s_map('S_front_inf_lh') = 'sfi__lh';
    aparc2009s_map('S_front_middle_lh') = 'sfm__lh';
    aparc2009s_map('S_front_sup_lh') = 'sfs__lh';
    aparc2009s_map('S_interm_prim-Jensen_lh') = 'sipj__lh';
    aparc2009s_map('S_intrapariet_and_P_trans_lh') = 'sipt__lh';
    aparc2009s_map('S_occipital_ant_lh') = 'soa__lh';
    aparc2009s_map('S_orbital-H_Shaped_lh') = 'sohs__lh';
    aparc2009s_map('S_orbital_lateral_lh') = 'sol__lh';
    aparc2009s_map('S_oc_middle_and_Lunatus_lh') = 'soml__lh';
    aparc2009s_map('S_orbital_med-olfact_lh') = 'somo__lh';
    aparc2009s_map('S_oc_sup_and_transversal_lh') = 'sost__lh';
    aparc2009s_map('S_oc-temp_lat_lh') = 'sotl__lh';
    aparc2009s_map('S_oc-temp_med_and_Lingual_lh') = 'sotml__lh';
    aparc2009s_map('S_pericallosal_lh') = 'spc__lh';
    aparc2009s_map('S_postcentral_lh') = 'spct__lh';
    aparc2009s_map('S_parieto_occipital_lh') = 'spo__lh';
    aparc2009s_map('S_precentral-inf-part_lh') = 'sprip__lh';
    aparc2009s_map('S_precentral-sup-part_lh') = 'sprsp__lh';
    aparc2009s_map('S_suborbital_lh') = 'sso__lh';
    aparc2009s_map('S_subparietal_lh') = 'ssp__lh';
    aparc2009s_map('S_temporal_inf_lh') = 'sti__lh';
    aparc2009s_map('S_temporal_sup_lh') = 'sts__lh';
    aparc2009s_map('S_temporal_transverse_lh') = 'stt__lh';
    
    % Right hemisphere structures (alphabetically by output abbreviation, matching left)
    aparc2009s_map('G_cuneus_rh') = 'gcn__rh';
    aparc2009s_map('G_cingul-Post-dorsal_rh') = 'gcpd__rh';
    aparc2009s_map('G_cingul-Post-ventral_rh') = 'gcpv__rh';
    aparc2009s_map('G_front_inf-Opercular_rh') = 'gfio__rh';
    aparc2009s_map('G_front_inf-Orbital_rh') = 'gfiob__rh';
    aparc2009s_map('G_front_inf-Triangul_rh') = 'gfit__rh';
    aparc2009s_map('G_front_middle_rh') = 'gfm__rh';
    aparc2009s_map('G_front_sup_rh') = 'gfs__rh';
    aparc2009s_map('G_Ins_lg_and_S_cent_ins_rh') = 'gilsci__rh';
    aparc2009s_map('G_insular_short_rh') = 'gis__rh';
    aparc2009s_map('G_orbital_rh') = 'go__rh';
    aparc2009s_map('G_occipital_middle_rh') = 'gom__rh';
    aparc2009s_map('G_occipital_sup_rh') = 'gos__rh';
    aparc2009s_map('G_oc-temp_lat-fusifor_rh') = 'gotlf__rh';
    aparc2009s_map('G_oc-temp_med-Lingual_rh') = 'gotml__rh';
    aparc2009s_map('G_oc-temp_med-Parahip_rh') = 'gotmp__rh';
    aparc2009s_map('G_precentral_rh') = 'gpc__rh';
    aparc2009s_map('G_pariet_inf-Angular_rh') = 'gpia__rh';
    aparc2009s_map('G_pariet_inf-Supramar_rh') = 'gpis__rh';
    aparc2009s_map('G_precuneus_rh') = 'gprcn__rh';
    aparc2009s_map('G_postcentral_rh') = 'gprct__rh';
    aparc2009s_map('G_parietal_sup_rh') = 'gps__rh';
    aparc2009s_map('G_rectus_rh') = 'gr__rh';
    aparc2009s_map('G_subcallosal_rh') = 'gs__rh';
    aparc2009s_map('G_and_S_cingul-Ant_rh') = 'gsca__rh';
    aparc2009s_map('G_and_S_cingul-Mid-Ant_rh') = 'gscma__rh';
    aparc2009s_map('G_and_S_cingul-Mid-Post_rh') = 'gscmp__rh';
    aparc2009s_map('G_and_S_frontomargin_rh') = 'gsfm__rh';
    aparc2009s_map('G_and_S_occipital_inf_rh') = 'gsoi__rh';
    aparc2009s_map('G_and_S_paracentral_rh') = 'gspc__rh';
    aparc2009s_map('G_and_S_subcentral_rh') = 'gssc__rh';
    aparc2009s_map('G_and_S_transv_frontopol_rh') = 'gstf__rh';
    aparc2009s_map('G_temporal_inf_rh') = 'gti__rh';
    aparc2009s_map('G_temporal_middle_rh') = 'gtm__rh';
    aparc2009s_map('G_temp_sup-G_T_transv_rh') = 'gtsgtt__rh';
    aparc2009s_map('G_temp_sup-Lateral_rh') = 'gtsl__rh';
    aparc2009s_map('G_temp_sup-Plan_polar_rh') = 'gtspp__rh';
    aparc2009s_map('G_temp_sup-Plan_tempo_rh') = 'gtspt__rh';
    aparc2009s_map('Lat_Fis-ant-Horizont_rh') = 'lfah__rh';
    aparc2009s_map('Lat_Fis-ant-Vertical_rh') = 'lfav__rh';
    aparc2009s_map('Lat_Fis-post_rh') = 'lfp__rh';
    aparc2009s_map('Medial_wall_rh') = 'mw__rh';
    aparc2009s_map('Pole_occipital_rh') = 'pocc__rh';
    aparc2009s_map('Pole_temporal_rh') = 'ptmp__rh';
    aparc2009s_map('S_calcarine_rh') = 'scc__rh';
    aparc2009s_map('S_central_rh') = 'sc__rh';
    aparc2009s_map('S_circular_insula_ant_rh') = 'scia__rh';
    aparc2009s_map('S_circular_insula_inf_rh') = 'scii__rh';
    aparc2009s_map('S_circular_insula_sup_rh') = 'scis__rh';
    aparc2009s_map('S_cingul-Marginalis_rh') = 'scm__rh';
    aparc2009s_map('S_collat_transv_ant_rh') = 'scta__rh';
    aparc2009s_map('S_collat_transv_post_rh') = 'sctp__rh';
    aparc2009s_map('S_front_inf_rh') = 'sfi__rh';
    aparc2009s_map('S_front_middle_rh') = 'sfm__rh';
    aparc2009s_map('S_front_sup_rh') = 'sfs__rh';
    aparc2009s_map('S_interm_prim-Jensen_rh') = 'sipj__rh';
    aparc2009s_map('S_intrapariet_and_P_trans_rh') = 'sipt__rh';
    aparc2009s_map('S_occipital_ant_rh') = 'soa__rh';
    aparc2009s_map('S_orbital-H_Shaped_rh') = 'sohs__rh';
    aparc2009s_map('S_orbital_lateral_rh') = 'sol__rh';
    aparc2009s_map('S_oc_middle_and_Lunatus_rh') = 'soml__rh';
    aparc2009s_map('S_orbital_med-olfact_rh') = 'somo__rh';
    aparc2009s_map('S_oc_sup_and_transversal_rh') = 'sost__rh';
    aparc2009s_map('S_oc-temp_lat_rh') = 'sotl__rh';
    aparc2009s_map('S_oc-temp_med_and_Lingual_rh') = 'sotml__rh';
    aparc2009s_map('S_pericallosal_rh') = 'spc__rh';
    aparc2009s_map('S_postcentral_rh') = 'spct__rh';
    aparc2009s_map('S_parieto_occipital_rh') = 'spo__rh';
    aparc2009s_map('S_precentral-inf-part_rh') = 'sprip__rh';
    aparc2009s_map('S_precentral-sup-part_rh') = 'sprsp__rh';
    aparc2009s_map('S_suborbital_rh') = 'sso__rh';
    aparc2009s_map('S_subparietal_rh') = 'ssp__rh';
    aparc2009s_map('S_temporal_inf_rh') = 'sti__rh';
    aparc2009s_map('S_temporal_sup_rh') = 'sts__rh';
    aparc2009s_map('S_temporal_transverse_rh') = 'stt__rh';
    
    % Special cases
    aparc2009s_map('Unknown') = 'unknown';
    % Whole-cortex DK metric (e.g. ...gm__dsk_mean) — not __dsk__<parcel>__...
    aparc2009s_map('all') = 'dst_mean';
    aparc2009s_map('lh') = 'lh';
    aparc2009s_map('rh') = 'rh';
    
    map_use = aparc2009s_map;
    map_tab2atlas = containers.Map(values(aparc2009s_map), keys(aparc2009s_map));
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
                % Unknown name - return empty
                warning('Unknown aparc.a2009s ROI name for direction %s: %s', direction, name);
                roinames_out{i} = '';
            end
        end
        atlas_name = local_atlasNameForCode(name, roinames_out{i}, direction, aparc2009s_map, map_tab2atlas);
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
    lut_path = fullfile(repo_root, 'support_files', 'aparc_a2009s_lut.csv');
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


