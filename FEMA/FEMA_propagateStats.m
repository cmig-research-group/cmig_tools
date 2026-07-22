function vol_stat = FEMA_propagateStats(roi_stat, parc_name, roinames, varargin) 
% FEMA_propagateStats - Propagate ROI-based statistics to voxel/vertex space
%
% Syntax:  vol_stat = FEMA_propagateStats(roi_stat, roi_atlas, parc_name)
%
% Required inputs:
%    roi_stat  - matrix containing ROI-based statistics
%    parc_name - tabulated data parcellation name, currently supporting 
%               'aseg' (Freesurfer subcortical segmentation) or 'at' (AtlasTrack) for voxelxwise data
%               'aparc' /'dsk' (Desikan-Killiany) or 'aparc.a2009s'/ 'dst' (Destrieux) for vertexwise data
%   ymat_names - cell array of roi names for each column in roi_stat
% Optional inputs:
%    'outdim' - output dimensions, if roi_atlas is volumetrix, whether to output vol_stat as 3D or 2D matrix
%               Default is '2D'
% Outputs:
%    vol_stat - matrix containing voxel/vertex-based statistics
%
    p = inputParser;
    addRequired(p, 'roi_stat', @isnumeric);
    addRequired(p, 'parc_name',    @(x) ischar(x) && ...
                                    ismember(x, {'aseg' 'at' 'dsk' 'dst' 'aparc' 'aparc_a2009s'}));
    addRequired(p, 'roinames', @(x) iscell(x) && all(cellfun(@ischar, x)));
    addParameter(p, 'outdim', '2D', @(x) ischar(x) && ismember(x, {'2D' '3D'}));

    parse(p, roi_stat, parc_name, roinames, varargin{:});
    roi_stat = p.Results.roi_stat;
    parc_name = p.Results.parc_name;
    roinames = p.Results.roinames;
    outdim   = p.Results.outdim;

    switch parc_name
        case {'aseg', 'at'}
            fname_roi2atlas = 'roi2ABCDAtlasMaps.mat';
        case {'dsk', 'dst', 'aparc', 'aparc_a2009s'}
            fname_roi2atlas = 'roi2SurfaceAtlasMaps.mat';
    end 
    % check if roi2atlas file exists
    if ~exist(fname_roi2atlas, 'file') 
        error('Mapping file %s not found. Check the path.', fname_roi2atlas);
    end 
    % load roi2atlas mapping
    load(fname_roi2atlas)

    % map tabulated data names to parcellation names
    parc_atlas = {'aseg' 'fiber' ...
                  'aparc' 'aparc_a2009s'};
    parc_tab = {'aseg' 'at' ...
                'dsk' 'dst'};
    parc_name_map = containers.Map(parc_tab, ...
                                   parc_atlas);
    % map tabulated parcellation name to atlas parcellation name
    parc_name_atlas = parc_name_map(parc_name);
    
    % check that parc_name exists in roi2atlas
    if ~isfield(roi2atlas, parc_name_atlas)
        error('Parcellation %s not found in mapping file %s.', parc_name_atlas, fname_roi2atlas);
    end
    roimat = roi2atlas.(parc_name_atlas).roimat;
    roinames_atlas = roi2atlas.(parc_name_atlas).roinames;
    roicodes_atlas = roi2atlas.(parc_name_atlas).roicodes;

    % extract tabulated roinames by taking everthing after parc_name and '__'
    roinames_tab = extractAfter(roinames, [parc_name '__']);
    % anything after __lh or __rh 
    %roinames_tab = regexprep(roinames_tab, '(lh|rh).*$', '$1');
    roinames_tab = regexprep(roinames_tab, '_wmean$|_mean$', '');
    % find empty roinames_tab and replace with 'all' 
    empty_idx = cellfun(@isempty, roinames_tab);
    roinames_tab(empty_idx) = {'all'};

    % atlas roinames 
    switch parc_name 
        case 'at' 
            roinames_atlas_out = fiber2tabroinames(roinames_atlas);
        case 'aseg'
            roinames_atlas_out = aseg2tabroinames(roinames_atlas);
        case {'dsk', 'aparc'}
            roinames_atlas_out = aparc2tabroinames(roinames_atlas);
        case {'dst', 'aparc_a2009s'}
            roinames_atlas_out = aparc2009s2tabroinames(roinames_atlas);
    end 

    % intersect roinames_tab and roinames_atlas_out to find common rois
    [roinames_common, ia, ib] = intersect(roinames_tab, roinames_atlas_out, 'stable');
    % check 
    if ~isequal(roinames_common, roinames_tab')
        error(('Not all roinames in roi_stat are found in the atlas parcellation %s.'), parc_name_atlas);
    end

    % map values in roi_stat to vol_stat using roimat
    vol_stat = roi_stat(:, ia) * roimat(ib, :);
end 



    
    

    
    
    
    

			



        
        

