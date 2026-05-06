function roi2atlas = roi2ABCDAtlas(fname_atlas, fname_out, varargin)
    % this only works for the ABCD3 atlas 
    % create a logical matrix to map roi data to volumetric ABCD atlas
    % 
    % parc_in can be a string or cell array of strings
    % outputs are saved as a structure with fields for each parcellation

    p = inputParser;
    addParameter(p, 'parc_in', {'aseg', 'aparcaseg', 'fiber'}, ...
                               @(x) ischar(x) || isstring(x) || iscell(x) && ...
                                ismember(x, {'aseg', 'aparcaseg', 'fiber'}));
    addParameter(p, 'save_flag', true, @(x) islogical(x));

    parse(p, varargin{:});
    parc_in = p.Results.parc_in;
    save_flag = p.Results.save_flag;

    % Convert single parcellation to cell array for uniform processing
    if ischar(parc_in) || isstring(parc_in)
        parc_in = {char(parc_in)};
    elseif iscell(parc_in)
        parc_in = cellfun(@char, parc_in, 'UniformOutput', false);
    end

    % check if fname_out is empty if so set save_flag to false
    if isempty(fname_out)
        save_flag = false;
    end

    % load atlas
    load(fname_atlas);

    % Process each parcellation
    roi2atlas = struct();
    vol_mask_aseg = atlas_dspace.vol_mask_aseg;
    vol_mask_sub = subsample_volume(vol_mask_aseg);
    ivec_mask = find(vol_mask_sub>=0.5);

    for pi = 1:length(parc_in)
        parc_name = parc_in{pi};
        parc_key = ['muvols_' parc_name];
        
        % Create logical matrix of roi voxels
        parcellation = subsample_volume(atlas_dspace.(parc_key));
        nroi = size(parcellation, 4);
        roimat = false(nroi, length(ivec_mask));
        
        for r = 1:nroi
            vol_roi = parcellation(:,:,:,r);
            vol_roi = vol_roi > 0.5;
            roimat(r,:) = vol_roi(ivec_mask);
        end
        
        % get roicodes, roinames, roirgb
        roicodes = atlas_dspace.(['indlist_' parc_name]);
        
        switch parc_key
            case {'muvols_aseg' 'muvols_aparcaseg'}
                [roicodes_lut, roinames_lut, rgbv_lut] = fs_colorlut();
                lut_rows = find(ismember(roicodes_lut, roicodes));
                roinames = roinames_lut(lut_rows);
                roirgb = rgbv_lut(lut_rows, 1:3) ./ 255;
            case 'muvols_fiber'
                fibers_csv = '/usr/pubsw/packages/MMPS/MMPS_DEV/documentation/AtlasTrack/DTI_Fiber_Legend.csv';
                fiber_legend = readtable(fibers_csv);
                idx_in = find(ismember(fiber_legend.FiberNumber, roicodes));
                roinames = fiber_legend.FiberName(idx_in);
                roirgb = zeros(length(idx_in), 3);
                for i = 1:length(idx_in)
                    roirgb(i,:) = str2num(fiber_legend.ColorRGB{idx_in(i)});
                end
                roirgb = roirgb / 255;
        end
        
        % Store in output structure with parcellation name as field
        roi2atlas.(parc_name).roimat = roimat;
        roi2atlas.(parc_name).roicodes = roicodes;
        roi2atlas.(parc_name).roinames = roinames;
        roi2atlas.(parc_name).roirgb = roirgb;
        roi2atlas.mask = vol_mask_sub;

        % make sure that all roinames, codes and rgb have the same length
        roi2atlas.(parc_name).roicodes = reshape(roi2atlas.(parc_name).roicodes, length(roi2atlas.(parc_name).roicodes), 1);
        roi2atlas.(parc_name).roinames = reshape(roi2atlas.(parc_name).roinames, length(roi2atlas.(parc_name).roinames), 1);
        roi2atlas.(parc_name).roirgb = reshape(roi2atlas.(parc_name).roirgb, length(roi2atlas.(parc_name).roirgb), 3);
    end