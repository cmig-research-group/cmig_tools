function fig = FEMA_imageMosaic(volData, cbarTxt, maskBinary, orientation, slices, ...
                                cmap,    climits, cLocation,  bgColor,     interp)
% Function that presents voxel-wise FEMA-derived statistics - fixed or 
% random effcts as a mosaic - sagittal, coronal, transverse, or multi-view 
% with all three orientations
%% Inputs:
% volData:      [R x C x P] 3D volume data containing voxel-wise statistics
%
% cbarTxt:      text to be displayed with the colorbar
%
% The following input arguments are optional:
%
% maskBinary:   [R x C x P] 3D binary volume used for masking
%
% orientation:  should be one of the following (see Notes):
%                   * 'sagittal'    [mosaic of sagittal orientation]
%                   * 'coronal'     [mosaic of coronal orientation]
%                   * 'transverse'  [mosaic of transverse orientation]
%                   * 'multi'       [mosaic of all three orientation]
% 
% slices:       vector of numbers indicating which slices to show 
%               (see Notes)
% 
% cmap:         name of the colormap to be used
%
% climits:      vector having the lower and upper limits for the colorbar
%
% cLocation:    location of colorbar; should be one of:
%                   * 'right' : colorbar on the right of the figure
%                   * 'bottom': colorbar at the bottom of the figure
%
% bgColor:      an RGB vector specifying the background color of the
%               figure; the values should be between 0 and 1 
%
% interp:       interpolation mode for imagesc which is used to generate
%               each slice; should be either of the following:
%                   * 'nearest'
%                   * 'bilinear'
%
%% Output(s):
% fig:          handle to the figure containing all the slices; the figure
%               handle can then be used for further customization or saving
%               as an image
%
%% Notes:
% The 'orientation' variable controls which of the three views are shown as
% a mosaic; if 'orientation' is specified as 'multi', the first row of the
% mosaic contains sagittal view, the second row contains coronal view, and
% the third row contains transverse view
%
% When the orientation is set to 'multi' and the user wants to specify the
% slices to be shown, then, 'slices' should be [3 x n] matrix where each 
% row contains the slices to be selected for each of the three orientations 
% (sagittal, coronal, and transverse - in this order); equal number of 
% slices should be specified in each orientation. If, instead of [3 x n],
% the user specifies [1 x n] slices, then the same slices are used for each
% of the orientations
%
% This script makes several assumptions - first, it assumes that input
% image and the generated/specified mask are all in alignment with each 
% other and have the same NIfTI header. Additionally, this script is 
% written for ABCD images and assumes that the first dimension of the 
% volume indexes the transverse slices, the second dimension of the volume 
% indexes the sagittal slices, and the third dimension of the volume 
% indexes the coronal slices. Further, it assumes that coronal slices need
% to be rotated clockwise by 90 degrees before visualization
%
% If the background mask and the volData size do not match, the maskBinary
% is either upsampled or downsampled - volData is untouched
%
% Note that the interpolation flag controls the internal interpolation of
% the image when slices are created using imagesc function; if the user
% selects 'bilinear', there might be a visible lag in drawing, depending on
% the number of slices, etc.
%
% In some situations, the resulting figure can have a large amount of empty
% space between each row of images (for example, when using multi-view for
% a large number of slices); this can be edited, as needed, using the
% returned figure handle
%
%% Examples:
% Display Z statistics for the 10th fixed effect
% FEMA_imageMosaic(vol_z(:, :, :, 10), '10th FFX');
% 
% Display 2nd random effect
% FEMA_imageMosaic(vol_sig2(:, :, :, 2), '2nd RFX');
%
% Display the 1st random effect, resclaed to original variance values
% FEMA_imageMosaic(vol_sig2(:, :, :, 1) .* vol_sig2t, '1st RFX');
%
% Save figure as a 900 DPI png image
% fig = FEMA_imageMosaic(vol_sig2(:, :, :, 3) .* vol_sig2t, '3rd RFX');
% print(fig, fullfile(outputDir, 'RFX_3.png'), '-dpng', '-r900');
% close(fig);
%
% Override the colorbar to parula
% FEMA_imageMosaic(vol_z(:, :, :, 5), '5th FFX', [], [], [], 'parula');
%
% Specify custom slices in coronal orientation 
% FEMA_imageMosaic(vol_sig2(:, :, :, 2), '2nd RFX', [], 'coronal', [30, 40, 50]);
%
% Specify custom slices in coronal orientation 
% FEMA_imageMosaic(vol_sig2(:, :, :, 2), '2nd RFX', [], 'multi', 45:5:60);
% 
%% Defaults:
% orientation:      'multi';
% maskBinary:       a whole brain mask based on ABCD atlas is created
% slices:           five linearly spaced slices between 25 and 75
% cmap:             'blueblackred' for data containing positive and
%                   negative values; 'fire' for data containing only 
%                   positive values
% climits:          [min(volData), max(volData)]
% cLocation:        'bottom'
% bgColor:          [0 0 0] i.e., black background
% interp:           'nearest'

%% Parse inputs and assign defaults
% Check volData
if ~exist('volData', 'var') || isempty(volData)
    error('Please specify the volumetric data to plot');
end

% Check orientation
if ~exist('orientation', 'var') || isempty(orientation)
    orientation = 'multi';
else
    if ~ismember(orientation, {'sagittal', 'coronal', 'transverse', 'multi'})
        error(['Unknown orientation specified: ', orientation, '; should be one of: sagittal, coronal, transverse, or multi']);
    end
end

% Check maskBinary
if ~exist('maskBinary', 'var') || isempty(maskBinary)
    disp('No mask specified; generating a whole brain mask based on the ABCD aseg atlas');

    % Load atlas information
    cfg = abcdConfig('showVol');
    aseg = load(fullfile(cfg.data.showVolData, 'Atlas', 'showVolAtlases_ABCD2_cor10.mat'), 'aseg');
    aseg = aseg.aseg;

    % Convert the probability values to full image
    asegProb                = zeros(aseg.prob{1},'single');
    asegProb(aseg.prob{2}) = single(aseg.prob{3});

    % Binarize
    asegBinarized = asegProb;
    for vols = 1:size(asegProb, 4)
        % Choosing an extremely liberal probability value for binarizing
        asegBinarized(:,:,:,vols) = single(asegBinarized(:,:,:,vols) > 0.2);
        % asegBinarized(:,:,:,vols) = single(asegBinarized(:,:,:,vols) > aseg.probabilityThreshold);
    end

    % Sum across regions to generate a full brain mask
    maskBinary = logical(sum(asegBinarized, 4));

    % Remove unwanted variables
    clear aseg*
end

% Check slices
if ~exist('slices', 'var') || isempty(slices)
    if size(volData,1) == 200
        % The data is 1mm isotropic (hopefully)
        slices = ceil(linspace(50, 150, 5));
    else
        % Either 2mm isotropic or some other data
        slices = ceil(linspace(30, 75, 5));
    end
    
    % Convert to matrix, if orientation is multi
    if strcmpi(orientation, 'multi')
        slices = repmat(slices, 3, 1);
    end
else
    tmpSZ = size(slices);
    if strcmpi(orientation, 'multi')
        % User has specified 1 x n slices
        if tmpSZ(1) == 1 || tmpSZ(2) == 1
            slices = reshape(slices, 1, length(slices));
            slices = repmat(slices, 3, 1);
        else
            % User has specified 3 x n slices
            if tmpSZ(1) == 3 || tmpSZ(2) == 3
                slices = reshape(slices, 3, length(slices));
            else
                error('Either specify slices as 1 x n vector or 3 x n matrix');
            end
        end
    else
        % Ensure that user has specified 1 x n slices
        if tmpSZ(1) == 1 || tmpSZ(2) == 1
            slices = reshape(slices, 1, length(slices));
        else
            error('Please specify slices as 1 x n vector');
        end
    end
end

% Check cmap
if ~exist('cmap', 'var') || isempty(cmap)
    minVal = min(volData, [], 'all');
    if minVal < 0
        cmap = 'blueblackred';
    else
        cmap = 'fire';
    end
end

% Check climits
if ~exist('climits', 'var') || isempty(climits)
    climits = [min(volData(:)), max(volData(:))];
else
    if length(climits) ~= 2
        error('Please specify a minimum and a maximum value for color bar limits');
    end
end

% Check cLocation
if ~exist('cLocation', 'var') || isempty(cLocation)
    cLocation = 'south';
else
    if strcmpi(cLocation, 'bottom')
        cLocation = 'south';
    else
        if strcmpi(cLocation, 'right')
            cLocation = 'eastoutside';
        else
            error(['Unknown cLocation specified: ', cLocation, '; should be either bottom or right']);
        end
    end
end

% Check bgColor
if ~exist('bgColor', 'var') || isempty(bgColor)
    bgColor = [0 0 0];
end

% Check interp
if ~exist('interp', 'var')|| isempty(interp)
    interp = 'nearest';
else
    if ~ismember(interp, {'nearest', 'bilinear'})
        error(['Unknown interp method specified: ', interp, '; should be either nearest or bilinear']);
    end
end

%% Determine if the brain mask needs to be resampled
ratio = unique(size(maskBinary)./size(volData));
if length(ratio) > 1
    warning(['The size of the data is: ', num2str(size(volData)), ' and the size of the generated mask is: ', num2str(size(maskBinary)), '; skipping masking of data']);
    maskBinary = false(size(volData));
else
    if ratio < 1
        warning(['The size of the data is: ', num2str(size(volData)), ' and the size of the generated mask is: ', num2str(size(maskBinary)), '; upsampling the mask']);
        maskBinary = logical(upsample_volume_nearest(maskBinary, 1/ratio, 1/ratio, 1/ratio));
    else
        if ratio > 1
            warning(['The size of the data is: ', num2str(size(volData)), ' and the size of the generated mask is: ', num2str(size(maskBinary)), '; downsampling the mask']);
            maskBinary = logical(subsample_volume(maskBinary, ratio, ratio, ratio));
        end
    end
end

%% Determine number of rows and columns in subplot
if strcmpi(orientation, 'multi')
    nrows = 3;
    ncols = size(slices, 2);
else
    % Need to find an optimum number of rows and columns
    [nrows, ncols] = calc_rows_cols_subplot(length(slices));
    
    % Prefer smaller number of rows
    if nrows > ncols
        tmp  = ncols;
        ncols = nrows;
        nrows = tmp;
    end
end

%% Determine the first and last slices which have data
% Add one more slice to the found values
switch orientation
    case 'sagittal'
        [tmp_upp_sag, tmp_left_sag] = find(squeeze(sum(logical(squeeze(maskBinary(:, slices,:))), 2)));
        left_sag                    = max(min(tmp_left_sag) - 1, 1);
        right_sag                   = min(max(tmp_left_sag) + 1, size(maskBinary, 2));
        low_sag                     = max(min(tmp_upp_sag)  - 1, 1);
        upp_sag                     = min(max(tmp_upp_sag)  + 1, size(maskBinary, 2));
        
    case 'coronal'
        [tmp_upp_cor, tmp_left_cor] = find(squeeze(sum(logical(squeeze(maskBinary(:, :, slices))), 3)));
        left_cor                    = max(min(tmp_left_cor) - 1, 1);
        right_cor                   = min(max(tmp_left_cor) + 1, size(maskBinary, 3));
        low_cor                     = max(min(tmp_upp_cor)  - 1, 1);
        upp_cor                     = min(max(tmp_upp_cor)  + 1, size(maskBinary, 3));

    case 'transverse'
        temp                        = rot90(squeeze(sum(logical(squeeze(maskBinary(slices,:,:))), 1)));
        [tmp_upp_tra, tmp_left_tra] = find(temp);
        left_tra                    = max(min(tmp_left_tra) - 1, 1);
        right_tra                   = min(max(tmp_left_tra) + 1, size(temp, 2));
        low_tra                     = max(min(tmp_upp_tra)  - 1, 1);
        upp_tra                     = min(max(tmp_upp_tra)  + 1, size(temp, 1));
        
    case 'multi'
        [tmp_upp_sag, tmp_left_sag] = find(squeeze(sum(logical(squeeze(maskBinary(:, slices,:))), 2)));
        left_sag                    = max(min(tmp_left_sag) - 1, 1);
        right_sag                   = min(max(tmp_left_sag) + 1, size(maskBinary, 2));
        low_sag                     = max(min(tmp_upp_sag)  - 1, 1);
        upp_sag                     = min(max(tmp_upp_sag)  + 1, size(maskBinary, 2));
        
        [tmp_upp_cor, tmp_left_cor] = find(squeeze(sum(logical(squeeze(maskBinary(:, :, slices))), 3)));
        left_cor                    = max(min(tmp_left_cor) - 1, 1);
        right_cor                   = min(max(tmp_left_cor) + 1, size(maskBinary, 3));
        low_cor                     = max(min(tmp_upp_cor)  - 1, 1);
        upp_cor                     = min(max(tmp_upp_cor)  + 1, size(maskBinary, 3));

        temp                        = rot90(squeeze(sum(logical(squeeze(maskBinary(slices,:,:))), 1)));
        [tmp_upp_tra, tmp_left_tra] = find(temp);
        left_tra                    = max(min(tmp_left_tra) - 1, 1);
        right_tra                   = min(max(tmp_left_tra) + 1, size(temp, 2));
        low_tra                     = max(min(tmp_upp_tra)  - 1, 1);
        upp_tra                     = min(max(tmp_upp_tra)  + 1, size(temp, 1));
end

%% Initialize figure
fig         = figure('Units', 'centimeters', 'Position', [10 10 16 12], 'Color', bgColor);
fgap        = 0.02;
leftMargin  = 0.01;
topMargin   = 0.01;
if strcmpi(cLocation, 'south')
    rightMargin  = 0.01;
    bottomMargin = 0.15;
else
    rightMargin  = 0.15;
    bottomMargin = 0.01;
end
allH = tight_subplot(nrows, ncols, [fgap fgap], [bottomMargin, topMargin], [leftMargin, rightMargin]);
count = 1;

% Actual plotting
switch orientation
    
    case 'sagittal'
        for slice = 1:length(slices)
            doSagittal(allH(count), squeeze(volData(:, slices(slice), :)), squeeze(maskBinary(:, slices(slice), :)), cmap, climits, [left_sag, right_sag], [low_sag, upp_sag], interp);
            count = count + 1;
        end
        
    case 'coronal'
        for slice = 1:length(slices)
            doCoronal(allH(count), squeeze(volData(:, :, slices(slice))), squeeze(maskBinary(:, :, slices(slice))), cmap, climits, [left_cor, right_cor], [low_cor, upp_cor], interp);
            count = count + 1;
        end
        
    case 'transverse'
        for slice = 1:length(slices)
            doTransverse(allH(count), rot90(squeeze(volData(slices(slice), :, :)), 1), rot90(squeeze(maskBinary(slices(slice), :, :))), cmap, climits, [left_tra, right_tra], [low_tra, upp_tra], interp);
            count = count + 1;
        end
        
    case 'multi'
        for slice = 1:size(slices, 2)
            doSagittal(allH(count), squeeze(volData(:, slices(1, slice), :)), squeeze(maskBinary(:, slices(1, slice), :)), cmap, climits, [left_sag, right_sag], [low_sag, upp_sag], interp);
            count = count + 1;
        end
        
        for slice = 1:size(slices, 2)
            doCoronal(allH(count), squeeze(volData(:, :, slices(2, slice))), squeeze(maskBinary(:, :, slices(2, slice))), cmap, climits, [left_cor, right_cor], [low_cor, upp_cor], interp);
            count = count + 1;
        end
        
        for slice = 1:size(slices, 2)
            doTransverse(allH(count), rot90(squeeze(volData(slices(3, slice), :, :)), 1), rot90(squeeze(maskBinary(slices(3, slice), :, :))), cmap, climits, [left_tra, right_tra], [low_tra, upp_tra], interp);
            count = count + 1;
        end
end

%% Delete any extra axes that were created
if (nrows * ncols) > (count - 1)
    delete(allH(count:end));
end

%% Add a colour bar 
if strcmpi(orientation, 'multi')
    
    % Axis number to the lower left plot where colorbar will be added
    hw = ncols * (nrows - 1) + 1;
    
    % Get some positions before the colorbar is shown
    pos_H1S = plotboxpos(allH(1));
    pos_H1C = plotboxpos(allH(ncols   + 1));
    pos_H1T = plotboxpos(allH(ncols*2 + 1));
    pos_HlS = plotboxpos(allH(ncols));
    pos_Hw  = plotboxpos(allH(hw));
    
    % Add colorbar
    cbar = colorbar(allH(hw), 'Location', cLocation);
    
    % Reposition colorbar
    if strcmpi(cLocation, 'south')
        
        % Estimate the horizontal gap between the coronal slices
        sGap = abs(allH(ncols + 1).Position(3) - pos_H1C(3));
        
        % Align colorbar
        cbar.Position(1) = pos_H1C(1);
        cbar.Position(2) = pos_Hw(2) - 0.145;
        cbar.Position(3) = pos_H1C(3) * ncols + sGap * (ncols - 1) + fgap * (ncols - 1);
        cbar.Position(4) = 0.04;
        
    else  
        % Estimate the vertical gap between slices
        sGap_C = abs(allH(ncols + 2).Position(4) - pos_H1C(4));
        sGap_T = abs(allH(2).Position(4) - pos_H1T(4));
        
        % Align colorbar
        cbar.Position(1) = pos_HlS(1) + pos_HlS(3) + 0.02;
        cbar.Position(2) = pos_H1T(2);
        cbar.Position(3) = 0.04;
        cbar.Position(4) = min(pos_H1S(4) + pos_H1C(4) + pos_H1T(4) + sGap_C + sGap_T + fgap * (nrows), 1 - fgap);
    end
else
    % Axis number to the lower left plot where colorbar will be added
    % HHandle the case of single row
    if nrows == 1
        hw          = ncols;
        pos_Hright  = plotboxpos(allH(end));
    else
        hw = ncols * (nrows - 1) + 1;
        pos_Hright = plotboxpos(allH(ncols * (nrows - 1)));
    end
    
    % Get some positions before displaying the colorbar
    pos_Hleft  = plotboxpos(allH(1));
    pos_Hw     = plotboxpos(allH(hw));
    
    % Add colorbar
    cbar = colorbar(allH(hw), 'Location', cLocation);
    
    if strcmpi(cLocation, 'south')

        % Estimate the horizontal gap between the slices
        sGap = abs(allH(1).Position(3) - pos_Hw(3));

        % Align colorbar
        cbar.Position(1) = pos_Hleft(1);
        cbar.Position(2) = pos_Hw(2) - 0.145;
        cbar.Position(3) = pos_Hleft(3) * ncols + sGap * (ncols - 1) + fgap * (ncols - 1);
        cbar.Position(4) = 0.04;

    else
        
        % Estimate the vertical gap between the slices
        if nrows == 1
            sGap_V = abs(allH(hw).Position(4) - pos_Hw(4));
        else
            sGap_V = abs(allH(hw + 1).Position(4) - pos_Hw(4));
        end
        
        % Align colorbar
        cbar.Position(1) = pos_Hright(1) + pos_Hright(3) + 0.02;
        cbar.Position(2) = pos_Hw(2);
        cbar.Position(3) = 0.04;
        cbar.Position(4) = pos_Hright(4) * nrows + sGap_V * (nrows - 1) + fgap * (nrows);
    end
end

% Customize colorbar
cbar.FontSize         = 10;
cbar.Label.String     = cbarTxt;
cbar.Label.FontSize   = 12;
cbar.Label.FontWeight = 'bold';

% In case of black background, make the font white
if sum(bgColor == [0 0 0]) == 3
    cbar.Color = [1 1 1];
end
end

function doSagittal(h, slicedData, slicedMask, cmap, climits, lrLims, udLims, interp)
if isempty(find(slicedData, 1))
    delete(h);
else
    im = imagesc(h, slicedData, 'AlphaData', slicedMask, climits);
    im.Interpolation = interp;
    colormap(h, cmap);
    axis(h, 'equal');
    box(h, 'off');
    xlim(h, lrLims);
    ylim(h, udLims);
    xticks(h, []);
    yticks(h, []);
    axis(h, 'off');
end
end

function doCoronal(h, slicedData, slicedMask, cmap, climits, lrLims, udLims, interp)
if isempty(find(slicedData, 1))
    delete(h);
else
    im = imagesc(h, slicedData, 'AlphaData', slicedMask, climits);
    im.Interpolation = interp;
    colormap(h, cmap);
    axis(h, 'equal');
    box(h, 'off');
    xlim(h, lrLims);
    ylim(h, udLims);
    xticks(h, []);
    yticks(h, []);
    axis(h, 'off');
end
end

function doTransverse(h, slicedData, slicedMask, cmap, climits, lrLims, udLims, interp)
if isempty(find(slicedData, 1))
    delete(h);
else
    im = imagesc(h, slicedData, 'AlphaData', slicedMask, climits);
    im.Interpolation = interp;
    colormap(h, cmap);
    axis(h, 'equal');
    box(h, 'off');
    xlim(h, lrLims);
    ylim(h, udLims);
    xticks(h, []);
    yticks(h, []);
    axis(h', 'off');
end
end

%% The following code section is work in progress for displaying an underlay image
% underlay:     [R x C x P] 3D volume data to be used as an underlay, or
%               either of the following:
%                   * 'abcd'
%                   * 'none'

% Check underlay
% if ~exist('underlay', 'var') || isempty(underlay)
%     underlay = [];
% else
%     if ischar(underlay)
%         if strcmpi(underlay, 'abcd')
%             underlay = atlas_T1('ABCD2_cor10', true);
%             underlay = underlay.imgs;
%         else
%             if strcmpi(underlay, 'none')
%                 underlay = [];
%             else
%                 error(['Unknown underlay value specified: ', underlay]);
%             end
%         end
%     end
% end

% %% Determine if underlay needs to be resampled
% if ~isempty(underlay)
%     ratio = unique(size(underlay)./size(volData));
%     if length(ratio) > 1
%         warning(['The size of the data is: ', num2str(size(volData)), ' and the size of the underlay is: ', num2str(size(underlay)), '; skipping displaying underlay']);
%         underlay = [];
%     else
%         if ratio < 1
%             warning(['The size of the data is: ', num2str(size(volData)), ' and the size of the underlay is: ', num2str(size(underlay)), '; upsampling the underlay']);
%             underlay = upsample_volume(underlay, 1/ratio, 1/ratio, 1/ratio);
%         else
%             if ratio > 1
%                 warning(['The size of the data is: ', num2str(size(volData)), ' and the size of the underlay is: ', num2str(size(underlay)), '; downsampling the underlay']);
%                 underlay = subsample_volume(underlay, ratio, ratio, ratio);
%             end
%         end
%     end
% end