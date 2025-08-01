function save_showVol_images(coords, outdir, fstem, varargin)
% save_showVol_images - saves screenshots from showVol with optional multi-parcellation ROI overlays.
%
% Usage:
%   save_showVol_images(coords, outdir, fstem, 'CoordSpace', 'RCS', ...
%       'drawROI', true, 'parcellation', {'aseg','fiber'}, ...
%       'roiNames', struct('aseg', {'Left-Caudate'}, 'fiber', {'all_fiber'}), ...
%       'lw', 2, 'fixedColor', [0 0.45 0.75])
%
% REQUIRED INPUTS:
%   coords      - Nx3 matrix of coordinates (LPH by default, or RCS if specified)
%   outdir      - directory to save images
%   fstem       - filename stem for output PNGs (string OR cell array)
%
% OPTIONAL INPUTS (varargin):
%   CoordSpace   - 'LPH' (default) or 'RCS'
%       Determines coordinate formatting and how coords are interpreted.
%       • RCS coords use 3-digit zero-padded numbers (e.g., X096_Y100_Z120).
%       • LPH coords use 2-digit zero-padded numbers with 'N' for negatives (e.g., X03_Y00_ZN17).
%
%   drawROI      - (logical, default: false)
%       Whether to draw ROI overlays for selected parcellations.
%
%   parcellation - (cell array)
%       List of parcellations to use (e.g., {'aseg','fiber'}).
%
%   roiNames     - (struct)
%       Fields are parcellation names, values are ROI name lists. Special keywords
%       like 'subcortical_all' or 'all_fiber' will expand to all ROIs for that parcellation.
%
%   lw           - (numeric, default: 1)
%       Line width for ROI outline drawing. Affects all ROI overlays drawn in the screenshot.
%
%   fixedColor   - (1x3 numeric RGB vector, default: [1 1 1])
%       Fixed color for ROI outlines. Default is white. Specify as [R G B],
%       e.g. [0 0.45 0.75] for MATLAB’s default blue.
%
% UI HIDING OPTIONS (all default to TRUE):
%   hideUI, hideCrosshairs, hideColorbar, hideROIText, hideButtons, hideVolName
%       These control hiding of GUI elements (crosshairs, colorbar, buttons, etc.)
%       in the saved screenshots.
%
% SHORTCUT ROI KEYWORDS (case-insensitive):
%   'subcortical_all', 'all_fiber', 'all_aparcaseg', 'all_thalamus', 'all_pauli'
%
% OUTPUT:
%   Saves PNG screenshots into a subfolder named after the orientation (coronal/sagittal/axial).
%   Filenames are built as:
%       coord[_fstem]_vol###
%   Example: X096_Y100_Z120_bmi_vol001.png

%% --------------------
% Parse inputs
%% --------------------
p = inputParser;

% Core functionality
p.addParameter('drawROI', false, @(x) islogical(x) || isnumeric(x));
p.addParameter('parcellation', {}, @(x) iscell(x) || isstring(x));
p.addParameter('roiNames', struct(), @isstruct);

% Coordinate space option
p.addParameter('CoordSpace', 'RCS', @(x) any(strcmpi(x, {'LPH','RCS'})));

% ROI drawing style
p.addParameter('lw', 1, @(x) isnumeric(x) && isscalar(x));  % default line width
p.addParameter('fixedColor', [1 1 1], @(x) isnumeric(x) && numel(x)==3);  % default = white [R G B]

% UI Hiding controls (all default to TRUE)
p.addParameter('hideUI', true, @(x) islogical(x) || isnumeric(x));          
p.addParameter('hideCrosshairs', true, @(x) islogical(x) || isnumeric(x));  
p.addParameter('hideColorbar', true, @(x) islogical(x) || isnumeric(x));    
p.addParameter('hideROIText', true, @(x) islogical(x) || isnumeric(x));     
p.addParameter('hideButtons', true, @(x) islogical(x) || isnumeric(x));     
p.addParameter('hideVolName', true, @(x) islogical(x) || isnumeric(x));     

p.parse(varargin{:});
opts = p.Results;

%% --------------------
% Setup
%% --------------------
scaleCoord = @(coord, scale) round(coord * scale);

if isempty(coords)
    coords = [99 99 129];
end
if ~exist(outdir, 'dir')
    mkdir(outdir);
end

% Normalize fstem into a cell array
if ischar(fstem) || isstring(fstem)
    fstem = {char(fstem)};
end

% Get the showVol window
showVolWins = findall(groot, 'Type', 'figure', 'Tag', 'showVol');
if isempty(showVolWins)
    error('No showVol window found.');
end
myWin = showVolWins(1);
handles = guidata(myWin);

nVols = numel(handles.vols);

% Expand fstem if needed
if numel(fstem) == 1
    fstem = repmat(fstem, 1, nVols);
elseif numel(fstem) ~= nVols
    error('fstem must be a single string or match number of vols (%d).', nVols);
end

% Orientation state
orientation = handles.ORIENTATION;
set(0, 'CurrentFigure', myWin);
axes(handles.axes1);
hold(handles.axes1, 'on');

%% --------------------
% Hide UI elements (one-time global hiding)
%% --------------------
if opts.hideUI
    if opts.hideCrosshairs
        set(handles.toggleLines1, 'ForegroundColor', 'k');
        set(handles.axes1_X, 'Visible', 'off');
    end
    if opts.hideColorbar && isfield(handles, 'colorbar') && isgraphics(handles.colorbar)
        set(handles.colorbar, 'Visible', 'off');
        set(allchild(handles.colorbar), 'Visible', 'off');
    end
    if opts.hideButtons
        delete(findall(myWin,'Style','pushbutton','String','O'))
        delete(findall(myWin,'Style','togglebutton','Tag','toggleLines1'))
        delete(findall(myWin,'Style','togglebutton','Tag','toggleLines2'))
        delete(findall(myWin,'Style','togglebutton','Tag','toggleLines3'))
        if isfield(handles,'overlay_opacity_slider')
            delete(handles.overlay_opacity_slider);
        end
    end
end

updateDisplay_newvol(handles, true);

%% --------------------
% Loop through volumes & coords
%% --------------------
for v = 1:nVols
    fprintf('Processing vol %d/%d…\n', v, nVols);

    % Set current volume state
    set(0, 'CurrentFigure', myWin);
    tmp = find(handles.hideVol == 0);
    handles.currentVol = tmp(v);
    handles.currentVolPB = handles.volPB(handles.currentVol);
    set(handles.currentVolPB, 'ForegroundColor', 'r');
    handles.Mvxl2lph = handles.vols{handles.currentVol}.Mvxl2lph;

    handles.ORIENTATION = orientation;
    handles = displayNewOrientation(handles);  
    guidata(myWin, handles);                    
    updateDisplay_newvol(handles, true);   

    V = handles.vols{handles.currentVol};

    % Track duplicate fstem count for numbering
    dupCount = sum(strcmp(fstem(1:v), fstem{v}));

    %% --------------------
    % Coord loop
    %% --------------------
    for ci = 1:size(coords,1)
        coord = coords(ci,:);    

        ud = struct();
        if strcmpi(opts.CoordSpace, 'LPH')
            ud.lph = coord;
        else
            Mvxl2lph = handles.Mvxl2lph_atlas;
            lph_coord = Mvxl2lph * [coord 1]';  % convert to LPH for showVol
            ud.lph = lph_coord(1:3)';
        end
        ud.M = handles.Mvxl2lph_atlas;
        ud.sender = myWin;
        set(myWin, 'UserData', ud);

        % Handle coordinate space conversion
        if strcmpi(opts.CoordSpace, 'LPH')
            % Convert LPH → voxel indices
            Mlph2vxl = inv(handles.Mvxl2lph_atlas);
            targetRCS = Mlph2vxl * [coord 1]';   
            handles.rr = targetRCS(1);
            handles.cc = targetRCS(2);
            handles.ss = targetRCS(3);
        else
            % Already in voxel indices
            handles.rr = coord(1);
            handles.cc = coord(2);
            handles.ss = coord(3);
        end

        % Refresh display for this coordinate
        handles = updateDisplay_newvol(handles, true);

        % Hide UI elements that need hiding every loop
        if opts.hideUI
            if opts.hideVolName && isfield(handles, 'volName_text') && isgraphics(handles.volName_text)
                set(handles.volName_text, 'Visible', 'off');
            end
            if opts.hideROIText
                roiFields = {'aseg','fiber','aparcaseg','thalamus','pauli'};
                for rf = 1:numel(roiFields)
                    if isfield(handles.anat, roiFields{rf}) && ...
                       isfield(handles.anat.(roiFields{rf}), 'roi_text')
                        try
                            set(handles.anat.(roiFields{rf}).roi_text, 'Visible', 'off');
                        catch
                        end
                    end
                end
            end
        end

        % ROI overlay logic
        if opts.drawROI && ~isempty(opts.parcellation)
            for a = 1:numel(opts.parcellation)
                atlas = opts.parcellation{a};
                if ~isfield(handles.anat, atlas)
                    warning('Parcellation "%s" not found in handles.anat. Skipping.', atlas);
                    continue;
                end

                atlasStruct = handles.anat.(atlas);

                if isfield(opts.roiNames, atlas)
                    roi_list = opts.roiNames.(atlas);
                else
                    roi_list = {};
                end
                if ischar(roi_list) || isstring(roi_list)
                    roi_list = {roi_list};
                end

                expanded_roi_list = {};
                for r = 1:numel(roi_list)
                    keyword = lower(roi_list{r});
                    switch keyword
                        case {'subcortical_all','all_fiber','all_aparcaseg','all_thalamus','all_pauli'}
                            expanded_roi_list = [expanded_roi_list; atlasStruct.uiNames(2:end)];
                        otherwise
                            expanded_roi_list = [expanded_roi_list; roi_list{r}];
                    end
                end

                keep_idx = find(ismember(atlasStruct.uiNames, expanded_roi_list));

                atlasStruct.uiRoiOverlaySelected = unique([atlasStruct.uiRoiOverlaySelected, keep_idx]);
                atlasStruct.uiRoiOverlayIdx = atlasStruct.uiRoiIdx(atlasStruct.uiRoiOverlaySelected);
                atlasStruct.uiRoiOverlayImg = atlasStruct.prob(:,:,:,atlasStruct.uiRoiOverlayIdx);

                if isfield(atlasStruct, 'params')
                    atlasStruct.params.showOverlay = true;
                    atlasStruct.params.showOutline = true;
                end

                handles.anat.(atlas) = atlasStruct;
            end

            handles.anat.params.showOverlay = true;
            handles.anat.params.showOutline = true;
            handles.anat.params.overlayAlpha = 0;
            set(handles.overlay_cb, 'Value', true);

            % Draw ROI outlines
            sc = calculateScale(handles);
            switch handles.ORIENTATION
                case 1 
                    sliceIndex = scaleCoord(handles.ss, sc.s);
                    ima = squeeze(V.imgs(:,:,sliceIndex));
                    rr = 1:size(ima,1); cc = 1:size(ima,2); ss = sliceIndex;
                    sliceLabel = 'CR';
                case 2
                    sliceIndex = scaleCoord(handles.cc, sc.c);
                    ima = squeeze(V.imgs(:,sliceIndex,:));
                    rr = 1:size(ima,1); ss = 1:size(ima,2); cc = sliceIndex;
                    sliceLabel = 'SR';
                case 3
                    sliceIndex = scaleCoord(handles.rr, sc.r);
                    ima = squeeze(V.imgs(sliceIndex,:,:));
                    cc = 1:size(ima,1); ss = 1:size(ima,2); rr = sliceIndex;
                    sliceLabel = 'CS';
            end
            [~, outline] = anatomyAddRoiOverlay_save_images(handles, ima, rr, cc, ss, 'outline');
            handles = anatomyDrawRoiOutline_save_images(handles, outline, sliceLabel, opts.lw, opts.fixedColor);
            guidata(myWin, handles);
        end

        drawnow; pause(0.3);

        % -------------------------
        % Save screenshots
        % -------------------------
        if strcmpi(opts.CoordSpace, 'RCS')
            % 3-digit padded for RCS (e.g., 096)
            coord_str = sprintf('X%03d_Y%03d_Z%03d', coord(1), coord(2), coord(3));
        else
            % 2-digit padded for LPH, with N for negatives (e.g., N03)
            coord_parts = cell(1,3);
            for ii = 1:3
                if coord(ii) < 0
                    coord_parts{ii} = sprintf('N%02d', abs(coord(ii)));
                else
                    coord_parts{ii} = sprintf('%02d', coord(ii));
                end
            end
            coord_str = sprintf('X%s_Y%s_Z%s', coord_parts{:});
        end

        switch handles.ORIENTATION
            case 1, orientation_name = 'coronal';
            case 2, orientation_name = 'sagittal';
            case 3, orientation_name = 'axial';
            otherwise, orientation_name = sprintf('orient%d', handles.ORIENTATION);
        end

        outdir_full = fullfile(outdir, orientation_name);
        if ~exist(outdir_full, 'dir'), mkdir(outdir_full); end

        vol_str = sprintf('vol%03d', v);

        if ~isempty(fstem{v})
            % coord_fstem_vol###
            fname_base = sprintf('%s_%s_%s', coord_str, fstem{v}, vol_str);
        else
            % coord_vol###
            fname_base = sprintf('%s_%s', coord_str, vol_str);
        end

        % Final PNG path
        fname_out = fullfile(outdir_full, [fname_base '.png']);

        shot = getframe(handles.axes1);
        imwrite(shot.cdata, fname_out);
        fprintf(' → Saved: %s\n', fname_out);
    end
end