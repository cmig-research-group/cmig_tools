function handles = anatomyDrawRoiOutline_save_images(handles, outline, slice, lw, fixedColor)
% Patched version of anatomyDrawRoiOutline
% Forces all ROI outlines and labels to draw into the correct axes
%    (so MATLAB doesn’t auto-create Figure 1).

if ~handles.hasABCDBrain
    return
end

% figure out which axis to draw into
% axname = ['axes' slice];
% if ~isfield(handles, axname)
%     warning('No axis found for slice %s', slice);
%     return
% end
% target_ax = handles.(axname);
target_ax = handles.axes1;

% make sure we're drawing into this axis
axes(target_ax);
hold(target_ax, 'on');

% clear previous outlines
fname = ['roiContour' slice];
if isfield(handles, fname)
    delete(handles.(fname));
    handles.(fname) = [];
end

delete(findobj(target_ax,'tag','roiOutline'));
delete(findobj(target_ax,'tag','roiLabel'));

if isempty(outline)
    return
end

doLabel = (target_ax == handles.axes1); %only show labels in main axis

%lw = 3;
[outlineCoord, outlineColor, outlineLabel] = deal(outline{:});
scalarVol = (isfield(handles.vols{handles.currentVol}, 'imgs') && size(handles.vols{handles.currentVol}.imgs,4)==1) ...
         || (isfield(handles.vols{handles.currentVol}, 'limits') && length(handles.vols{handles.currentVol}.limits) > 2);

for iR = 1:length(outlineCoord)
    % explicitly pass parent axis to line/text to prevent Figure 1
    if scalarVol && ~handles.anat.params.isTracking
        h1 = line(outlineCoord{iR}(1,:), outlineCoord{iR}(2,:), ...
                  'color', 'k', 'linewidth', lw*2, ...
                  'tag', 'roiOutline', 'Parent', target_ax);
    else
        h1 = [];
    end

    %fixedColor = [0, 0.45, 0.75];

    h2 = line(outlineCoord{iR}(1,:), outlineCoord{iR}(2,:), ...
              'color', fixedColor, 'linewidth', lw, ...
              'tag', 'roiOutline', 'Parent', target_ax);

    % h2 = line(outlineCoord{iR}(1,:), outlineCoord{iR}(2,:), ...
    %           'color', outlineColor(iR,:), 'linewidth', lw, ...
    %           'tag', 'roiOutline', 'Parent', target_ax);

    label = outlineLabel{iR,3};
    if ~isempty(label) && ~handles.anat.params.isTracking
        h3 = text(outlineLabel{iR,1}, outlineLabel{iR,2}, label, ...
                  'color', 'w', 'fontsize', 7, 'fontweight', 'bold', ...
                  'horizontalalignment', 'center', ...
                  'verticalalignment', 'middle', ...
                  'tag', 'roiLabel', 'Parent', target_ax);
    else
        h3 = [];
    end

    handles.(fname) = [handles.(fname) h1 h2 h3];
end

guidata(handles.figure, handles);