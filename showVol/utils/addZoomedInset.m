
%% ------------------------------------------------------------
%% --- add zoomed inset for FOD images
%% ------------------------------------------------------------

function [ima,handles] = addZoomedInset(handles, ima, r0, c0, scale, location)
% ima = addZoomedInset(handles, ima, r0, c0, scale, location)
%   show zoomed FOD image in main axis
%
% disables inset if current volume is not FOD
%
% addZoomedInset(handles, true) called after changing the zoom to
%   turn off when the FOD image itself is sufficiently zoomed in

handles.zoomedInsetAxes.Visible = 'off';

if gca ~= handles.axes1, return, end %only show in main axis; presumes caller has set gca to current axis (which updateDisplay_newVol does)

V = handles.vols{handles.currentVol};
isFOD = isfield(V,'imgs1');
zoomLimit = 2.8;

%turn off axes if this is not FOD, or if showing axis is turned off by user (` key toggles)
if ~isFOD || ~handles.showZoomedInset
  set(handles.zoomedInsetAxes.Children, 'Visible', 'off');
  return
end

% Usage 2: when zooming, if zoomed in, turn off above a certain zoom
if (numel(ima)==1 && ima==true)
  if handles.zoom > zoomLimit
    set(handles.zoomedInsetAxes.Children, 'Visible', 'off');
  else
    set(handles.zoomedInsetAxes.Children, 'Visible', 'on');
  end
  return
end

srcSize = round(5 * scale); %+/- 5mm around cursor
shift = floor(scale/2); %FIXME for 1/2mm shift of FOD images relative to other volumes
try
  subImage = ima(r0*scale + (-srcSize:srcSize)+shift, c0*scale + (-srcSize:srcSize)+shift , :);
  if strcmp(location,'ur') %axial slice, flipped up/down
    subImage = flipud(subImage);
  end
  if ~ishandle(handles.zoomedInsetImage)
    cla(handles.zoomedInsetAxes)
    handles.zoomedInsetImage = image(handles.zoomedInsetAxes, subImage);
    cx = mean(xlim(handles.zoomedInsetAxes));
    cy = mean(ylim(handles.zoomedInsetAxes));
    hold(handles.zoomedInsetAxes, 'on')
    plot(handles.zoomedInsetAxes, cx,cy,'w+')
  else
    handles.zoomedInsetImage.CData = subImage;
  end
  
  if handles.zoom > zoomLimit
    set(handles.zoomedInsetAxes.Children, 'Visible', 'off');
  else
    set(handles.zoomedInsetAxes.Children, 'Visible', 'on');
  end
  
catch
  
end

%% end of function ima = addZoomedInset(ima, r0, c0, show, scale)
%% ------------------------------------------------------------
