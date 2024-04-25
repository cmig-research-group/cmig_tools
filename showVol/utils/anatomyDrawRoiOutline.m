%% ------------------------------------------------------------
%% --- Draw roi outlines
%% ------------------------------------------------------------
function handles = anatomyDrawRoiOutline(handles,outline,slice)

if ~handles.hasABCDBrain
  return
end

%axis is the axis slice we're drawing into
fname = ['roiContour' slice];
if ~isfield(handles, fname), return; end
delete(handles.(fname))
handles.(fname) = [];

axname = ['axes' slice];
delete(findobj(handles.(axname),'tag','roiOutline'))
delete(findobj(handles.(axname),'tag','roiLabel'))

if isempty(outline), return, end

doLabel = (handles.(axname) == handles.axes1); %only show labels in main axis

lw = 1.25;
[outlineCoord, outlineColor, outlineLabel] = deal(outline{:});
scalarVol = isfield(handles.vols{handles.currentVol}, 'imgs') && size(handles.vols{handles.currentVol}.imgs,4)==1 ...
  || (isfield(handles.vols{handles.currentVol}, 'limits') && length(handles.vols{handles.currentVol}.limits) > 2); %

for iR = 1:length(outlineCoord)
  if  scalarVol && ~handles.anat.params.isTracking
    h1 = line(outlineCoord{iR}(1,:), outlineCoord{iR}(2,:), 'color', 'k', 'linewidth', lw*2,'tag','roiOutline'); %black background for scalar volumes only; don't draw when tracking, for speed
  else
    h1 = [];
  end
  h2 = line(outlineCoord{iR}(1,:), outlineCoord{iR}(2,:), 'color', outlineColor(iR,:), 'linewidth', lw,'tag','roiOutline');
  label = outlineLabel{iR,3};
  if ~isempty(label) && ~handles.anat.params.isTracking %&& doLabel %don't draw labels when tracking
    h3 = text(outlineLabel{iR,1}, outlineLabel{iR,2}, label, 'color','w', 'fontsize',7,'fontweight','bold',...
      'horizontalalignment','center','verticalalignment','middle','tag','roiLabel');
  else
    h3=[];
  end
  handles.(fname) = [handles.(fname) h1 h2 h3];
end
guidata(handles.figure, handles);

% this old method was efficient, but can't handle lines with different colors
% if isempty(handles.(fname)) || any(~ishandle(handles.(fname)))
%     delete(handles.(fname))
%     handles.(fname)(1) = line(1,1,'color','k','linewidth',lw*2,'visible','off');
%     handles.(fname)(2) = line(1,1,'color','w','linewidth',lw,'visible','off');
%     handles.(fname)(3) = line(1,1,'color','w','linewidth',lw/2,'linestyle',':','visible','off');
% end
% if isempty(outline)
%     set(handles.(fname),'Visible','off')
% else
%     if iscell(outline)
%         set(handles.(fname)(1), 'XData', outline{1}(1,:), 'YData', outline{1}(2,:),'visible','on')
%         set(handles.(fname)(2), 'XData', outline{1}(1,:), 'YData', outline{1}(2,:),'visible','on')
%         set(handles.(fname)(3), 'XData', outline{2}(1,:), 'YData', outline{2}(2,:),'visible','on')
%     else
%         set(handles.(fname)(1), 'XData', outline(1,:), 'YData', outline(2,:),'visible','on')
%         set(handles.(fname)(2), 'XData', outline(1,:), 'YData', outline(2,:),'visible','on')
%         set(handles.(fname)(3), 'visible','off')
%     end
% end

%% ------------------------------------------------------------
