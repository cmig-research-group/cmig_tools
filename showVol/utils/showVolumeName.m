%% ------------------------------------------------------------
%% --- Display volume name on figure and title
%% ------------------------------------------------------------
function handles = showVolumeName(handles)
%this is called from updateDisplay_zoom

%get name, if present, and add to figure title
try
  name = handles.vols{handles.currentVol}.name;
  set(handles.figure,'name',['showVol - ' name])
catch
  if isfield(handles,'volName_text') && ishandle(handles.volName_text)
    set(handles.volName_text,'String','', 'visible','off')
  end
  set(handles.figure,'name','showVol')
  return
end

%show on figure when overlays are displayed (needed?)
if handles.anat.params.showOverlay && handles.showZoomedInset
  
  tmpax = gca;
  tmpf = gcf;
  set(handles.figure, 'CurrentAxes', handles.axes1);
  xl = get(handles.axes1,'XLim');
  yl = get(handles.axes1,'YLim');
  revy = strcmp(get(handles.axes1,'YDir'), 'reverse');
  %show in upper left of main axis
  xx = xl(1) + 0.05*range(xl);
  if revy
    yy = yl(1);
  else
    yy = yl(2);
  end
  
  if ~isfield(handles,'volName_text') || ~ishandle(handles.volName_text)
    handles.volName_text = text(xx,yy,name,'fontsize',12,'VerticalAlignment','top','color',[1 1 1],'interpreter','none','fontweight','bold');
  end
  set(handles.volName_text, 'String',name,'Position',[xx yy 0.1],'visible','on');
  
else
  if isfield(handles,'volName_text') && ishandle(handles.volName_text)
    set(handles.volName_text,'visible','off')
  end
end

%% ------------------------------------------------------------