

%% ------------------------------------------------------------
%% --- Anatomical ROI Text Display
%% ------------------------------------------------------------
function handles = anatomyShowROIString(handles)
%this is called from updateDisplay_zoom

if ~handles.hasABCDBrain
  return
end

if ~handles.anat.params.showOverlay || ~( isfield(handles.anat,'aseg') || isfield(handles.anat,'fiber') || isfield(handles.anat,'aparcaseg'))  ...
    || ~handles.showZoomedInset
  try
    set(handles.anat.aseg.roi_text,'visible','off')
    set(handles.anat.fiber.roi_text,'visible','off')
    set(handles.anat.aparcaseg.roi_text,'visible','off')
  catch
  end
  return
end

tmpax = gca;
tmpf = gcf;
set(handles.figure, 'CurrentAxes', handles.axes1);
xl = get(handles.axes1,'XLim');
yl = get(handles.axes1,'YLim');
revy = strcmp(get(handles.axes1,'YDir'), 'reverse');
%show in upper right of main axis
xx = xl(1) + 0.6*range(xl);
xxc = xl(1) + 0.05*range(xl);
space = 0.1 * range(yl);
if revy
  yya = yl(1);
  yyf = yl(1) + space;
  yyc = yl(1) + 0.5*space;
else
  yya = yl(2);
  yyf = yl(2) - space;
  yyc = yl(2) - 0.5*space;
end

if isfield(handles.anat,'aseg') && handles.anat.aseg.showNames
  roi = anatomyFromRCS(handles, 'aseg');
  roiStr = join(roi.str,newline);
  roiStr = ['-ASEG-' newline roiStr{1}];
  if isempty(handles.anat.aseg.roi_text) || ~ishandle(handles.anat.aseg.roi_text)
    handles.anat.aseg.roi_text  = text(xx,yya,'aseg ROI','fontsize',10,'VerticalAlignment','top','color',[1 1 1],'interpreter','none','fontweight','bold');
  end
  set(handles.anat.aseg.roi_text, 'String',roiStr,'Position',[xx yya 0.1],'visible','on');
else
  try
    set(handles.anat.aseg.roi_text,'visible','off')
  catch
  end
end

if isfield(handles.anat,'fiber') && handles.anat.fiber.showNames
  roi = anatomyFromRCS(handles, 'fiber');
  roiStr = join(roi.str,newline);
  roiStr = ['-FIBER-' newline roiStr{1}];
  if isempty(handles.anat.fiber.roi_text) || ~ishandle(handles.anat.fiber.roi_text)
    handles.anat.fiber.roi_text = text(xx,yyf,'fiber ROI','fontsize',10,'VerticalAlignment','top','color',[1 1 1],'interpreter','none','fontweight','bold');
  end
  set(handles.anat.fiber.roi_text, 'String',roiStr,'Position',[xx yyf 0.1],'visible','on');
else
  try
    set(handles.anat.fiber.roi_text,'visible','off')
  catch
  end
end

if isfield(handles.anat,'aparcaseg') && handles.anat.aparcaseg.showNames
  roi = anatomyFromRCS(handles, 'aparcaseg');
  if roi.prob(1)>0
    isCtx = contains(roi.str,'ctx-'); %limit to cortical ROIs
    if ~any(isCtx)
      roi.str = {'none'};
    else
      roi.str = roi.str(isCtx);
    end
  end
  roiStr = join(roi.str, newline);
  roiStr = ['-APARC-' newline roiStr{1}];
  if isempty(handles.anat.aparcaseg.roi_text) || ~ishandle(handles.anat.aparcaseg.roi_text)
    handles.anat.aparcaseg.roi_text  = text(xxc,yyc,'aparc ROI','fontsize',10,'VerticalAlignment','top','color',[1 1 1],'interpreter','none','fontweight','bold');
  end
  set(handles.anat.aparcaseg.roi_text, 'String',roiStr,'Position',[xxc yyc 0.1],'visible','on');
else
  try
    set(handles.anat.aparcaseg.roi_text,'visible','off')
  catch
  end
end

guidata(handles.figure, handles)
set(tmpf, 'CurrentAxes', tmpax);

%% ------------------------------------------------------------
