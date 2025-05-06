%% ------------------------------------------------------------
%% --- Compute new zoom values and update display
%% ------------------------------------------------------------

function handles = updateDisplay_zoom(handles) % Should be called upon initial display and change of volume

if handles.zoom < sqrt(eps)
  
  handles.rrLim = handles.rrLimOrig;
  handles.ccLim = handles.ccLimOrig;
  handles.ssLim = handles.ssLimOrig;
  
else
  
  zoom = 2^handles.zoom;
  halfRange_mm = handles.halfRange_mm/zoom;
  tmp = [-halfRange_mm +halfRange_mm];
  handles.rrLim = handles.rrZoom + tmp/handles.rcsScale(1);
  handles.ccLim = handles.ccZoom + tmp/handles.rcsScale(2);
  handles.ssLim = handles.ssZoom + tmp/handles.rcsScale(3);
  
end

%sf = calculateScale(handles);
sf.r=1; %we've made the axis always in mm (in updateDisplay_newVol) so scaling no longer needed
sf.c=1;
sf.s=1;

set(handles.axesCR, 'Xlim', sf.c * handles.ccLim);
set(handles.axesCR, 'Ylim', sf.r * handles.rrLim);
set(handles.axesSR, 'Xlim', sf.s * handles.ssLim);
set(handles.axesSR, 'Ylim', sf.r * handles.rrLim);
set(handles.axesCS, 'Xlim', sf.c * handles.ccLim);
set(handles.axesCS, 'Ylim', sf.s * handles.ssLim);

set(handles.CtextCR, 'Position', [sf.c*handles.cc       sf.r*handles.rrLim(1)]);
set(handles.RtextCR, 'Position', [sf.c*handles.ccLim(1) sf.r*handles.rr]);
set(handles.StextSR, 'Position', [sf.s*handles.ss       sf.r*handles.rrLim(1)]);
set(handles.RtextSR, 'Position', [sf.s*handles.ssLim(1) sf.r*handles.rr]);
set(handles.CtextCS, 'Position', [sf.c*handles.cc       sf.s*handles.ssLim(2)]);
set(handles.StextCS, 'Position', [sf.c*handles.ccLim(1) sf.s*handles.ss]);

handles.ClineCR.YData = sf.r*handles.rrLim;
handles.SlineSR.YData = sf.r*handles.rrLim;

addZoomedInset(handles, true); %turns off inset at higher zooms

handles = anatomyShowROIString(handles);

handles = showVolumeName(handles);

handles = annotationDraw(handles);

%% end of function handles = updateDisplay_zoom(handles)