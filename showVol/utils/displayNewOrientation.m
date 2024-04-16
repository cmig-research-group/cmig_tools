%% ------------------------------------------------------------
%% --- Change the display upon setting a new orientation
%% ------------------------------------------------------------
function handles = displayNewOrientation(handles)
switch handles.ORIENTATION
  
  case 1 %coronal main
    handles.axesCR = handles.axes1;
    handles.axesSR = handles.axes2;
    handles.axesCS = handles.axes3;
    handles = newDisplay(handles);
    handles.axes1_X = [handles.ClineCR handles.RlineCR ...
      handles.CtextCR handles.RtextCR];
    handles.axes2_X = [handles.SlineSR handles.RlineSR ...
      handles.StextSR handles.RtextSR];
    handles.axes3_X = [handles.ClineCS handles.SlineCS ...
      handles.CtextCS handles.StextCS];
    handles.ccDecrChar = 'j';
    handles.ccIncrChar = 'l';
    handles.rrDecrChar = 'i';
    handles.rrIncrChar = 'k';
    handles.ssDecrChar = ',';
    handles.ssIncrChar = '.';
    
  case 2 %sagittal main
    handles.axesSR = handles.axes1;
    handles.axesCS = handles.axes2;
    handles.axesCR = handles.axes3;
    handles = newDisplay(handles);
    handles.axes1_X = [handles.SlineSR handles.RlineSR ...
      handles.StextSR handles.RtextSR];
    handles.axes2_X = [handles.ClineCS handles.SlineCS ...
      handles.CtextCS handles.StextCS];
    handles.axes3_X = [handles.ClineCR handles.RlineCR ...
      handles.CtextCR handles.RtextCR];
    handles.ssDecrChar = 'j';
    handles.ssIncrChar = 'l';
    handles.rrDecrChar = 'i';
    handles.rrIncrChar = 'k';
    handles.ccDecrChar = ',';
    handles.ccIncrChar = '.';
    
  case 3 %axial main
    handles.axesCS = handles.axes1;
    handles.axesCR = handles.axes2;
    handles.axesSR = handles.axes3;
    handles = newDisplay(handles);
    handles.axes1_X = [handles.ClineCS handles.SlineCS ...
      handles.CtextCS handles.StextCS];
    handles.axes2_X = [handles.ClineCR handles.RlineCR ...
      handles.CtextCR handles.RtextCR];
    handles.axes3_X = [handles.SlineSR handles.RlineSR ...
      handles.StextSR handles.RtextSR];
    handles.ccDecrChar = 'j';
    handles.ccIncrChar = 'l';
    handles.ssDecrChar = 'k';
    handles.ssIncrChar = 'i';
    handles.rrDecrChar = ',';
    handles.rrIncrChar = '.';
    
end

if handles.HIDE1
  set(handles.axes1_X, 'Visible', 'off');
end
if handles.HIDE2
  set(handles.axes2_X, 'Visible', 'off');
end
if handles.HIDE3
  set(handles.axes3_X, 'Visible', 'off');
end

%handles = updateDisplay_newvol(handles);

%handles = updateDisplay_zoom(handles);

%% end of function handles = displayNewOrientation(handles)
%% ------------------------------------------------------------