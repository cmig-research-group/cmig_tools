%% ------------------------------------------------------------
%% --- Executes on button press in the volume pushbuttons
%% ------------------------------------------------------------

function volPB_Callback(hObject, eventdata, handles)
if isnumeric(hObject)
  [minVal, chosenVol] = min(abs(handles.volPB - hObject));
else
  [minVal, chosenVol] = min(cellfun(@(x)norm(x-hObject.Position),{handles.volPB(:).Position}'));
end
if chosenVol == handles.currentVol
  if sum(handles.hideVol) < handles.numVols-1
    handles.hideVol(chosenVol) = 1;
    set(hObject, 'ForegroundColor', [0.5 0.5 0.5]);
    set(hObject, 'TooltipString', ['Show ' handles.volnames{handles.currentVol}]);
    tmp = find(handles.hideVol == 0);
    tmpInd = find(tmp > handles.currentVol);
    if isempty(tmpInd)
      tmpInd = 1;
    end
    handles.currentVol = tmp(tmpInd(1));
    handles.currentVolPB = handles.volPB(handles.currentVol);
  end
else
  set(handles.currentVolPB, 'ForegroundColor', 'k');
  set(handles.currentVolPB, 'TooltipString', ['Show ' handles.volnames{handles.currentVol}]);
  handles.hideVol(chosenVol) = 0;
  handles.currentVol = chosenVol;
  handles.currentVolPB = hObject;
end
set(handles.currentVolPB, 'ForegroundColor', 'r');
set(handles.currentVolPB, 'TooltipString', ['Hide ' handles.volnames{handles.currentVol}]);
handles.Mvxl2lph = handles.vols{handles.currentVol}.Mvxl2lph;
handles = updateDisplay_newvol(handles);
guidata(hObject, handles);

%% ------------------------------------------------------------
