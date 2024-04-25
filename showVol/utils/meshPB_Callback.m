%% ------------------------------------------------------------
%% --- Executes on button press in the mesh pushbuttons
%% ------------------------------------------------------------

function meshPB_Callback(hObject, eventdata, handles)
[minVal, chosenSurf] = min(abs(handles.meshPB - hObject));
handles.showSurf(chosenSurf) = 1 - handles.showSurf(chosenSurf);
if handles.showSurf(chosenSurf)
  set(handles.meshCrossSectionCR(chosenSurf).lh, 'Visible', 'on');
  set(handles.meshCrossSectionSR(chosenSurf).lh, 'Visible', 'on');
  set(handles.meshCrossSectionCS(chosenSurf).lh, 'Visible', 'on');
  set(hObject, 'ForegroundColor', handles.surfs(chosenSurf).color);
  set(hObject, 'TooltipString', 'Hide this mesh');
else
  set(handles.meshCrossSectionCR(chosenSurf).lh, 'Visible', 'off');
  set(handles.meshCrossSectionSR(chosenSurf).lh, 'Visible', 'off');
  set(handles.meshCrossSectionCS(chosenSurf).lh, 'Visible', 'off');
  set(hObject, 'ForegroundColor', [0.5 0.5 0.5]);
  set(hObject, 'TooltipString', 'Show this mesh');
end
guidata(hObject, handles);

%% ------------------------------------------------------------