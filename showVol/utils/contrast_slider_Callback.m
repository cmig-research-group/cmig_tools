%% ------------------------------------------------------------
%% --- Executes on contrast slider movement.
%% ------------------------------------------------------------

function contrast_slider_Callback(hObject, eventdata, handles)
cVal = get(handles.contrast_slider, 'Value');
handles.cVal(handles.currentVol) = cVal;
handles = updateDisplay_clims(handles);
guidata(hObject, handles);

%% ------------------------------------------------------------