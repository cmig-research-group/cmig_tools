%% ------------------------------------------------------------
%% --- Executes on button press in contrast_reset.
%% ------------------------------------------------------------

function contrast_reset_Callback(hObject, eventdata, handles)
set(handles.contrast_slider, 'Value', 1);
handles.cVal(handles.currentVol) = 1;
handles = updateDisplay_clims(handles);
guidata(hObject, handles);

%% ------------------------------------------------------------