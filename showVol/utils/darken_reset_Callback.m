%% ------------------------------------------------------------
%% --- Executes on button press in darken_reset.
%% ------------------------------------------------------------

function darken_reset_Callback(hObject, eventdata, handles)
set(handles.darken_slider, 'Value', 0.5);
handles.dVal(handles.currentVol) = 0.5;
handles = updateDisplay_clims(handles);
guidata(hObject, handles);

%% ------------------------------------------------------------