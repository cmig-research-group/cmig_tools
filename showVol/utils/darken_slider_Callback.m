%% ------------------------------------------------------------
%% --- Executes on darken slider movement.
%% ------------------------------------------------------------

function darken_slider_Callback(hObject, eventdata, handles)
dVal = get(handles.darken_slider, 'Value');
handles.dVal(handles.currentVol) = dVal;
handles = updateDisplay_clims(handles);
guidata(hObject, handles);

%% ------------------------------------------------------------