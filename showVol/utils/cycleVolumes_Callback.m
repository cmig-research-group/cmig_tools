%% ------------------------------------------------------------
%% --- Executes on button press in cycleVolumes.
%% ------------------------------------------------------------

function cycleVolumes_Callback(hObject, eventdata, handles)
handles = cycleVolumes(handles);
guidata(hObject, handles);
%% ------------------------------------------------------------