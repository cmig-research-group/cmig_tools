%% ------------------------------------------------------------
%% --- Executes on button press in storePoint.
%% ------------------------------------------------------------

function storePoint_Callback(hObject, eventdata, handles)
handles = storeCurrentPoint(handles);
guidata(hObject, handles);

%% ------------------------------------------------------------