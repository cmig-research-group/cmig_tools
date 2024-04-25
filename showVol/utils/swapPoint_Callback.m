%% ------------------------------------------------------------
%% --- Executes on button press in swapPoint.
%% ------------------------------------------------------------

function swapPoint_Callback(hObject, eventdata, handles)
handles = swapCurrentAndStoredPoints(handles);
guidata(hObject, handles);

%% ------------------------------------------------------------