%% ------------------------------------------------------------
%% --- Executes on button press in incrementColumnNumber.
%% ------------------------------------------------------------

function incrementColumnNumber_Callback(hObject, eventdata, handles)
handles = incrementColumnNumber(handles);
guidata(hObject, handles);

%% ------------------------------------------------------------