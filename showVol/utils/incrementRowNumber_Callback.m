%% ------------------------------------------------------------
%% --- Executes on button press in incrementRowNumber.
%% ------------------------------------------------------------

function incrementRowNumber_Callback(hObject, eventdata, handles)
handles = incrementRowNumber(handles);
guidata(hObject, handles);

%% ------------------------------------------------------------
