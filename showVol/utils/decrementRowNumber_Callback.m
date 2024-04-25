%% ------------------------------------------------------------
%% --- Executes on button press in decrementRowNumber.
%% ------------------------------------------------------------

function decrementRowNumber_Callback(hObject, eventdata, handles)
handles = decrementRowNumber(handles);
guidata(hObject, handles);

%% ------------------------------------------------------------