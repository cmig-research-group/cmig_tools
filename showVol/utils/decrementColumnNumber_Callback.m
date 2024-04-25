%% ------------------------------------------------------------
%% --- Executes on button press in decrementColumnNumber.
%% ------------------------------------------------------------

function decrementColumnNumber_Callback(hObject, eventdata, handles)
handles = decrementColumnNumber(handles);
guidata(hObject, handles);

%% ------------------------------------------------------------