%% ------------------------------------------------------------
%% --- Executes on button press in decrementSliceNumber.
%% ------------------------------------------------------------

function decrementSliceNumber_Callback(hObject, eventdata, handles)
handles = decrementSliceNumber(handles);
guidata(hObject, handles);

%% ------------------------------------------------------------
