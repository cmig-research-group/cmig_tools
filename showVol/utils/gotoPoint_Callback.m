
%% ------------------------------------------------------------
%% --- Executes on button press in gotoPoint.
%% ------------------------------------------------------------

function gotoPoint_Callback(hObject, eventdata, handles)
handles = gotoStoredPoint(handles);
guidata(hObject, handles);

%% ------------------------------------------------------------
