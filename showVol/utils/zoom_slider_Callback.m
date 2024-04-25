%% ------------------------------------------------------------
%% --- Executes on zoom slider movement.
%% ------------------------------------------------------------

function zoom_slider_Callback(hObject, eventdata, handles)
handles.zoom = get(hObject,'Value');
handles.rrZoom = handles.rr;
handles.ccZoom = handles.cc;
handles.ssZoom = handles.ss;
handles = updateDisplay_zoom(handles);

guidata(hObject, handles);

%% ------------------------------------------------------------