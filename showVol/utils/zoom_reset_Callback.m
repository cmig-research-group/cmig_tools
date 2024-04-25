%% ------------------------------------------------------------
%% --- Executes on button press in zoom_reset.
%% ------------------------------------------------------------

function zoom_reset_Callback(hObject, eventdata, handles)
defaultZoom = 0.48; %was 0, but this fills out the window better
handles.zoom = defaultZoom;
set(handles.zoom_slider, 'Value', defaultZoom);
%recenter as well
handles.rrZoom = ceil(handles.rrMax/2 + sqrt(eps));
handles.ccZoom = ceil(handles.ccMax/2 + sqrt(eps));
handles.ssZoom = ceil(handles.ssMax/2 + sqrt(eps));
handles = updateDisplay_zoom(handles);
guidata(hObject, handles);

%% ------------------------------------------------------------