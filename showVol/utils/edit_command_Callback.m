%% ------------------------------------------------------------
%% --- Executes on entering into edit_command edit box
%% ------------------------------------------------------------

function edit_command_Callback(hObject, eventdata, handles)
str = get(hObject, 'String');
set(handles.text_last_echo, 'String', str);
set(hObject, 'String', '');
handles = parseCommand(handles, str);
guidata(hObject, handles);

%% ------------------------------------------------------------