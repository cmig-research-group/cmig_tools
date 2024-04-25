%% ------------------------------------------------------------
%% --- Executes on entering into edit_ss edit box
%% ------------------------------------------------------------

function edit_ss_Callback(hObject, eventdata, handles)
ss = round(str2double(get(hObject,'String')));
if isfinite(ss)
  if (ss >= 1) & (ss <= handles.ssMax)
    handles.ss = ss;
    handles = updateDisplay_ss(handles);
    printInfo(handles);
    guidata(hObject, handles);
  end
end
set(hObject,'String',num2str(handles.ss));

%% ------------------------------------------------------------