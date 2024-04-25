%% ------------------------------------------------------------
%% --- Executes on entering into edit_cc edit box
%% ------------------------------------------------------------

function edit_cc_Callback(hObject, eventdata, handles)
cc = round(str2double(get(hObject,'String')));
if isfinite(cc)
  if (cc >= 1) & (cc <= handles.ccMax)
    handles.cc = cc;
    handles = updateDisplay_cc(handles);
    printInfo(handles);
    guidata(hObject, handles);
  end
end
set(hObject,'String',num2str(handles.cc));

%% ------------------------------------------------------------