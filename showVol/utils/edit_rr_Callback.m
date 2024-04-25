%% ------------------------------------------------------------
%% --- Executes on entering into edit_rr edit box
%% ------------------------------------------------------------

function edit_rr_Callback(hObject, eventdata, handles)
rr = round(str2double(get(hObject,'String')));
if isfinite(rr)
  if (rr >= 1) & (rr <= handles.rrMax)
    handles.rr = rr;
    handles = updateDisplay_rr(handles);
    printInfo(handles);
    guidata(hObject, handles);
  end
end
set(hObject,'String',num2str(handles.rr));

%% ------------------------------------------------------------