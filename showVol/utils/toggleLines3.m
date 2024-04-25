%% ------------------------------------------------------------
%% --- Toggle the display of crosshair in axes3
%% ------------------------------------------------------------
function handles = toggleLines3(handles)
handles.HIDE3 = 1 - handles.HIDE3;
set(handles.toggleLines3, 'Value', handles.HIDE3);
if handles.HIDE3
  set(handles.toggleLines3, 'ForegroundColor', 'k');
  set(handles.axes3_X, 'Visible', 'off');
else
  set(handles.toggleLines3, 'ForegroundColor', 'r');
  set(handles.axes3_X, 'Visible', 'on');
end
%% ------------------------------------------------------------