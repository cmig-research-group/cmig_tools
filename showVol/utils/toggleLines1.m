%% ------------------------------------------------------------
%% --- Toggle the display of crosshair in axes1
%% ------------------------------------------------------------
function handles = toggleLines1(handles)
handles.HIDE1 = 1 - handles.HIDE1;
set(handles.toggleLines1, 'Value', handles.HIDE1);
if handles.HIDE1
  set(handles.toggleLines1, 'ForegroundColor', 'k');
  set(handles.axes1_X, 'Visible', 'off');
else
  set(handles.toggleLines1, 'ForegroundColor', 'r');
  set(handles.axes1_X, 'Visible', 'on');
end
%% ------------------------------------------------------------
