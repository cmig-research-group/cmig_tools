%% ------------------------------------------------------------
%% --- Toggle the display of crosshair in axes2
%% ------------------------------------------------------------
function handles = toggleLines2(handles)
handles.HIDE2 = 1 - handles.HIDE2;
set(handles.toggleLines2, 'Value', handles.HIDE2);
if handles.HIDE2
  set(handles.toggleLines2, 'ForegroundColor', 'k');
  set(handles.axes2_X, 'Visible', 'off');
else
  set(handles.toggleLines2, 'ForegroundColor', 'r');
  set(handles.axes2_X, 'Visible', 'on');
end
%% ------------------------------------------------------------
