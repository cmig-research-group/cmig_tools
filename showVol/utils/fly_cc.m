%% ------------------------------------------------------------
%% --- Fly to column number stop_cc
%% ------------------------------------------------------------

function handles = fly_cc(handles, stop_cc, pauseTime)

% Start
start_cc = handles.cc;

% Stop
stop_cc = max(stop_cc, 1);
stop_cc = min(stop_cc, handles.ccMax);

% Step
if stop_cc > handles.cc
  step = +1;
else
  step = -1;
end

% Begin loop
set(handles.toggleLines1, 'Value', 1);
set(handles.toggleLines1, 'ForegroundColor', 'k');
set(handles.axes1_X, 'Visible', 'off');
for cc = start_cc : step : stop_cc
  handles.cc = cc;
  handles = updateDisplay_cc(handles);
  printInfo(handles); % !@#
  pause(pauseTime);
end
if handles.HIDE1 == 0
  set(handles.toggleLines1, 'Value', 0);
  set(handles.toggleLines1, 'ForegroundColor', 'r');
  set(handles.axes1_X, 'Visible', 'on');
end

%% end of function handles = fly_cc(handles, stop_cc, pauseTime)
%% ------------------------------------------------------------