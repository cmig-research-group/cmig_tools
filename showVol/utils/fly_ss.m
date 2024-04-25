%% ------------------------------------------------------------
%% --- Fly to slice number stop_ss
%% ------------------------------------------------------------

function handles = fly_ss(handles, stop_ss, pauseTime)

% Start
start_ss = handles.ss;

% Stop
stop_ss = max(stop_ss, 1);
stop_ss = min(stop_ss, handles.ssMax);

% Step
if stop_ss > handles.ss
  step = +1;
else
  step = -1;
end

% Begin loop
set(handles.toggleLines1, 'Value', 1);
set(handles.toggleLines1, 'ForegroundColor', 'k');
set(handles.axes1_X, 'Visible', 'off');
for ss = start_ss : step : stop_ss
  handles.ss = ss;
  handles = updateDisplay_ss(handles);
  printInfo(handles); % !@#
  pause(pauseTime);
end
if handles.HIDE1 == 0
  set(handles.toggleLines1, 'Value', 0);
  set(handles.toggleLines1, 'ForegroundColor', 'r');
  set(handles.axes1_X, 'Visible', 'on');
end
printInfo(handles);

%% end of function handles = fly_ss(handles, stop_ss, pauseTime)
%% ------------------------------------------------------------