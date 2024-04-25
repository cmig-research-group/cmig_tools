%% ------------------------------------------------------------
%% --- Fly to row number stop_rr
%% ------------------------------------------------------------

function handles = fly_rr(handles, stop_rr, pauseTime)

% Start
start_rr = handles.rr;

% Stop
stop_rr = max(stop_rr, 1);
stop_rr = min(stop_rr, handles.rrMax);

% Step
if stop_rr > handles.rr
  step = +1;
else
  step = -1;
end

% Begin loop
set(handles.toggleLines1, 'Value', 1);
set(handles.toggleLines1, 'ForegroundColor', 'k');
set(handles.axes1_X, 'Visible', 'off');
for rr = start_rr : step : stop_rr
  handles.rr = rr;
  handles = updateDisplay_rr(handles);
  printInfo(handles); % !@#
  pause(pauseTime);
end
if handles.HIDE1 == 0
  set(handles.toggleLines1, 'Value', 0);
  set(handles.toggleLines1, 'ForegroundColor', 'r');
  set(handles.axes1_X, 'Visible', 'on');
end

%% end of function function handles = fly_rr(handles, stop_rr, pauseTime)

%% ------------------------------------------------------------