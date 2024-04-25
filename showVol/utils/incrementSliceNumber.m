%% ------------------------------------------------------------
%% --- Increment slice number
%% ------------------------------------------------------------

function handles = incrementSliceNumber(handles)
if handles.ss < handles.ssMax
  handles.ss = handles.ss + 1;
end
handles = updateDisplay_ss(handles);
broadcastLPH(handles)
printInfo(handles);

%% ------------------------------------------------------------