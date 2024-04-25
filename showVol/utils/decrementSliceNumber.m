%% ------------------------------------------------------------
%% --- Decrement slice number
%% ------------------------------------------------------------

function handles = decrementSliceNumber(handles)
if handles.ss > 1
  handles.ss = handles.ss - 1;
end
handles = updateDisplay_ss(handles);
broadcastLPH(handles)
printInfo(handles);

%% ------------------------------------------------------------