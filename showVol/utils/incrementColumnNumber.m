%% ------------------------------------------------------------
%% --- Increment column number
%% ------------------------------------------------------------

function handles = incrementColumnNumber(handles)
if handles.cc < handles.ccMax
  handles.cc = handles.cc + 1;
end
handles = updateDisplay_cc(handles);
broadcastLPH(handles)
printInfo(handles);

%% ------------------------------------------------------------