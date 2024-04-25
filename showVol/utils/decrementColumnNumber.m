%% ------------------------------------------------------------
%% --- Decrement column number
%% ------------------------------------------------------------

function handles = decrementColumnNumber(handles)
if handles.cc > 1
  handles.cc = handles.cc - 1;
end
handles = updateDisplay_cc(handles);
broadcastLPH(handles)
printInfo(handles);

%% ------------------------------------------------------------