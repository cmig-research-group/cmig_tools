%% ------------------------------------------------------------
%% --- Forward-cycle the orientations shown in each axis
%% ------------------------------------------------------------

function handles = cycleOrientationsForward(handles)
handles.ORIENTATION = handles.ORIENTATION + 1;
if handles.ORIENTATION > 3
  handles.ORIENTATION = 1;
end
handles = displayNewOrientation(handles);

%% ------------------------------------------------------------