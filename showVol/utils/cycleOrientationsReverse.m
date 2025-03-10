%% ------------------------------------------------------------
%% --- Reverse-cycle the orientations shown in each axis
%% ------------------------------------------------------------

function handles = cycleOrientationsReverse(handles)
    handles.ORIENTATION = handles.ORIENTATION - 1;
    if handles.ORIENTATION < 1
      handles.ORIENTATION = 3;
    end
    handles = displayNewOrientation(handles);
    
    %% ------------------------------------------------------------