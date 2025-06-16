%% ------------------------------------------------------------
%% -- Swap the current RCS point with the stored RCS point
%% ------------------------------------------------------------

function handles = swapCurrentAndStoredPoints(handles)
    set(handles.swapPoint, 'Value', 0);
    tmp = [handles.rr handles.cc handles.ss];
    handles.rr = handles.rrStored;
    handles.cc = handles.ccStored;
    handles.ss = handles.ssStored;
    handles.rrStored = tmp(1);
    handles.ccStored = tmp(2);
    handles.ssStored = tmp(3);
    str = sprintf('Go to voxel (%d,%d,%d) -- keyboard shortcut: ;', ...
      handles.rrStored, handles.ccStored, handles.ssStored);
    set(handles.gotoPoint,'TooltipString', str);
    str = sprintf('Swap current voxel with voxel (%d,%d,%d) -- keyboard shortcut: /', ...
      handles.rrStored, handles.ccStored, handles.ssStored);
    set(handles.swapPoint,'TooltipString', str);
    handles = updateDisplay_newvol(handles);
    %handles = updateDisplay_rr(handles);
    %handles = updateDisplay_cc(handles);
    %handles = updateDisplay_ss(handles);
    printInfo(handles);
    
    %% ------------------------------------------------------------