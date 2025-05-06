%% ------------------------------------------------------------
%% -- Store the current RCS point
%% ------------------------------------------------------------

function handles = storeCurrentPoint(handles)
    set(handles.storePoint, 'Value', 0);
    handles.rrStored = handles.rr;
    handles.ccStored = handles.cc;
    handles.ssStored = handles.ss;
    str = sprintf('Go to voxel (%d,%d,%d) -- keyboard shortcut: ;', ...
      handles.rrStored, handles.ccStored, handles.ssStored);
    set(handles.gotoPoint,'TooltipString', str);
    str = sprintf('Swap current voxel with voxel (%d,%d,%d) -- keyboard shortcut: /', ...
      handles.rrStored, handles.ccStored, handles.ssStored);
    set(handles.swapPoint,'TooltipString', str);
    
    %% ------------------------------------------------------------