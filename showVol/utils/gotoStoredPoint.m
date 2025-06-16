%% ------------------------------------------------------------
%% -- Go to the stored RCS point
%% ------------------------------------------------------------

function handles = gotoStoredPoint(handles)
    set(handles.gotoPoint, 'Value', 0);
    handles.rr = handles.rrStored;
    handles.cc = handles.ccStored;
    handles.ss = handles.ssStored;
    handles = updateDisplay_newvol(handles);
    %handles = updateDisplay_rr(handles);
    %handles = updateDisplay_cc(handles);
    %handles = updateDisplay_ss(handles);
    printInfo(handles);
    
    %% ------------------------------------------------------------