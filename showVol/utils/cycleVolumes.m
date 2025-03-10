%% ------------------------------------------------------------
%% --- Cycle through the different input volumes
%% ------------------------------------------------------------

function handles = cycleVolumes(handles)
    set(handles.currentVolPB, 'ForegroundColor', 'k');
    tmp = find(handles.hideVol == 0);
    tmpInd = find(tmp > handles.currentVol);
    if isempty(tmpInd)
      tmpInd = 1;
    end
    handles.currentVol = tmp(tmpInd(1));
    handles.currentVolPB = handles.volPB(handles.currentVol);
    set(handles.currentVolPB, 'ForegroundColor', 'r');
    handles.Mvxl2lph = handles.vols{handles.currentVol}.Mvxl2lph;
    updateDisplay_newvol(handles);
    
    %% ------------------------------------------------------------