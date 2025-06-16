%% ------------------------------------------------------------
%% --- Movie-cycle through the different input volumes
%% ------------------------------------------------------------

function handles = cycleVolumesMovie(handles, numLoops, pauseTime)
    numVisible = length(find(handles.hideVol == 0));
    if numVisible > 1
      for counter = 1 : numLoops * numVisible
        handles = cycleVolumes(handles);
        if pauseTime > -eps
          pause(pauseTime);
        end
      end
    end
    
    %% ------------------------------------------------------------