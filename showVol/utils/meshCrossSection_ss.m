

%% ------------------------------------------------------------
%% --- Compute meshCrossSection for slice = constant
%% ------------------------------------------------------------
function handles = meshCrossSection_ss(handles)

for mmm = 1:handles.numSurfs
  
  %% Draw in axesCR, which draws imgs(:,:,k)
  axes(handles.axesCR);
  [startPts, endPts] = findContour(handles.surfs(mmm), ...
    handles.surfs(mmm).vertices(:,3), ...
    handles.ss);
  
  startR = startPts(:,1); endR = endPts(:,1);
  startC = startPts(:,2); endC = endPts(:,2);
  
  hold on;
  lh = quiver(startC, startR, endC-startC, endR-startR, 0, '.');
  hold off
  
  set(lh, 'Color', handles.surfs(mmm).color);
  handles.meshCrossSectionCR(mmm).lh = lh;
  
  if handles.showSurf(mmm) == 0
    set(handles.meshCrossSectionCR(mmm).lh, 'Visible', 'off');
  end
  
end

%% end of function handles = meshCrossSection_ss(handles)
%% ------------------------------------------------------------
