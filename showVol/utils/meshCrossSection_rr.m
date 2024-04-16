%% ------------------------------------------------------------
%% --- Compute meshCrossSection for row = constant
%% ------------------------------------------------------------
function handles = meshCrossSection_rr(handles)

for mmm = 1:handles.numSurfs
  
  %% Draw in axis CS, which draws imgs(i,:,:)
  axes(handles.axesCS);
  [startPts, endPts] = findContour(handles.surfs(mmm), ...
    handles.surfs(mmm).vertices(:,1), ...
    handles.rr);
  
  startC = startPts(:,2); endC = endPts(:,2);
  startS = startPts(:,3); endS = endPts(:,3);
  
  hold on;
  lh = quiver(startC, startS, endC-startC, endS-startS, 0, '.');
  hold off
  
  set(lh, 'Color', handles.surfs(mmm).color);
  handles.meshCrossSectionCS(mmm).lh = lh;
  
  if handles.showSurf(mmm) == 0
    set(handles.meshCrossSectionCS(mmm).lh, 'Visible', 'off');
  end
  
end

%% end of function handles = meshCrossSection_ss(handles)
%% ------------------------------------------------------------