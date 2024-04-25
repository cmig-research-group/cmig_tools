

%% ------------------------------------------------------------
%% --- Compute meshCrossSection for col = constant
%% ------------------------------------------------------------
function handles = meshCrossSection_cc(handles)

for mmm = 1:handles.numSurfs
  
  %% Draw in axis SR, which draws imgs(:,j,:)
  axes(handles.axesSR);
  [startPts, endPts] = findContour(handles.surfs(mmm), ...
    handles.surfs(mmm).vertices(:,2), ...
    handles.cc);
  
  startR = startPts(:,1); endR = endPts(:,1);
  startS = startPts(:,3); endS = endPts(:,3);
  
  hold on;
  lh = quiver(startS, startR, endS-startS, endR-startR, 0, '.');
  hold off
  
  set(lh, 'Color', handles.surfs(mmm).color);
  handles.meshCrossSectionSR(mmm).lh = lh;
  
  if handles.showSurf(mmm) == 0
    set(handles.meshCrossSectionSR(mmm).lh, 'Visible', 'off');
  end
  
end

%% end of function handles = meshCrossSection_ss(handles)
%% ------------------------------------------------------------

