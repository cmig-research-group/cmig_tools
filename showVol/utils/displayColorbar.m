
%% ------------------------------------------------------------

function displayColorbar(handles)
V = handles.vols{handles.currentVol};
hasColorbar = isfield(V,'colormap');

if ~hasColorbar
  set(findobj(handles.colorbar),'Visible','off')
else
  if isfield(V,'imgs') && size(V.imgs,4)==1
    V.limits = [handles.CLim(handles.currentVol,1) handles.CLim(handles.currentVol,2)];
  end
  drawColorbar(handles.colorbar, V.colormap, V.limits)
end