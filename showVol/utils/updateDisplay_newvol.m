function handles = updateDisplay_newvol(handles,doNew)
% this has become the common drawing routine, called not only for new
% volume (upon switching orientations, which creates new images), but also for newDisplay (doNew=true:
% initialize image and other objects),
% as well as updateDisplay_rr, _cc, _ss.
% updateDisplay_zoom and _clims are unique

if nargin < 2
  doNew = false;
end

frontFig = gcf;
set(0,'CurrentFigure',handles.figure)

sf = calculateScale(handles);
V = handles.vols{handles.currentVol};

isFOD = isfield(V,'imgs1');
if ~isFOD
  addZoomedInset(handles,false);
  shift = [0 0];
else
  shift = [-1 -1]; %needed to bring center of FOD in line with cursor, saggital is slightly off in 'S' direction, however
  
  %needed map R C S [imgs 1, 2, 3 = CS, SR. CR = axial, sagittal, coronal] to proper slice (test on [113, 100, 141], a single beautifully isolated FOD)
  % also verified at edges [100,100,89->90], and [100, 102->103, 90], and [30,97,130]. The shapes and presence of FOD matches perfectly.
  FODSliceFudge = [1 0 1]; 
  %however, the following points have a FOD present/absent in sagittal and not the others
  %[100, 103, 90] or [99,105,96] [108,103,96]. [112, 99, 144]
  % it seems to be caused by a left-right flip in the sagittal images! 
end