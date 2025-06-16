%% ------------------------------------------------------------
%% --- Update the display to a new volume (change in currentVol)
%% ------------------------------------------------------------
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

%% Draw in axesCR, which draws imgs(:,:,k) % Coronal in std view
set(handles.figure,'CurrentAxes',handles.axesCR)
if isfield(V,'imgs')
  if size(V.imgs,4)==3
    ima = squeeze(V.imgs(:,:,scaleCoord(handles.ss, sf.s),:));
    [ima, outline] = anatomyAddRoiOverlay(handles, ima, 1:size(ima,1), 1:size(ima,2), scaleCoord(handles.ss, sf.s), 'fill');
    if doNew,  handles.imCR = image(ima); end
  elseif size(V.imgs,4)==1
    ima = squeeze(V.imgs(:,:,scaleCoord(handles.ss, sf.s)));
    [ima, outline] = anatomyAddRoiOverlay(handles, ima, 1:size(ima,1), 1:size(ima,2), scaleCoord(handles.ss, sf.s), 'fill');
    if doNew,  handles.imCR = imagesc(ima); colormap(handles.cmap); end
  end
elseif isFOD
  sliceFactor = handles.rcsScale(3) / V.sliceSpacing;
  ima = V.imgs3{scaleCoord(handles.ss,sliceFactor)+FODSliceFudge(3)};
  if V.compressed
    ima = expandFOD(ima);
  end
  [ima, handles] = addZoomedInset(handles, ima, handles.rr, handles.cc, sf.scale, 'lr');
  [ima, outline] = anatomyAddRoiOverlay(handles, ima, 1:size(ima,1), 1:size(ima,2), scaleCoord(handles.ss, sf.s), 'outline');
  %ima = imresize(ima,0.5);
  if doNew
    %shift = floor(sf.scale/2); %shift FOD images by 1/2 mm towards origin to properly center on 1mm grid
    handles.imCR = image(ima);
  end
end
if ~doNew
  set(handles.imCR, 'CData', ima)
  if isfield(V,'imgs') && size(V.imgs,4)==1
    set(handles.imCR, 'CDataMapping','scaled')
  end
end
handles = anatomyDrawRoiOutline(handles, outline, 'CR');
handles.imCR.XData = [1 size(ima,2)/sf.c] + shift(1);
handles.imCR.YData = [1 size(ima,1)/sf.r] + shift(2);
%axis([0.5 handles.ccMax+0.5 0.5 handles.rrMax+0.5]) %unneded--is set in updateDisplay_zoom
axis(handles.axesCR, 'image');
set(handles.axesCR, 'DataAspectRatio', [handles.voxSize(1) handles.voxSize(2) 1]);

if doNew
  handles.CtextCR = text(handles.cc, 1,          'C (R -> L)');
  handles.RtextCR = text(1,          handles.rr, 'R');
  handles.ClineCR = line([handles.cc handles.cc],  [0.5  0.5+handles.rrMax]);
  handles.RlineCR = line([0.5  0.5+handles.ccMax], [handles.rr handles.rr]);
  handles.roiContourCR = [];
else
  handles.CtextCR.Position = [handles.cc 1];
  handles.RtextCR.Position = [1          handles.rr];
  handles.ClineCR.XData = [handles.cc handles.cc];   handles.ClineCR.YData = [0.5  0.5+handles.rrMax];
  handles.RlineCR.XData = [0.5  0.5+handles.ccMax];  handles.RlineCR.YData = [handles.rr handles.rr];
end

%% Draw in axis SR, which draws imgs(:,j,:) % Sagittal in std view
set(handles.figure,'CurrentAxes',handles.axesSR)
if isfield(V,'imgs')
  if size(V.imgs,4)==3
    ima = squeeze(V.imgs(:,scaleCoord(handles.cc,sf.c),:,:));
    [ima, outline] = anatomyAddRoiOverlay(handles, ima, 1:size(ima,1), scaleCoord(handles.cc,sf.c), 1:size(ima,2), 'fill');
    if doNew, handles.imSR = image(ima); end
  else
    ima = squeeze(V.imgs(:,scaleCoord(handles.cc,sf.c),:));
    [ima, outline] = anatomyAddRoiOverlay(handles, ima, 1:size(ima,1), scaleCoord(handles.cc,sf.c), 1:size(ima,2), 'fill');
    if doNew, handles.imSR = imagesc(ima); end
  end
elseif isFOD
  sliceFactor = handles.rcsScale(2) / V.sliceSpacing;
  cIdx = handles.ccMax - scaleCoord(handles.cc,sliceFactor)+FODSliceFudge(2);
  cIdx = max(1,cIdx); %catch any 0 index
  ima = V.imgs2{cIdx}; %FIXME: swap left and right in FOD series
  if V.compressed
    ima = expandFOD(ima);
  end
  [ima, handles] = addZoomedInset(handles, ima, handles.rr, handles.ss, sf.scale, 'lr');
  [ima, outline] = anatomyAddRoiOverlay(handles, ima, 1:size(ima,1), scaleCoord(handles.cc,sf.c), 1:size(ima,2), 'outline');
  if doNew
    handles.imSR = image(ima);
  end
end
if ~doNew
  set(handles.imSR, 'CData', ima)
  if isfield(V,'imgs') && size(V.imgs,4)==1
    set(handles.imSR, 'CDataMapping','scaled')
  end
end
handles = anatomyDrawRoiOutline(handles, outline, 'SR');
handles.imSR.XData = [1 size(ima,2)/sf.s] + shift(1);
handles.imSR.YData = [1 size(ima,1)/sf.r] + shift(2);
%axis([0.5 handles.ssMax+0.5 0.5 handles.rrMax+0.5])
axis(handles.axesSR, 'image');
set(handles.axesSR, 'DataAspectRatio', [handles.voxSize(1) handles.voxSize(3) 1]);

if doNew
  handles.StextSR = text(handles.ss, 1,          'S');
  handles.RtextSR = text(1,          handles.rr, 'R');
  handles.SlineSR = line([handles.ss handles.ss], [0.5 0.5+handles.rrMax]);
  handles.RlineSR = line([0.5 0.5+handles.ssMax], [handles.rr handles.rr]);
  handles.roiContourSR = [];
else
  handles.StextSR.Position = [handles.ss 1];
  handles.RtextSR.Position = [1          handles.rr];
  handles.SlineSR.XData = [handles.ss handles.ss]; handles.SlineSR.YData = [0.5 0.5+handles.rrMax];
  handles.RlineSR.XData = [0.5 0.5+handles.ssMax]; handles.RlineSR.YData = [handles.rr handles.rr];
end

%% Draw in axis CS, which draws imgs(i,:,:) % Axial in std view
set(handles.figure,'CurrentAxes',handles.axesCS)
if isfield(V,'imgs')
  if size(V.imgs,4)==3
    ima = squeeze(V.imgs(scaleCoord(handles.rr,sf.r),:,:,:));
    ima = permute(ima,[2 1 3]); %why just this one?
    [ima, outline] = anatomyAddRoiOverlay(handles, ima, scaleCoord(handles.rr,sf.r), 1:size(ima,2), 1:size(ima,1),'fill');
    if doNew, handles.imCS = image(ima); end
  else
    ima = squeeze(V.imgs(scaleCoord(handles.rr,sf.r),:,:))';
    [ima, outline] = anatomyAddRoiOverlay(handles, ima, scaleCoord(handles.rr,sf.r), 1:size(ima,2), 1:size(ima,1),'fill');
    if doNew, handles.imCS = imagesc(ima); end
  end
elseif isFOD
  sliceFactor = handles.rcsScale(1) / V.sliceSpacing;
  ima = V.imgs1{scaleCoord(handles.rr,sliceFactor)+FODSliceFudge(1)};
  if V.compressed
    ima = expandFOD(ima);
  end
  [ima, handles] = addZoomedInset(handles, ima, handles.ss, handles.cc, sf.scale, 'ur');
  [ima, outline] = anatomyAddRoiOverlay(handles, ima, scaleCoord(handles.rr,sf.r), 1:size(ima,2), 1:size(ima,1),'outline');
  if doNew
    handles.imCS = image(ima);
  end
end
if ~doNew
  set(handles.imCS, 'CData', ima)
  if isfield(V,'imgs') && size(V.imgs,4)==1
    set(handles.imCS, 'CDataMapping','scaled')
  end
end
handles = anatomyDrawRoiOutline(handles, outline, 'CS');
handles.imCS.XData = [1 size(ima,2)/sf.c] + shift(1);
handles.imCS.YData = [1 size(ima,1)/sf.s] + shift(2);
%axis([0.5 handles.ccMax+0.5 0.5 handles.ssMax+0.5])
axis(handles.axesCS, 'image');
axis(handles.axesCS, 'xy'); %FIXME: axis xy flips up/down anteror/posterior - it's needed to put anterior at top of axis
set(handles.axesCS, 'DataAspectRatio', [handles.voxSize(3) handles.voxSize(2) 1]);

if doNew
  handles.CtextCS = text(handles.cc, handles.ssMax, 'C (R -> L)');
  handles.StextCS = text(1,          handles.ss,    'S');
  handles.ClineCS = line([handles.cc handles.cc], [0.5 0.5+handles.ssMax]);
  handles.SlineCS = line([0.5 0.5+handles.ccMax], [handles.ss handles.ss]);
  handles.roiContourCS = [];
else
  handles.CtextCS.Position = [sf.c*handles.cc sf.s*handles.ssMax];
  handles.StextCS.Position = [sf.c*1          sf.s*handles.ss];
  handles.ClineCS.XData = [handles.cc handles.cc]; handles.ClineCS.YData = [0.5 0.5+handles.ssMax];
  handles.SlineCS.XData = [0.5 0.5+handles.ccMax]; handles.SlineCS.YData = [handles.ss handles.ss];
end

%% colormap (for indexed FOD)
if isfield(V,'colormap')
  colormap(V.colormap(:,1:3))
else
  colormap(gray(128))
end

%show colorbar for grayscale or stat maps
displayColorbar(handles)

%% For all axes
axesHandles = [handles.axesCR handles.axesSR handles.axesCS];
set(axesHandles, 'Color', 'k');
set(axesHandles, 'CLim', handles.CLim(handles.currentVol, :));

textHandles = [handles.CtextCR handles.RtextCR ...
  handles.StextSR handles.RtextSR ...
  handles.CtextCS handles.StextCS];
set(textHandles, 'VerticalAlignment', 'top');
set(textHandles, 'HorizontalAlignment', 'left');
set(textHandles, 'Color', 'r');

lineHandles = [handles.ClineCR handles.RlineCR ...
  handles.SlineSR handles.RlineSR ...
  handles.ClineCS handles.SlineCS];
set(lineHandles, 'Color', [0.5 0.5 0], 'LineStyle', ':','LineWidth',handles.LineWidth);

%% Draw mesh cross sections TODO: handle scaling
handles = meshCrossSection_ss(handles);
handles = meshCrossSection_cc(handles);
handles = meshCrossSection_rr(handles);

printInfo(handles);

handles = updateDisplay_zoom(handles);
guidata(handles.figure, handles)

%restore
set(0,'CurrentFigure',frontFig)