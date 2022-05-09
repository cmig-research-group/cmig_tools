function varargout = showVol(varargin)
% Viewer for 3-D images (voxel volumes). 
%
% showVol(vol_1, [vol_2], ..., [vol_N], ...
%         [initialLPH])
%
% NOTE 1: vol_* is a struct with (at least) the fields 
%        'imgs' and 'Mvxl2lph' 
%
% NOTE 2: initialLPH is a vector of length 3 that specifies the 
%         LPH coordinates of the starting point of the crosshair
%
% Commands that can be entered in the "Command" edit box:
% ------------------------------------------------------
% 1) To do a flythrough in the main image:
%    fly stopIndex [pauseTime]
%     -- stopIndex: the row, column or slice to fly to in the main
%        (large) image. 
%     -- pauseTime: the time spent on each image in the flythrough. 
%        (The default is 0.)
%
% 2) To cycle through the input volumes as an animation:
%    cycle [numLoops] [pauseTime_1] ... [pauseTime_N]
%     -- numLoops: the number of loops in the animation.
%        (The default is 10.)
%     -- pauseTime_n: the time spent on each image in the animation.
%        (The default is 0.1. If negative, then vol_n is not shown.)
% 
% 3) To go to a particular point in LPH coordinates
%    lph L_coord P_coord H_coord
%
% 4) To go to a particular point in RCS coordinates
%    lph R_coord C_coord S_coord
%
% List of keyboard shortcuts:
% --------------------------
% NOTE: sometimes the UI pushbuttons, sliders, etc. do not release
% their focus (control), so you have click on the axes or on the
% figure background in order to use the keyboard shortcuts.
%
% -- 'v': cycle to the next volume
% -- 'c': cycle between volumes as an animation
% -- 'o': cycle orientations
% -- 't': toggle crosshair display in the main image
% -- 'i': move crosshair up (w.r.t. the main image)
% -- 'j': move crosshair left (w.r.t. the main image)
% -- 'k': move crosshair right (w.r.t. the main image)
% -- 'm': move crosshair down (w.r.t. the main image)
% -- ',': move crosshair "out of the screen" (w.r.t. the main image)
% -- '.': move crosshair "into the screen" (w.r.t. the main image)
% -- 'p': store the current voxel/point
% -- ';': move to the stored voxel
% -- '/': swap the current voxel with the stored voxel
% -- 'q': constrast down
% -- 'w': contrast up
% -- 'a': darken down
% -- 's': darken up
% -- 'x': zoom out
% -- 'z': zoom in

% (mukundb)

%% ============================================================
%% 
%%               SECTION 1: INITIALIZATION
%% 
%% ============================================================


%% ------------------------------------------------------------ 
%% Begin initialization code - DO NOT EDIT
%% ------------------------------------------------------------ 
gui_Singleton = 0;
gui_State = struct('gui_Name',       mfilename, ...
                   'gui_Singleton',  gui_Singleton, ...
                   'gui_OpeningFcn', @showVol_OpeningFcn, ...
                   'gui_OutputFcn',  @showVol_OutputFcn, ...
                   'gui_LayoutFcn',  [] , ...
                   'gui_Callback',   []);
if nargin & isstr(varargin{1})
    gui_State.gui_Callback = str2func(varargin{1});
end

if nargout
    [varargout{1:nargout}] = gui_mainfcn(gui_State, varargin{:});
else
    gui_mainfcn(gui_State, varargin{:});
end
%% End initialization code - DO NOT EDIT
%% ------------------------------------------------------------ 


%% ------------------------------------------------------------ 
%% --- Executes just before showVol is made visible.
%% ------------------------------------------------------------ 
function showVol_OpeningFcn(hObject, eventdata, handles, varargin)

% Initialize figure
if strcmp(version('-release'), '14')
  cameratoolbar('Close'); % this takes 1 second!
else
  cameratoolbar('hide');
end
colormap(gray);
set(handles.figure, 'Renderer', 'zbuffer');

if isempty(varargin)
  fprintf('--------------------------------\n')
  fprintf('%s: no inputs given\n', mfilename);
  fprintf('--------------------------------\n')
  error('Exiting ...');
end

% Read in inputs
handles.numVols = 0;
initialLPH = [];
for iii = 1:nargin-3

  if (isfield(varargin{iii}, 'imgs')) && (isfield(varargin{iii}, 'Mvxl2lph'))
    % volStruct
    handles.numVols = handles.numVols + 1;
    handles.vols{handles.numVols} = varargin{iii};
  elseif length(varargin{iii}) == 3
    % initial LPH coords
    initialLPH = varargin{iii};
    initialLPH = initialLPH(:);
    initialLPH(4) = 1;
  end
    
end

% Check that we have at least one vol
if handles.numVols == 0
  fprintf('--------------------------------\n')
  fprintf('%s: must have at least one vol struct input\n', ...
	  mfilename);
  fprintf('--------------------------------\n')
  error('Exiting ...');
end

% Load up from the first vol
handles.Mvxl2lph = handles.vols{1}.Mvxl2lph;

% Error checking I: all vols must have the same voxel number
voxNumCheck = [];
for vvv = 1:handles.numVols
  voxNumCheck = [voxNumCheck; size(handles.vols{vvv}.imgs)];
end
if (handles.numVols == 1) || (sum(range(voxNumCheck)) < eps)
  handles.rrMax = voxNumCheck(1,1);
  handles.ccMax = voxNumCheck(1,2);
  handles.ssMax = voxNumCheck(1,3);
else
  fprintf('--------------------------------\n')
  fprintf('%s: all the vols must have the same voxel number\n', ...
	  mfilename);
  fprintf('--------------------------------\n')
  error('Exiting ...');
end

% Error checking II: all vols must have the same voxel size
voxSizeCheck = [];
for vvv = 1:handles.numVols
  voxSizeCheck = [voxSizeCheck; 
		  [norm(handles.vols{vvv}.Mvxl2lph(:,1)) ...
		   norm(handles.vols{vvv}.Mvxl2lph(:,2)) ...
		   norm(handles.vols{vvv}.Mvxl2lph(:,3))];
		 ];
end
if (handles.numVols == 1) || (sum(range(voxSizeCheck)) < 1e-3)
  handles.voxSize = voxSizeCheck(1,:);
  if ( (handles.voxSize(1) < eps) | ...
       (handles.voxSize(2) < eps) | ...
       (handles.voxSize(3) < eps) )
    fprintf('--------------------------------\n')
    fprintf('%s: bad info in Mvxl2lph\n', mfilename);
    fprintf('--------------------------------\n')
    error('Exiting ...');
  end
else
  fprintf('--------------------------------\n')
  fprintf('%s: all the vols must have the same voxel size\n', ...
	  mfilename);
  fprintf('--------------------------------\n')
  error('Exiting ...');
end

% Contrast info
maxI = zeros(handles.numVols,1);
minI = zeros(handles.numVols,1);
for vvv = 1:handles.numVols
  if (isfield(handles.vols{vvv}, 'minI') & ...
      isfield(handles.vols{vvv}, 'maxI'))
    maxI(vvv) = handles.vols{vvv}.maxI;
    minI(vvv) = handles.vols{vvv}.minI;
  else
    if exist('maxmin')
      [maxI(vvv), minI(vvv)] = maxmin(handles.vols{vvv}.imgs);
    else
      maxI(vvv) = max(max(max(handles.vols{vvv}.imgs)));
      minI(vvv) = min(min(min(handles.vols{vvv}.imgs)));
    end
  end
end
handles.maxI = maxI;
handles.minI = minI;
handles.CLim = [minI maxI+eps];

handles.cVal = zeros(handles.numVols,1);
set(handles.contrast_slider, 'Value', 0);
set(handles.contrast_slider,   'Min', 0);
set(handles.contrast_slider,   'Max', 1);

handles.dVal = zeros(handles.numVols,1);
set(handles.darken_slider, 'Value', 0);
set(handles.darken_slider,   'Min', 0);
set(handles.darken_slider,   'Max', 1);

% Starting points
handles.rr = ceil(handles.rrMax/2 + sqrt(eps));
handles.cc = ceil(handles.ccMax/2 + sqrt(eps));
handles.ss = ceil(handles.ssMax/2 + sqrt(eps));
handles.rrStored = handles.rr;
handles.ccStored = handles.cc;
handles.ssStored = handles.ss;

str = sprintf('Store this voxel -- keyboard shortcut: p');
set(handles.storePoint,'TooltipString', str);
str = sprintf('Go to voxel (%d,%d,%d) -- keyboard shortcut: ;', ...
	      handles.rrStored, handles.ccStored, handles.ssStored);
set(handles.gotoPoint,'TooltipString', str);
str = sprintf('Swap current voxel with voxel (%d,%d,%d) -- keyboard shortcut: /', ...
	      handles.rrStored, handles.ccStored, handles.ssStored);
set(handles.swapPoint,'TooltipString', str);

% Zoom info
handles.zoom = 0;
set(handles.zoom_slider, 'Value', 0);
set(handles.zoom_slider, 'Min', 0);
set(handles.zoom_slider, 'Max', 6);

tmp_rr = (handles.rr - 0.5) * handles.voxSize(1);
tmp_cc = (handles.cc - 0.5) * handles.voxSize(2);
tmp_ss = (handles.ss - 0.5) * handles.voxSize(3);
handles.halfRange_mm = max(max(tmp_rr, tmp_cc), tmp_ss);

tmp = [-handles.halfRange_mm +handles.halfRange_mm];
handles.rrLimOrig = handles.rr + tmp/handles.voxSize(1);
handles.ccLimOrig = handles.cc + tmp/handles.voxSize(2);
handles.ssLimOrig = handles.ss + tmp/handles.voxSize(3);

handles.rrZoom = handles.rr;
handles.ccZoom = handles.cc;
handles.ssZoom = handles.ss;

%% ------------------------------------------------------------ 
% Initial display 

str = 'showVol(''volPB_Callback'',gcbo,[],guidata(gcbo))';
for vvv = 1:handles.numVols
  handles.volPB(vvv) = uicontrol('Tag', 'volButton', ...
				 'Style', 'pushbutton', ...
				 'Unit', 'characters', ...
				 'Position', [10 10 10 10], ...
				 'TooltipString', 'Show this volume', ...
				 'String', ['vol' num2str(vvv)]);
  set(handles.volPB(vvv), 'Callback', str);
end
handles.hideVol = zeros(handles.numVols, 1);
handles.currentVol = 1;
handles.currentVolPB = handles.volPB(1);
set(handles.currentVolPB, 'ForegroundColor', 'r');
set(handles.currentVolPB, 'TooltipString', 'Hide this volume');
handles.ORIENTATION = 1;
handles.HIDE1 = 0;
handles.HIDE2 = 0;
handles.HIDE3 = 0;
handles = displayNewOrientation(handles);
printInfo(handles);
handles.output = hObject;
guidata(hObject, handles);
figure_ResizeFcn(hObject, eventdata, handles);
pos = [0.5 4.5 152 59];
set(handles.figure, 'Position', pos);

if ~isempty(initialLPH)
  handles = gotoPointLPH(handles, initialLPH);
  printInfo(handles);
  guidata(hObject, handles);
end

%% end of function showVol_OpeningFcn(hObject, eventdata, handles, varargin)
%% ------------------------------------------------------------ 


%% ------------------------------------------------------------ 
%% --- Outputs from this function are returned to the command line.
%% ------------------------------------------------------------ 
function varargout = showVol_OutputFcn(hObject, eventdata, handles)
varargout{1} = handles.output;
%% ------------------------------------------------------------ 


%% ============================================================
%% 
%%               SECTION 2: BASIC FUNCTIONS
%% 
%% ============================================================


%% ------------------------------------------------------------ 
%% --- Print info into text boxes
%% ------------------------------------------------------------ 
function printInfo(handles)

rr = handles.rr; cc = handles.cc; ss = handles.ss;
voxVal = handles.vols{handles.currentVol}.imgs(rr,cc,ss);
CLim1 = handles.CLim(handles.currentVol, 1);
CLim2 = handles.CLim(handles.currentVol, 2);

% Convert RCS to LPH
tmp = [rr; cc; ss; 1];
tmp = handles.Mvxl2lph * tmp;
ll = tmp(1); pp = tmp(2); hh = tmp(3);

set(handles.edit_rr, 'String', sprintf('%d', rr));
set(handles.edit_cc, 'String', sprintf('%d', cc));
set(handles.edit_ss, 'String', sprintf('%d', ss));

set(handles.text_ll, 'String', sprintf('%.3f', ll));
set(handles.text_pp, 'String', sprintf('%.3f', pp));
set(handles.text_hh, 'String', sprintf('%.3f', hh));

set(handles.text_voxval, 'String', sprintf('%.2f', voxVal));

scrVal = (voxVal - CLim1)/(CLim2 - CLim1);
scrVal = min(scrVal, 1);
scrVal = max(scrVal, 0);
set(handles.text_scrval, 'String', sprintf('%.2f', 100*scrVal));

%% end of function printInfo(handles)
%% ------------------------------------------------------------ 


%% ------------------------------------------------------------ 
%% --- Display the images and lines from scratch
%% ------------------------------------------------------------ 
function handles = newDisplay(handles)

%% Draw in axesCR, which draws imgs(:,:,k)
axes(handles.axesCR);
ima = handles.vols{handles.currentVol}.imgs(:,:,handles.ss);
handles.imCR = imagesc(ima); 
axis image
set(handles.axesCR, 'DataAspectRatio', [handles.voxSize(1) handles.voxSize(2) 1]);
handles.CtextCR = text(handles.cc, 1, 'C');
handles.ClineCR = line([handles.cc handles.cc],[0.5 0.5+handles.rrMax]);
handles.RtextCR = text(1, handles.rr, 'R');
handles.RlineCR = line([0.5 0.5+handles.ccMax],[handles.rr handles.rr]);

%% Draw in axis SR, which draws imgs(:,j,:)
axes(handles.axesSR);
ima = squeeze(handles.vols{handles.currentVol}.imgs(:,handles.cc,:));
handles.imSR = imagesc(ima); 
axis image
set(handles.axesSR, 'DataAspectRatio', [handles.voxSize(1) handles.voxSize(3) 1]);
handles.StextSR = text(handles.ss, 1, 'S');
handles.SlineSR = line([handles.ss handles.ss],[0.5 0.5+handles.rrMax]);
handles.RtextSR = text(1, handles.rr, 'R');
handles.RlineSR = line([0.5 0.5+handles.ssMax],[handles.rr handles.rr]);

%% Draw in axis CS, which draws imgs(i,:,:)
axes(handles.axesCS);
ima = squeeze(handles.vols{handles.currentVol}.imgs(handles.rr,:,:))';
handles.imCS = imagesc(ima); 
axis image; axis xy
set(handles.axesCS, 'DataAspectRatio', [handles.voxSize(3) handles.voxSize(2) 1]);
handles.CtextCS = text(handles.cc, handles.ssMax, 'C');
handles.ClineCS = line([handles.cc handles.cc],[0.5 0.5+handles.ssMax]);
handles.StextCS = text(1, handles.ss, 'S');
handles.SlineCS = line([0.5 0.5+handles.ccMax],[handles.ss handles.ss]);

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
set(lineHandles, 'Color', 'r');
set(lineHandles, 'LineStyle', ':');

%% end of function handles = newDisplay(handles)
%% ------------------------------------------------------------ 


%% ------------------------------------------------------------ 
%% --- Update display upon change in slice number
%% ------------------------------------------------------------ 
function handles = updateDisplay_ss(handles)
ima = handles.vols{handles.currentVol}.imgs(:,:,handles.ss);
set(handles.imCR, 'CData', ima);
set(handles.StextSR, 'Position', [handles.ss handles.rrLim(1)]);
set(handles.SlineSR, 'XData', [handles.ss handles.ss]);
set(handles.StextCS, 'Position', [handles.ccLim(1) handles.ss]);
set(handles.SlineCS, 'YData', [handles.ss handles.ss]);
%% ------------------------------------------------------------ 


%% ------------------------------------------------------------ 
%% --- Update display upon change in column number
%% ------------------------------------------------------------ 
function handles = updateDisplay_cc(handles)
ima = squeeze(handles.vols{handles.currentVol}.imgs(:,handles.cc,:));
set(handles.imSR, 'CData', ima);
set(handles.CtextCR, 'Position', [handles.cc handles.rrLim(1)]);
set(handles.ClineCR, 'XData', [handles.cc handles.cc]);
set(handles.CtextCS, 'Position', [handles.cc handles.ssLim(2)]);
set(handles.ClineCS, 'XData', [handles.cc handles.cc]);
%% ------------------------------------------------------------ 


%% ------------------------------------------------------------ 
%% --- Update display upon change in row number
%% ------------------------------------------------------------ 
function handles = updateDisplay_rr(handles)
ima = squeeze(handles.vols{handles.currentVol}.imgs(handles.rr,:,:))';
set(handles.imCS, 'CData', ima);
set(handles.RtextCR, 'Position', [handles.ccLim(1) handles.rr]);
set(handles.RlineCR, 'YData', [handles.rr handles.rr]);
set(handles.RtextSR, 'Position', [handles.ssLim(1) handles.rr]);
set(handles.RlineSR, 'YData', [handles.rr handles.rr]);
%% ------------------------------------------------------------ 


%% ------------------------------------------------------------ 
%% --- Compute new zoom values and update display
%% ------------------------------------------------------------ 
function handles = updateDisplay_zoom(handles)

if handles.zoom < sqrt(eps)

  handles.rrLim = handles.rrLimOrig;
  handles.ccLim = handles.ccLimOrig;
  handles.ssLim = handles.ssLimOrig;
  
else
  
  zoom = 2^handles.zoom;
  halfRange_mm = handles.halfRange_mm/zoom;
  tmp = [-halfRange_mm +halfRange_mm];
  handles.rrLim = handles.rrZoom + tmp/handles.voxSize(1);
  handles.ccLim = handles.ccZoom + tmp/handles.voxSize(2);
  handles.ssLim = handles.ssZoom + tmp/handles.voxSize(3);

end

set(handles.axesCR, 'Xlim', handles.ccLim);
set(handles.axesCR, 'Ylim', handles.rrLim);
set(handles.axesSR, 'Xlim', handles.ssLim);
set(handles.axesSR, 'Ylim', handles.rrLim);
set(handles.axesCS, 'Xlim', handles.ccLim);
set(handles.axesCS, 'Ylim', handles.ssLim);

set(handles.CtextCR, 'Position', [handles.cc handles.rrLim(1)]);
set(handles.RtextCR, 'Position', [handles.ccLim(1) handles.rr]);
set(handles.StextSR, 'Position', [handles.ss handles.rrLim(1)]);
set(handles.RtextSR, 'Position', [handles.ssLim(1) handles.rr]);
set(handles.CtextCS, 'Position', [handles.cc handles.ssLim(2)]);
set(handles.StextCS, 'Position', [handles.ccLim(1) handles.ss]);

%% end of function handles = updateDisplay_zoom(handles)
%% ------------------------------------------------------------ 


%% ------------------------------------------------------------ 
%% --- Compute new CLim values and update display
%% ------------------------------------------------------------ 
function handles = updateDisplay_clims(handles)

cVal = handles.cVal(handles.currentVol);
dVal = handles.dVal(handles.currentVol);
minI = handles.minI(handles.currentVol);
maxI = handles.maxI(handles.currentVol);

X = cot((pi/4)*(cVal+1));
CLim1_scaled = dVal*(1-X);
CLim2_scaled = CLim1_scaled + X;
CLim1 = minI + CLim1_scaled*(maxI-minI);
CLim2 = minI + CLim2_scaled*(maxI-minI) + eps;

set(handles.axes1, 'CLim', [CLim1 CLim2]);
set(handles.axes2, 'CLim', [CLim1 CLim2]);
set(handles.axes3, 'CLim', [CLim1 CLim2]);

handles.CLim(handles.currentVol, :) = [CLim1 CLim2];

rr = handles.rr; cc = handles.cc; ss = handles.ss;
voxVal = handles.vols{handles.currentVol}.imgs(rr,cc,ss);
scrVal = (voxVal - CLim1)/(CLim2 - CLim1);
scrVal = min(scrVal, 1);
scrVal = max(scrVal, 0);
set(handles.text_scrval, 'String', sprintf('%.2f', 100*scrVal));

%% end of function handles = updateDisplay_clims(handles)
%% ------------------------------------------------------------ 


%% ------------------------------------------------------------ 
%% --- Update the display to a new volume (change in currentVol)
%% ------------------------------------------------------------ 
function updateDisplay_newvol(handles)

set(handles.contrast_slider, 'Value', ...
		  handles.cVal(handles.currentVol));
set(handles.darken_slider, 'Value', ...
		  handles.dVal(handles.currentVol));
set(handles.axes1, 'CLim', handles.CLim(handles.currentVol,:));
set(handles.axes2, 'CLim', handles.CLim(handles.currentVol,:));
set(handles.axes3, 'CLim', handles.CLim(handles.currentVol,:));

ima = handles.vols{handles.currentVol}.imgs(:,:,handles.ss);
set(handles.imCR, 'CData', ima);
ima = squeeze(handles.vols{handles.currentVol}.imgs(:,handles.cc,:));
set(handles.imSR, 'CData', ima);
ima = squeeze(handles.vols{handles.currentVol}.imgs(handles.rr,:,:))';
set(handles.imCS, 'CData', ima);

printInfo(handles); 

%% ------------------------------------------------------------ 


%% ------------------------------------------------------------ 
%% --- Toggle the display of crosshair in axes1
%% ------------------------------------------------------------ 
function handles = toggleLines1(handles)
handles.HIDE1 = 1 - handles.HIDE1;
set(handles.toggleLines1, 'Value', handles.HIDE1);
if handles.HIDE1
  set(handles.toggleLines1, 'ForegroundColor', 'k');
  set(handles.axes1_X, 'Visible', 'off');
else
  set(handles.toggleLines1, 'ForegroundColor', 'r');
  set(handles.axes1_X, 'Visible', 'on');
end
%% ------------------------------------------------------------ 


%% ------------------------------------------------------------ 
%% --- Toggle the display of crosshair in axes2
%% ------------------------------------------------------------ 
function handles = toggleLines2(handles)
handles.HIDE2 = 1 - handles.HIDE2;
set(handles.toggleLines2, 'Value', handles.HIDE2);
if handles.HIDE2
  set(handles.toggleLines2, 'ForegroundColor', 'k');
  set(handles.axes2_X, 'Visible', 'off');
else
  set(handles.toggleLines2, 'ForegroundColor', 'r');
  set(handles.axes2_X, 'Visible', 'on');
end
%% ------------------------------------------------------------ 


%% ------------------------------------------------------------ 
%% --- Toggle the display of crosshair in axes3
%% ------------------------------------------------------------ 
function handles = toggleLines3(handles)
handles.HIDE3 = 1 - handles.HIDE3;
set(handles.toggleLines3, 'Value', handles.HIDE3);
if handles.HIDE3
  set(handles.toggleLines3, 'ForegroundColor', 'k');
  set(handles.axes3_X, 'Visible', 'off');
else
  set(handles.toggleLines3, 'ForegroundColor', 'r');
  set(handles.axes3_X, 'Visible', 'on');
end
%% ------------------------------------------------------------ 


%% ============================================================
%% 
%%               SECTION 3: WRAPPER FUNCTIONS
%% 
%% ============================================================


%% ------------------------------------------------------------ 
%% --- Change the display upon setting a new orientation
%% ------------------------------------------------------------ 
function handles = displayNewOrientation(handles)

switch handles.ORIENTATION
  
 case 1
  handles.axesCR = handles.axes1;
  handles.axesSR = handles.axes2;
  handles.axesCS = handles.axes3;
  handles = newDisplay(handles);
  handles.axes1_X = [handles.ClineCR handles.RlineCR ...
		     handles.CtextCR handles.RtextCR];
  handles.axes2_X = [handles.SlineSR handles.RlineSR ...
		     handles.StextSR handles.RtextSR];
  handles.axes3_X = [handles.ClineCS handles.SlineCS ...
		     handles.CtextCS handles.StextCS];
  handles.ccDecrChar = 'j';
  handles.ccIncrChar = 'l';
  handles.rrDecrChar = 'i';
  handles.rrIncrChar = 'k';
  handles.ssDecrChar = ',';
  handles.ssIncrChar = '.';
  
 case 2
  handles.axesSR = handles.axes1;
  handles.axesCS = handles.axes2;
  handles.axesCR = handles.axes3;
  handles = newDisplay(handles);
  handles.axes1_X = [handles.SlineSR handles.RlineSR ...
		     handles.StextSR handles.RtextSR];
  handles.axes2_X = [handles.ClineCS handles.SlineCS ...
		     handles.CtextCS handles.StextCS];
  handles.axes3_X = [handles.ClineCR handles.RlineCR ...
		     handles.CtextCR handles.RtextCR];
  handles.ssDecrChar = 'j';
  handles.ssIncrChar = 'l';
  handles.rrDecrChar = 'i';
  handles.rrIncrChar = 'k';
  handles.ccDecrChar = ',';
  handles.ccIncrChar = '.';
  
 case 3
  handles.axesCS = handles.axes1;
  handles.axesCR = handles.axes2;
  handles.axesSR = handles.axes3;
  handles = newDisplay(handles);
  handles.axes1_X = [handles.ClineCS handles.SlineCS ...
		     handles.CtextCS handles.StextCS];
  handles.axes2_X = [handles.ClineCR handles.RlineCR ...
		     handles.CtextCR handles.RtextCR];
  handles.axes3_X = [handles.SlineSR handles.RlineSR ...
		     handles.StextSR handles.RtextSR];
  handles.ccDecrChar = 'j';
  handles.ccIncrChar = 'l';
  handles.ssDecrChar = 'k';
  handles.ssIncrChar = 'i';
  handles.rrDecrChar = ',';
  handles.rrIncrChar = '.';

end

if handles.HIDE1
  set(handles.axes1_X, 'Visible', 'off');
end
if handles.HIDE2
  set(handles.axes2_X, 'Visible', 'off');
end
if handles.HIDE3
  set(handles.axes3_X, 'Visible', 'off');
end

handles = updateDisplay_zoom(handles);

%% end of function handles = displayNewOrientation(handles)
%% ------------------------------------------------------------ 


%% ------------------------------------------------------------ 
%% --- Reverse-cycle the orientations shown in each axis 
%% ------------------------------------------------------------ 
function handles = cycleOrientationsReverse(handles)
handles.ORIENTATION = handles.ORIENTATION - 1;
if handles.ORIENTATION < 1
  handles.ORIENTATION = 3;
end
handles = displayNewOrientation(handles);
%% ------------------------------------------------------------ 


%% ------------------------------------------------------------ 
%% --- Forward-cycle the orientations shown in each axis 
%% ------------------------------------------------------------ 
function handles = cycleOrientationsForward(handles)
handles.ORIENTATION = handles.ORIENTATION + 1;
if handles.ORIENTATION > 3
  handles.ORIENTATION = 1;
end
handles = displayNewOrientation(handles);
%% ------------------------------------------------------------ 


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


%% ------------------------------------------------------------ 
%% -- Go to the stored RCS point
%% ------------------------------------------------------------ 
function handles = gotoStoredPoint(handles)
set(handles.gotoPoint, 'Value', 0);
handles.rr = handles.rrStored;
handles.cc = handles.ccStored;
handles.ss = handles.ssStored;
handles = updateDisplay_rr(handles);
handles = updateDisplay_cc(handles);
handles = updateDisplay_ss(handles);
printInfo(handles);
%% ------------------------------------------------------------ 


%% ------------------------------------------------------------ 
%% -- Swap the current RCS point with the stored RCS point
%% ------------------------------------------------------------ 
function handles = swapCurrentAndStoredPoints(handles)
set(handles.swapPoint, 'Value', 0);
tmp = [handles.rr handles.cc handles.ss];
handles.rr = handles.rrStored;
handles.cc = handles.ccStored;
handles.ss = handles.ssStored;
handles.rrStored = tmp(1);
handles.ccStored = tmp(2);
handles.ssStored = tmp(3);
str = sprintf('Go to voxel (%d,%d,%d) -- keyboard shortcut: ;', ...
	      handles.rrStored, handles.ccStored, handles.ssStored);
set(handles.gotoPoint,'TooltipString', str);
str = sprintf('Swap current voxel with voxel (%d,%d,%d) -- keyboard shortcut: /', ...
	      handles.rrStored, handles.ccStored, handles.ssStored);
set(handles.swapPoint,'TooltipString', str);
handles = updateDisplay_rr(handles);
handles = updateDisplay_cc(handles);
handles = updateDisplay_ss(handles);
printInfo(handles);
%% ------------------------------------------------------------ 


%% ------------------------------------------------------------ 
%% -- Go to the given RCS point
%% ------------------------------------------------------------ 
function handles = gotoPointRCS(handles, targetRCS)
rr = round(targetRCS(1));
cc = round(targetRCS(2));
ss = round(targetRCS(3));
if ((rr >= 1) & (rr <= handles.rrMax) & ...
    (cc >= 1) & (cc <= handles.ccMax) & ...
    (ss >= 1) & (ss <= handles.ssMax))
  handles.rr = rr;
  handles.cc = cc;
  handles.ss = ss;
  handles = updateDisplay_rr(handles);
  handles = updateDisplay_cc(handles);
  handles = updateDisplay_ss(handles);
  printInfo(handles);
else
  set(handles.edit_command, 'String', ...
		    'Outside the volume!');
end
%% ------------------------------------------------------------ 


%% ------------------------------------------------------------ 
%% -- Go to the given LPH point
%% ------------------------------------------------------------ 
function handles = gotoPointLPH(handles, targetLPH)
Mlph2vxl = inv(handles.Mvxl2lph);
targetRCS = Mlph2vxl * targetLPH;
rr = round(targetRCS(1));
cc = round(targetRCS(2));
ss = round(targetRCS(3));
if ((rr >= 1) & (rr <= handles.rrMax) & ...
    (cc >= 1) & (cc <= handles.ccMax) & ...
    (ss >= 1) & (ss <= handles.ssMax))
  handles.rr = rr;
  handles.cc = cc;
  handles.ss = ss;
  handles = updateDisplay_rr(handles);
  handles = updateDisplay_cc(handles);
  handles = updateDisplay_ss(handles);
  printInfo(handles);
else
  set(handles.edit_command, 'String', ...
		    'Outside the volume!');
end
%% ------------------------------------------------------------ 


%% ------------------------------------------------------------ 
%% --- Decrement row number 
%% ------------------------------------------------------------ 
function handles = decrementRowNumber(handles)
if handles.rr > 1
  handles.rr = handles.rr - 1;
end
handles = updateDisplay_rr(handles);
printInfo(handles);
%% ------------------------------------------------------------ 


%% ------------------------------------------------------------ 
%% --- Increment row number 
%% ------------------------------------------------------------ 
function handles = incrementRowNumber(handles)
if handles.rr < handles.rrMax
  handles.rr = handles.rr + 1;
end
handles = updateDisplay_rr(handles);
printInfo(handles);
%% ------------------------------------------------------------ 


%% ------------------------------------------------------------ 
%% --- Decrement column number 
%% ------------------------------------------------------------ 
function handles = decrementColumnNumber(handles)
if handles.cc > 1
  handles.cc = handles.cc - 1;
end
handles = updateDisplay_cc(handles);
printInfo(handles);
%% ------------------------------------------------------------ 


%% ------------------------------------------------------------ 
%% --- Increment column number 
%% ------------------------------------------------------------ 
function handles = incrementColumnNumber(handles)
if handles.cc < handles.ccMax
  handles.cc = handles.cc + 1;
end
handles = updateDisplay_cc(handles);
printInfo(handles);
%% ------------------------------------------------------------ 


%% ------------------------------------------------------------ 
%% --- Decrement slice number 
%% ------------------------------------------------------------ 
function handles = decrementSliceNumber(handles)
if handles.ss > 1
  handles.ss = handles.ss - 1;
end
handles = updateDisplay_ss(handles);
printInfo(handles);
%% ------------------------------------------------------------ 


%% ------------------------------------------------------------ 
%% --- Increment slice number 
%% ------------------------------------------------------------ 
function handles = incrementSliceNumber(handles)
if handles.ss < handles.ssMax
  handles.ss = handles.ss + 1;
end
handles = updateDisplay_ss(handles);
printInfo(handles);
%% ------------------------------------------------------------ 


%% ------------------------------------------------------------ 
%% --- Fly to row/col/slice number stopIndex
%% ------------------------------------------------------------ 
function handles = fly(handles, stopIndex, varargin)

if nargin < 3
  pauseTime = 0;
else
  pauseTime = varargin{1};
end

switch handles.ORIENTATION
 case 1
  handles = fly_ss(handles, stopIndex, pauseTime);
 case 2
  handles = fly_cc(handles, stopIndex, pauseTime);
 case 3
  handles = fly_rr(handles, stopIndex, pauseTime);
end

%% end of function handles = fly(handles, stopIndex, varargin)
%% ------------------------------------------------------------ 


%% ------------------------------------------------------------ 
%% --- Fly to row number stop_rr
%% ------------------------------------------------------------ 
function handles = fly_rr(handles, stop_rr, pauseTime)

% Start
start_rr = handles.rr;

% Stop
stop_rr = max(stop_rr, 1);
stop_rr = min(stop_rr, handles.rrMax);

% Step
if stop_rr > handles.rr
  step = +1;
else
  step = -1;
end

% Begin loop
set(handles.toggleLines1, 'Value', 1);
set(handles.toggleLines1, 'ForegroundColor', 'k');
set(handles.axes1_X, 'Visible', 'off');
for rr = start_rr : step : stop_rr
  handles.rr = rr;
  handles = updateDisplay_rr(handles);
  printInfo(handles); % !@#
  pause(pauseTime);
end
if handles.HIDE1 == 0
  set(handles.toggleLines1, 'Value', 0);
  set(handles.toggleLines1, 'ForegroundColor', 'r');
  set(handles.axes1_X, 'Visible', 'on');
end

%% end of function function handles = fly_rr(handles, stop_rr, pauseTime)
%% ------------------------------------------------------------ 


%% ------------------------------------------------------------ 
%% --- Fly to column number stop_cc
%% ------------------------------------------------------------ 
function handles = fly_cc(handles, stop_cc, pauseTime)

% Start
start_cc = handles.cc;

% Stop
stop_cc = max(stop_cc, 1);
stop_cc = min(stop_cc, handles.ccMax);

% Step
if stop_cc > handles.cc
  step = +1;
else
  step = -1;
end

% Begin loop
set(handles.toggleLines1, 'Value', 1);
set(handles.toggleLines1, 'ForegroundColor', 'k');
set(handles.axes1_X, 'Visible', 'off');
for cc = start_cc : step : stop_cc
  handles.cc = cc;
  handles = updateDisplay_cc(handles);
  printInfo(handles); % !@#
  pause(pauseTime);
end
if handles.HIDE1 == 0
  set(handles.toggleLines1, 'Value', 0);
  set(handles.toggleLines1, 'ForegroundColor', 'r');
  set(handles.axes1_X, 'Visible', 'on');
end

%% end of function handles = fly_cc(handles, stop_cc, pauseTime)
%% ------------------------------------------------------------ 


%% ------------------------------------------------------------ 
%% --- Fly to slice number stop_ss
%% ------------------------------------------------------------ 
function handles = fly_ss(handles, stop_ss, pauseTime)

% Start
start_ss = handles.ss;

% Stop
stop_ss = max(stop_ss, 1);
stop_ss = min(stop_ss, handles.ssMax);

% Step
if stop_ss > handles.ss
  step = +1;
else
  step = -1;
end

% Begin loop
set(handles.toggleLines1, 'Value', 1);
set(handles.toggleLines1, 'ForegroundColor', 'k');
set(handles.axes1_X, 'Visible', 'off');
for ss = start_ss : step : stop_ss
  handles.ss = ss;
  handles = updateDisplay_ss(handles);
  printInfo(handles); % !@#
  pause(pauseTime);
end
if handles.HIDE1 == 0
  set(handles.toggleLines1, 'Value', 0);
  set(handles.toggleLines1, 'ForegroundColor', 'r');
  set(handles.axes1_X, 'Visible', 'on');
end
printInfo(handles);

%% end of function handles = fly_ss(handles, stop_ss, pauseTime)
%% ------------------------------------------------------------ 


%% ------------------------------------------------------------ 
%% --- Parse string given in command edit box
%% ------------------------------------------------------------ 
function handles = parseCommand(handles, str)

%% Empty -- do nothing
str = deblank(str);
if isempty(str)
  return
end

%% Jump to slice/col/row?
num = str2double(str);
if isfinite(num)
  switch handles.ORIENTATION
   case 1,
    ss = round(num);
    if (ss >= 1) & (ss <= handles.ssMax)
      handles.ss = ss;
      handles = updateDisplay_ss(handles);
      printInfo(handles);
    end
   case 2,
    cc = round(num);
    if (cc >= 1) & (cc <= handles.ccMax)
      handles.cc = cc;
      handles = updateDisplay_cc(handles);
      printInfo(handles);
    end
   case 3,
    rr = round(num);
    if (rr >= 1) & (rr <= handles.rrMax)
      handles.rr = rr;
      handles = updateDisplay_rr(handles);
      printInfo(handles);
    end
  end
  return;
end

%% Possible command in str
[strFirst, strRest] = strtok(str);
switch strFirst

  %% ------------------------------------------------------------ 
 case 'o'
  %% ------------------------------------------------------------ 
  
  handles = cycleOrientationsForward(handles);

  %% ------------------------------------------------------------ 
 case {'f', 'fly'}
  %% ------------------------------------------------------------ 

  WRONG_FORMAT = 1;
  if ~isempty(strRest)
    [strStopIndex, strRest] = strtok(strRest);
    stopIndex = round(str2double(strStopIndex));
    if isfinite(stopIndex)
      if isempty(strRest)
	handles = fly(handles, stopIndex);
	WRONG_FORMAT = 0;
      else
	[strPauseTime, strRest] = strtok(strRest);
	pauseTime = str2double(strPauseTime);
	if isfinite(pauseTime)
	  if pauseTime >= 0
	    handles = fly(handles, stopIndex, pauseTime);
	    WRONG_FORMAT = 0;
	  end
	end
      end
    end
  end
  if WRONG_FORMAT
    set(handles.edit_command, 'String', ...
		      'fly stopIndex [pauseTime]');
  end

  %% ------------------------------------------------------------ 
 case {'c', 'cycle'}
  %% ------------------------------------------------------------ 
  
  WRONG_FORMAT = 1;
  numLoops = 10;
  pauseTime = 0.1;
  if isempty(strRest)
    WRONG_FORMAT = 0;
    handles = cycleVolumesMovie(handles, numLoops, pauseTime);    
  else
    [strNumLoops, strRest] = strtok(strRest);
    numLoops = round(str2double(strNumLoops));
    if numLoops > 0
      if isempty(strRest)    
	WRONG_FORMAT = 0;
	handles = cycleVolumesMovie(handles, numLoops, pauseTime);    
      else
	[strPauseTime, strRest] = strtok(strRest);
	pauseTime = str2double(strPauseTime);
	if pauseTime >= 0
	  WRONG_FORMAT = 0;
	  handles = cycleVolumesMovie(handles, numLoops, pauseTime);    
	end
      end
    end
  end
  if WRONG_FORMAT
    set(handles.edit_command, 'String', ...
		      'cycle [numLoops] [pauseTime]');
  end

  %% ------------------------------------------------------------ 
 case 'lph'
  %% ------------------------------------------------------------ 
  
  WRONG_FORMAT = 1;
  if ~isempty(strRest)
    [strL, strRest] = strtok(strRest);
    targetL = str2double(strL);
    if isfinite(targetL)
      if ~isempty(strRest)
	[strP, strRest] = strtok(strRest);
	targetP = str2double(strP);
	if isfinite(targetP)
	  [strH, strRest] = strtok(strRest);
	  targetH = str2double(strH);
	  if isfinite(targetH)
	    targetLPH = [targetL; targetP; targetH; 1];
	    handles = gotoPointLPH(handles, targetLPH);
	    WRONG_FORMAT = 0;
	  end
	end
      end
    end
  end
  if WRONG_FORMAT
    set(handles.edit_command, 'String', ...
		      'lph gotoL gotoP gotoH');
  end

  %% ------------------------------------------------------------ 
 case 'rcs'
  %% ------------------------------------------------------------ 
  
  WRONG_FORMAT = 1;
  if ~isempty(strRest)
    [strR, strRest] = strtok(strRest);
    targetR = str2double(strR);
    if isfinite(targetR)
      if ~isempty(strRest)
	[strC, strRest] = strtok(strRest);
	targetC = str2double(strC);
	if isfinite(targetC)
	  [strS, strRest] = strtok(strRest);
	  targetS = str2double(strS);
	  if isfinite(targetS)
	    targetRCS = [targetR; targetC; targetS; 1];
	    handles = gotoPointRCS(handles, targetRCS);
	    WRONG_FORMAT = 0;
	  end
	end
      end
    end
  end
  if WRONG_FORMAT
    set(handles.edit_command, 'String', ...
		      'rcs gotoR gotoC gotoS');
  end

end % switch strFirst

%% end of function handles = parseCommand(handles, str)
%% ------------------------------------------------------------ 


%% ============================================================
%% 
%%               SECTION 4: UI CALLBACK FUNCTIONS
%% 
%% ============================================================


%% ------------------------------------------------------------ 
%% --- Executes on button press in cycleVolumes.
%% ------------------------------------------------------------ 
function cycleVolumes_Callback(hObject, eventdata, handles)
handles = cycleVolumes(handles);
guidata(hObject, handles);
%% ------------------------------------------------------------ 


%% ------------------------------------------------------------ 
%% --- Executes on button press in the volume pushbuttons
%% ------------------------------------------------------------ 
function volPB_Callback(hObject, eventdata, handles)
[minVal, chosenVol] = min(abs(handles.volPB - hObject));
if chosenVol == handles.currentVol
  if sum(handles.hideVol) < handles.numVols-1
    handles.hideVol(chosenVol) = 1;
    set(hObject, 'ForegroundColor', [0.5 0.5 0.5]);
    set(hObject, 'TooltipString', 'Show this volume');
    tmp = find(handles.hideVol == 0);
    tmpInd = find(tmp > handles.currentVol);
    if isempty(tmpInd)
      tmpInd = 1;
    end
    handles.currentVol = tmp(tmpInd(1));
    handles.currentVolPB = handles.volPB(handles.currentVol);
  end
else
  set(handles.currentVolPB, 'ForegroundColor', 'k');
  set(handles.currentVolPB, 'TooltipString', 'Show this volume');
  handles.hideVol(chosenVol) = 0;
  handles.currentVol = chosenVol;
  handles.currentVolPB = hObject;
end
set(handles.currentVolPB, 'ForegroundColor', 'r');
set(handles.currentVolPB, 'TooltipString', 'Hide this volume');
handles.Mvxl2lph = handles.vols{handles.currentVol}.Mvxl2lph;
updateDisplay_newvol(handles);
guidata(hObject, handles);
%% ------------------------------------------------------------ 


%% ------------------------------------------------------------ 
%% --- Executes on button press in zoom_reset.
%% ------------------------------------------------------------ 
function zoom_reset_Callback(hObject, eventdata, handles)
handles.zoom = 0;
set(handles.zoom_slider, 'Value', 0);
handles = updateDisplay_zoom(handles);
guidata(hObject, handles);
%% ------------------------------------------------------------ 


%% ------------------------------------------------------------ 
%% --- Executes on zoom slider movement.
%% ------------------------------------------------------------ 
function zoom_slider_Callback(hObject, eventdata, handles)
handles.zoom = get(hObject,'Value');
handles.rrZoom = handles.rr;
handles.ccZoom = handles.cc;
handles.ssZoom = handles.ss;
handles = updateDisplay_zoom(handles);
guidata(hObject, handles);
%% ------------------------------------------------------------ 


%% ------------------------------------------------------------ 
%% --- Executes on button press in contrast_reset.
%% ------------------------------------------------------------ 
function contrast_reset_Callback(hObject, eventdata, handles)
set(handles.contrast_slider, 'Value', 0);
handles.cVal(handles.currentVol) = 0;
handles = updateDisplay_clims(handles);
guidata(hObject, handles);
%% ------------------------------------------------------------ 


%% ------------------------------------------------------------ 
%% --- Executes on contrast slider movement.
%% ------------------------------------------------------------ 
function contrast_slider_Callback(hObject, eventdata, handles)
cVal = get(handles.contrast_slider, 'Value');
handles.cVal(handles.currentVol) = cVal;
handles = updateDisplay_clims(handles);
guidata(hObject, handles);
%% ------------------------------------------------------------ 


%% ------------------------------------------------------------ 
%% --- Executes on button press in darken_reset.
%% ------------------------------------------------------------ 
function darken_reset_Callback(hObject, eventdata, handles)
set(handles.darken_slider, 'Value', 0);
handles.dVal(handles.currentVol) = 0;
handles = updateDisplay_clims(handles);
guidata(hObject, handles);
%% ------------------------------------------------------------ 


%% ------------------------------------------------------------ 
%% --- Executes on darken slider movement.
%% ------------------------------------------------------------ 
function darken_slider_Callback(hObject, eventdata, handles)
dVal = get(handles.darken_slider, 'Value');
handles.dVal(handles.currentVol) = dVal;
handles = updateDisplay_clims(handles);
guidata(hObject, handles);
%% ------------------------------------------------------------ 


%% ------------------------------------------------------------ 
%% --- Executes on entering into edit_command edit box
%% ------------------------------------------------------------ 
function edit_command_Callback(hObject, eventdata, handles)
str = get(hObject, 'String');
set(handles.text_last_echo, 'String', str);
set(hObject, 'String', '');
handles = parseCommand(handles, str);
guidata(hObject, handles);
%% ------------------------------------------------------------ 


%% ------------------------------------------------------------ 
%% --- Executes on button press in storePoint.
%% ------------------------------------------------------------ 
function storePoint_Callback(hObject, eventdata, handles)
handles = storeCurrentPoint(handles);
guidata(hObject, handles);
%% ------------------------------------------------------------ 


%% ------------------------------------------------------------ 
%% --- Executes on button press in gotoPoint.
%% ------------------------------------------------------------ 
function gotoPoint_Callback(hObject, eventdata, handles)
handles = gotoStoredPoint(handles);
guidata(hObject, handles);
%% ------------------------------------------------------------ 


%% ------------------------------------------------------------ 
%% --- Executes on button press in swapPoint.
%% ------------------------------------------------------------ 
function swapPoint_Callback(hObject, eventdata, handles)
handles = swapCurrentAndStoredPoints(handles);
guidata(hObject, handles);
%% ------------------------------------------------------------ 


%% ------------------------------------------------------------ 
%% --- Executes on button press in decrementRowNumber.
%% ------------------------------------------------------------ 
function decrementRowNumber_Callback(hObject, eventdata, handles)
handles = decrementRowNumber(handles);
guidata(hObject, handles);
%% ------------------------------------------------------------ 


%% ------------------------------------------------------------ 
%% --- Executes on button press in decrementColumnNumber.
%% ------------------------------------------------------------ 
function decrementColumnNumber_Callback(hObject, eventdata, handles)
handles = decrementColumnNumber(handles);
guidata(hObject, handles);
%% ------------------------------------------------------------ 


%% ------------------------------------------------------------ 
%% --- Executes on button press in decrementSliceNumber.
%% ------------------------------------------------------------ 
function decrementSliceNumber_Callback(hObject, eventdata, handles)
handles = decrementSliceNumber(handles);
guidata(hObject, handles);
%% ------------------------------------------------------------ 


%% ------------------------------------------------------------ 
%% --- Executes on entering into edit_rr edit box
%% ------------------------------------------------------------ 
function edit_rr_Callback(hObject, eventdata, handles)
rr = round(str2double(get(hObject,'String')));
if isfinite(rr)
  if (rr >= 1) & (rr <= handles.rrMax)
    handles.rr = rr;
    handles = updateDisplay_rr(handles);
    printInfo(handles);
    guidata(hObject, handles);
  end
end
set(hObject,'String',num2str(handles.rr));
%% ------------------------------------------------------------ 


%% ------------------------------------------------------------ 
%% --- Executes on entering into edit_cc edit box
%% ------------------------------------------------------------ 
function edit_cc_Callback(hObject, eventdata, handles)
cc = round(str2double(get(hObject,'String')));
if isfinite(cc)
  if (cc >= 1) & (cc <= handles.ccMax)
    handles.cc = cc;
    handles = updateDisplay_cc(handles);
    printInfo(handles);
    guidata(hObject, handles);
  end
end
set(hObject,'String',num2str(handles.cc));
%% ------------------------------------------------------------ 


%% ------------------------------------------------------------ 
%% --- Executes on entering into edit_ss edit box
%% ------------------------------------------------------------ 
function edit_ss_Callback(hObject, eventdata, handles)
ss = round(str2double(get(hObject,'String')));
if isfinite(ss)
  if (ss >= 1) & (ss <= handles.ssMax)
    handles.ss = ss;
    handles = updateDisplay_ss(handles);
    printInfo(handles);
    guidata(hObject, handles);
  end
end
set(hObject,'String',num2str(handles.ss));
%% ------------------------------------------------------------ 


%% ------------------------------------------------------------ 
%% --- Executes on button press in incrementRowNumber.
%% ------------------------------------------------------------ 
function incrementRowNumber_Callback(hObject, eventdata, handles)
handles = incrementRowNumber(handles);
guidata(hObject, handles);
%% ------------------------------------------------------------ 


%% ------------------------------------------------------------ 
%% --- Executes on button press in incrementColumnNumber.
%% ------------------------------------------------------------ 
function incrementColumnNumber_Callback(hObject, eventdata, handles)
handles = incrementColumnNumber(handles);
guidata(hObject, handles);
%% ------------------------------------------------------------ 


%% ------------------------------------------------------------ 
%% --- Executes on button press in incrementSliceNumber.
%% ------------------------------------------------------------ 
function incrementSliceNumber_Callback(hObject, eventdata, handles)
handles = incrementSliceNumber(handles);
guidata(hObject, handles);
%% ------------------------------------------------------------ 


%% ------------------------------------------------------------ 
%% --- Executes on button press in toggleLines1.
function toggleLines1_Callback(hObject, eventdata, handles)
handles = toggleLines1(handles);
guidata(hObject, handles);
%% ------------------------------------------------------------ 


%% ------------------------------------------------------------ 
%% --- Executes on button press in toggleLines2.
%% ------------------------------------------------------------ 
function toggleLines2_Callback(hObject, eventdata, handles)
handles = toggleLines2(handles);
guidata(hObject, handles);
%% ------------------------------------------------------------ 


%% ------------------------------------------------------------ 
%% --- Executes on button press in toggleLines3.
%% ------------------------------------------------------------ 
function toggleLines3_Callback(hObject, eventdata, handles)
handles = toggleLines3(handles);
guidata(hObject, handles);
%% ------------------------------------------------------------ 


%% ------------------------------------------------------------ 
%% --- Executes on button press in cycleOrientations1.
%% ------------------------------------------------------------ 
function cycleOrientations1_Callback(hObject, eventdata, handles)
handles = cycleOrientationsForward(handles);
guidata(hObject, handles);
%% ------------------------------------------------------------ 


%% ------------------------------------------------------------ 
%% --- Executes on button press in cycleOrientations2.
%% ------------------------------------------------------------ 
function cycleOrientations2_Callback(hObject, eventdata, handles)
handles = cycleOrientationsReverse(handles);
guidata(hObject, handles);
%% ------------------------------------------------------------ 


%% ------------------------------------------------------------ 
%% --- Executes on mouse button press in figure window
%% ------------------------------------------------------------ 
function figure_WindowButtonDownFcn(hObject, eventdata, handles)

% Here, we hack around the fact that MATLAB R13 v6.5 does not seem to
% trap the axes button down callback -- instead, this is caught by the
% figure button down callback. So the code belows figures out which
% axes you clicked in, and calls the appropriate axes callback. 

pos_axCR = get(handles.axesCR,'CurrentPoint');
XLim_axCR = get(handles.axesCR, 'XLim');
YLim_axCR = get(handles.axesCR, 'YLim');

pos_axSR = get(handles.axesSR,'CurrentPoint');
XLim_axSR = get(handles.axesSR, 'XLim');
YLim_axSR = get(handles.axesSR, 'YLim');

pos_axCS = get(handles.axesCS,'CurrentPoint');
XLim_axCS = get(handles.axesCS, 'XLim');
YLim_axCS = get(handles.axesCS, 'YLim');

if     ((pos_axCR(1,1) >= XLim_axCR(1)) & ...
	(pos_axCR(1,1) <= XLim_axCR(2)) & ...
	(pos_axCR(1,2) >= YLim_axCR(1)) & ...
	(pos_axCR(1,2) <= YLim_axCR(2)))   
  axesCR_ButtonDownFcn(handles.axesCR, eventdata, handles);
elseif ((pos_axSR(1,1) >= XLim_axSR(1)) & ...
	(pos_axSR(1,1) <= XLim_axSR(2)) & ...
	(pos_axSR(1,2) >= YLim_axSR(1)) & ...
	(pos_axSR(1,2) <= YLim_axSR(2)))   
  axesSR_ButtonDownFcn(handles.axesSR, eventdata, handles);
elseif ((pos_axCS(1,1) >= XLim_axCS(1)) & ...
	(pos_axCS(1,1) <= XLim_axCS(2)) & ...
	(pos_axCS(1,2) >= YLim_axCS(1)) & ...
	(pos_axCS(1,2) <= YLim_axCS(2)))   
  axesCS_ButtonDownFcn(handles.axesCS, eventdata, handles);
end
%% end of function figure_WindowButtonDownFcn(hObject, eventdata, handles)
%% ------------------------------------------------------------ 


%% ------------------------------------------------------------ 
%% --- Executes on mouse press over axesCR background.
%% ------------------------------------------------------------ 
function axesCR_ButtonDownFcn(hObject, eventdata, handles)

pos = get(hObject, 'CurrentPoint');
handles.cc = round(pos(1,1));
handles.cc = max(handles.cc, 1);
handles.cc = min(handles.cc, handles.ccMax);
handles.rr = round(pos(1,2));
handles.rr = max(handles.rr, 1);
handles.rr = min(handles.rr, handles.rrMax);
handles = updateDisplay_cc(handles);
handles = updateDisplay_rr(handles);
printInfo(handles);
guidata(hObject, handles);

%% end of function axesCR_ButtonDownFcn(hObject, eventdata, handles)
%% ------------------------------------------------------------ 


%% ------------------------------------------------------------ 
%% --- Executes on mouse press over axes3 background.
%% ------------------------------------------------------------ 
function axesSR_ButtonDownFcn(hObject, eventdata, handles)

pos = get(hObject, 'CurrentPoint');
handles.ss = round(pos(1,1));
handles.ss = max(handles.ss, 1);
handles.ss = min(handles.ss, handles.ssMax);
handles.rr = round(pos(1,2));
handles.rr = max(handles.rr, 1);
handles.rr = min(handles.rr, handles.rrMax);
handles = updateDisplay_ss(handles);
handles = updateDisplay_rr(handles);
printInfo(handles);
guidata(hObject, handles);

%% end of function axesSR_ButtonDownFcn(hObject, eventdata, handles)
%% ------------------------------------------------------------ 


%% ------------------------------------------------------------ 
%% --- Executes on mouse press over axes background.
%% ------------------------------------------------------------ 
function axesCS_ButtonDownFcn(hObject, eventdata, handles)

pos = get(hObject, 'CurrentPoint');
handles.cc = round(pos(1,1));
handles.cc = max(handles.cc, 1);
handles.cc = min(handles.cc, handles.ccMax);
handles.ss = round(pos(1,2));
handles.ss = max(handles.ss, 1);
handles.ss = min(handles.ss, handles.ssMax);
handles = updateDisplay_cc(handles);
handles = updateDisplay_ss(handles);
printInfo(handles);
guidata(hObject, handles);

%% function axesCS_ButtonDownFcn(hObject, eventdata, handles)
%% ------------------------------------------------------------ 


%% ------------------------------------------------------------ 
%% --- Executes on key press anywhwere in figure window
%% ------------------------------------------------------------ 
function figure_KeyPressFcn(hObject, eventdata, handles)

key = lower(get(hObject, 'CurrentCharacter'));
switch key
  
 case 'c'
  handles = cycleVolumesMovie(handles, 10, 0.1);
 
 case 'v'
  handles = cycleVolumes(handles);

 case 'o'
  handles = cycleOrientationsForward(handles);
  
 case 't'
  handles = toggleLines1(handles);

 case handles.rrDecrChar
  handles = decrementRowNumber(handles);
  
 case handles.rrIncrChar
  handles = incrementRowNumber(handles);
  
 case handles.ccDecrChar
  handles = decrementColumnNumber(handles);
  
 case handles.ccIncrChar
  handles = incrementColumnNumber(handles);
  
 case handles.ssDecrChar
  handles = decrementSliceNumber(handles);  
  
 case handles.ssIncrChar
  handles = incrementSliceNumber(handles);
  
 case {'p', ' '}
  handles = storeCurrentPoint(handles);
  set(handles.storePoint, 'Value', 1);
  pause(0.05);
  set(handles.storePoint, 'Value', 0);

 case ';'
  handles = gotoStoredPoint(handles);
  set(handles.gotoPoint, 'Value', 1);
  pause(0.05);
  set(handles.gotoPoint, 'Value', 0);

 case '/'
  handles = swapCurrentAndStoredPoints(handles);
  set(handles.swapPoint, 'Value', 1);
  pause(0.05);
  set(handles.swapPoint, 'Value', 0);

 case 'q'
  cVal = handles.cVal(handles.currentVol) - 0.1;
  cVal = max(cVal, get(handles.contrast_slider, 'Min'));
  handles.cVal(handles.currentVol) = cVal;
  set(handles.contrast_slider, 'Value', cVal);
  handles = updateDisplay_clims(handles);
  
 case 'w'
  cVal = handles.cVal(handles.currentVol) + 0.1;
  cVal = min(cVal, get(handles.contrast_slider, 'Max'));
  handles.cVal(handles.currentVol) = cVal;
  set(handles.contrast_slider, 'Value', cVal);
  handles = updateDisplay_clims(handles);
  
 case 'a'
  dVal = handles.dVal(handles.currentVol) - 0.1;
  dVal = max(dVal, get(handles.darken_slider, 'Min'));
  handles.dVal(handles.currentVol) = dVal;
  set(handles.darken_slider, 'Value', dVal);
  handles = updateDisplay_clims(handles);
  
 case 's'
  dVal = handles.dVal(handles.currentVol) + 0.1;
  dVal = min(dVal, get(handles.darken_slider, 'Max'));
  handles.dVal(handles.currentVol) = dVal;
  set(handles.darken_slider, 'Value', dVal);
  handles = updateDisplay_clims(handles);
  
 case 'z'
  zoom = handles.zoom - 0.6;
  zoom = max(zoom, get(handles.zoom_slider, 'Min'));
  handles.zoom = zoom;
  set(handles.zoom_slider, 'Value', zoom);
  handles.rrZoom = handles.rr;
  handles.ccZoom = handles.cc;
  handles.ssZoom = handles.ss;
  handles = updateDisplay_zoom(handles);
  
 case 'x'
  zoom = handles.zoom + 0.6;
  zoom = min(zoom,  get(handles.zoom_slider, 'Max'));
  handles.zoom = zoom;
  set(handles.zoom_slider, 'Value', zoom);
  handles.rrZoom = handles.rr;
  handles.ccZoom = handles.cc;
  handles.ssZoom = handles.ss;
  handles = updateDisplay_zoom(handles);
  
end

guidata(hObject, handles);

%% end of function figure_KeyPressFcn(hObject, eventdata, handles)
%% ------------------------------------------------------------ 


%% ------------------------------------------------------------ 
%% --- Executes when figure window is resized.
%% ------------------------------------------------------------ 
function figure_ResizeFcn(hObject, eventdata, handles)

if isempty(handles)
  return;
end

% Set up resize function
set(handles.figure, 'Units', 'centimeters');
figpos_cm = get(handles.figure, 'Position');
set(handles.figure, 'Units', 'characters');
figpos_ch = get(handles.figure, 'Position');
handles.chFactor = figpos_ch(4)/figpos_ch(3)*figpos_cm(3)/figpos_cm(4);

%% ---------------------------------------
lB = 6;       % left border
tB = 1;       % top border
bB = 10;      % bottom border
rB = 6;       % right border
	      
tAx123 = 1;   % vertical space above axes1, axes2 and axes3
rAx1 = 7;     % horizontal space right of axes1
rAx23 = 0;    % horizontal space right of axes2 and axes3
sAx2Ax3 = 2;  % vertical space between axes2 and axes3
	      
wPBt = 7;     % width of topmost panel pushbuttons
hPBt = 2;     % height of topmost panel pushbuttons
sPBt = 1;     % horizontal space between topmost panel pushbuttons

wPBr = 0;     % width of rightmost panel pushbuttons
hPBr = 2;     % height of rightmost panel pushbuttons
sPBr = 0;     % vertical space between rightmost panel pushbuttons

ax1MinH = 45; % don't let ax1 height get smaller than this

figPos_orig = get(handles.figure, 'Position');
W_orig = figPos_orig(3);
H_orig = figPos_orig(4);
%% ---------------------------------------

%% Compute the axes sizes for a window of the current size
H_padding = bB + tAx123 + hPBt + tB;
W_padding = lB + rAx1 + rAx23 + wPBr + rB;
ax1_h = max(H_orig - H_padding, ax1MinH);
ax1_w = ax1_h/handles.chFactor;
ax23_h = (ax1_h - sAx2Ax3)/2;
ax23_w = ax23_h/handles.chFactor;

%% ------------------------------------------------------------ 
set(handles.axes1, 'Position', [lB bB ax1_w ax1_h]);
tmp1 = lB + ax1_w + rAx1;
set(handles.axes3, 'Position', [tmp1 bB ax23_w ax23_h]);
tmp2 = bB + ax23_h + sAx2Ax3;
set(handles.axes2, 'Position', [tmp1 tmp2 ax23_w ax23_h]);

%% ------------------------------------------------------------ 
bXO_w = 3.5;
bXO_h = 1.5;
tmp1 = lB + ax1_w - bXO_w;
tmp2 = bB + ax1_h - bXO_h;
set(handles.toggleLines1, 'Position', [tmp1 tmp2 bXO_w bXO_h]);
tmp1 = lB + ax1_w + rAx1 + ax23_w - bXO_w;
tmp2 = bB + ax1_h - bXO_h;
set(handles.toggleLines2, 'Position', [tmp1 tmp2 bXO_w bXO_h]);
tmp2 = bB + ax23_h - bXO_h;
set(handles.toggleLines3, 'Position', [tmp1 tmp2 bXO_w bXO_h]);
tmp2 = bB + ax23_h + sAx2Ax3;
set(handles.cycleOrientations1, 'Position', [tmp1 tmp2 bXO_w bXO_h]);
tmp2 = bB;
set(handles.cycleOrientations2, 'Position', [tmp1 tmp2 bXO_w bXO_h]);

%% ------------------------------------------------------------ 
tmp1 = lB;
tmp2 = bB + ax1_h + tAx123;
set(handles.cycleVolumes, 'Position', [tmp1 tmp2 wPBt hPBt]);
if isfield(handles, 'volPB')
  for vvv = 1:length(handles.volPB)
    tmp1 = tmp1 + wPBt + sPBt;
    set(handles.volPB(vvv), 'Position', [tmp1 tmp2 wPBt hPBt]);
  end
  tmp1 = lB + wPBt + sPBt/2;
  tmp2 = bB + ax1_h + tAx123/2;
  tmp3 = length(handles.volPB)*(wPBt + sPBt) + sPBt/6;
  tmp4 = hPBt + tAx123;
  set(handles.volFrame, 'Position', [tmp1 tmp2 tmp3 tmp4]);
end

%% ------------------------------------------------------------ 
switch version('-release')
 case '13', 
  backGray = [0.7 0.7 0.7];
  liteGray = [0.85 0.85 0.85];
  fontSize = 10;
 case '14', 
  backGray = [0.69 0.71 0.70]*0.97;
  liteGray = [0.84 0.86 0.85]*0.97;
  fontSize = 9;
  otherwise, 
  backGray = [0.69 0.71 0.70]*0.97;
  liteGray = [0.84 0.86 0.85]*0.97;
  fontSize = 9;
end

set(handles.figure, 'Color', [0.7 0.7 0.7]);

set(handles.toggleLines1, 'BackgroundColor', backGray);
set(handles.toggleLines2, 'BackgroundColor', backGray);
set(handles.toggleLines3, 'BackgroundColor', backGray);
set(handles.cycleOrientations1, 'BackgroundColor', backGray);
set(handles.cycleOrientations2, 'BackgroundColor', backGray);
set(handles.cycleVolumes, 'BackgroundColor', backGray);

if (isfield(handles, 'volPB'))
    set(handles.volPB, 'BackgroundColor', backGray);
end
set(handles.volFrame, 'BackgroundColor', backGray);
set(handles.frame5, 'BackgroundColor', backGray);
set(handles.contrast_reset, 'BackgroundColor', backGray);
set(handles.darken_reset, 'BackgroundColor', backGray);
set(handles.zoom_reset, 'BackgroundColor', backGray);
set(handles.contrast_slider, 'BackgroundColor', backGray);
set(handles.darken_slider, 'BackgroundColor', backGray);
set(handles.zoom_slider, 'BackgroundColor', backGray);
set(handles.frame4, 'BackgroundColor', backGray);
set(handles.text_last, 'BackgroundColor', backGray);
set(handles.text_last_echo, 'BackgroundColor', backGray);
set(handles.text_command, 'BackgroundColor', backGray);
set(handles.frame3, 'BackgroundColor', backGray);
set(handles.text_voxvalstr, 'BackgroundColor', backGray);
set(handles.text_voxval, 'BackgroundColor', backGray);
set(handles.text_scrvalstr, 'BackgroundColor', backGray);
set(handles.text_scrval, 'BackgroundColor', backGray);
set(handles.frame6, 'BackgroundColor', backGray);
set(handles.storePoint, 'BackgroundColor', backGray);
set(handles.gotoPoint, 'BackgroundColor', backGray);
set(handles.swapPoint, 'BackgroundColor', backGray);
set(handles.frame2, 'BackgroundColor', backGray);
set(handles.text_L, 'BackgroundColor', backGray);
set(handles.text_ll, 'BackgroundColor', backGray);
set(handles.text_P, 'BackgroundColor', backGray);
set(handles.text_pp, 'BackgroundColor', backGray);
set(handles.text_H, 'BackgroundColor', backGray);
set(handles.text_hh, 'BackgroundColor', backGray);
set(handles.frame1, 'BackgroundColor', backGray);
set(handles.text_row, 'BackgroundColor', backGray);
set(handles.decrementRowNumber, 'BackgroundColor', backGray);
set(handles.incrementRowNumber, 'BackgroundColor', backGray);
set(handles.text_col, 'BackgroundColor', backGray);
set(handles.decrementColumnNumber, 'BackgroundColor', backGray);
set(handles.incrementColumnNumber, 'BackgroundColor', backGray);
set(handles.text_slice, 'BackgroundColor', backGray);
set(handles.decrementSliceNumber, 'BackgroundColor', backGray);
set(handles.incrementSliceNumber, 'BackgroundColor', backGray);

set(handles.edit_command, 'BackgroundColor', liteGray); %%%
set(handles.edit_rr, 'BackgroundColor', liteGray);
set(handles.edit_cc, 'BackgroundColor', liteGray);
set(handles.edit_ss, 'BackgroundColor', liteGray);

set(handles.cycleVolumes, 'FontSize', fontSize);

if (isfield(handles, 'volPB'))
    set(handles.volPB, 'FontSize', fontSize);
end
set(handles.contrast_reset, 'FontSize', fontSize);
set(handles.darken_reset, 'FontSize', fontSize);
set(handles.zoom_reset, 'FontSize', fontSize);
set(handles.text_last, 'FontSize', fontSize);
set(handles.text_last_echo, 'FontSize', fontSize);
set(handles.text_command, 'FontSize', fontSize);
set(handles.text_voxvalstr, 'FontSize', fontSize);
set(handles.text_voxval, 'FontSize', fontSize);
set(handles.text_scrvalstr, 'FontSize', fontSize);
set(handles.text_scrval, 'FontSize', fontSize);
set(handles.storePoint, 'FontSize', fontSize);
set(handles.gotoPoint, 'FontSize', fontSize);
set(handles.swapPoint, 'FontSize', fontSize);
set(handles.text_L, 'FontSize', fontSize);
set(handles.text_ll, 'FontSize', fontSize);
set(handles.text_P, 'FontSize', fontSize);
set(handles.text_pp, 'FontSize', fontSize);
set(handles.text_H, 'FontSize', fontSize);
set(handles.text_hh, 'FontSize', fontSize);
set(handles.text_row, 'FontSize', fontSize);
set(handles.decrementRowNumber, 'FontSize', fontSize);
set(handles.incrementRowNumber, 'FontSize', fontSize);
set(handles.text_col, 'FontSize', fontSize);
set(handles.decrementColumnNumber, 'FontSize', fontSize);
set(handles.incrementColumnNumber, 'FontSize', fontSize);
set(handles.text_slice, 'FontSize', fontSize);
set(handles.decrementSliceNumber, 'FontSize', fontSize);
set(handles.incrementSliceNumber, 'FontSize', fontSize);
 
set(handles.edit_command, 'FontSize', fontSize); %%%
set(handles.edit_rr, 'FontSize', fontSize);
set(handles.edit_cc, 'FontSize', fontSize);
set(handles.edit_ss, 'FontSize', fontSize);
