function varargout = showVol(varargin)
% Viewer for 3-D images with mesh cross-sections (optionally) overlaid.
%
% If used with ABCD brain volumes, additional features are available:
%   atlas ROI display; atlas lookup; descriptive anatomical annotations
%
% showVol(vol_1, [vol_2], ..., [vol_N], ...
%             [mesh_1], [mesh_2], ..., [mesh_N], ...
%             [initialLPH | configStruct])
%
% NOTE 1: vol_* is a struct with (at least) the fields
%        'imgs' and 'Mvxl2lph', a struct or cell array of same, or an RCS x N matrix
%
% NOTE 2: mesh_* is a struct which has (at least) the fields
%         'vertices' and 'faces'. If, in addition, the field 'color'
%         is given (as an RGB triple or a single char e.g. 'r'), then
%         this specifies the color with which that particular mesh is
%         drawn.
%
% NOTE 3: initialLPH is a vector of length 3 that specifies the
%         LPH coordinates of the starting point of the crosshair
%
% NOTE 4: Final agument can be a  configuration structure with fields:
%           RCS         initial RCS coords
%           crange      color axis range
%           cmap        colormap
%           winpos      initial window position (in characters)
%                         default is [0.5 1.5 167 59] % minimal size
%                         another option: [0.5 0.5 200 100];  % big coronal view
%           linewidth   line width of crosshairs (default 2)
%           link        link coordinates across figures (default: true)
%
%         ABCD-specific features
%           roiatlas    REQUIRED: atlas version for ROIs: ABCD1, ABCD2 OR 5.0_ABCD3
%           roiset      path to ROI settings file (a set of displayed ROIs, saved previously--very useful)
%           annotations true/false(default): Enable anatomical annotations, a 'yelp' for the brain
%
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
% -- '!': save screenshot of main axis
% -- '`': toggle zoomed inset for FOD images (default on)
% -- '\': toggle linking of RCS coordinates across multiple showVol figures (default on)
%
%  ABCD annotations-specific keys
% -- '+': add an annotation to the current point.
% -- '-': edit (or delete) an annotation at the current point
% -- 'shift+cmd+ '-': undo last add/edit/delete
% -- ']': list all annotations
%
% This software is Copyright (c) 2022 The Regents of the University of California. All Rights Reserved.
% See LICENSE.


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
  'gui_OpeningFcn', @showVolMesh_OpeningFcn, ...
  'gui_OutputFcn',  @showVolMesh_OutputFcn, ...
  'gui_LayoutFcn',  [] , ...
  'gui_Callback',   []);
if nargin && isstr(varargin{1})
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
%% --- Executes just before showVolMesh is made visible.
%% ------------------------------------------------------------
function showVolMesh_OpeningFcn(hObject, eventdata, handles, varargin)

% Initialize figure
if strcmp(version('-release'), '14')
  cameratoolbar('Close'); % this takes 1 second!
else
  cameratoolbar('hide');
end
colormap(gray);
set(handles.figure, 'Renderer', 'zbuffer');

set(handles.figure, 'Tag','showVol')

%check dependencies (make sure utils is on path)
handles = checkDependencies(handles);

if isempty(varargin)
  fprintf('--------------------------------\n')
  fprintf('%s: no inputs given\n', mfilename);
  fprintf('--------------------------------\n')
  error('Exiting ...');
end

% Read in inputs
handles.LineWidth = 2;
handles.numVols = 0;
handles.numSurfs = 0;
initialLPH = [];
initialRCS = [];
handles.crange = [];
handles.cmap = gray;
winposArg = [];
handles.doLinkCoordinates = true;
handles.hasABCDBrain = false; %keep track if any volumes are ABCD Brain
handles.doAnnotations = false; %only enable if explicitly selected
handles.anat.roiatlas = [];
handles.anat.roifile = [];
handles.anat.roiset = [];

for iii = 1:nargin-3
  if iscell(varargin{iii}) % AMD: Handle cell arrays of vols and/or meshes -- recursively call function for cell arrays (avoid duplication)
    v = varargin{iii};
    for vi = 1:length(v)
      if (isfield(v{vi}, 'imgs')) && (isfield(v{vi}, 'Mvxl2lph'))
        % volStruct
        handles.numVols = handles.numVols + 1;
        handles.vols{handles.numVols} = v{vi};
        canonicalSize = [200 200 260];
        sz = size(handles.vols{handles.numVols}.imgs);
        isABCDBrain = isequal(sz(1:3), canonicalSize);
        handles.vols{handles.numVols}.isABCDBrain = isABCDBrain;
        handles.hasABCDBrain = isABCDBrain;
      elseif isfield(v{vi}, 'vertices')
        % meshStruct
        handles.numSurfs = handles.numSurfs + 1;
        handles.surfs(handles.numSurfs).vertices = v{vi}.vertices;
        handles.surfs(handles.numSurfs).faces = v{vi}.faces;
        handles.surfs(handles.numSurfs).color = [1 1 0];
        if isfield(v{vi}, 'color')
          if ~isempty(v{vi}.color)
            handles.surfs(handles.numSurfs).color = v{vi}.color;
          end
        end
      elseif length(v{vi}) == 3
        % initial LPH coords
        initialLPH = v{vi};
        initialLPH = initialLPH(:);
        initialLPH(4) = 1;
      end
    end
  elseif isstruct(varargin{iii}) %Handle image volumes, meshes, or configuration structure
    if (isfield(varargin{iii}, 'imgs')) && (isfield(varargin{iii}, 'Mvxl2lph'))
      % volStruct
      handles.numVols = handles.numVols + 1;
      handles.vols{handles.numVols} = varargin{iii};
      canonicalSize = [200 200 260];
      sz = size(handles.vols{handles.numVols}.imgs);
      isABCDBrain = isequal(sz(1:3), canonicalSize);
      handles.vols{handles.numVols}.isABCDBrain = isABCDBrain;
      handles.hasABCDBrain = handles.hasABCDBrain | isABCDBrain;
    elseif (isfield(varargin{iii}, 'imgs1')) && (isfield(varargin{iii}, 'imgs2')) && (isfield(varargin{iii}, 'imgs3'))
      % pre-rendered images in each plane
      handles.numVols = handles.numVols + 1;
      handles.vols{handles.numVols} = varargin{iii};
      isABCDBrain = true; %FIXME just going to assume it for now, as this option only used for FOD, could also check for same aspect ratio    
      handles.vols{handles.numVols}.isABCDBrain = isABCDBrain;
      handles.hasABCDBrain = isABCDBrain;
    elseif isfield(varargin{iii}, 'vertices')
      % meshStruct
      handles.numSurfs = handles.numSurfs + 1;
      handles.surfs(handles.numSurfs).vertices = varargin{iii}.vertices;
      handles.surfs(handles.numSurfs).faces = varargin{iii}.faces;
      handles.surfs(handles.numSurfs).color = [1 1 0];
      if isfield(varargin{iii}, 'color')
        if ~isempty(varargin{iii}.color)
          handles.surfs(handles.numSurfs).color = varargin{iii}.color;
        end
      end
    end
    if isfield(varargin{iii}, 'RCS') % Should also handle RAS, LPH
      initialRCS = varargin{iii}.RCS; initialRCS(4) = 1;
    end
    if isfield(varargin{iii}, 'crange')
      handles.crange = varargin{iii}.crange;
    end
    if isfield(varargin{iii}, 'cmap')
      handles.cmap = varargin{iii}.cmap;
    end
    if isfield(varargin{iii}, 'winpos')
      winposArg = varargin{iii}.winpos;
    end
    if isfield(varargin{iii}, 'linewidth')
      handles.LineWidth = varargin{iii}.linewidth;
    end
    if isfield(varargin{iii}, 'link')
      handles.doLinkCoordinates = logical(varargin{iii}.link(1));
    end
    
    % %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % Handle ABCD-specific configuration options
    % %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%       
    if isfield(varargin{iii}, 'roifile')
      handles.anat.roifile = varargin{iii}.roifile;
    end
    if isfield(varargin{iii}, 'roiset')
      fname = varargin{iii}.roiset;
      if exist(fname,'file')
        handles.anat.roiset = fname;
      else
        disp('Your ROI settings file (roiset argument) does not exist.')
      end
    end
    if isfield(varargin{iii}, 'roiatlas')
      handles.anat.roiatlas = varargin{iii}.roiatlas;
    end
    if isfield(varargin{iii}, 'annotations')
      handles.doAnnotations = varargin{iii}.annotations;
    end
    % %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
  else
    if ndims(varargin{iii}) >= 3 % AMD: handle volumes specified as matrices -- convert to struct
      v = varargin{iii};
      assumeAtlas = '6.0_ABCD3'; %FIXME: need to know atlas version, assume 6.0 (Nov 24)
      
      %check if these might be brain volumes in ABCD atlas space; if so, must add correct transform for ROIs to display properly
      canonicalSize = [200 200 260];
      sz = size(v);
      isABCDBrain = isequal(sz(1:3), canonicalSize);
      if size(v,4)==3 %assume should be an RGB volume FIXME: better solution: ask folks to create a struct
        handles.numVols = handles.numVols + 1;
        if isABCDBrain
          fprintf('%s: Volume %d [%d x %d x %d x 3] initialized as RGB ABCD brain volume. Assuming %s atlas.\n',...
            mfilename, iii, sz(1), sz(2), sz(3), assumeAtlas)
          Mvxl2lph = Mvxl2lph_atlas(assumeAtlas);
          handles.vols{handles.numVols} = ctx_mgh2ctx(v, M_LPH_TO_RAS*Mvxl2lph);
          handles.vols{handles.numVols}.isABCDBrain = true;
          handles.hasABCDBrain = true;
        else
          handles.vols{handles.numVols} = ctx_mgh2ctx(v, eye(4));
          handles.vols{handles.numVols}.isABCDBrain = false;
        end
      else
        for vi = 1:size(v,4)
          handles.numVols = handles.numVols + 1;
          if isABCDBrain
            if vi == 1
              if length(sz)==3, sz(4)=1; end
              fprintf('%s: Volume %d [%d x %d x %d x %d] initialized as %d scalar ABCD brain volumes. Assuming %s atlas.\n',...
                mfilename, iii, sz(1), sz(2), sz(3), sz(4), sz(4), assumeAtlas)
            end
            Mvxl2lph = Mvxl2lph_atlas(assumeAtlas);
            handles.vols{handles.numVols} = ctx_mgh2ctx(v(:,:,:,vi), M_LPH_TO_RAS*Mvxl2lph);
            handles.vols{handles.numVols}.isABCDBrain = true;
            handles.hasABCDBrain = true;
          else
            handles.vols{handles.numVols} = ctx_mgh2ctx(v(:,:,:,vi), eye(4));
            handles.vols{handles.numVols}.isABCDBrain = false;
          end
        end
      end
      
      %JRI: auto-hide empty vols
      if sum(handles.vols{handles.numVols}.imgs(:))==0
        handles.hideVol(handles.numVols) = 1;
      end
      
    elseif length(varargin{iii}) == 3
      % initial LPH coords
      initialLPH = varargin{iii};
      initialLPH = initialLPH(:);
      initialLPH(4) = 1;
    end
  end
end

if ~handles.hasABCDBrain
  fprintf('%s: The volumes you''ve passed in don''t seem to be standard 1mm ABCD brain volumes, so things like ROIs will not work\n',mfilename)
  if handles.doAnnotations
    fprintf('%s: ...and annotations will be disabled\n', mfilename)
    handles.doAnnotations = false;
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

%% -------------------------------

%handle atlas defaults if a roifile or roiatlas was not specified
%FIXME: needs updating as we add future ABCD atlas versions
if handles.hasABCDBrain
  
  cfg = abcdConfig('showVol');
  
  if isempty(cfg.data.showVolData)
    error('%s: You must download showVolData first. See https://github.com/cmig-research-group/showVol for instructions.',mfilename) %TODO: UPDATE repository
  end
  
  if isempty(handles.anat.roiatlas)
    error('%s: You must explicitly specify an ABCD atlas version to use ROIs: ''ABCD1'', ''ABCD2'', ''5.0_ABCD3'' are valid as of Feb 2024.\n\n To do so, the final argument to showVol must be: struct(''roiatlas'',<atlasVersion>)',mfilename)
  end
  
  if isempty(handles.anat.roifile)
    handles.anat.roifile = showVolAtlasFile(handles.anat.roiatlas);
    disp([mfilename ': Using default ' handles.anat.roiatlas ' ROI file: ' handles.anat.roifile])
  else
    disp([mfilename ': Using custom ROI file: ' handles.anat.roifile])
  end

end

%% -------------------------------

%FIXME, this stuff is pretty obsolete. use defalt axis volume size
%find a 1mm volume; if none, fall back on old methods using first volume
voxSizes = [];
for vvv = 1:handles.numVols
  voxSizes = [voxSizes; norm(handles.vols{vvv}.Mvxl2lph(:,1)) ...
    norm(handles.vols{vvv}.Mvxl2lph(:,2)) ...
    norm(handles.vols{vvv}.Mvxl2lph(:,3))];
end

vol1mm = find(mean(voxSizes,2)-1 < eps);
if isempty(vol1mm)
  defaultVolume = 1;
else
  defaultVolume = vol1mm(1);
end
DV = handles.vols{defaultVolume};
rcs = size(DV.imgs(:,:,:,1));

%set up RCS coordinates
handles.rrMax = rcs(1);
handles.ccMax = rcs(2);
handles.ssMax = rcs(3);
handles.rcsScale = voxSizes(defaultVolume,:);

%initialze
handles.voxSize = voxSizes(defaultVolume,:);
if handles.hasABCDBrain
  handles.Mvxl2lph_atlas = Mvxl2lph_atlas(handles.anat.roiatlas);
else
  handles.Mvxl2lph_atlas = eye(4);
end
handles.Mvxl2lph = DV.Mvxl2lph; %FIXME: IS only used for transforming meshes. Will break if no volume is passed in

% Error checking I: all vols must have the same aspect ratio
voxNumCheck = [];
for vvv = 1:handles.numVols
  if isfield(handles.vols{vvv},'imgs')
    voxNumCheck = [voxNumCheck; size(handles.vols{vvv}.imgs(:,:,:,1))];
  elseif isfield(handles.vols{vvv},'imgs1')
    voxNumCheck = [voxNumCheck; handles.vols{vvv}.dimr handles.vols{vvv}.dimc handles.vols{vvv}.dimd];
  end
end
aspectRatio = voxNumCheck .* voxSizes;
if any(abs(sum(aspectRatio - aspectRatio(defaultVolume,:),2)) > 1e-3)
  fprintf('--------------------------------\n')
  fprintf('%s: all the vols must have the same aspect ratio \n', ...
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
    if ~isempty(handles.crange)
      minI(vvv) = handles.crange(1);
      maxI(vvv) = handles.crange(2);
    elseif isfield(handles.vols{vvv},'imgs')
      if exist('maxmin','file')
        [maxI(vvv), minI(vvv)] = maxmin(handles.vols{vvv}.imgs);
      else
        maxI(vvv) = max(handles.vols{vvv}.imgs(:));
        minI(vvv) = min(handles.vols{vvv}.imgs(:));
      end
%       if isfield(handles.vols{vvv},'limits')
%         minI(vvv)= handles.vols{vvv}.limits(1); %stat maps
%         maxI(vvv)= handles.vols{vvv}.limits(2);
%       end
    elseif isfield(handles.vols{vvv},'colormap') && isfield(handles.vols{vvv},'imgs1') %indexed FOD
      maxI(vvv)=255-eps;
      minI(vvv)=0;
    end
  end
end
handles.maxI = maxI;
handles.minI = minI;
handles.CLim = [minI maxI+eps];


% AMD: should change contrast controls to standard window/level?

handles.cVal = ones(handles.numVols,1);
set(handles.contrast_slider, 'Value', 1.0);
set(handles.contrast_slider,   'Min', 0.5);
set(handles.contrast_slider,   'Max', 1.5);

handles.dVal = 0.5*ones(handles.numVols,1);
set(handles.darken_slider, 'Value', 0.5);
set(handles.darken_slider,   'Min', 0);
set(handles.darken_slider,   'Max', 1);

% Cursor starting points
handles.rr = ceil(handles.rrMax/2 ); %NB 1mm images in ABCD1 and ABCD2 are centered on 100,100,130 not 101,101,131 as before.
handles.cc = ceil(handles.ccMax/2 );
handles.ss = ceil(handles.ssMax/2 );

%UM ABCD3 seems to be centered at 99,99,129!!
if contains('ABCD3', handles.anat.roiatlas)
  handles.rr = ceil(handles.rrMax/2 )-1;
  handles.cc = ceil(handles.ccMax/2 )-1;
  handles.ss = ceil(handles.ssMax/2 )-1;
end

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
defaultZoom = 0.24; %new default to maximize space
handles.zoom = defaultZoom;
set(handles.zoom_slider, 'Value', defaultZoom);
set(handles.zoom_slider, 'Min', 0);
set(handles.zoom_slider, 'Max', 6);

tmp_rr = (handles.rr - 0.5) * handles.rcsScale(1);
tmp_cc = (handles.cc - 0.5) * handles.rcsScale(2);
tmp_ss = (handles.ss - 0.5) * handles.rcsScale(3);
handles.halfRange_mm = max(max(tmp_rr, tmp_cc), tmp_ss);

tmp = [-handles.halfRange_mm +handles.halfRange_mm];
handles.rrLimOrig = handles.rr + tmp/handles.rcsScale(1);
handles.ccLimOrig = handles.cc + tmp/handles.rcsScale(2);
handles.ssLimOrig = handles.ss + tmp/handles.rcsScale(3);

handles.rrZoom = handles.rr;
handles.ccZoom = handles.cc;
handles.ssZoom = handles.ss;

% show a zoomed FOD around cursor as inset (hacky side effect, this also will toggle display of volume's name)
handles.showZoomedInset = true; % UI (key) command to toggle off is '`'

% brodcast/receive coords to/from other showVol figures
if handles.doLinkCoordinates
  handles = registerCoordinateListener(handles);
end

% initialize atlas lookup (ABCD only, nop if not)
handles = anatomySetup(handles);

%initialize annotations (ABCD only, nop if not)
handles = annotationSetup(handles);


%set up FOD zoomed inset
set(handles.zoomedInsetAxes, 'Visible','off')
handles.zoomedInsetImage = inf;

%there seems to be some invisible toolbars that take a lot of time processing
%mouse move events
delete(handles.CameraToolBar)
delete(handles.uitoolbar1)

%% ------------------------------------------------------------
% Initialize meshes

handles.showSurf = ones(handles.numSurfs, 1);

% For each given mesh
Mlph2vxl = inv(handles.Mvxl2lph);
%str = 'showVolMesh(''meshPB_Callback'',gcbo,[],guidata(gcbo))';
str = sprintf('%s(''meshPB_Callback'',gcbo,[],guidata(gcbo))',mfilename);
for mmm = 1:handles.numSurfs
  
  % Convert surfaces vertex coords from LPH to RCS
  tmp = handles.surfs(mmm).vertices;
  tmp(:,4) = 1;
  tmp = tmp*transpose(Mlph2vxl);
  handles.surfs(mmm).vertices = tmp(:,1:3);
  
  % Assign pushbutton for each mesh, if not given
  handles.meshPB(mmm) = uicontrol('Tag', 'meshButton', ...
    'Style', 'pushbutton', ...
    'Unit', 'characters', ...
    'Position', [10 10 10 10], ...
    'ForegroundColor', handles.surfs(mmm).color, ...
    'TooltipString', 'Hide this mesh', ...
    'String', ['mesh' num2str(mmm)]);
  set(handles.meshPB(mmm), 'Callback', str);
  
end

%% ------------------------------------------------------------
% Initial display

%str = 'showVolMesh(''volPB_Callback'',gcbo,[],guidata(gcbo))';
str = sprintf('%s(''volPB_Callback'',gcbo,[],guidata(gcbo))',mfilename);
for vvv = 1:handles.numVols
  if isfield(handles.vols{vvv}, 'name')
    handles.volnames{vvv} = handles.vols{vvv}.name;
  else
    handles.volnames{vvv} = 'this volume'; %default if no name
  end
  handles.volPB(vvv) = uicontrol('Tag', 'volButton', ...
    'Style', 'pushbutton', ...
    'Unit', 'characters', ...
    'Position', [10 10 10 10], ...
    'TooltipString', ['Show ' handles.volnames{vvv}], ...
    'String', ['vol' num2str(vvv)]);
  set(handles.volPB(vvv), 'Callback', str);
end
if isfield(handles,'hideVol')
  if length(handles.hideVol)<handles.numVols
    handles.hideVol(handles.numVols)=0;
  end
else
  handles.hideVol = zeros(handles.numVols, 1);
end
handles.currentVol = 1;
handles.currentVolPB = handles.volPB(1);
set(handles.currentVolPB, 'ForegroundColor', 'r');
set(handles.currentVolPB, 'TooltipString', ['Hide ' handles.volnames{handles.currentVol}]);
handles.ORIENTATION = 1;
handles.HIDE1 = 0;
handles.HIDE2 = 0;
handles.HIDE3 = 0;
handles = displayNewOrientation(handles);
printInfo(handles);
handles.output = hObject;

guidata(hObject, handles);
figure_ResizeFcn(hObject, eventdata, handles);

%set window size
if handles.numSurfs == 0
  %pos = [90 4.5 174 59];
  %  pos = [0.5 4.5 152 59];
  %  pos = [0.5 0.5 200 100];  % Anders's new fig size (May 2020)
  pos = [0.5 1.5 167 59];   % minimal size (fixes right side images clipped)
else
  %pos = [79 4.5 185 59];
  %  pos = [0.5 4.5 163 59];
  %  pos = [0.5 0.5 200 100];
  pos = [0.5 1.5 178 59];   % minimal size
end

%override default pos with passed-in argument
if ~isempty(winposArg)
  pos = winposArg;
end

set(handles.figure, 'Position', pos);

if ~isempty(initialRCS)
  initialLPH = inv(Mlph2vxl) * initialRCS(:); % colvec(initialRCS);
end

if ~isempty(initialLPH)
  handles = gotoPointLPH(handles, initialLPH);
  printInfo(handles);
  guidata(hObject, handles);
end

%% end of function showVolMesh_OpeningFcn(hObject, eventdata, handles, varargin)
%% ------------------------------------------------------------


%% ------------------------------------------------------------
%% --- Outputs from this function are returned to the command line.
%% ------------------------------------------------------------
function varargout = showVolMesh_OutputFcn(hObject, eventdata, handles)
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

if isfield(handles.vols{handles.currentVol},'imgs')
  sf = calculateScale(handles);
  voxVal = handles.vols{handles.currentVol}.imgs(scaleCoord(rr, sf.r), scaleCoord(cc,sf.c), scaleCoord(ss,sf.s), 1);
  CLim1 = handles.CLim(handles.currentVol, 1);
  CLim2 = handles.CLim(handles.currentVol, 2);
  scrVal = (voxVal - CLim1)/(CLim2 - CLim1);
  scrVal = min(scrVal, 1);
  scrVal = max(scrVal, 0);
else %pre-rendered, e.g. FOD, vox val unintersting
  voxVal = 0;
  scrVal = 0;
end

% Convert RCS to LPH
tmp = [rr; cc; ss; 1];
tmp = handles.Mvxl2lph_atlas * tmp;
ll = tmp(1); pp = tmp(2); hh = tmp(3);

set(handles.edit_rr, 'String', sprintf('%d', rr));
set(handles.edit_cc, 'String', sprintf('%d', cc));
set(handles.edit_ss, 'String', sprintf('%d', ss));

set(handles.text_ll, 'String', sprintf('%.3f', ll));
set(handles.text_pp, 'String', sprintf('%.3f', pp));
set(handles.text_hh, 'String', sprintf('%.3f', hh));

set(handles.text_voxval, 'String', sprintf('%.2f', voxVal));
set(handles.text_scrval, 'String', sprintf('%.2f', 100*scrVal));


%% end of function printInfo(handles)
%% ------------------------------------------------------------


%% ------------------------------------------------------------
%% --- Display the images and lines from scratch
%% ------------------------------------------------------------

function handles = newDisplay(handles)
% refactored to use single function, second arg=true means 'newDisplay'

handles = updateDisplay_newvol(handles, true);

%% end of function handles = newDisplay(handles)

%% ------------------------------------------------------------


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


%% ------------------------------------------------------------
%% --- Update display upon change in slice number
%% ------------------------------------------------------------

function handles = updateDisplay_ss(handles)

handles = updateDisplay_newvol(handles);


%% ------------------------------------------------------------


%% ------------------------------------------------------------
%% --- Update display upon change in column number
%% ------------------------------------------------------------

function handles = updateDisplay_cc(handles)

handles = updateDisplay_newvol(handles);


%% ------------------------------------------------------------


%% ------------------------------------------------------------
%% --- Update display upon change in row number
%% ------------------------------------------------------------

function handles = updateDisplay_rr(handles)

handles = updateDisplay_newvol(handles);


%% ------------------------------------------------------------

%% ------------------------------------------------------------
%% --- Compute scaling factors
%% ------------------------------------------------------------

function sf = calculateScale(handles, vol)

%calculate scaling factors, using correct voxel sizes from vol struct
% FIXME as currently written, image size is taken from first image, which is how the r,c,s axes are defined, so
% make scales relative to that.
vol1 = handles.vols{1};
if nargin < 2
  vol = handles.vols{handles.currentVol};
end
try
    sf.r = vol1.vx/vol.vx;
    sf.c = vol1.vy/vol.vy;
    sf.s = vol1.vz/vol.vz;
catch
    sf.r = 1;
    sf.c = 1;
    sf.s = 1;
end
sf.scale = mean([sf.r sf.c sf.s]);

return

% FIXME: what sets rrMax? rrMax gets unreliable with multi-scaled voxels
if isfield(handles.vols{handles.currentVol},'imgs')
  ima = handles.vols{handles.currentVol}.imgs;
  sf.r = size(ima,1)/handles.rrMax;
  sf.c = size(ima,2)/handles.ccMax;
  sf.s = size(ima,3)/handles.ssMax;
elseif isfield(handles.vols{handles.currentVol},'imgs1')
  ima = handles.vols{handles.currentVol}.imgs3{handles.ss};
  sf.r = size(ima,1)/handles.rrMax;
  sf.c = size(ima,2)/handles.ccMax;
  ima = handles.vols{handles.currentVol}.imgs2{handles.cc};
  sf.s = size(ima,2)/handles.ssMax;
else %simply an image volume
  ima = handles;
  sf.r = size(ima,1)/handles.rrMax;
  sf.c = size(ima,2)/handles.ccMax;
  sf.s = size(ima,3)/handles.ssMax;
end
sf.scale = mean([sf.r sf.c sf.s]);

%% ------------------------------------------------------------


%% ------------------------------------------------------------
%% --- Expand compressed FOD image
%% ------------------------------------------------------------
function img = expandFOD(a)
if isempty(a{3})
  img = zeros(a{1},a{2},3,'uint8');
else
  img = zeros(a{1}*a{2},3,'uint8');
  img(a{3},:)= a{4};
  img = reshape(img,a{1},a{2},3);
end

%% ------------------------------------------------------------

%% ------------------------------------------------------------
%% --- Compute new zoom values and update display
%% ------------------------------------------------------------

function handles = updateDisplay_zoom(handles) % Should be called upon initial display and change of volume

if handles.zoom < sqrt(eps)
  
  handles.rrLim = handles.rrLimOrig;
  handles.ccLim = handles.ccLimOrig;
  handles.ssLim = handles.ssLimOrig;
  
else
  
  zoom = 2^handles.zoom;
  halfRange_mm = handles.halfRange_mm/zoom;
  tmp = [-halfRange_mm +halfRange_mm];
  handles.rrLim = handles.rrZoom + tmp/handles.rcsScale(1);
  handles.ccLim = handles.ccZoom + tmp/handles.rcsScale(2);
  handles.ssLim = handles.ssZoom + tmp/handles.rcsScale(3);
  
end

%sf = calculateScale(handles);
sf.r=1; %we've made the axis always in mm (in updateDisplay_newVol) so scaling no longer needed
sf.c=1;
sf.s=1;

set(handles.axesCR, 'Xlim', sf.c * handles.ccLim);
set(handles.axesCR, 'Ylim', sf.r * handles.rrLim);
set(handles.axesSR, 'Xlim', sf.s * handles.ssLim);
set(handles.axesSR, 'Ylim', sf.r * handles.rrLim);
set(handles.axesCS, 'Xlim', sf.c * handles.ccLim);
set(handles.axesCS, 'Ylim', sf.s * handles.ssLim);

set(handles.CtextCR, 'Position', [sf.c*handles.cc       sf.r*handles.rrLim(1)]);
set(handles.RtextCR, 'Position', [sf.c*handles.ccLim(1) sf.r*handles.rr]);
set(handles.StextSR, 'Position', [sf.s*handles.ss       sf.r*handles.rrLim(1)]);
set(handles.RtextSR, 'Position', [sf.s*handles.ssLim(1) sf.r*handles.rr]);
set(handles.CtextCS, 'Position', [sf.c*handles.cc       sf.s*handles.ssLim(2)]);
set(handles.StextCS, 'Position', [sf.c*handles.ccLim(1) sf.s*handles.ss]);

handles.ClineCR.YData = sf.r*handles.rrLim;
handles.SlineSR.YData = sf.r*handles.rrLim;

addZoomedInset(handles, true); %turns off inset at higher zooms

handles = anatomyShowROIString(handles);

handles = showVolumeName(handles);

handles = annotationDraw(handles);

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

CLim1 = minI + (dVal-0.5/cVal)*(maxI-minI);
CLim2 = minI + (dVal+0.5/cVal)*(maxI-minI) + eps;

set(handles.axes1, 'CLim', [CLim1 CLim2]);
set(handles.axes2, 'CLim', [CLim1 CLim2]);
set(handles.axes3, 'CLim', [CLim1 CLim2]);

handles.CLim(handles.currentVol, :) = [CLim1 CLim2];

if isfield(handles.vols{handles.currentVol}, 'imgs')
  sf = calculateScale(handles);
  rr = scaleCoord(handles.rr,sf.r); cc = scaleCoord(handles.cc,sf.c); ss = scaleCoord(handles.ss,sf.s);
  voxVal = handles.vols{handles.currentVol}.imgs(rr,cc,ss,1);
  scrVal = (voxVal - CLim1)/(CLim2 - CLim1);
  scrVal = min(scrVal, 1);
  scrVal = max(scrVal, 0);
  set(handles.text_scrval, 'String', sprintf('%.2f', 100*scrVal));
end

%% end of function handles = updateDisplay_clims(handles)
%% ------------------------------------------------------------


%% ------------------------------------------------------------
%% --- add zoomed inset for FOD images
%% ------------------------------------------------------------

function [ima,handles] = addZoomedInset(handles, ima, r0, c0, scale, location)
% ima = addZoomedInset(handles, ima, r0, c0, scale, location)
%   show zoomed FOD image in main axis
%
% disables inset if current volume is not FOD
%
% addZoomedInset(handles, true) called after changing the zoom to
%   turn off when the FOD image itself is sufficiently zoomed in

handles.zoomedInsetAxes.Visible = 'off';

if gca ~= handles.axes1, return, end %only show in main axis; presumes caller has set gca to current axis (which updateDisplay_newVol does)

V = handles.vols{handles.currentVol};
isFOD = isfield(V,'imgs1');
zoomLimit = 2.8;

%turn off axes if this is not FOD, or if showing axis is turned off by user (` key toggles)
if ~isFOD || ~handles.showZoomedInset
  set(handles.zoomedInsetAxes.Children, 'Visible', 'off');
  return
end

% Usage 2: when zooming, if zoomed in, turn off above a certain zoom
if (numel(ima)==1 && ima==true)
  if handles.zoom > zoomLimit
    set(handles.zoomedInsetAxes.Children, 'Visible', 'off');
  else
    set(handles.zoomedInsetAxes.Children, 'Visible', 'on');
  end
  return
end

srcSize = round(5 * scale); %+/- 5mm around cursor
shift = floor(scale/2); %FIXME for 1/2mm shift of FOD images relative to other volumes
try
  subImage = ima(r0*scale + (-srcSize:srcSize)+shift, c0*scale + (-srcSize:srcSize)+shift , :);
  if strcmp(location,'ur') %axial slice, flipped up/down
    subImage = flipud(subImage);
  end
  if ~ishandle(handles.zoomedInsetImage)
    cla(handles.zoomedInsetAxes)
    handles.zoomedInsetImage = image(handles.zoomedInsetAxes, subImage);
    cx = mean(xlim(handles.zoomedInsetAxes));
    cy = mean(ylim(handles.zoomedInsetAxes));
    hold(handles.zoomedInsetAxes, 'on')
    plot(handles.zoomedInsetAxes, cx,cy,'w+')
  else
    handles.zoomedInsetImage.CData = subImage;
  end
  
  if handles.zoom > zoomLimit
    set(handles.zoomedInsetAxes.Children, 'Visible', 'off');
  else
    set(handles.zoomedInsetAxes.Children, 'Visible', 'on');
  end
  
catch
  
end

%% end of function ima = addZoomedInset(ima, r0, c0, show, scale)
%% ------------------------------------------------------------


%% ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
%% +++++++ Linking coordinates across windows +++++++
%% ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

%% ------------------------------------------------------------
%% --- notify other viewers of our LPH coords
%% ------------------------------------------------------------
function broadcastLPH(handles)

if ~isfield(handles,'doLinkCoordinates') || ~handles.doLinkCoordinates, return, end

showVolWins = findobj(get(0,'children'),'flat','Tag','showVol','-or','Tag','showVol_listener');
for win = showVolWins'
  if isequal(win,handles.figure), continue, end
  rcs = [handles.rr handles.cc handles.ss 1];
  M = handles.Mvxl2lph_atlas;
  lph = M * rcs(:);
  ud.lph = lph(1:3);
  ud.M = M;
  ud.sender = handles.figure;
  set(win, 'UserData',ud)
end

%% ------------------------------------------------------------
%% --- receive broadcasted LPH coords
%% ------------------------------------------------------------
function receiveLPH(~,event)

destFig = event.AffectedObject;
ud = destFig.UserData;
handles = guidata(destFig);
srcFig = ud.sender;

if ~isfield(handles,'doLinkCoordinates') || ~handles.doLinkCoordinates, return, end
if ~isfield(ud,'lph'), return, end %|| isequal(ud.rcs, [handles.rr handles.cc handles.ss])
M = handles.Mvxl2lph_atlas;
rcs = floor(inv(M) * [ud.lph(:); 1]);

handles.rr = min(max(1,rcs(1)), handles.rrMax);
handles.cc = min(max(1,rcs(2)), handles.ccMax);
handles.ss = min(max(1,rcs(3)), handles.ssMax);

handles.doLinkCoordinates = false; %temporarily turn off so we don't re-broadcast
set(0,'CurrentFigure',destFig)
handles = updateDisplay_newvol(handles);
handles.doLinkCoordinates = true;
guidata(destFig,handles)
set(0,'CurrentFigure',srcFig)

%% ------------------------------------------------------------
%% --- register to act on any changes to my UserData Property
%% ------------------------------------------------------------
function handles = registerCoordinateListener(handles)

handles.CoordinateListener = addlistener(handles.figure,'UserData', 'PostSet', @receiveLPH);
guidata(handles.figure, handles)

%% +++++++ end of functions related to broadcasting coords across viewers
%% ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++


%% ------------------------------------------------------------
%% --- Check dependencies
%% ------------------------------------------------------------

function handles = checkDependencies(handles)
%make sure utils are in path
handles.showVolPath = fileparts(mfilename('fullpath'));
addpath(fullfile(handles.showVolPath,'utils','')); %make sure utils in path

%make sure cmig_utils are present (public release will include, so should not be an issue)
hasCmigUtils = exist('atlas_T1.m','file'); %pick a file found only in cmig_utils
if ~hasCmigUtils
  error('You must also get cmig_utils from github (https://github.com/cmig-research-group/cmig_utils) and add to your path.') %FIXME: verify this is release version or is it included in cmig_tools?
end

guidata(handles.figure, handles)


%% ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
%% +++++++ ANNOTATION FUNCTIONS +++++++
%% ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

%% ------------------------------------------------------------
%% --- Annotation Initialization
%% ------------------------------------------------------------

function handles = annotationSetup(handles)

% annotations are ABCD-specific; Disable if using generic volumes
if ~handles.hasABCDBrain
  handles.annotation = [];
  guidata(handles.figure, handles)
  setAnnotationGUIVisibility(handles, false)
  return
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% If you're having problems with git, please set this variable to true to
% skip it so you can continue working until we figure it out
disableGit = false;

%When you are ready to debug, set this  true to print more info about git commands
handles.annotation.doDebugGit = false;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

annotationDir = fullfile(userhome,'cmig-annotations','');
if ~exist(annotationDir,'dir') %First time using: git hasn't been set up. Try to clone it.
  if ~disableGit
    %test if git is setup or not
    disp('Testing for github setup.')
    [s,o]=jsystem('ssh -T git@github.com');
  end
  if disableGit || s ~=1 %failed to ssh to github (or git checking is disabled)
    if ~disableGit
      fprintf(2,'Github SSH access may not be set up properly.\nPlease follow the instructions on github under the cmig-annotaitons repository here:\nhttps://github.com/cmig-research-group/cmig-annotations\n')
      fprintf('Suggested key generation commands:\n > ssh-keygen -t rsa -b 4096 -C "your_github_email@example.com" -f "%s/.ssh/github_id_rsa"\n\n',userhome);
      fprintf(' NOTES: suggested to not use a passphrase, and skip the ssh-agent and ssh-add steps. You''ll have to restart matlab as well.\n\n')
    else
      disp('Git is disabled')
    end
    disp('For now, it will store annotations locally for testing purposes, but these will not be merged into the group list.')
    annotationDir = fullfile(userhome,'cmig-annotations-local','');
    if ~exist(annotationDir,'dir'), mkdir(annotationDir); end
    handles.annotation.gitRepoDir = '';
  else
    disp('Cloning annotation repo from github')
    handles.annotation.gitRepoDir = userhome;
    git(handles, 'clone git@github.com:cmig-research-group/cmig-annotations');
    handles.annotation.gitRepoDir = annotationDir;
  end
else
  handles.annotation.gitRepoDir = annotationDir;
end
handles.annotation.useGithub = ~isempty(handles.annotation.gitRepoDir);

varnames = {'label', 'R',     'A',     'S',     'r',     'c',     's',     'abbrev','note',  'vol', 'category','extent', 'parent', 'author','date', 'uuid'};
vartypes = {'string','double','double','double','double','double','double','string','string','string','string','double','string', 'string','datetime','string'};

handles.annotation.filename = fullfile(annotationDir,'showVolAnnotations');

if exist([handles.annotation.filename '.csv'], 'file')
  [handles.annotation.A , handles.annotation.C ] = annotationLoad(handles);
else % FIXME: for debugging only
  handles.annotation.A = table('size',[0 length(varnames)],'VariableTypes',vartypes, 'variableNames', varnames);
  handles.annotation.C = cell2table({'' 1 '' ''},'VariableNames',{'category','colorIdxUnique','parent','uuid'});
end
handles.annotation.visible = false; %toggle these with 'Show' button 
handles.annotation.showAuthor = false;

%per-user hidden settings (these are not tracked)
handles.annotation.uuidStatusFilename = fullfile(annotationDir,['showVolAnnotations_UUID_status_' username '.mat' ]);
if exist(handles.annotation.uuidStatusFilename,'file')
  handles.annotation.uuid = load(handles.annotation.uuidStatusFilename);
else
  handles.annotation.uuid.hidden = {};  % uuid of annotations to hide
  handles.annotation.uuid.hiddenCategories = {};
  handles.annotation.uuid.unsynced = {};  % uuid of annotations added/changed since last sync
  handles.annotation.uuid.changed = {}; % our annotations changed on github
  handles.annotation.uuid.deleted = {}; % our annotions deleted on github
  handles.annotation.uuid.inView = {}; %keep track of on-screen visible handles
  uuid = handles.annotation.uuid;
  save(handles.annotation.uuidStatusFilename, '-struct', 'uuid')
end

handles.annotation.annotationUI = [];
handles.annotation.listener = addlistener(handles.figure,'UserData', 'PostSet', @annotationReceive);
handles.annotation.updateUI = true;

guidata(handles.figure, handles)

%% ------------------------------------------------------------


%% ------------------------------------------------------------
%% --- receive notification of annotationsUI update
%% ------------------------------------------------------------
function annotationReceive(~,event)

showVolFig = event.AffectedObject;
annotationFig = gcf;
handles = guidata(showVolFig);
ud = showVolFig.UserData;
if ~isfield(ud,'annotation') || ~handles.annotation.visible, return, end

switch ud.annotation.action
  case 'updateVisible'
    handles.annotation.uuid.hidden = ud.annotation.hidden;
    handles.annotation.uuid.hiddenCategories = ud.annotation.hiddenCategories;
    
    set(0,'CurrentFigure',showVolFig)
    handles = updateDisplay_newvol(handles);
    set(0,'CurrentFigure',annotationFig)
    
    uuid = handles.annotation.uuid;
    save(handles.annotation.uuidStatusFilename, '-struct', 'uuid')
    
    updateAnnotationUI(handles); %push change back to annotation UI
    
  case 'sync'
    handles = annotationGithubSync(handles);
    
    %     case 'edit'
    %         A = handles.annotation.A;
    %         sel = contains(A.uuid, ud.annotation.editUuid);
    %         if ~any(sel), return, end
    %         handles = annotationEdit(handles, A(sel,:));
    
  case 'delete'
    A = handles.annotation.A;
    Ad = ud.annotation.deleteRow;
    sel = contains(A.uuid, Ad.uuid) & contains(A.author, Ad.author);
    if ~any(sel), warning('update annotation not found: should be impossible'), return, end
    set(0,'CurrentFigure',showVolFig)
    handles = annotationDelete(handles, sel);
    set(0,'CurrentFigure',annotationFig)
    
    updateAnnotationUI(handles); %push change back to annotation UI
    
  case 'updateAnnotation'
    A = handles.annotation.A;
    C = handles.annotation.C;
    Au = ud.annotation.updateRow; %updated annotation, with category name as string
    sel = contains(A.uuid, Au.uuid) & contains(A.author, Au.author);
    if ~(sum(sel)==1), warning('update annotation not found: should be impossible. No change made.'), return, end
    
    %Au.author = {username}; %don't change author, but add a note if
    %modified by different author
    rightnow = datetime('now');
    if ~strcmp(Au.author, username)
      Au.note = {sprintf('%s; ; {modified by %s on %s}', Au.note{1}, username, rightnow)};
    end
    Au.date = rightnow;
    
    [handles.annotation.C, catuuid] = category2uuid(handles.annotation.C, Au.category);
    Au.category = catuuid;
    
    A(sel,:) = Au;
    handles.annotation.uuid.unsynced = union(handles.annotation.uuid.unsynced, Au.uuid);
    
    handles.annotation.A = A;
    set(0,'CurrentFigure',showVolFig)
    handles = updateDisplay_newvol(handles);
    handles = annotationSave(handles);
    set(0,'CurrentFigure',annotationFig)
    
    updateAnnotationUI(handles); %push change back to annotation UI
    
  case 'goto'
    rcs = ud.annotation.gotoRcs;
    handles.rr = rcs(1); handles.cc = rcs(2); handles.ss = rcs(3);
    %reset zoom center coordinates to new coordinate
    handles.rrZoom = handles.rr;
    handles.ccZoom = handles.cc;
    handles.ssZoom = handles.ss;
    figure(showVolFig) %bring, keep showVol in front
    handles = updateDisplay_newvol(handles);
    
  case 'undo'
    handles = annotationUndo(handles);
    
end

guidata(showVolFig,handles)

% ---------- update annotation UI (if it exists) ----------
function updateAnnotationUI(handles)

if ~isempty(handles.annotation) && handles.annotation.updateUI && ~isempty(handles.annotation.annotationUI) && ishandle(handles.annotation.annotationUI) % && strcmp(handles.annotation.annotationUI.Visible,'on')
  handles.annotation.annotationUI.UserData = handles;
end

% ---------- helper function for categories ----------
% get uuid for category, adding if necessary
function [C, catuuid] = category2uuid(C, category)
if isempty(category), catuuid = ''; return, end %handle 'no category' case
[~,int] = intersect(C.category, category);
if isempty(int) %new category?
  catuuid = {getuuid};
  cidx = max(C.colorIdxUnique) + 1; %keep a unique, increasing color index that stays with the category
  C = [C ; ...
    struct2table(struct('category',category, 'colorIdxUnique', cidx, 'parent','','uuid',catuuid ),'AsArray',true)];
else
  catuuid = C.uuid(int);
end


%% ------------------------------------------------------------
%% --- Add an Annotation
%% ------------------------------------------------------------
function handles = annotationAdd(handles, r, c, s, str)
% annotationAdd(handles, r, c, s, [str])
% str is {'label' 'abbreviation' 'notes'}. if no str, opens dialog to request
if isempty(handles.annotation) || ~handles.annotation.visible, return, end

if nargin<5, str = ''; end

%open a dialog for label and description if needed
if isempty(str)
  %grab a list of categories
  catStr = join(handles.annotation.C.category(2:end), '; ');
  catStr = ['Category: (e.g. ' catStr{1} ')'];
  a = inputdlg({'Label','Abbrev', 'Extent (+/- #voxels to display annotation) N or R A S (R/L; A/P; S/I)',...
    catStr, 'Notes'},'Add Annotation',[1 64; 1 6; 1 38; 1 38; 7 80],...
    [{''} {''}, num2str([1 1 1]), {''} {''}]);
  if isempty(a) || isempty(a{1}), return, end
  label = char(string(a{1}));
  abbrev = char(string(a{2}));
  if isempty(abbrev), abbrev = label(1:min(length(label),4)); end
  extent = str2num(a{3});
  if any(isnan(extent)) || ~(length(extent) == 1 || length(extent) == 3)
    fprintf(2,'bad value for extent (%s), using [1 1 1]\n',sprintf('%d ', extent))
    extent = [1 1 1];
  else
    if length(extent) == 1
      extent = extent * [1 1 1];
    end
  end
  category = char(string(a{4}));
  note = strjoin(cellstr(a{5}),'; ');
else
  label = str{1};
  abbrev = str{2};
  note = str{3};
end

[Rr,Aa,Ss] = rcs2ras(handles.Mvxl2lph_atlas, r, c, s);

vol = handles.volnames{handles.currentVol};
if strcmp(vol,'this volume'), vol = ''; end %don't include generic name

[handles.annotation.C, catuuid] = category2uuid(handles.annotation.C, category);

uuid = getuuid;
handles.annotation.A = [handles.annotation.A ; ...
  struct2table(struct('label',label, 'R',Rr, 'A',Aa, 'S',Ss,...
  'r',r, 'c',c, 's',s, 'abbrev',abbrev, 'note',note,...
  'vol',vol,'category',catuuid,'extent',extent, 'parent','',...
  'author',username, 'date',datetime('now'), 'uuid', uuid),'AsArray',true) ];

handles.annotation.uuid.unsynced = union(handles.annotation.uuid.unsynced, uuid);

handles.annotation.A = sortrows(handles.annotation.A,{'A','S','R'});
handles = annotationSave(handles);

handles = updateDisplay_newvol(handles);
updateAnnotationUI(handles)
guidata(handles.figure,handles)

%% ------------------------------------------------------------

%% ------------------------------------------------------------
%% --- Delete an Annotation
%% ------------------------------------------------------------
function handles = annotationDelete(handles, row)
%delete annotation at row
Ad = handles.annotation.A(row,:);
handles.annotation.uuid.unsynced = setdiff(handles.annotation.uuid.unsynced, Ad.uuid);
handles.annotation.uuid.changed = setdiff(handles.annotation.uuid.changed, Ad.uuid);
handles.annotation.uuid.deleted = setdiff(handles.annotation.uuid.deleted, Ad.uuid);
handles.annotation.uuid.hidden = setdiff(handles.annotation.uuid.hidden, Ad.uuid);
handles.annotation.A(row,:) = [];

handles = annotationSave(handles);

handles = updateDisplay_newvol(handles);
updateAnnotationUI(handles)
guidata(handles.figure,handles)

%% ------------------------------------------------------------
%% --- UI to Pick annotations for various actions -- simple UI
%% ------------------------------------------------------------

function handles = annotationPick(handles,action,A)
if nargin < 3
  A = handles.annotation.A;
end
C = handles.annotation.C;
%create listing string
A = sortrows(A, {'category', 'label'});
list = {};
for iA = 1:height(A)
  x = A(iA,:);
  category = C.category{strcmp(x.category{1}, C.uuid)};
  list{iA} = sprintf('%-20s:  %s (%s) [%s], [%d %d %d]',category,x.label{1},x.abbrev{1},x.author{1}(1:2),x.r,x.c,x.s);
end
sel = 'single';
[sel, ok] = listdlg('ListString',list, ...
  'Name', ['Annotations'], 'PromptString', ['Select annotation to ' action],...
  'ListSize', [750 750],'OKString', action,'selectionMode',sel);
if ~ok, return, end

switch lower(action)
  case 'go'
    x=A(sel(1),:);
    handles.rr = x.r; handles.cc = x.c; handles.ss = x.s;
    handles = updateDisplay_newvol(handles);
  case 'edit'
    x=A(sel(1),:);
    handles = annotationEdit(handles, x);
end
guidata(handles.figure, handles)

%% ------------------------------------------------------------

%% ------------------------------------------------------------
%% --- Edit annotations at point
%% ------------------------------------------------------------

function handles = annotationEdit(handles, r, c, s)
if isempty(handles.annotation) || ~handles.annotation.visible, return, end

if nargin == 2 %2 arg calling case, second arg is single table row
  A = r;
  for iA = 1:height(handles.annotation.A)
    if isequal(A, handles.annotation.A(iA,:))
      atCoord = iA;
      break
    end
  end
else
  A = handles.annotation.A;
  
  atCoord = find(annotationsAtPoint(A,r,c,s,0));
  
  %if edit at location with no annotations, put up the full list
  if isempty(atCoord)
    annotationPick(handles, 'Edit');
    return
  end
  % if > 1 at this location, put up a list to pick
  if length(atCoord)>1
    annotationPick(handles,'Edit',A(atCoord,:));
    return
  end
  
  A=handles.annotation.A(atCoord,:);
end

%format default answers
origCategory = handles.annotation.C.category{strcmp(A.category, handles.annotation.C.uuid)};

if ~isempty(A.note)
  note = strrep(A.note,'; ',newline);
else
  note = {''};
end

a = inputdlg({'Label  [If blank will delete this annotation]','Abbrev', 'Extent (+/- voxels) N or R A S (R/L; A/P; S/I)', 'Category', 'Notes'},...
  'Edit/Delete Annotation',[1 64; 1 6; 1 38; 1 38; 7 80],...
  [A.label, A.abbrev, num2str(A.extent), {origCategory}, note]);
if isempty(a), return, end
if isempty(a{1}) %empty label means delete this
  handles = annotationDelete(handles, atCoord);
else
  label = char(string(a{1}));
  abbrev = char(string(a{2}));
  if isempty(abbrev), abbrev = label(1:min(length(label),4)); end
  extent = str2num(a{3});
  if any(isnan(extent)) || ~(length(extent) == 1 || length(extent) == 3)
    fprintf(2,'bad value for extent (%s), using [1 1 1]\n',sprintf('%d ', extent))
    extent = [1 1 1];
  else
    if length(extent) == 1
      extent = extent * [1 1 1];
    end
  end
  category = char(string(a{4}));
  note = strjoin(cellstr(a{5}),'; ');
  
  %changed/new category?
  if ~isempty(category) && ~strcmp(category, origCategory)
    [handles.annotation.C, catuuid] = category2uuid(handles.annotation.C, category);
  else
    catuuid = A.category;
  end
  
  %note: edit keeps original user's name, keeps original uuid
  %Aadd a note if modified by different author
  if ~strcmp(A.author, username)
    note = sprintf('%s; ; {modified by %s}', note, username);
  end
  handles.annotation.A(atCoord,:) = struct2table(struct('label',label, 'R',A.R, 'A',A.A, 'S',A.S,...
    'r',A.r, 'c',A.c, 's',A.s, 'abbrev',abbrev, 'note',note,...
    'vol',A.vol,'category',catuuid,'extent',extent, 'parent','',...
    'author',A.author, 'date',datetime('now'), 'uuid', A.uuid),'AsArray',true);
  handles.annotation.uuid.unsynced = union(handles.annotation.uuid.unsynced, A.uuid);
end

handles = annotationSave(handles);

frontFig = gcf;
set(0,'CurrentFigure',handles.figure)
handles = updateDisplay_newvol(handles);
set(0,'CurrentFigure',frontFig)

guidata(handles.figure,handles)
updateAnnotationUI(handles);

%% ------------------------------------------------------------

function tf = annotationsAtPoint(A,rr,cc,ss,extentSwitch)
% find annotations in neighborhood of rcs points, taking into account each extent
% rr, cc, ss can be single points, or ranges defining a slice or volume
if nargin < 5
  extentSwitch = 1; %pass in 0 to turn off extents
end
tf = false(height(A),1);
for iA = 1:height(A)
  extent = A.extent(iA,:) * extentSwitch;
  re = [rr(1) rr(end)] + extent(3) * [-1 1]; % S/A
  ce = [cc(1) cc(end)] + extent(1) * [-1 1]; % R/L
  se = [ss(1) ss(end)] + extent(2) * [-1 1]; % A/P
  r = A.r(iA); c = A.c(iA); s = A.s(iA);
  tf(iA) = (r >= re(1) & r <= re(end) & ...
    c >= ce(1) & c <= ce(end) & ...
    s >= se(1) & s <= se(end) );
end

%% ------------------------------------------------------------
%% --- Get or display an annotations at point
%% ------------------------------------------------------------
function A = annotationGet(handles,r,c,s)

if isempty(handles.annotation) || ~handles.hasABCDBrain
  return
end

A = handles.annotation.A;
A = A(annotationsAtPoint(A,r,c,s),:);

if ~nargout && ~isempty(A)
  clc
  fprintf('\n===== Annotations at [%d, %d, %d]:\n', r, c, s)
  for iA = 1:height(A)
    if any(contains(handles.annotation.uuid.changed, A.uuid{iA}))
      if strcmp(A.author{iA}, username)
        str = ' [changed on github (my original)]';
      else
        str = ' [changed on github (new edit)]';
      end
      fid = 2;
    elseif any(contains(handles.annotation.uuid.deleted, A.uuid{iA}))
      str = ' [deleted on Github]';
      fid = 2;
    elseif any(contains(handles.annotation.uuid.hidden, A.uuid{iA}))
      str = ' [hidden]';
      fid = 2;
    else
      str = '';
      fid = 1;
    end
    note = strrep(A.note{iA},'; ', [newline ': ']);
    category = handles.annotation.C.category{strcmp(A.category{iA}, handles.annotation.C.uuid)};
    
    fprintf(fid, '------------%s\n:', str);
    fprintf(fid,' %s (''%s'')  [%s] %s\n:  %s (vol=%s)\n:\n: %s\n------------\n\n',...
      A.label{iA}, A.abbrev{iA}, A.author{iA}(1:2), A.date(iA), category,  A.vol{iA}, note);
  end
end

%% ------------------------------------------------------------


%% ------------------------------------------------------------
%% --- Draw all Annotations
%% ------------------------------------------------------------
function handles = annotationDraw(handles,scale)

if ~isfield(handles,'ORIENTATION'), return, end %this gets called once before init is complete

% draw labels on big axis, dits on smaller ones
delete(findobj(handles.figure,'tag','annotation'))

if isempty(handles.annotation) || ~handles.annotation.visible || isempty(handles.annotation.A), return, end

%find unhidden annotations
A = handles.annotation.A;
handles.annotation.uuid.inView = {}; %keep track of on-screen visible handles

for iView = 1:3 %loop through three views
  switch iView
    case 1
      drawax = handles.axesCR;
      rr = [1 handles.rrMax];
      cc = [1 handles.ccMax];
      ss = handles.ss;
      scale = size(handles.imCR.CData,1)/handles.rrMax;
    case 2
      drawax = handles.axesSR;
      rr = [1 handles.rrMax];
      cc = handles.cc;
      ss = [1 handles.ssMax];
      scale = size(handles.imSR.CData,1)/handles.rrMax;
    case 3
      drawax = handles.axesCS;
      rr = handles.rr;
      cc = [1 handles.ccMax];
      ss = [1 handles.ssMax];
      scale = size(handles.imCS.CData,2)/handles.ccMax;
  end
  doLabel = iView == handles.ORIENTATION;
  
  inView = annotationsAtPoint(A,rr,cc,ss);
  if ~any(inView), continue, end
  As = A(inView,:);
  
  handles.annotation.uuid.inView = union(handles.annotation.uuid.inView, As.uuid);
  
  [~,iVisible] = setdiff(As.uuid, handles.annotation.uuid.hidden);
  As = As(iVisible,:);
  
  authentic = annotationsAtPoint(As,rr,cc,ss,0); %points actually in this slice (not extended)
  
  %draw the annotations
  len = 3 * scale / 2^handles.zoom;
  inset = .25;
  nA = height(As);
  displayed = false(nA,1);
  
  % loop over annotations in view
  for iA = 1:nA
    
    if displayed(iA), continue, end %skip if label already was displayed
    
    r0 = As.r(iA); c0 = As.c(iA); s0 = As.s(iA);
    switch iView
      case 1 %CR
        x = scale*c0 + [inset  len];
        y = scale*r0 + [inset  len];
      case 2 %SR
        x = scale*s0 + [inset  len];
        y = scale*r0 + [inset  len];
      case 3 %CS
        x = scale*c0 + [inset  len];
        y = scale*s0 + [inset  -len];
    end
    
    isFOD = isfield(handles.vols{handles.currentVol},'imgs1');
    
    if doLabel && ~isFOD% DRAW LABEL in main axis
      %prepare label
      label = As.abbrev{iA};
      if isempty(label)
        label = As.label{iA}(1:min(length(label),4));
      end
      if handles.annotation.showAuthor, label = sprintf('{%s}_{%s}', label, As.author{iA}(1:2)); end
      %add any additional annotations at this same point to the label
      atCoord = find(annotationsAtPoint(As, r0, c0, s0) );
      atCoord(atCoord == iA) = [];
      if ~isempty(atCoord)
        for iC = atCoord'
          tmp = As.abbrev{iC};
          if isempty(tmp)
            tmp = As.label{iC}(1:min(length(label),4));
          end
          if handles.annotation.showAuthor, tmp = sprintf('{%s}_{%s}', tmp, As.author{iC}(1:2)); end
          label = [label ';' tmp];
          displayed(iC)=true;
        end
      end
      
      if authentic(iA)
        color = 'w';
      else
        color = [.65 .65 .65];
      end
      
      h = [ outlinetext(handles.axes1, x(2), y(2), label, color, 'fontsize',11,'fontweight','bold','horizontalalignment','left','tag','annotation') ...
        line(handles.axes1, x, y, 'color','k','linewidth',3,'tag','annotation')...
        line(handles.axes1, x, y, 'color',color,'linewidth',1,'tag','annotation') ];
      
      % TRY: add invisible ui element with tooltip containing description of label
      % cool: replace button with screengrab of image below it, making it invisible
      % But it is very slow, kills interactivity
      %     tipstr = sprintf('%s (%s)  [%s]\n\n%s',As.label, As.abbrev, As.author(1:2), As.note);
      %     pixpos = get(h(5),'position'); %text label
      %     ax0 = getpixelposition(handles.axes1,true);
      %     axpos = [pixpos(1) pixpos(2)-10 25 20];
      %     figpos = axpos + [ax0(1) ax0(2) 0 0];
      %     shot = getframe(handles.axes1, axpos);
      %     h2 = uicontrol(handles.axes1.Parent,'style','pushbutton','String','','foregroundColor',[1 0 0],...
      %         'Tooltip',tipstr,'visible','on','unit','pixels','position',figpos,'tag','annotation',...
      %         'CData',shot.cdata);
      
    else % DRAW A DOT in small axes
      atCoord = find( annotationsAtPoint(As, r0, c0, s0) );
      atCoord(atCoord == iA) = [];
      displayed(atCoord) = true; %so only draw the first marker at a point with multiple
      
      if authentic(iA)
        color = 'w';
      else
        color = [.75 .75 .75];
      end
      
      h =  [ line(drawax, [x(1) x(1)]-inset, [y(1) y(1)]-inset,'marker','o', ...
        'markerfacecolor','k','markeredgecolor','k', 'markersize',6,'tag','annotation'), ...
        line(drawax, [x(1) x(1)]-inset, [y(1) y(1)]-inset,'marker','o', ...
        'markerfacecolor',color,'markeredgecolor','none', 'markersize',3,'tag','annotation') ];
    end
  end %annotations in slice
end %view loop

% notify annotationUI of changes
updateAnnotationUI(handles)

%% ------------------------------------------------------------


%% ------------------------------------------------------------
%% --- Load Annotations
%% ------------------------------------------------------------
function [A,C,handles] = annotationLoad(handles)
A = readtable([handles.annotation.filename '.csv'], ...
  'Format','%q%f%f%f%f%f%f%q%q%q%q%f%f%f%q%q%D%q', ...  % from s, but color split into 3 vars
  'TextType','char','DatetimeType','datetime'); % AMD: this craps out in 2020a (Invalid format. error)
%recombine color
%A = mergevars(A,{'extent_1','extent_2','extent_3'},'NewVariableName','extent'); % AMD this doesn't work in 2017a (function mergevars missing)
A.extent = cat(2,A.extent_1,A.extent_2,A.extent_3); % Get rid of dependence on mergevars
C = readtable([handles.annotation.filename '_category.csv'],...
  'Format','%q%f%q%q', ...
  'TextType','char','DatetimeType','datetime');
if nargout > 2
  handles.annotation.uuid = load(handles.annotation.uuidStatusFilename);
end

%% ------------------------------------------------------------


%% ------------------------------------------------------------
%% --- Save Annotations
%% ------------------------------------------------------------

function handles = annotationSave(handles)
%make a backup for restore 'undo'
backupfiles(handles)
fbase =  handles.annotation.filename;
handles.annotation.A = sortrows(handles.annotation.A,{'A','S','R'});
writetable(handles.annotation.A, [fbase '.csv']);
handles.annotation.C = sortrows(handles.annotation.C,{'category'});
writetable(handles.annotation.C, [fbase '_category.csv']);
uuid = handles.annotation.uuid;
save(handles.annotation.uuidStatusFilename, '-struct', 'uuid')

%% ------------------------------------------------------------

%% ------------------------------------------------------------
%% --- Undo last change
%% ------------------------------------------------------------
function handles = annotationUndo(handles)
fbase =  handles.annotation.filename;
if ~exist([fbase '.csv.backup'],'file'), return, end
swapfile( [fbase '.csv.backup'],                                [fbase '.csv'] )
swapfile( [fbase '_category.csv.backup'],                       [fbase '_category.csv'] )
swapfile( [handles.annotation.uuidStatusFilename '.backup'],    handles.annotation.uuidStatusFilename )

[A , C, handles] = annotationLoad(handles);
handles.annotation.A = A;
handles.annotation.C = C;

handles = updateDisplay_newvol(handles);
updateAnnotationUI(handles);
guidata(handles.figure,handles)

%% ------------------------------------------------------------

function backupfiles(handles)
fbase =  handles.annotation.filename;
cmd{1} = sprintf('cp ''%s.csv'' ''%s.csv.backup''', fbase, fbase);
cmd{2} = sprintf('cp ''%s_category.csv'' ''%s_category.csv.backup''', fbase, fbase);
cmd{3} = sprintf('cp ''%s'' ''%s.backup''', handles.annotation.uuidStatusFilename, handles.annotation.uuidStatusFilename);
for i = 1:3
  [s,o] = jsystem(cmd{i});
  if s, warning('backup failed: %s\nerror: %s\n',cmd{i},o), end
end

function swapfile(A,B)
cmd{1} = sprintf('mv ''%s'' ''%s''', A,[B '___temp___']);
cmd{2} = sprintf('mv ''%s'' ''%s''', B,A );
cmd{3} = sprintf('mv ''%s'' ''%s''', [B '___temp___'],B);
for i = 1:3
  [s,o] = jsystem(cmd{i});
  if s, warning('swap failed: %s\nerror: %s\n',cmd{i},o), end
end


%% ------------------------------------------------------------
%% --- Sync Annotations with github
%% ------------------------------------------------------------
% considerations: multiuser is a little tricky: we want to merge edits across people in a way
%   that is different from the usual methods for code: If someone edits a
%   label I created, give me a choice to resolve it.
% similarly, if someone deletes a label I made, give me a chance to
% accept/deny that deletion. This assumes that any changes since last sync
% have not yet been comitted and all changes are unstaged.
function handles = annotationGithubSync(handles)
if ~handles.annotation.useGithub, return, end

disp('Sync with github...')

% check if any uncommited changes; commit
[s,o] = git(handles, 'diff --name-only --no-color');
if s, error('some problem with git diff: %s\n',o); end
if ~isempty(o) %we have changes to commit
  s = annotationCommitGit(handles);
  if s, fprintf(2,'Commit cancelled, so we can''t sync. Stopping.\n'); return, end
  handles.annotation.uuid.unsynced = {};
end

%store current commit, in case we need to recover
[s,recoverHash] = git(handles, 'rev-parse --verify HEAD');
if s, fprintf(2,'couldn''t look up current head, bailing.\n'); return, end

%check upstream
git(handles,'fetch')
[s,o] = git(handles,'status')
if contains(o,'have diverged') %we both have new commits, so pull&rebase then push
  try
    disp(' github is ahead of us. Getting updates before merging our recent changes and pushing.')
    [s,o] = git(handles, 'pull --rebase -s stop');
    if contains(o,'Successfully rebased and updated') %simple case, no merge needed
      
    else
      if (~ispc && s~=1) || ~contains(o,'git-merge-stop'), error('problem rebasing: %s\n',o), end %return codes not being set on PC
      disp(' merging')
      handles = annotationGithubMerge(handles);
      disp(' pushing merged annotations')
      git(handles, 'add .');
      o='';
      if ~contains(o,'Successfully rebased')
        [s,o] = git(handles, 'rebase --continue');
      end
      git(handles, 'push');
      handles.annotation.uuid.unsynced = {};
    end
    handles = annotationSave(handles);
  catch
    warning('There has been an error (%s). Restoring git to state prior to sync',lasterr)
    git(handles, 'rebase --abort');
    git(handles, 'reset --mixed', recoverHash);
    return
  end
elseif contains(o,'Your branch is behind') %github is ahead, but we have nothing new, so get it
  if ~isempty(handles.annotation.uuid.unsynced), error('unsynced changes'), end %FIXME: sanity check
  disp(' We have no changes, but there are new annotations in github. Merge')
  [s,o] = git(handles, 'pull');
  if s, error('some problem with git pull: %s\n',o); end
  try
    handles = annotationGithubMerge(handles);
  catch
    warning('There has been an error (%s). Restoring git to state prior to sync.',lasterr)
    git(handles, 'reset --hard', recoverHash);
    return
  end
  handles = annotationSave(handles);
elseif contains(o,'Your branch is ahead') %we have changes, but none on github: just push
  disp(' Pushing our new annotations to github')
  [s,o] = git(handles, 'push');
  if s || contains(o, 'error: failed to push'), error('problem pushing'), end
  
elseif contains(o,'up to date') || contains(o,'up-to-date')
  disp('  Nothing new. Annotations are up-to-date.')
end

disp('Sync complete.')

handles = updateDisplay_newvol(handles);
updateAnnotationUI(handles);
guidata(handles.figure, handles)


%% ------------------------------------------------------------ }


%% ------------------------------------------------------------
%% --- Merge our annotations with new version from Github with others' changes
%% ------------------------------------------------------------

function handles = annotationGithubMerge(handles)
%this is called at point where our changes are in handles.annotation.A,C
% and their modifications are now in the files
ourA = handles.annotation.A; ourC = handles.annotation.C;
[theirA, theirC] = annotationLoad(handles);

%annotations - somewhat complex merge keeping our original in case of their
%change or deletion so we may decide which one we prefer -- i.e. lets
%author override another's change or deletion of their own annotations
%unchangedA = intersect(ourA, theirA);
ourUnique = setdiff(ourA, theirA);
theirUnique = setdiff(theirA, ourA);

%keep older version of annotations we created and others changed so we can
%evaluate. OTOH if others changed others annotations, accept
[~, iOurChange, iTheirChange] = intersect(ourUnique.uuid, theirUnique.uuid);
%any changes originally not ours, don't keep our olde version
originallyOursChanged  =  contains(ourUnique.author(iOurChange),username); %we don't care about changes to annotations authored by others (that is for them)
changedByUs = contains(theirUnique.author(iTheirChange),username); %nor annotations made by us (can happen if working on two systems).
conflicts = originallyOursChanged & ~changedByUs;
ourChangedA = ourUnique(iOurChange(originallyOursChanged),:);

if any(conflicts)
  %add to labels, indicating change
  ourUnique.label(iOurChange(conflicts)) = ...
    cellfun(@(x) ['*M ' x], ourUnique.label(iOurChange(conflicts)), 'UniformOutput',false);
  ourUnique.abbrev(iOurChange(conflicts)) = ...
    cellfun(@(x) ['*M-' x], ourUnique.abbrev(iOurChange(conflicts)), 'UniformOutput',false);
  
  ourChangedA = ourUnique(iOurChange(conflicts),:);
  changingAuthors = theirUnique.author(iTheirChange(conflicts));
  fprintf('-- Merge with Github found conflicting changes --\n-   Changes: %d of your annotations were changed by another user:\n', sum(conflicts));
  disp(ourChangedA)
  fprintf('    changed by: %s\n\n',strjoin(changingAuthors, ', '));
  fprintf(2,'- Please check these and delete your original if you like the changes\n   and modify the new version as desired.\n After that, you''ll have to sync again to finalize the resolution.\n\n')
end
ourUnique(iOurChange(~conflicts),:) = []; %delete old version of annotations changed by their author

[uniqueMinusUnsynced, iSynced] = setdiff(ourUnique.uuid, handles.annotation.uuid.unsynced);
[~, iTheyDeleted] = setdiff(uniqueMinusUnsynced, theirUnique.uuid);
originallyOursDeleted = contains(ourUnique.author(iSynced(iTheyDeleted)), username); %any deletions that were originally ours, keep
ourDeletedA = ourUnique(iSynced(iTheyDeleted(originallyOursDeleted)),:);

if any(originallyOursDeleted)
  ourUnique.label(iSynced(iTheyDeleted(originallyOursDeleted))) = ...
    cellfun( @(x) ['*D ' x], ourUnique.label(iSynced(iTheyDeleted(originallyOursDeleted))), 'UniformOutput',false );
  ourUnique.abbrev(iSynced(iTheyDeleted(originallyOursDeleted))) = ...
    cellfun( @(x) ['*D-' x], ourUnique.abbrev(iSynced(iTheyDeleted(originallyOursDeleted))), 'UniformOutput',false );
  
  ourDeletedA = ourUnique(iSynced(iTheyDeleted(originallyOursDeleted)),:);
  fprintf('-   Deletions: %d of your annotations were deleted by another user:\n', sum(originallyOursDeleted));
  disp(ourDeletedA)
  fprintf(2,'- Please check these and delete your versions if you approve of their deletion,  and resync.\n\n')
end
ourUnique(iTheyDeleted(~originallyOursDeleted),:) = [];

handles.annotation.A = [theirA; ourUnique];
handles.annotation.A = sortrows(handles.annotation.A,{'A','S','R'});
handles.annotation.uuid.changed = ourChangedA.uuid;
handles.annotation.uuid.deleted = ourDeletedA.uuid;
handles.annotation.uuid.hidden = intersect(handles.annotation.uuid.hidden, handles.annotation.A.uuid);

%as a little extra backup, save changed/deleted
writetable([ourChangedA; ourDeletedA],fullfile(handles.annotation.gitRepoDir,'changes.csv'))

%categories
%one edge case--two people each create a category with the same name. add
%author qualifier to each to distinguish. It's a lot of complication.
[sameCat, iSameOur, iSameTheir] = intersect(ourC.category, theirC.category);
diffUuid = ~strcmp(ourC.uuid(iSameOur), theirC.uuid(iSameTheir));
sameCat = sameCat(diffUuid);
iSameOur = iSameOur(diffUuid);
iSameTheir = iSameTheir(diffUuid);

for iC = 1:length(iSameOur)
  author = username;
  ourC.category{iSameOur(iC)} = [sameCat{iC} ' (' author(1:2) ')'];
  sel = find(strcmp(theirC.uuid(iSameTheir(iC)), theirA.category));
  if ~isempty(sel)
    author = theirA.author{sel(1)};
  else
    author = '??';
  end
  theirC.category{iSameTheir(iC)} = [sameCat{iC} ' (' author(1:2) ')'];
end
if ~isempty(sameCat)
  fprintf('-   Category collisions: %d categories were added by both you and another user:\n', length(sameCat));
  disp( join( cat(2, ourC.category(iSameOur), theirC.category(iSameTheir)), ' <--> ', 2 ) )
  fprintf(2,'- You can ''merge'' these by hand by changing the categories for all affected annotations. If this ever comes upm, we can do something easier.\n\n')
end

%adjust colorIndex values to make unique
[theirNewC, iTheirNewC] = setdiff(theirC, ourC); %any new additions by others
[~,~,their] = intersect(ourC.colorIdxUnique, theirNewC.colorIdxUnique); %find conflicts in color index
usedIdx = unique(cat(1,ourC.colorIdxUnique, theirNewC.colorIdxUnique));
maxIdx = max(usedIdx);
unusedIdx = [ setdiff(1:maxIdx, usedIdx),  maxIdx+(1:length(their)) ];
for i = 1:length(their)
  theirC.colorIdxUnique(iTheirNewC(their(i))) = unusedIdx(i);
end

handles.annotation.C = union(ourC, theirC);
handles.annotation.C = sortrows(handles.annotation.C,{'category'});

%% ------------------------------------------------------------
%% --- Commit Annotations to local git repo
%% ------------------------------------------------------------

function s = annotationCommitGit(handles)
% commiting local changes is necessary prior to pushing or pulling
if ~handles.annotation.useGithub, return, end
a = inputdlg({'Message'},'Git commit message',[1 100]);
if isempty(a), s = 1; return, end
message = char(string(a{1}));
[s,o] = git(handles, 'commit -am', ['"' message '"']);
if s || contains(o,'fatal') , fprintf(2,'Commit failed due to: %s\n',o); end

%% ------------------------------------------------------------


%% ------------------------------------------------------------
%% --- Annotation Helper Functions
%% ------------------------------------------------------------

% FIXME: temporary: add correct transform for standard ABCD 2mm
% base it on Mvxl2lph for 1mm images (assuming L/R are same!)
function [Rr,Aa,Ss] = rcs2ras(Mvxl2lph, r, c, s)

lph = Mvxl2lph * [r(:)'; c(:)'; s(:)'; ones(1,length(r))];
RAS = M_LPH_TO_RAS * lph;
Rr = RAS(1,:)';
Aa = RAS(2,:)';
Ss = RAS(3,:)';


%% ------------------------------------------------------------

function uuid = getuuid
uuid = char(java.util.UUID.randomUUID);

%% ------------------------------------------------------------

function h = outlinetext(ax, x,y,str,color,varargin)
%plot white text with a black outline
d=2; %shadow displacement in pixels
h1 = text(ax, x,y,str,varargin{:},'visible','off');
set(h1,'units','pixels')
pixpos = get(h1,'position');
delete(h1);
sx = pixpos(1) + d*[-1 1 -1  1];
sy = pixpos(2) + d*[-1 1  1 -1];
h = [ text(ax, sx,sy,str, varargin{:}, 'color','k','units','pixels')' ...
  text(ax, pixpos(1),pixpos(2),str, varargin{:}, 'color',color,'units','pixels') ];

%% ------------------------------------------------------------

function [status, output] = git(handles, varargin)
% run git commands
if ~ispc
  [s,gitbin] = jsystem('which git');
  shell = '/bin/bash -c';
  cmd = strjoin( [ {'PATH=.:$PATH'}, {'TERM=ansi'}, {gitbin}, varargin ] );
else
  %[s,gitbin] = jsystem('where.exe git'); %FIXME totally untested on PC
  %gitbin = ['''' gitbin ''''];
  gitbin = 'git'; %cygwin installs it
  s=0;
  cmd = strjoin( [ {'set PATH=%cd%;C:\cygwin64\bin;%PATH%&&'}, {['START "" /B "' gitbin '"']}, varargin] );
  %this monstrosity gives a command like: set PATH=%cd%;C:\cygwin64\bin;%PATH%&& START "" /B "git" diff --name-only --no-color
  shell = 'cmd.exe /c';
end
if s, disp(gitbin), error('Couldn''t find git. Please make sure it is installed.'), end

[status, output] = jsystem(cmd, shell, handles.annotation.gitRepoDir);
output = deblank(regexprep(output, '(\x9B|\x1B\[)[0-?]*[ -/]*[@-~]', '')); %remove ansi codes
if handles.annotation.doDebugGit
  fprintf(2,'GIT: %s\nStatus: %d\nOutput: %s\n',cmd,status,output); % FIXME: for debugging only
end

if status && ~nargout %handle error
  error('Problem with git command: %s\nstatus = %d\noutput=\n%s\n\n',cmd, status,output)
end

%% end of annotation helper functions
%% ------------------------------------------------------------

%% +++++++ end of Annotation functions
%% ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

% convert rcs coords of one volume to another. Will need to round to be used as image indices
function [r2, c2, s2] = rcs2rcs(M1, M2, M2_size, r, c, s)

if isequal(M1, M2)
  r2 = r; c2 = c; s2 = s;
  return
end

%FIXME: is this really the simplest way?
vox1 = {r c s};
vox2 = {};
for i = 1:3
  vox = zeros(4,length(vox1{i}));
  vox(i,:) = vox1{i};
  vox(4,:) = 1;
  
  lph = M1 * vox;
  lph(4,:) = 1;
  
  vox = inv(M2)*lph;
  vox2{i} = vox(i,:);
end
[r2, c2, s2] = deal(vox2{:});

%clamp to edges of vol2
r2 = min(M2_size(1), max(1, r2));
c2 = min(M2_size(2), max(1, c2));
s2 = min(M2_size(3), max(1, s2));


%% ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
%% +++++++ ANATOMICAL ROI FUNCTIONS +++++++
%% ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

%% ------------------------------------------------------------
%% --- Anatomical ROI atlas Initialization
%% ------------------------------------------------------------


function handles = anatomySetup(handles)

% if no ABCD brain volulmes or invalid roifile, turn off ROI display
if ~handles.hasABCDBrain || isempty(handles.anat.roifile)
  handles.anat = [];
  handles.anat.params.showOverlay = true; %slightly hacky--this also controls showing the volume name
  set(handles.overlay_cb,'Value',false)
  guidata(handles.figure, handles)
  setAnatomyGUIVisibility(handles, false)
  return
end

%file pre-processed by prepareAtlases*.m. Specific atlases are expanded on
%demand in atlas_popup_Callback
atlasfile = handles.anat.roifile;
roiset = handles.anat.roiset;
roiatlas = handles.anat.roiatlas;
if isempty(roiatlas)
  [~,roifile] = fileparts(atlasfile); 
  atlasnum = regexp(roifile,'ABCD(\d)', 'tokens');
  atlasnum = str2double(atlasnum{1}{1});
  roiatlas = sprintf('ABCD%d',atlasnum);
else
  [~,roifile] = fileparts(atlasfile);
  atlasnum = regexp(roifile,'ABCD(\d)', 'tokens');
  atlasnum = str2double(atlasnum{1}{1});
end
  
handles.anat = load(atlasfile);
handles.anat.params.atlasfile = atlasfile;
handles.anat.params.roiset = roiset;
handles.anat.params.roiatlas = roiatlas;
handles.anat.params.roiatlasversion = atlasnum;

if ~isempty(roiset) && exist(roiset,'file')
  fprintf('%s: Loading ROI settings from %s\n',mfilename,roiset)
  handles = roi_load_button_Callback([], [], handles);
end

%set up the UI
alpha = handles.anat.params.overlayAlpha;
if alpha < handles.overlay_opacity_slider.Min, alpha = handles.overlay_opacity_slider.Min; end
set(handles.overlay_opacity_slider,'Value',alpha);
set(handles.image_fade_slider,'Value',handles.anat.params.imageFade);

atlasNames = cellfun(@(x) sprintf('%s[%d]',x,handles.anat.params.roiatlasversion), handles.anat.atlases(:)','UniformOutput',false);

set(handles.atlas_popup,'String',['-pick atlas' atlasNames],'Value',1);
handles.anat.params.showOverlay = true; %slightly hacky--this also controls showing the volume name
set(handles.overlay_cb,'Value',true)
handles.anat.params.isTracking = false; %use to speed up drawing when tracking mouse

%FIXME -- early atlas version may have forgot to subtract 10000 from fiber roicodes!
try
  if any(handles.anat.fiber.roicodes > 10000)
    handles.anat.fiber.roicodes = handles.anat.fiber.roicodes - 10000;
  end
catch
end

handles.anat.aseg.showNames = true;
handles.anat.aseg.prob = expandROI(handles.anat.aseg.prob);
handles.anat.fiber.showNames = true;
handles.anat.fiber.prob = expandROI(handles.anat.fiber.prob);
if isfield(handles.anat,'aparcaseg')
  handles.anat.aparcaseg.showNames = true;
  handles.anat.aparcaseg.prob = expandROI(handles.anat.aparcaseg.prob);
end

%% ------------------------------------------------------------

function prob = expandROI(prob)
if iscell(prob)
  c = prob;
  prob = zeros(c{1},'single');
  prob(c{2}) = single(c{3});
end

%% ------------------------------------------------------------
%% --- Anatomical ROI Lookup from RCS coordinates
%% ------------------------------------------------------------

function roi = anatomyFromRCS(handles, type)

if ~handles.hasABCDBrain
  return
end

sf = calculateScale(handles,handles.anat.(type));

prob = squeeze( handles.anat.(type).prob( scaleCoord(handles.rr, sf.r), scaleCoord(handles.cc, sf.c), scaleCoord(handles.ss, sf.s), :) );
[sp,si] = sort(prob,'descend');

switch type
  case 'aseg'
    probThreshold = 0.5;
  case 'fiber'
    probThreshold = 0.1; %fiber probs in general much lower
  case 'aparcaseg'
    probThreshold = 0.5;
end
show = find(sp >= probThreshold);
if isempty(show), show = find(sp >= probThreshold/5); end %FIXME - show at least one ROI as long as it has non-negligible probability
sp = sp(show);
si = si(show);

str = '';
if ~isempty(sp)
  for i = 1:length(sp)
    %FIXME: this has changed a lot
    try
      roiIdx   = find(handles.anat.(type).indlist(si(i)) == handles.anat.(type).roicodes); %as of May 2021, ABCD2
    catch
      try
      roiIdx   = find(handles.anat.(type).filist(si(i)) == handles.anat.(type).roicodes);
      catch
        roiIdx = si(i); %as of Sept 2021
      end
    end
    if isempty(roiIdx), continue, end
    roiProb(i)  = sp(i);
    roinames{i} = handles.anat.(type).roinames{roiIdx};
    roicode(i)  = handles.anat.(type).roicodes(roiIdx);
    %rgb(i,:)    = handles.anat.(type).roicolors(roiIdx,1:3);
    str{i}      = sprintf('  %.2f: %s',roiProb(i), roinames{i});
  end
end
if isempty(str)
  roiProb = 0;
  roinames{1} = 'none';
  roicode = nan;
  %rgb = [0 0 0];
  str = {'  none'};
end
roi.roinames = roinames;
roi.roicode = roicode;
roi.prob = roiProb;
%roi.rgb = rgb;
roi.str = str;

%% ------------------------------------------------------------
%% --- Display volume name on figure and title
%% ------------------------------------------------------------
function handles = showVolumeName(handles)
%this is called from updateDisplay_zoom

%get name, if present, and add to figure title
try
  name = handles.vols{handles.currentVol}.name;
  set(handles.figure,'name',['showVol - ' name])
catch
  if isfield(handles,'volName_text') && ishandle(handles.volName_text)
    set(handles.volName_text,'String','', 'visible','off')
  end
  set(handles.figure,'name','showVol')
  return
end

%show on figure when overlays are displayed (needed?)
if handles.anat.params.showOverlay && handles.showZoomedInset
  
  tmpax = gca;
  tmpf = gcf;
  set(handles.figure, 'CurrentAxes', handles.axes1);
  xl = get(handles.axes1,'XLim');
  yl = get(handles.axes1,'YLim');
  revy = strcmp(get(handles.axes1,'YDir'), 'reverse');
  %show in upper left of main axis
  xx = xl(1) + 0.05*range(xl);
  if revy
    yy = yl(1);
  else
    yy = yl(2);
  end
  
  if ~isfield(handles,'volName_text') || ~ishandle(handles.volName_text)
    handles.volName_text = text(xx,yy,name,'fontsize',12,'VerticalAlignment','top','color',[1 1 1],'interpreter','none','fontweight','bold');
  end
  set(handles.volName_text, 'String',name,'Position',[xx yy 0.1],'visible','on');
  
else
  if isfield(handles,'volName_text') && ishandle(handles.volName_text)
    set(handles.volName_text,'visible','off')
  end
end

%% ------------------------------------------------------------



%% ------------------------------------------------------------
%% --- Anatomical ROI Text Display
%% ------------------------------------------------------------
function handles = anatomyShowROIString(handles)
%this is called from updateDisplay_zoom

if ~handles.hasABCDBrain
  return
end

if ~handles.anat.params.showOverlay || ~( isfield(handles.anat,'aseg') || isfield(handles.anat,'fiber') || isfield(handles.anat,'aparcaseg'))  ...
    || ~handles.showZoomedInset
  try
    set(handles.anat.aseg.roi_text,'visible','off')
    set(handles.anat.fiber.roi_text,'visible','off')
    set(handles.anat.aparcaseg.roi_text,'visible','off')
  catch
  end
  return
end

tmpax = gca;
tmpf = gcf;
set(handles.figure, 'CurrentAxes', handles.axes1);
xl = get(handles.axes1,'XLim');
yl = get(handles.axes1,'YLim');
revy = strcmp(get(handles.axes1,'YDir'), 'reverse');
%show in upper right of main axis
xx = xl(1) + 0.6*range(xl);
xxc = xl(1) + 0.05*range(xl);
space = 0.1 * range(yl);
if revy
  yya = yl(1);
  yyf = yl(1) + space;
  yyc = yl(1) + 0.5*space;
else
  yya = yl(2);
  yyf = yl(2) - space;
  yyc = yl(2) - 0.5*space;
end

if isfield(handles.anat,'aseg') && handles.anat.aseg.showNames
  roi = anatomyFromRCS(handles, 'aseg');
  roiStr = join(roi.str,newline);
  roiStr = ['-ASEG-' newline roiStr{1}];
  if isempty(handles.anat.aseg.roi_text) || ~ishandle(handles.anat.aseg.roi_text)
    handles.anat.aseg.roi_text  = text(xx,yya,'aseg ROI','fontsize',10,'VerticalAlignment','top','color',[1 1 1],'interpreter','none','fontweight','bold');
  end
  set(handles.anat.aseg.roi_text, 'String',roiStr,'Position',[xx yya 0.1],'visible','on');
else
  try
    set(handles.anat.aseg.roi_text,'visible','off')
  catch
  end
end

if isfield(handles.anat,'fiber') && handles.anat.fiber.showNames
  roi = anatomyFromRCS(handles, 'fiber');
  roiStr = join(roi.str,newline);
  roiStr = ['-FIBER-' newline roiStr{1}];
  if isempty(handles.anat.fiber.roi_text) || ~ishandle(handles.anat.fiber.roi_text)
    handles.anat.fiber.roi_text = text(xx,yyf,'fiber ROI','fontsize',10,'VerticalAlignment','top','color',[1 1 1],'interpreter','none','fontweight','bold');
  end
  set(handles.anat.fiber.roi_text, 'String',roiStr,'Position',[xx yyf 0.1],'visible','on');
else
  try
    set(handles.anat.fiber.roi_text,'visible','off')
  catch
  end
end

if isfield(handles.anat,'aparcaseg') && handles.anat.aparcaseg.showNames
  roi = anatomyFromRCS(handles, 'aparcaseg');
  if roi.prob(1)>0
    isCtx = contains(roi.str,'ctx-'); %limit to cortical ROIs
    if ~any(isCtx)
      roi.str = {'none'};
    else
      roi.str = roi.str(isCtx);
    end
  end
  roiStr = join(roi.str, newline);
  roiStr = ['-APARC-' newline roiStr{1}];
  if isempty(handles.anat.aparcaseg.roi_text) || ~ishandle(handles.anat.aparcaseg.roi_text)
    handles.anat.aparcaseg.roi_text  = text(xxc,yyc,'aparc ROI','fontsize',10,'VerticalAlignment','top','color',[1 1 1],'interpreter','none','fontweight','bold');
  end
  set(handles.anat.aparcaseg.roi_text, 'String',roiStr,'Position',[xxc yyc 0.1],'visible','on');
else
  try
    set(handles.anat.aparcaseg.roi_text,'visible','off')
  catch
  end
end

guidata(handles.figure, handles)
set(tmpf, 'CurrentAxes', tmpax);

%% ------------------------------------------------------------


%% ------------------------------------------------------------
%% --- UI to Select ROIs for overlay
%% ------------------------------------------------------------

function handles = anatomyRoiOverlaySelectUI(handles,atlas)

if ~handles.hasABCDBrain
  return
end

[sel, ok] = listdlg('ListString',handles.anat.(atlas).uiNames, ...
  'Name', [atlas ' ROI'], 'PromptString', 'Select ROI(s) to overlay',...
  'InitialValue', handles.anat.(atlas).uiRoiOverlaySelected,...
  'ListSize', [250 750]);
if ~ok, return, end
if isempty(sel) || (length(sel)==1 && sel == 1) % 1 is the first 'None' entry
  handles.anat.(atlas).uiRoiOverlayIdx = [];
  handles.anat.(atlas).uiRoiOverlayImg = [];
  handles.anat.(atlas).uiRoiOverlaySelected = [];
else
  sel(sel==1)=[]; %remove 'None' entry
  handles.anat.(atlas).uiRoiOverlayIdx = handles.anat.(atlas).uiRoiIdx(sel);
  img = handles.anat.(atlas).prob(:,:,:,handles.anat.(atlas).uiRoiOverlayIdx);
  handles.anat.(atlas).uiRoiOverlayImg = img;
  handles.anat.(atlas).uiRoiOverlaySelected = sel;
  handles.anat.params.showOverlay = true;
  set(handles.overlay_cb,'Value',true)
end
guidata(handles.figure, handles)

%% ------------------------------------------------------------


%% ------------------------------------------------------------
%% --- add ROI overlays to image and/or compute outlines
%% ------------------------------------------------------------
function [ima, outline] = anatomyAddRoiOverlay(handles, ima, rr, cc, ss, mode)
% overlay probability for selected ROIs & calculate contour outlines
%mode is 'fill' (default) or 'outline' or 'filloutline' for both
% default is fill for lowres images and outline for FOD.

%return if nothing to do
if ~handles.hasABCDBrain || ~handles.anat.params.showOverlay
  outline = [];
  return
end

atlases = handles.anat.atlases;
shown = false;
for iA = 1:length(atlases)
  atlas = atlases{iA};
  if isfield(handles.anat,atlas) && ~isempty(handles.anat.(atlas).uiRoiOverlayImg)
    shown = true;
    break
  end
end
if ~shown
  outline = [];
  return
end

%add outline based on checkbox
if strcmp(mode,'fill') && handles.anat.params.showOutline
  mode = 'filloutline';
end

%for scalar volumes, convert to RGB if filling ROI, but preserve ability to adjust brightness and contrast
% FIXME: one inconvenience is that the contrast/darkness sliders don't change the image look when
% ROIs are displayed. Need to switch ROI off then back on to update it
if size(ima,3) == 1 && contains(mode,'fill')
  handles = updateDisplay_clims(handles);
  CLim = get(gca,'CLim');
  ima = (ima - CLim(1))/(CLim(2) - CLim(1));
  ima(isnan(ima)) = 0;
  ima = repmat(ima, 1, 1, 3);
end

isCoronal = length(ss)==1;
isSagittal = length(cc)==1;
isAxial = length(rr)==1; %x/y are transposed in axial view

outlineCoord = cell(0);
outlineColor = [];
outlineLabel = cell(0,3);
mask = [];
roiImg = [];

M_img = handles.vols{handles.currentVol}.Mvxl2lph;

%loop over atlases
for iA = 1:length(atlases)
  atlas = atlases{iA};
  if ~isfield(handles.anat, atlas), continue; end
  
  A = handles.anat.(atlas);
  if isempty(A.uiRoiOverlayImg); continue; end
  
  aimg = A.uiRoiOverlayImg;
  M_atlas = A.Mvxl2lph;
  
  %resample atlas into current image coordinates
  [ra, ca, sa] = rcs2rcs(M_img, M_atlas, A.size, rr, cc, ss); %image rcs to atlas rcs
  ra = floor(ra); ca = floor(ca); sa = floor(sa); % convert to matrix indices
  
  %collapse rois into single rgb image & mask
  if contains(mode, 'fill')
    
    overlay = squeeze(aimg(ra, ca, sa, :));
    if isAxial, overlay = permute(overlay, [2 1 3]); end
        
    if isempty(mask)
      mask = ( sum(overlay,3) > 0 );
      roiImg = zeros(size(overlay,1), size(overlay,2));
    else
      mask = (mask | sum(overlay,3) > 0);
    end
    for iR = 1:size(overlay,3)
      o = overlay(:,:,iR);
      if any(o(:))
        if handles.anat.params.binarizeOverlay
          o(o>A.probabilityThreshold) = 1;
        end
        c = A.roicolors(A.uiRoiOverlayIdx(iR),:);
        %c = A.roicolors(iR,:); %somehow this recreates colors in Anders 5/23/21 email ??
        m = max(c(:)); %FIXME: temporary fix, will be fixed permanentyly in prepareAtlases.m
        if m > 1
          c = c/255;
        end
        roiImg = roiImg + repmat(o,[1 1 3]) .* permute(c, [3 1 2]);
      end
    end
  end
  
  %make outlines
  if contains(mode, 'outline')
    
    %keep full atlas slice for outline calculation, at original scales
    if isCoronal
      overlayUnthresholded = squeeze(aimg(:, :, sa, :));
    elseif isSagittal
      overlayUnthresholded = squeeze(aimg(:, ca, :, :));
    elseif isAxial
      overlayUnthresholded = squeeze(aimg(ra, :, :, :));
      overlayUnthresholded = permute(overlayUnthresholded, [2 1 3]);
    end
    
    level = A.probabilityThreshold;
    for iR = 1:size(overlayUnthresholded,3)
      o = overlayUnthresholded(:,:,iR);
      if ~any(o(:)), continue, end
      cmat = contourc(double(o),[level level]);
      if isempty(cmat), continue, end
      cmat(:,1) = [];
      breaks = find(cmat(1,:)==level);
      cmat(:,breaks) = nan;
      
      % NOTE: We are using atlas coords for axis now, even for FODs so no
      % need to scale up outline coords any longer
      %convert cmat to coords of current image, from atlas coords
%       if isCoronal
%         [r2, c2, ~] = rcs2rcs(M_atlas, M_img, [size(ima) ss], cmat(1,:), cmat(2,:),ss);
%         cmat = [r2; c2];
%       elseif isSagittal
%         [r2, ~, s2] = rcs2rcs(M_atlas, M_img, [size(ima,1) cc size(ima,2)], cmat(1,:), cc, cmat(2,:));
%         cmat = [r2; s2];
%       elseif isAxial
%         [~, c2, s2] = rcs2rcs(M_atlas, M_img, [rr size(ima)], rr, cmat(1,:), cmat(2,:));
%         cmat = [c2; s2];
%       end
%       cmat(:,breaks) = nan;
      
      if isempty(breaks), breaks = size(cmat,2);end
      center = nanmean(cmat(:,1:breaks(1)-1),2);
      outlineCoord{end+1} = cmat;
      if contains(mode,'fill') && handles.anat.params.overlayAlpha>0
        c = [1 1 1]; %when showing colored roiprob, use white outline FIXME: make this UI-selectable
      else
        c = A.roicolors(A.uiRoiOverlayIdx(iR),:);
        m = max(c(:)); %FIXME: temporary fix, will be fixed permanentyly in prepareAtlases.m
        if m > 1
          c = c/255;
        end
        c = min(1, [.75 .75 .75] + c); %desaturate towards white
      end
      outlineColor = cat(1, outlineColor, c);
      
      %if we have an abbreviation, add it as a label if showNames is true
      if A.showNames && ~isempty(A.roiabbreviations)
        label = A.roiabbreviations(A.uiRoiOverlayIdx(iR));
      else
        label = '';
      end
      outlineLabel(end+1,:) = {center(1) center(2) label};
      
    end %loop over ROIs
  end %if outline
  
end %loop over atlas

%this was nice, but not sure needed any more
%         %outline _only_ add a second contour at 50%, useful in FOD display,
%         outline50 = [];
%         if ~contains(mode, 'fill') && handles.anat.fiberProbabilityThreshold < 0.4
%             level = 0.5;
%             cmat = [];
%             for iR = 1:size(overlay,3)
%                 cmat = cat(2, cmat, contourc(double(overlayUnthresholded(:,:,iR)),[level level]) );
%             end
%             cmat(:,cmat(1,:)==level) = nan;
%             outline50 = cat(2, outline50, cmat*a2i);
%         end
%         if ~isempty(outline50), outline = {outline outline50};end


%combine original image (possibly faded) with alpha blended roi image
if contains(mode,'fill')
  fade = handles.anat.params.imageFade;
  if isinteger(ima), ima = cast(ima,'single')/255; end
  % experimental: want low probs to be more transparent than higher ones, so that fringes turn more transparent
  fadeProb = 0.15;
  lum = mean(roiImg,3);
  inFringe = (lum<fadeProb);
  shadedMask = cast(mask,'single');
  shadedMask(inFringe) = lum(inFringe)/fadeProb;
  
  %ima = fade * ima.*~mask + (1-handles.anat.params.overlayAlpha) * ima.*mask + handles.anat.params.overlayAlpha * roiImg.*shadedMask;
  ima = ima.*(1-handles.anat.params.overlayAlpha*shadedMask/4) + handles.anat.params.overlayAlpha * roiImg.*shadedMask; %try new method, get rid of 'fade' param
  %ima = ima + handles.anat.params.overlayAlpha * roiImg.*shadedMask; %try new method, get rid of 'fade' param
  ima = min(ima,1); %clamp to 1
end
if contains(mode,'outline')
  outline = {outlineCoord outlineColor outlineLabel};
else
  outline = [];
end

%% ------------------------------------------------------------


%% ------------------------------------------------------------
%% --- Draw roi outlines
%% ------------------------------------------------------------
function handles = anatomyDrawRoiOutline(handles,outline,slice)

if ~handles.hasABCDBrain
  return
end

%axis is the axis slice we're drawing into
fname = ['roiContour' slice];
if ~isfield(handles, fname), return; end
delete(handles.(fname))
handles.(fname) = [];

axname = ['axes' slice];
delete(findobj(handles.(axname),'tag','roiOutline'))
delete(findobj(handles.(axname),'tag','roiLabel'))

if isempty(outline), return, end

doLabel = (handles.(axname) == handles.axes1); %only show labels in main axis

lw = 1.25;
[outlineCoord, outlineColor, outlineLabel] = deal(outline{:});
scalarVol = isfield(handles.vols{handles.currentVol}, 'imgs') && size(handles.vols{handles.currentVol}.imgs,4)==1 ...
  || (isfield(handles.vols{handles.currentVol}, 'limits') && length(handles.vols{handles.currentVol}.limits) > 2); %

for iR = 1:length(outlineCoord)
  if  scalarVol && ~handles.anat.params.isTracking
    h1 = line(outlineCoord{iR}(1,:), outlineCoord{iR}(2,:), 'color', 'k', 'linewidth', lw*2,'tag','roiOutline'); %black background for scalar volumes only; don't draw when tracking, for speed
  else
    h1 = [];
  end
  h2 = line(outlineCoord{iR}(1,:), outlineCoord{iR}(2,:), 'color', outlineColor(iR,:), 'linewidth', lw,'tag','roiOutline');
  label = outlineLabel{iR,3};
  if ~isempty(label) && ~handles.anat.params.isTracking %&& doLabel %don't draw labels when tracking
    h3 = text(outlineLabel{iR,1}, outlineLabel{iR,2}, label, 'color','w', 'fontsize',7,'fontweight','bold',...
      'horizontalalignment','center','verticalalignment','middle','tag','roiLabel');
  else
    h3=[];
  end
  handles.(fname) = [handles.(fname) h1 h2 h3];
end
guidata(handles.figure, handles);

% this old method was efficient, but can't handle lines with different colors
% if isempty(handles.(fname)) || any(~ishandle(handles.(fname)))
%     delete(handles.(fname))
%     handles.(fname)(1) = line(1,1,'color','k','linewidth',lw*2,'visible','off');
%     handles.(fname)(2) = line(1,1,'color','w','linewidth',lw,'visible','off');
%     handles.(fname)(3) = line(1,1,'color','w','linewidth',lw/2,'linestyle',':','visible','off');
% end
% if isempty(outline)
%     set(handles.(fname),'Visible','off')
% else
%     if iscell(outline)
%         set(handles.(fname)(1), 'XData', outline{1}(1,:), 'YData', outline{1}(2,:),'visible','on')
%         set(handles.(fname)(2), 'XData', outline{1}(1,:), 'YData', outline{1}(2,:),'visible','on')
%         set(handles.(fname)(3), 'XData', outline{2}(1,:), 'YData', outline{2}(2,:),'visible','on')
%     else
%         set(handles.(fname)(1), 'XData', outline(1,:), 'YData', outline(2,:),'visible','on')
%         set(handles.(fname)(2), 'XData', outline(1,:), 'YData', outline(2,:),'visible','on')
%         set(handles.(fname)(3), 'visible','off')
%     end
% end

%% ------------------------------------------------------------

%% +++++++ end of Anatomical ROI functions
%% ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

% r,c,s coords of UI remain in scale of first volume passed in; use this to convert to
% current volume's coords
function x = scaleCoord(x,sf)
%x = floor(sf*(x-1)+1); %this is for one-based indexing
x = floor(sf*x); %but new 1mm space is centered as if 0-based indexing


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
  
  case 1 %coronal main
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
    
  case 2 %sagittal main
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
    
  case 3 %axial main
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

%handles = updateDisplay_newvol(handles);

%handles = updateDisplay_zoom(handles);

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
handles = updateDisplay_newvol(handles);
%handles = updateDisplay_rr(handles);
%handles = updateDisplay_cc(handles);
%handles = updateDisplay_ss(handles);
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
handles = updateDisplay_newvol(handles);
%handles = updateDisplay_rr(handles);
%handles = updateDisplay_cc(handles);
%handles = updateDisplay_ss(handles);
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
  handles = updateDisplay_newvol(handles);
  broadcastLPH(handles)
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
Mlph2vxl = inv(handles.Mvxl2lph_atlas);
targetRCS = Mlph2vxl * targetLPH;
handles = gotoPointRCS(handles, targetRCS);

%% ------------------------------------------------------------


%% ------------------------------------------------------------
%% --- Decrement row number
%% ------------------------------------------------------------

function handles = decrementRowNumber(handles)
if handles.rr > 1
  handles.rr = handles.rr - 1;
end
handles = updateDisplay_rr(handles);
broadcastLPH(handles)
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
broadcastLPH(handles)
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
broadcastLPH(handles)
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
broadcastLPH(handles)
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
broadcastLPH(handles)
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
broadcastLPH(handles)
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
if isnumeric(hObject)
  [minVal, chosenVol] = min(abs(handles.volPB - hObject));
else
  [minVal, chosenVol] = min(cellfun(@(x)norm(x-hObject.Position),{handles.volPB(:).Position}'));
end
if chosenVol == handles.currentVol
  if sum(handles.hideVol) < handles.numVols-1
    handles.hideVol(chosenVol) = 1;
    set(hObject, 'ForegroundColor', [0.5 0.5 0.5]);
    set(hObject, 'TooltipString', ['Show ' handles.volnames{handles.currentVol}]);
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
  set(handles.currentVolPB, 'TooltipString', ['Show ' handles.volnames{handles.currentVol}]);
  handles.hideVol(chosenVol) = 0;
  handles.currentVol = chosenVol;
  handles.currentVolPB = hObject;
end
set(handles.currentVolPB, 'ForegroundColor', 'r');
set(handles.currentVolPB, 'TooltipString', ['Hide ' handles.volnames{handles.currentVol}]);
handles.Mvxl2lph = handles.vols{handles.currentVol}.Mvxl2lph;
handles = updateDisplay_newvol(handles);
guidata(hObject, handles);

%% ------------------------------------------------------------


%% ------------------------------------------------------------
%% --- Executes on button press in the mesh pushbuttons
%% ------------------------------------------------------------

function meshPB_Callback(hObject, eventdata, handles)
[minVal, chosenSurf] = min(abs(handles.meshPB - hObject));
handles.showSurf(chosenSurf) = 1 - handles.showSurf(chosenSurf);
if handles.showSurf(chosenSurf)
  set(handles.meshCrossSectionCR(chosenSurf).lh, 'Visible', 'on');
  set(handles.meshCrossSectionSR(chosenSurf).lh, 'Visible', 'on');
  set(handles.meshCrossSectionCS(chosenSurf).lh, 'Visible', 'on');
  set(hObject, 'ForegroundColor', handles.surfs(chosenSurf).color);
  set(hObject, 'TooltipString', 'Hide this mesh');
else
  set(handles.meshCrossSectionCR(chosenSurf).lh, 'Visible', 'off');
  set(handles.meshCrossSectionSR(chosenSurf).lh, 'Visible', 'off');
  set(handles.meshCrossSectionCS(chosenSurf).lh, 'Visible', 'off');
  set(hObject, 'ForegroundColor', [0.5 0.5 0.5]);
  set(hObject, 'TooltipString', 'Show this mesh');
end
guidata(hObject, handles);

%% ------------------------------------------------------------


%% ------------------------------------------------------------
%% --- Executes on button press in zoom_reset.
%% ------------------------------------------------------------

function zoom_reset_Callback(hObject, eventdata, handles)
defaultZoom = 0.48; %was 0, but this fills out the window better
handles.zoom = defaultZoom;
set(handles.zoom_slider, 'Value', defaultZoom);
%recenter as well
handles.rrZoom = ceil(handles.rrMax/2 + sqrt(eps));
handles.ccZoom = ceil(handles.ccMax/2 + sqrt(eps));
handles.ssZoom = ceil(handles.ssMax/2 + sqrt(eps));
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
set(handles.contrast_slider, 'Value', 1);
handles.cVal(handles.currentVol) = 1;
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
set(handles.darken_slider, 'Value', 0.5);
handles.dVal(handles.currentVol) = 0.5;
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

tracking = ~isempty(hObject.WindowButtonMotionFcn); %not in tracking currently

% ABCD-specific optimizations
if handles.hasABCDBrain
  if tracking
    handles.annotation.updateUI = false; %disable updates of annotation table for speed, reduce drawing and don't transmot coordinates to other windows
    handles.anat.params.isTracking = true;
  else
    handles.annotation.updateUI = true;
    handles.anat.params.isTracking = false;
  end
end

pos_axCR = get(handles.axesCR,'CurrentPoint');
XLim_axCR = get(handles.axesCR, 'XLim');
YLim_axCR = get(handles.axesCR, 'YLim');

pos_axSR = get(handles.axesSR,'CurrentPoint');
XLim_axSR = get(handles.axesSR, 'XLim');
YLim_axSR = get(handles.axesSR, 'YLim');

pos_axCS = get(handles.axesCS,'CurrentPoint');
XLim_axCS = get(handles.axesCS, 'XLim');
YLim_axCS = get(handles.axesCS, 'YLim');

inAxis = false; %was the cursor within an image axis?
if ((pos_axCR(1,1) >= XLim_axCR(1)) && ...
    (pos_axCR(1,1) <= XLim_axCR(2)) && ...
    (pos_axCR(1,2) >= YLim_axCR(1)) && ...
    (pos_axCR(1,2) <= YLim_axCR(2)))
  axesCR_ButtonDownFcn(handles.axesCR, eventdata, handles);
  inAxis = true;
elseif ((pos_axSR(1,1) >= XLim_axSR(1)) && ...
    (pos_axSR(1,1) <= XLim_axSR(2)) && ...
    (pos_axSR(1,2) >= YLim_axSR(1)) && ...
    (pos_axSR(1,2) <= YLim_axSR(2)))
  axesSR_ButtonDownFcn(handles.axesSR, eventdata, handles);
  inAxis = true;
elseif ((pos_axCS(1,1) >= XLim_axCS(1)) && ...
    (pos_axCS(1,1) <= XLim_axCS(2)) && ...
    (pos_axCS(1,2) >= YLim_axCS(1)) && ...
    (pos_axCS(1,2) <= YLim_axCS(2)))
  axesCS_ButtonDownFcn(handles.axesCS, eventdata, handles);
  inAxis = true;
end
handles = guidata(handles.figure);
%broadcastLPH(handles) %don't broadcast coords while tracking

if inAxis %enable mouse tracking
  if ~tracking %not already tracking, so begin tracking
    if handles.hasABCDBrain
      handles.annotation.updateUI = true;
      handles.anat.params.isTracking = false;
    end
    hObject.Pointer = 'circle';
    hObject.WindowButtonMotionFcn = {@figure_WindowButtonMotionFcn, handles};
    hObject.WindowButtonUpFcn = {@figure_WindowButtonUpFcn, handles};
    
    annotationGet(handles, handles.rr, handles.cc, handles.ss); %show any annotations present at clicked location
    
  end
else %out of axis: cancel tracking (as if button released)
  figure_WindowButtonUpFcn(hObject, eventdata, handles)
end

%% end of function figure_WindowButtonDownFcn(hObject, eventdata, handles)
%% ------------------------------------------------------------

%% ------------------------------------------------------------
%% --- Executes on mouse button move within figure window
%% ------------------------------------------------------------

% During axis drag use button down function to continually update the display
function figure_WindowButtonMotionFcn(hObject, eventdata, handles)
figure_WindowButtonDownFcn(hObject, eventdata, handles)

%% ------------------------------------------------------------


%% ------------------------------------------------------------
%% --- Executes on mouse button release in figure window
%% ------------------------------------------------------------

% Stop tracking on button up, or if cursor exits axis
function figure_WindowButtonUpFcn(hObject, eventdata, handles)
hObject.WindowButtonMotionFcn = '';
hObject.WindowButtonUpFcn = '';
hObject.Pointer = 'arrow';
handles = guidata(handles.figure);
if handles.hasABCDBrain
  handles.annotation.updateUI = true;
  handles.anat.params.isTracking = false;
end
handles = updateDisplay_newvol(handles);
broadcastLPH(handles)
guidata(hObject, handles);

%% ------------------------------------------------------------

%% ------------------------------------------------------------
%% --- Executes on mouse press over axesCR background.
%% ------------------------------------------------------------

function axesCR_ButtonDownFcn(hObject, eventdata, handles)
pos = get(hObject, 'CurrentPoint');
ima = get(handles.imCR,'CData');
%handles.cc = round(interp1([0 size(ima,2)+1],[0 handles.ccMax+1],pos(1,1),'linear','extrap'));
handles.cc = round(pos(1,1)); %Note: now using image XData to place all on same 1mm scale, so no scaling from image to world needed here any longer
handles.cc = max(handles.cc, 1);
handles.cc = min(handles.cc, handles.ccMax);
%handles.rr = round(interp1([0 size(ima,1)+1],[0 handles.rrMax+1],pos(1,2),'linear','extrap'));
handles.rr = round(pos(1,2));
handles.rr = max(handles.rr, 1);
handles.rr = min(handles.rr, handles.rrMax);
handles = updateDisplay_newvol(handles);
printInfo(handles);
guidata(hObject, handles);

%% end of function axesCR_ButtonDownFcn(hObject, eventdata, handles)
%% ------------------------------------------------------------


%% ------------------------------------------------------------
%% --- Executes on mouse press over axes3 background.
%% ------------------------------------------------------------

function axesSR_ButtonDownFcn(hObject, eventdata, handles)
pos = get(hObject, 'CurrentPoint');
ima = get(handles.imSR,'CData');
%handles.ss = round(interp1([0 size(ima,2)+1],[0 handles.ssMax+1],pos(1,1),'linear','extrap'));
handles.ss = round(pos(1,1));
handles.ss = max(handles.ss, 1);
handles.ss = min(handles.ss, handles.ssMax);
%handles.rr = round(interp1([0 size(ima,1)+1],[0 handles.rrMax+1],pos(1,2),'linear','extrap'));
handles.rr = round(pos(1,2));
handles.rr = max(handles.rr, 1);
handles.rr = min(handles.rr, handles.rrMax);
handles = updateDisplay_newvol(handles);
printInfo(handles);
guidata(hObject, handles);

%% end of function axesSR_ButtonDownFcn(hObject, eventdata, handles)
%% ------------------------------------------------------------


%% ------------------------------------------------------------
%% --- Executes on mouse press over axes background.
%% ------------------------------------------------------------

function axesCS_ButtonDownFcn(hObject, eventdata, handles)
pos = get(hObject, 'CurrentPoint');
ima = get(handles.imCS,'CData');
%handles.cc = round(interp1([0 size(ima,2)+1],[0 handles.ccMax+1],pos(1,1),'linear','extrap'));
handles.cc = round(pos(1,1));
handles.cc = max(handles.cc, 1);
handles.cc = min(handles.cc, handles.ccMax);
%handles.ss = round(interp1([0 size(ima,1)+1],[0 handles.ssMax+1],pos(1,2),'linear','extrap'));
handles.ss = round(pos(1,2));
handles.ss = max(handles.ss, 1);
handles.ss = min(handles.ss, handles.ssMax);
handles = updateDisplay_newvol(handles);
printInfo(handles);
guidata(hObject, handles);

%% function axesCS_ButtonDownFcn(hObject, eventdata, handles)
%% ------------------------------------------------------------


%% ------------------------------------------------------------
%% --- Executes on key press anywhwere in figure window
%% ------------------------------------------------------------

function figure_KeyPressFcn(hObject, eventdata, handles)
key = lower(get(hObject, 'CurrentCharacter'));
if isempty(key), return, end
mod = get(hObject, 'CurrentModifier');
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
    
  case '!' % AMD: Save snapshot of main window
    shot = getframe(handles.axes1);
    if ~isfield(handles,'screenshotnum')
      handles.screenshotnum = 0;
    end
    handles.screenshotnum = handles.screenshotnum+1;
    fname_out = sprintf('ScreenShot%03d.png',handles.screenshotnum);
    imwrite(shot.cdata,fname_out); 
    fprintf(1,'File %s written.\n',fname_out);

  case {'`','~'} %toggle zoomed FOD inset (& other UI acoutrements)
    handles.showZoomedInset = ~handles.showZoomedInset;
    handles = updateDisplay_newvol(handles);
    
  case '\' %toggle coordinate linking between showVol windows
    handles.doLinkCoordinates = ~handles.doLinkCoordinates;
    if(handles.doLinkCoordinates), disp('Coordinate Linking ON'), else, disp('Coordinate Linking OFF'),end
end

% Test for additional ABCD-specific annotation optimization
if handles.doAnnotations
    switch key
        case {'=' '+'} %add an annotation at current cursor
            handles = annotationAdd(handles, handles.rr, handles.cc, handles.ss);
            
        case '-' %edit annotation
            handles = annotationEdit(handles, handles.rr, handles.cc, handles.ss);
            
        case char(31) %revert last add/delete (control shift '-')
            if any(contains(mod,'control'))
            disp('Restoring saved annotations.')
            handles = annotationUndo(handles);
            handles = updateDisplay_newvol(handles);
            end
            
        case ']' %FIXME Debug-only to dump annotation tables
            clc
            disp('================= Annotations =======================')
            disp(handles.annotation.A)
            disp('================= Categories ========================')
            disp(handles.annotation.C)
            disp('================= UUID ==============================')
            disp(handles.annotation.uuid)
            disp('=====================================================')
            
        case '}' % FIXME debug only to reload annotations (e.g. after hand-editing files)
            disp('reloading annotations')
            [A , C, handles] = annotationLoad(handles);
            
            handles.annotation.A = A; %TODO clean up this calling
            handles.annotation.C = C;
            handles = updateDisplay_newvol(handles);
            
        case '"' %FIXME git Debug-only
            clc
            disp('===== git status =====')
            [s,o] = git(handles, 'status');
            fprintf('git status:\nStatus: %d\nOutput: %s\n',s,o);
            
            [s,o] = git(handles, 'diff');
            fprintf('git diff:\nStatus: %d\nOutput: %s\n',s,o);
    end
end % ABCD annotation-specific keys

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
rB = 4;       % right border

tAx123 = 1;   % vertical space above axes1, axes2 and axes3
rAx1 = 7;     % horizontal space right of axes1
rAx23 = 3;    % horizontal space right of axes2 and axes3
sAx2Ax3 = 2;  % vertical space between axes2 and axes3

wPBt = 7;     % width of topmost panel pushbuttons
hPBt = 2;     % height of topmost panel pushbuttons
sPBt = 1;     % horizontal space between topmost panel pushbuttons

wPBr = 10;    % width of rightmost panel pushbuttons
hPBr = 2;     % height of rightmost panel pushbuttons
sPBr = 0;     % vertical space between rightmost panel pushbuttons

ax1MinH = 45; % don't let ax1 height get smaller than this

if (isfield(handles, 'numSurfs'))
  if handles.numSurfs == 0
    rB = 6;
    rAx23 = 0;
    wPBr = 0;
  end
end

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
tmp2 = tmp2;
tmp3 = bB + ax23_h + sAx2Ax3;
set(handles.colorbar, 'Position', [tmp1 tmp3 bXO_w tmp2-tmp3]); % colorbar, main axis under toggle lines button
zoomedInsetAxisSize = 0.33;  %proportion of main axis
set(handles.zoomedInsetAxes, 'Position', [lB+ax1_w*(1-zoomedInsetAxisSize-0.01) bB+ax1_h*.01 ax1_w*zoomedInsetAxisSize ax1_h*zoomedInsetAxisSize], 'Visible','off','Xtick',[],'YTick',[]);

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
tmp1 = lB + ax1_w + rAx1 + ax23_w + rAx23;
tmp2 = bB + ax1_h - hPBr;
if isfield(handles, 'meshPB')
  for mmm = 1:length(handles.meshPB)
    set(handles.meshPB(mmm), 'Position', [tmp1 tmp2 wPBr hPBr]);
    tmp2 = tmp2 - hPBr - sPBr;
  end
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
%set(handles.volPB, 'BackgroundColor', backGray);
set(handles.volFrame, 'BackgroundColor', backGray);
if isfield(handles, 'meshPB')
  set(handles.meshPB, 'BackgroundColor', backGray);
end
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
set(handles.frame8, 'BackgroundColor', backGray);
set(handles.overlay_cb, 'BackgroundColor', backGray);
set(handles.binary_cb, 'BackgroundColor', backGray);
set(handles.outline_cb, 'BackgroundColor', backGray);
set(handles.atlas_popup, 'BackgroundColor', backGray);
set(handles.roi_label_cb, 'BackgroundColor', backGray);
set(handles.overlay_opacity_slider, 'BackgroundColor', backGray);
set(handles.image_fade_slider, 'BackgroundColor', backGray);

set(handles.frame10, 'BackgroundColor', backGray);
set(handles.annotationAdd_button, 'BackgroundColor', backGray);
set(handles.annotationEdit_button, 'BackgroundColor', backGray);
set(handles.annotationShow_button, 'BackgroundColor', backGray);
set(handles.annotationList_button, 'BackgroundColor', backGray);

set(handles.edit_command, 'BackgroundColor', liteGray); %%%
set(handles.edit_rr, 'BackgroundColor', liteGray);
set(handles.edit_cc, 'BackgroundColor', liteGray);
set(handles.edit_ss, 'BackgroundColor', liteGray);
set(handles.atlas_pthresh_edit, 'BackgroundColor', liteGray);

set(handles.cycleVolumes, 'FontSize', fontSize);
if (isfield(handles, 'volPB'))
  set(handles.volPB, 'FontSize', fontSize);
end

if isfield(handles, 'meshPB')
  set(handles.meshPB, 'FontSize', fontSize);
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
set(handles.overlay_cb, 'FontSize', fontSize);
set(handles.binary_cb, 'FontSize', fontSize);
set(handles.outline_cb, 'FontSize', fontSize);
set(handles.atlas_popup, 'FontSize', fontSize);
set(handles.annotationAdd_button, 'FontSize', fontSize);
set(handles.annotationEdit_button, 'FontSize', fontSize);
set(handles.annotationShow_button, 'FontSize', fontSize);
set(handles.annotationList_button, 'FontSize', fontSize);


set(handles.edit_command, 'FontSize', fontSize); %%%
set(handles.edit_rr, 'FontSize', fontSize);
set(handles.edit_cc, 'FontSize', fontSize);
set(handles.edit_ss, 'FontSize', fontSize);
set(handles.atlas_pthresh_edit, 'FontSize', fontSize);

handles = annotationDraw(handles);

% --- Executes during object deletion, before destroying properties.
function figure_DeleteFcn(hObject, eventdata, handles)
%close annotation window
try
  if ~isempty(handles.annotation.annotationUI) && ishandle(handles.annotation.annotationUI)
    delete(handles.annotation.annotationUI)
  end
end

%% +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
%% +++ Callback functions for ROI overlay
%% +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

%% ------------------------------------------------------------
%% --- Executes on toggle of overlay_cb
%% ------------------------------------------------------------

function overlay_cb_Callback(hObject, eventdata, handles)

handles.anat.params.showOverlay = get(hObject,'Value');
handles = updateDisplay_newvol(handles);
guidata(handles.figure, handles);

%% ------------------------------------------------------------

%% ------------------------------------------------------------
%% --- Executes on toggle of binary_cb: change overlay from binary to prob-weighted
%% ------------------------------------------------------------

function binary_cb_Callback(hObject, eventdata, handles)

handles.anat.params.binarizeOverlay = get(hObject,'Value');
handles = updateDisplay_newvol(handles);
guidata(hObject, handles);

%% ------------------------------------------------------------


%% ------------------------------------------------------------
%% --- Executes on button press in outline_cb: toggle outline
%% ------------------------------------------------------------

function outline_cb_Callback(hObject, eventdata, handles)

handles.anat.params.showOutline = get(hObject,'Value');
handles = updateDisplay_newvol(handles);
guidata(hObject, handles);

%% ------------------------------------------------------------

%% ------------------------------------------------------------
%% --- Executes on selection change in atlas_popup.
%% ------------------------------------------------------------

function atlas_popup_Callback(hObject, eventdata, handles)

choices = cellstr(get(hObject,'String'));
selected = get(hObject,'Value');
if selected==1
  atlas = '';
else
  atlas = choices{selected};
  if any(atlas=='[') %in multi-atlas case, we add [n] to each atlas name for the popup to make it clear which atlas version
    atlas = atlas(1:end-3);
  end
end
handles.anat.params.currentAtlas = atlas;
if isempty(atlas)
  set(handles.atlas_pthresh_edit,'String','');
  set(handles.roi_label_cb,'Value',false);
else
  %uncompress if needed
  c = handles.anat.(atlas).prob;
  if iscell(c)
    prob = zeros(c{1}, 'single');
    prob(c{2}) = c{3};
    handles.anat.(atlas).prob = prob;
    handles.anat.(atlas).showNames = true;
  end
  set(handles.atlas_pthresh_edit,'String',num2str(handles.anat.(atlas).probabilityThreshold,2));
  set(handles.roi_label_cb,'Value',handles.anat.(atlas).showNames);
end

%call the picker
handles = anatomyRoiOverlaySelectUI(handles,atlas);
set(handles.roi_label_cb,'Value',true);
handles.anat.params.showRoiLabel = true;

handles = updateDisplay_newvol(handles);
guidata(handles.figure, handles)

%% ------------------------------------------------------------
%% --- Executes on change in in edit box for current atlas pthresh
%% ------------------------------------------------------------
function atlas_pthresh_edit_Callback(hObject, eventdata, handles)

atlas = handles.anat.params.currentAtlas;
if isempty(atlas)
  set(hObject, 'String','');
else
  thresh = str2double(get(hObject,'String'));
  if isnan(thresh)
    thresh = handles.anat.(atlas).probabilityThreshold;
  else
    thresh = max(0,min(1,abs(thresh)));
    handles.anat.(atlas).probabilityThreshold = thresh;
    guidata(handles.figure, handles);
  end
  set(hObject,'String', num2str(thresh,2));
  handles = updateDisplay_newvol(handles);
  guidata(handles.figure, handles)
end


%% ------------------------------------------------------------
%% --- Executes on button press in roi_label_cb.
%% ------------------------------------------------------------

function roi_label_cb_Callback(hObject, eventdata, handles)
% handles.annotation.showAsegLabel = get(hObject,'Value');
% handles = updateDisplay_newvol(handles);
% guidata(handles.figure, handles);

atlas = handles.anat.params.currentAtlas;
if isempty(atlas)
  set(hObject, 'Value',false);
else
  handles.anat.(atlas).showNames = get(hObject,'Value');
  handles.anat.params.showRoiLabel = get(hObject,'Value');
  handles = updateDisplay_newvol(handles);
  guidata(handles.figure, handles)
end


%% ------------------------------------------------------------


%% ------------------------------------------------------------
%% --- Executes on overlay opacity slider movement
%% ------------------------------------------------------------
function overlay_opacity_slider_Callback(hObject, eventdata, handles)

handles.anat.params.overlayAlpha = get(hObject,'Value');
if handles.anat.params.overlayAlpha < 0.15, handles.anat.params.overlayAlpha=0; end
handles = updateDisplay_newvol(handles);
guidata(handles.figure, handles);

%% ------------------------------------------------------------

%% ------------------------------------------------------------
%% --- Executes on image fade slider movement
%% ------------------------------------------------------------
function image_fade_slider_Callback(hObject, eventdata, handles)

handles.anat.params.imageFade = get(hObject,'Value');
handles = updateDisplay_newvol(handles);
guidata(handles.figure, handles);

%% ------------------------------------------------------------

% --- set visibility of ROI overlay UI
function setAnatomyGUIVisibility(handles, showit)
if showit
  v='on';
else
  v='off';
end
h=handles;
h = [h.frame8 h.overlay_cb h.binary_cb h.outline_cb h.atlas_popup h.atlas_pthresh_edit h.roi_label_cb ...
  h.roi_save_button h.roi_load_button h.overlay_opacity_slider h.image_fade_slider];

set(h,'visible',v)

%% end of callback functions for ROI overlay UI
%% +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++


%% +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
%% +++ Callback functions for annotations
%% +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

% --- Executes on button press in annotationAdd_button.
function annotationAdd_button_Callback(hObject, eventdata, handles)
showVolWins = findobj(get(0,'children'),'flat','Name','showVol');
if length(showVolWins) >1, fprintf(2,'For now, adding annotations only works if a single showVol window is open, sorry.\n'), return, end
handles = annotationAdd(handles, handles.rr, handles.cc, handles.ss);
guidata(handles.figure, handles)

% --- Executes on button press in annotationEdit_button.
function annotationEdit_button_Callback(hObject, eventdata, handles)
handles = annotationEdit(handles, handles.rr, handles.cc, handles.ss);

% --- Executes on button press in annotationShow_button.
function annotationShow_button_Callback(hObject, eventdata, handles)

if isempty(handles.annotation), return, end
if ~handles.annotation.visible
  handles.annotation.visible = true;
  handles.annotation.showAuthor = false;
elseif ~handles.annotation.showAuthor
  handles.annotation.showAuthor = true;
else
  handles.annotation.visible = false;
end
handles = updateDisplay_newvol(handles);
guidata(handles.figure, handles)

% --- Executes on button press in annotationList_button.
function annotationList_button_Callback(~, eventdata, handles)

if isempty(handles.annotation), return, end

%check if multiple showVol open--that'll lead to problems
showVolWins = findobj(get(0,'children'),'flat','Name','showVol');
if length(showVolWins) > 1, fprintf(2,'For now, the annotation list UI only works if a single showVol window is open, sorry.\n'), return, end
% create UI if needed, or show/hide if already created
if isempty(handles.annotation.annotationUI) || ~ishandle(handles.annotation.annotationUI)
  handles.annotation.annotationUI = annotationsUI('UserData',handles,'visible','on');
else
  if strcmp(handles.annotation.annotationUI.Visible,'off')
    handles.annotation.annotationUI.Visible = 'on';
  else
    handles.annotation.annotationUI.Visible = 'off';
  end
end
guidata(handles.figure, handles)

% --- set visibility of annotation UI
function setAnnotationGUIVisibility(handles, showit)
if showit
  v='on';
else
  v='off';
end
h=handles;
h = [h.frame10 h.annotationAdd_button h.annotationEdit_button h.annotationShow_button h.annotationList_button];

set(h,'visible',v)

%% end of callback functions for annotations
%% +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++


% --- Executes on button press in roi_save_button.
function roi_save_button_Callback(hObject, eventdata, handles)
% hObject    handle to roi_save_button (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
if isempty(handles.anat), return, end

anat.params = handles.anat.params;
atlases = handles.anat.atlases;
for i = 1:length(handles.anat.atlases)
  atlas = atlases{i};
  anat.(atlas).uiRoiOverlaySelected = handles.anat.(atlas).uiRoiOverlaySelected;
  anat.(atlas).uiRoiOverlayImg = handles.anat.(atlas).uiRoiOverlayImg;
  anat.(atlas).uiRoiOverlayIdx = handles.anat.(atlas).uiRoiOverlayIdx;
  anat.(atlas).probabilityThreshold = handles.anat.(atlas).probabilityThreshold;
  anat.(atlas).showNames = handles.anat.(atlas).showNames;
end

[fname, fpath] = uiputfile(fullfile(userhome,'*.mat'),'Save a ROI settings file');
fname = fullfile(fpath,fname);
save(fname, 'anat','-v7.3')
handles.anat.params.roiset = fname;
guidata(handles.figure, handles)

% --- Executes on button press in roi_load_button.
function handles = roi_load_button_Callback(hObject, eventdata, handles)
% hObject    handle to roi_load_button (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

if isempty(handles.anat), return, end

if isempty(hObject) && ~isempty(handles.anat.params.roiset) %special case called during initialization with a roi settings in configuration
  fname = handles.anat.params.roiset;
else
  [fname, fpath] = uigetfile(fullfile(userhome,'*.mat'),'Pick a ROI settings file');
  if ~fpath, return, end
  fname = fullfile(fpath,fname);
  handles.anat.params.roiset = fname;
end
load(fname,'anat')

%check if roiset uses same atlas version
if ~strcmp(handles.anat.params.atlasfile, anat.params.atlasfile)
  disp('This roiset file is incompatible with the current atlas version.')
  return
end

handles.anat.params = anat.params;
handles.anat.params.roiset = fname;

atlases = handles.anat.atlases;
for i = 1:length(handles.anat.atlases)
  atlas = atlases{i};
  handles.anat.(atlas).uiRoiOverlaySelected = anat.(atlas).uiRoiOverlaySelected;
  handles.anat.(atlas).uiRoiOverlayImg = anat.(atlas).uiRoiOverlayImg;
  handles.anat.(atlas).uiRoiOverlayIdx = anat.(atlas).uiRoiOverlayIdx;
  handles.anat.(atlas).probabilityThreshold = anat.(atlas).probabilityThreshold;
  handles.anat.(atlas).showNames = anat.(atlas).showNames;
  if ~isempty(handles.anat.(atlas).uiRoiOverlaySelected)
    %we don't save the prob, so uncompress as needed
    c = handles.anat.(atlas).prob;
    if iscell(c)
      prob = zeros(c{1}, 'single');
      prob(c{2}) = c{3};
      handles.anat.(atlas).prob = prob;
    end
  end
end
if ~nargout
  handles = updateDisplay_newvol(handles);
  guidata(handles.figure, handles)
end
