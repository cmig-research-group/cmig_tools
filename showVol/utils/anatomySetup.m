
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

%file pre-processed by prepareAtlases.m. Specific atlases are expanded on
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
  atlasnum = str2double(roiatlas(end));
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

%FIXME -- forgot to subtract 10000 from fiber roicodes!
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