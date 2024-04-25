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