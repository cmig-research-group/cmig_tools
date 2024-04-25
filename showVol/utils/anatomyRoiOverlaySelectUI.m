

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