function varargout = annotationsUI(varargin)
% ANNOTATIONS MATLAB code for annotations.fig
%      ANNOTATIONS, by itself, creates a new ANNOTATIONS or raises the existing
%      singleton*.
%
%      H = ANNOTATIONS returns the handle to a new ANNOTATIONS or the handle to
%      the existing singleton*.
%
%      ANNOTATIONS('CALLBACK',hObject,eventData,handles,...) calls the local
%      function named CALLBACK in ANNOTATIONS.M with the given input arguments.
%
%      ANNOTATIONS('Property','Value',...) creates a new ANNOTATIONS or raises the
%      existing singleton*.  Starting from the left, property value pairs are
%      applied to the GUI before annotations_OpeningFcn gets called.  An
%      unrecognized property name or invalid value makes property application
%      stop.  All inputs are2 passed to annotations_OpeningFcn via varargin.
%
%      *See GUI Options on GUIDE's Tools menu.  Choose "GUI allows only one
%      instance to run (singleton)".
%
% See also: GUIDE, GUIDATA, GUIHANDLES

% TO DO: selecting row 'highlights' annotation on showVol figure
% 

% Edit the above text to modify the response to help annotations

% Last Modified by GUIDE v2.5 01-Oct-2020 08:37:05

%% ============================================================
%% 
%%               SECTION 1: INITIALIZATION
%% 
%% ============================================================

% Begin initialization code - DO NOT EDIT
gui_Singleton = 1;
gui_State = struct('gui_Name',       mfilename, ...
                   'gui_Singleton',  gui_Singleton, ...
                   'gui_OpeningFcn', @annotations_OpeningFcn, ...
                   'gui_OutputFcn',  @annotations_OutputFcn, ...
                   'gui_LayoutFcn',  [] , ...
                   'gui_Callback',   []);
if nargin && ischar(varargin{1})
    gui_State.gui_Callback = str2func(varargin{1});
end

if nargout
    [varargout{1:nargout}] = gui_mainfcn(gui_State, varargin{:});
else
    gui_mainfcn(gui_State, varargin{:});
end
% End initialization code - DO NOT EDIT

% --- Executes just before annotations is made visible.
function annotations_OpeningFcn(hObject, eventdata, handles, varargin)
% This function has no output args, see OutputFcn.
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
% varargin   command line arguments to annotations (see VARARGIN)

% Choose default command line output for annotations
handles.output = hObject;
handles.parent = hObject.UserData;

if ~isfield(handles,'UDListener') || isempty(handles.UDListener)
    handles.UDListener = addlistener(gcf,'UserData', 'PostSet', @receiveUD);
end
handles.filterMatchUuid = {};
handles.cmap = [1 1 1; pastelmap; pastelmap];

guidata(gcf, handles)

handles = update(hObject, eventdata, handles);

displayAnnotation(handles,'')
handles.undo_button.Enable = 'off';

%position in upper right of screen
ss = get(0,'ScreenSize');
origUnit = hObject.Units;
hObject.Units = 'pixels';
p = hObject.Position;
p(1) = ss(3) - p(3);
p(2) = ss(4) - p(4)-64;
hObject.Position = p;
hObject.Units = origUnit;

%set up annotation editing fields with cut/copy/paste as well as open URL
%contextual menu
copy_paste(handles.label_edit, handles.abbrev_edit, handles.extent_edit, handles.category_edit, handles.note_edit)

% Update handles structure
guidata(hObject, handles);

% --- Outputs from this function are returned to the command line.
function varargout = annotations_OutputFcn(hObject, eventdata, handles) 
% Get default command line output from handles structure
varargout{1} = handles.output;


%% ============================================================
%% 
%%               SECTION 2: TABLE CALLBACKS
%% 
%% ============================================================

%% ANNOTATION TABLE

% --- Executes when selected cell(s) is changed in annotationTable.
function annotationTable_CellSelectionCallback(hObject, eventdata, handles)

if isempty(eventdata.Indices)
    eventdata.Source.UserData = [];
    displayAnnotation(handles, [])
    return
else
    selrow = eventdata.Indices(1);
end
selUuid = handles.annotationTableRowUuid(selrow);
selAuthor = handles.annotationTableRowAuthor(selrow);
A = handles.parent.annotation.A;
Asel = A(strcmp(selUuid, A.uuid) & strcmp(selAuthor, A.author) ,:);
eventdata.Source.UserData = {selrow Asel};

displayAnnotation(handles, Asel)

%% ------------------------------------------------------------ 

% --- Executes when entered data in editable cell(s) in annotationTable.
function annotationTable_CellEditCallback(hObject, eventdata, handles)

%only one column is editable
selrow = eventdata.Indices(1);
selUuid = handles.annotationTableRowUuid(selrow);

uuid = handles.parent.annotation.uuid;
if eventdata.NewData
    uuid.hidden = setdiff(uuid.hidden, selUuid);
else
    uuid.hidden = union(uuid.hidden, selUuid);
end
handles.parent.annotation.uuid.hidden = uuid.hidden;
guidata(handles.figure1, handles)

ud.annotation.action = 'updateVisible';
ud.annotation.hidden = uuid.hidden;
ud.annotation.hiddenCategories = uuid.hiddenCategories;
set(handles.parent.figure, 'UserData', ud);

%% ------------------------------------------------------------ 

%% CATEGORY TABLE

% --- Executes when entered data in editable cell(s) in categoryTable.
function categoryTable_CellEditCallback(hObject, eventdata, handles)

% show/hide visibility
selrow = eventdata.Indices(1);
%selCategory = handles.categoryTableRowCategory{selrow};
selCategoryUuid = handles.categoryTableRowUuid{selrow};
%mod = get(get(hObject,'parent'), 'CurrentModifier');%if shift, hide all others
% if contains(mod,'shift')
%     toHide = setdiff(2:length(selCategory),selrow);
% end
updateVisibleCategories(handles, selCategoryUuid, eventdata.NewData)

%% ------------------------------------------------------------ 

% --- Executes when selected cell(s) is changed in categoryTable.
function categoryTable_CellSelectionCallback(hObject, eventdata, handles)

selrows = eventdata.Indices(:,1);
selcols = eventdata.Indices(:,2);
if any(selcols==1), return, end %ignore selection events on checkboxes, which occur when toggling
eventdata.Source.UserData = selrows;

%% ------------------------------------------------------------ 

% --- Executes on button press in categoryShowOnly_pb.
function categoryShowOnly_pb_Callback(hObject, eventdata, handles)
% hObject    handle to categoryShowOnly_pb (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% get selected rows
selrows =  handles.categoryTable.UserData;
showUuid = handles.categoryTableRowUuid(selrows);
hideUuid = setdiff(handles.categoryTableRowUuid, showUuid);

updateVisibleCategories(handles, hideUuid, 0, 1);

%% ------------------------------------------------------------ 

% --- Executes on button press in categoryShowAll_pb.
function categoryShowAll_pb_Callback(hObject, eventdata, handles)
% hObject    handle to categoryShowAll_pb (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

updateVisibleCategories(handles, handles.categoryTableRowUuid, 1, 1); %show all

%% ------------------------------------------------------------ 
% --- set category visibility based on uuid.hiddenCategories
function updateVisibleCategories(handles, selCategoryUuid, showHide, force)
% updateVisibleCategories(handles, categoryUuid, showHide, force)
%  showHide: show=1, hide=0
%  force: true sets shown/hidden to categoryUuid; else toggles
if ~iscell(selCategoryUuid); selCategoryUuid = {selCategoryUuid}; end
if nargin < 4, force = false; end

%show/hide is done at annotation level, look up which annotations have
%these categories
categoryUuid = annotationsWithCategory(handles, selCategoryUuid);
if isempty(categoryUuid), return, end

uuid = handles.parent.annotation.uuid;
if showHide %show
  if force
    uuid.hidden = setdiff(handles.parent.annotation.A.uuid, categoryUuid);
    uuid.hiddenCategories = setdiff(handles.categoryTableRowUuid, selCategoryUuid);
  else
    uuid.hidden = setdiff(uuid.hidden, categoryUuid);
    uuid.hiddenCategories = setdiff(uuid.hiddenCategories, selCategoryUuid);
  end
else %hide
  if force
    uuid.hidden = categoryUuid;
    uuid.hiddenCategories = selCategoryUuid;
  else
    uuid.hidden = union(uuid.hidden, categoryUuid);
    uuid.hiddenCategories = union(uuid.hiddenCategories, selCategoryUuid);
  end
end
handles.parent.annotation.uuid = uuid;

%update annotation list
guidata(handles.figure1, handles)
update(handles.figure1, '', handles);

%update showVol viewer
ud.annotation.action = 'updateVisible';
ud.annotation.hidden = uuid.hidden;
ud.annotation.hiddenCategories = uuid.hiddenCategories;
handles.parent.figure.UserData = ud;

%% ------------------------------------------------------------ 
% return all annotation uuids with category in selCategoryUuid
function categoryUuid = annotationsWithCategory(handles, selCategoryUuid)

%'none' category needs special attention, since its uuid is empty, which
%matches every other uuid
isNone = cellfun(@(x) isempty(x), selCategoryUuid);

A = handles.parent.annotation.A;
categoryUuid = A.uuid(contains(A.category, selCategoryUuid(~isNone)));
if any(isNone)
  AisNone = cellfun(@(x) isempty(x), A.category);
  if all(isNone)
    categoryUuid = A.uuid(AisNone);
  else
    categoryUuid = union(categoryUuid, A.uuid(AisNone));
  end
end

%% ============================================================
%% 
%%               SECTION 3: UI CALLBACKS
%% 
%% ============================================================

% --- Executes on button press in githubSync_button.
function githubSync_button_Callback(hObject, eventdata, handles)
ud.annotation.action = 'sync';
handles.parent.figure.UserData = ud;

%% ------------------------------------------------------------  

% --- Executes on button press in goto_button.
function goto_button_Callback(hObject, eventdata, handles)

if isempty(handles.annotationTable.UserData), return; end
Ag = handles.annotationTable.UserData{2};
ud.annotation.action = 'goto';
ud.annotation.gotoRcs = [Ag.r(1) Ag.c(1) Ag.s(1)];
handles.parent.figure.UserData = ud;


%% ------------------------------------------------------------ 

% --- Executes on button press in edit_button.
function edit_button_Callback(hObject, eventdata, handles)
% if isempty(handles.annotationTable.UserData), return; end
% Ad = handles.annotationTable.UserData{2};
% ud.annotation.action = 'edit';
% ud.annotation.editUuid = Ad.uuid;
% handles.parent.figure.UserData = ud;

%% ------------------------------------------------------------ 

% --- Executes on button press in savechange_button.
function savechange_button_Callback(hObject, eventdata, handles)
%grab new values from boxes
if isempty(handles.annotationTable.UserData), return, end
Au = handles.annotationTable.UserData{2};

newLabel = handles.label_edit.String;
oldLabel = Au.label{1};
newAbbrev = handles.abbrev_edit.String;
oldAbbrev = Au.abbrev{1};
newExtent = str2num(handles.extent_edit.String);
if any(isnan(newExtent)) || ~(length(newExtent) == 1 || length(newExtent) == 3)
    fprintf(2,'bad value for extent (%s), using [1 1 1]\n',sprintf('%d ', newExtent))
    newExtent = [1 1 1];
else
    if length(newExtent) == 1
        newExtent = newExtent * [1 1 1];
    end
end
oldExtent = Au.extent;
if strcmp(handles.category_edit.Visible, 'on')
    newCategory = handles.category_edit.String;
else
    cats = handles.category_popup.String;
    sel = handles.category_popup.Value;
    if sel == 1
        newCategory = {''};
    else
        newCategory = cats(sel);
    end
end
C = handles.parent.annotation.C;
oldCategory = C.category(strcmp(Au.category{1}, C.uuid));
newNote = join(handles.note_edit.String,'; ');
oldNote = Au.note{1};

changed = ~strcmp(newLabel, oldLabel) || ~strcmp(newAbbrev, oldAbbrev) || ...
    ~strcmp(newCategory, oldCategory) || ~strcmp(newNote, oldNote) || ~isequal(newExtent, oldExtent);

if changed
    Au.label = newLabel;
    Au.abbrev = newAbbrev;
    Au.extent = newExtent;
    Au.category = newCategory; %string; convered to UUID in showvol
    Au.note = newNote;

    ud.annotation.action = 'updateAnnotation';
    ud.annotation.updateRow = Au;
    handles.parent.figure.UserData = ud;
    
    handles.undo_button.Enable = 'on';
end

%% ------------------------------------------------------------ 

% --- Executes on button press in delete_button.
function delete_button_Callback(hObject, eventdata, handles)

if isempty(handles.annotationTable.UserData), return; end
Ad = handles.annotationTable.UserData{2};
if strcmp(Ad.label,'README') && Ad.r == 51 && Ad.c == 61 && Ad.s == 61
    disp('Don''t delete the README')
    return
end
ud.annotation.action = 'delete';
ud.annotation.deleteRow = Ad;
handles.parent.figure.UserData = ud;
handles.undo_button.Enable = 'on';

%% ------------------------------------------------------------ 

% --- Executes on button press in undo_button.
function undo_button_Callback(hObject, eventdata, handles)
ud.annotation.action = 'undo';
handles.parent.figure.UserData = ud;
handles.undo_button.Enable = 'off';


%% ------------------------------------------------------------ 


%% FILTERING annotations shown in table

% --- Executes on button press in onlyVisible_cb.
function onlyVisible_cb_Callback(hObject, eventdata, handles)
update(handles.figure1, eventdata, handles);

%% ------------------------------------------------------------ 

% --- Executes on button press in showhidden_cb.
function showhidden_cb_Callback(hObject, eventdata, handles)
update(handles.figure1, eventdata, handles);

%% ------------------------------------------------------------ 
% --- Executes on any text entered into filtering fields
function filter_button_Callback(hObject, eventdata, handles)

label = handles.labelFilter_box.String;
abbr = handles.abbrFilter_box.String;
author = handles.authorFilter_box.String;
description = handles.descriptionFilter_box.String;
A = handles.parent.annotation.A;

matchFilter = contains(A.label,label,'IgnoreCase',true) & contains(A.abbrev,abbr,'IgnoreCase',true) & contains(A.author, author,'IgnoreCase',true) & contains(A.note, description,'IgnoreCase',true);
handles.filterMatchUuid = A.uuid(matchFilter);
handles = update(hObject, eventdata, handles);
guidata(handles.figure1, handles)

%% ------------------------------------------------------------ 

% --- Executes on button press in clearFilter_button.
function clearFilter_button_Callback(hObject, eventdata, handles)

%clear all text from filter boxes and filter
handles.labelFilter_box.String='';
handles.abbrFilter_box.String='';
handles.authorFilter_box.String='';
handles.descriptionFilter_box.String='';

filter_button_Callback(hObject, eventdata, handles)

%% ------------------------------------------------------------

%% ============================================================
%% 
%%               SECTION 4: EVENTS from showVol window
%% 
%% ============================================================

%% LISTENER
function receiveUD(hObject,event)
thisFig = event.AffectedObject;
showVolHandles = thisFig.UserData; %showVol simply passes its full handles
handles = guidata(thisFig);
handles.parent = showVolHandles;

try
    handles = update(hObject,[], handles);
    guidata(thisFig, handles)
catch
    lasterr
    disp('no action')
end

%% ============================================================
%% 
%%               SECTION 5: CONTENT
%% 
%% ============================================================

%% TABLE CONTENT
function handles = update(hObject, eventdata, handles)

an = handles.parent.annotation;

%% CATEGORY TABLE
C = an.C;
D = {};
rowCategory = {};
rowUuid = {};
categoryColors = [];
isNone = cellfun(@(x) isempty(x), an.uuid.hiddenCategories);

handles.categoryVisible = ~contains(C.uuid, an.uuid.hiddenCategories(~isNone));
if any(isNone)
  handles.categoryVisible(1) = false;
else
  handles.categoryVisible(1) = true;
end

for iC = 1:height(C)
    x = C(iC,:);
    D{iC,1} = handles.categoryVisible(iC);
    rowCategory(iC) = x.category;
    rowUuid(iC) = x.uuid;
    
    if isempty(rowCategory{iC})
        cat = '[None]';
    else
        cat = rowCategory{iC};
    end
    
    %add # of anotations
    ann = annotationsWithCategory(handles, x.uuid);
    cat = sprintf('%s     (%d)', cat, length(ann));
    
    D{iC,2}= cat;
    cidx = 1+mod(x.colorIdxUnique-1, size(handles.cmap,1));
    categoryColors(iC,:) = handles.cmap(cidx,:);
end

handles.categoryTable.Data = D;
handles.categoryTable.BackgroundColor = categoryColors;
handles.categoryTableRowCategory = rowCategory;
handles.categoryTableRowUuid = rowUuid;

%% ANNOTATION TABLE
A = an.A;
if handles.onlyVisible_cb.Value
    [~, inViewIdx] = intersect(A.uuid, an.uuid.inView);
    A = A(inViewIdx,:);
end
if ~isempty(handles.filterMatchUuid)
    [~, showIdx] = intersect(A.uuid, handles.filterMatchUuid);
    A = A(showIdx,:);
end
if ~handles.showhidden_cb.Value
    [~,hiddenIdx] = intersect(A.uuid, handles.parent.annotation.uuid.hidden);
    A(hiddenIdx,:) = [];
end

categoryStr = cellfun(@(x) C.category(strcmp(x,C.uuid)), A.category);
A.categoryStr = categoryStr;
A = sortrows(A, {'categoryStr', 'label'});

D = {};
rowUuid = {};
rowAuthor = {}; %because of way we handle merge, need author to fully identify a row
rowColors = [];
for iA=1:height(A)
    x = A(iA,:);
    D{iA,1}=~any(contains(handles.parent.annotation.uuid.hidden, x.uuid));
    D{iA,2}=sprintf('<html><font size=-2>[%d %d %d]</font></html>',x.R,x.A,x.S);
    D{iA,3}=sprintf('<html><font size=-2>[%d %d %d]</font></html>',x.r,x.c,x.s);
    isConflict = ~isempty(regexp(x.abbrev{1},'^\*[DM]-','once'));
    if isConflict
        abbrev = ['<html><span style="color:red"><b>' x.abbrev{1}];
        label = ['<html><span style="color:red">' x.label{1}];
    else
        abbrev = ['<html><b>' x.abbrev{1}];
        label = x.label{1};
    end
    D{iA,4} = abbrev;
    D{iA,5} = label;
    D{iA,6} = x.author{1}(1:2);
    note = x.note{1};
    %if the note has Q: or Q(initials):, but not A:, it has an 'unanswered question'
    nQuestion = length(regexp(note, 'Q\(?\w*\)?:'));
    nAnswer   = length(regexp(note, 'A\(?\w*\)?:'));
    if nQuestion>nAnswer
        note = ['<html><span style="color:red">' note];
    else
        if nQuestion>0
            note = ['<html><span style="color:green">' note];
        end
    end
    D{iA,7} = note;
    rowUuid(iA) = x.uuid;
    rowAuthor(iA) = x.author;
    %which category
    if isempty(x.category)
        rowColors(iA,:) = [1 1 1];
    else
        catIdx = find(contains(C.uuid, x.category));
        if ~isempty(catIdx)
            rowColors(iA,:) = categoryColors(catIdx(1),:);
        else
            rowColors(iA,:) = [1 .5 0]; %should not happen
        end
    end
end

handles.annotationTable.Data = D;
if ~isempty(D)
    handles.annotationTable.BackgroundColor = rowColors;
end

handles.annotationTableRowUuid = rowUuid;
handles.annotationTableRowAuthor = rowAuthor;

guidata(handles.figure1,handles)

%% ------------------------------------------------------------ 

%% EDIT PANEL
function displayAnnotation(handles, A)

if ~isempty(A)
    C = handles.parent.annotation.C;
    if ~isempty(A.category{1})
        cat = C.category{strcmp(C.uuid, A.category)};
    else
        cat = {'[None]'};
    end
    if ~isempty(A.note)
        note = strrep(A.note,'; ',newline);
    else
        note = {''};
    end
    
    handles.label_edit.String    = A.label;
    handles.abbrev_edit.String   = A.abbrev;
    handles.extent_edit.String   = num2str(A.extent);
    handles.category_edit.String = cat;
    handles.author_text.String   = A.author;
    handles.date_text.String     = char(A.date);
    handles.volume_text.String   = A.vol;
    handles.note_edit.String     = note;
    handles.savechange_button.Enable = 'on';
    handles.delete_button.Enable     = 'on';
    handles.category_edit.Visible  = 'off';
    cats = [{'[None]'}; C.category(2:end); {'[New Category]'}];
    handles.category_popup.String = cats;
    handles.category_popup.Value = find(strcmp(cat, cats));
    handles.category_popup.Visible = 'on';
    handles.category_popup.Enable    = 'on';
else
    handles.label_edit.String    = '';
    handles.abbrev_edit.String   = '';
    handles.extent_edit.String   = '';
    handles.category_edit.String = '';
    handles.author_text.String   = '';
    handles.date_text.String     = '';
    handles.volume_text.String   = '';
    handles.note_edit.String     = '';
    handles.savechange_button.Enable = 'off';
    handles.delete_button.Enable     = 'off';
    handles.category_popup.Enable    = 'off';
    handles.category_popup.Visible = 'on';
    handles.category_popup.String = {''};
    handles.category_popup.Value = 1;
    handles.category_edit.String = '';
    handles.category_edit.Visible  = 'off';
end

% --- Executes on selection change in category_popup.
function category_popup_Callback(hObject, eventdata, handles)
% hObject    handle to category_popup (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: contents = cellstr(get(hObject,'String')) returns category_popup contents as cell array
%        contents{get(hObject,'Value')} returns selected item from category_popup
cats = cellstr(get(hObject,'String'));
category = cats(hObject.Value);
if strcmp(category,'[New Category]')
    handles.category_popup.Visible = 'off';
    handles.category_edit.String = '';
    handles.category_edit.Visible  = 'on';
end

%% ------------------------------------------------------------ 

%% HELPER
function cmap = pastelmap
% Modified Brewer pastel 1 colormap for category backgrounds
cmap = [[251,180,174];[179,205,227];[204,235,197];[222,203,228];[254,217,166];...
  [255,255,204];[229,216,189];[253,218,236];[242,242,242]]/255;

cmaph = rgb2hsv(cmap);
cmaph = [cmaph.*[1 .5 1.05]; cmaph]; %add a faded version to start, followed by original
cmaph = min(cmaph,1); %clamp
cmap = hsv2rgb(cmaph);
cmap(9,:) = cmap(18,:); %fix up the grays
cmap(18,:) = [210 210 210]/255;
