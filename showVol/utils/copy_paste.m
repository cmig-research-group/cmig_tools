function copy_paste(varargin)

% Copy Paste ver. 1.0
% Amitabh Verma
% amtukv@gmail.com
% 
% This is a very simple cut/copy/paste context menu using java handles. This code uses findjobj by Yair Altman 
% (http://www.mathworks.com/matlabcentral/fileexchange/14317-findjobj-find-java-handles-of-matlab-graphic-objects) to locate the java handles.
% 
% Usgae:
% copy_paste(hObject)
% 
% Example:
% 
% f1 = figure('Name','Copy-Paste','Position', [480 480 340 100],'NumberTitle','off');
% h1 = uicontrol(f1, 'Style', 'edit',...
%        'String', 'Cut',...
%        'BackgroundColor', [1 1 1],...
%        'Position', [10 10 100 50]);
% h2 = uicontrol(f1, 'Style', 'edit',...
%        'String', 'Copy',...
%        'BackgroundColor', [1 1 1],...
%        'Position', [120 10 100 50]);   
% h3 = uicontrol(f1, 'Style', 'edit',...
%        'String', 'Paste',...
%        'BackgroundColor', [1 1 1],...
%        'Position', [230 10 100 50]);
% 
% 
% For single edit box object
% copy_paste(h2);
% 
% For selected edit box object
% copy_paste(h1,h3);
% 
% For all the edit boxes in a figure
% copy_paste(gcf) or copy_paste(f1)
% 
% Please report bugs, feature request, comments and suggestions to Amitabh Verma <amtukv@gmail.com>
% END
% From file exchange. JRI edited


% Define a context menu; it is not attached to anything
copypastemenu = uicontextmenu;

% Define the context menu items and install their callbacks
uimenu(copypastemenu, 'Text', 'Cut', 'Accelerator','X', 'MenuSelectedFcn', 'cut(jObject)');
uimenu(copypastemenu, 'Text', 'Copy', 'Accelerator','C', 'MenuSelectedFcn', 'copy(jObject)');
uimenu(copypastemenu, 'Text', 'Paste', 'Accelerator','V', 'MenuSelectedFcn', 'paste(jObject)');
uimenu(copypastemenu, 'Text', 'Select All', 'Accelerator','A', 'MenuSelectedFcn', 'selectAll(jObject)');
uimenu(copypastemenu, 'Text', 'Open URL', 'MenuSelectedFcn', 'web(jObject.SelectedText,''-new'',''-browser'')');

% Find the Java handle using findjobj by Yair Altman and is called on
% right-click in/on the object. (Plase findjobj in the same dir as copy_paste)
set(copypastemenu, 'Callback', '[jObject,hObject] = javahandle();')

% Loop to check all the passed handles for a figure object
for loop1 = 1:size(varargin,2)
% Check to see if the call was made for the complete figure or single object
    if strcmp(get(varargin{loop1},'Type'),'figure')
    
% Locate edit box objects
    hedit_box = findall(varargin{loop1},'Style','edit');

% Attach the context menu to each edit box
    for no_ebox = 1:length(hedit_box)
        set(hedit_box(no_ebox),'uicontextmenu',copypastemenu)
    end

    else
    for loop2 = 1:size(varargin,2)
% Attach the context menu (to specific edit box objects)
        if strcmp(get(varargin{loop2},'Style'),'edit')
            set(varargin{loop2},'uicontextmenu',copypastemenu)
        end
    end
    end
end
end