function output = FEMA_info(action)
% Function that returns FEMA version and citation information

if ~exist('action', 'var') || isempty(action)
    action = 'version';
else
    action = lower(action);
    if ~ismember(action, {'version', 'ver', 'v', 'citation', 'cite', 'full'})
        action = 'full';
    end
end

dirDocs = fullfile(fileparts(fileparts(which('caller_FEMA'))), 'docs');
version = 'FEMA 3.0.0';

switch action
    case {'version', 'ver', 'v'}
        output = version;

    case {'citation', 'cite'}
        if ~exist(dirDocs, 'dir')
            output = 'Could not find the docs folder';
        else
            output = fileread(fullfile(dirDocs, 'citationInfo.txt'));
        end

    case {'full'}
        txt1 = version;
        if ~exist(dirDocs, 'dir')
            txt2 = 'Could not find the docs folder';
        else
            txt2 = fileread(fullfile(dirDocs, 'citationInfo.txt'));
        end

        output = [txt1, sprintf('%s\n', ' '), txt2];
end