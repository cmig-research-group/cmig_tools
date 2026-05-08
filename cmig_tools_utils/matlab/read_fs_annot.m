function out_struct = read_fs_annot(dirFreesurfer, icoNum, fname, splitLR, offset)
% Function that reads an annotation file from Freesurfer and returns
% concatenated list of cData

if ~exist('splitLR', 'var') || isempty(splitLR)
    splitLR = false;
end

if ~exist('icoNum', 'var') || isempty(icoNum)
    icoNum = 5;
else
    if icoNum < 5
        error('Supported ico numbers for atlas files are: 5, 6, and 7');
    end
end

if ~exist('offset', 'var') || isempty(offset)
    offset = false;
else
    if ~islogical(offset)
        error('offset should be either true or false');
    end
end

% Determine fsaverage folder name
switch icoNum
    case 3
        fsaverage_name = 'fsaverage3';
    case 4
        fsaverage_name = 'fsaverage4';
    case 5
        fsaverage_name = 'fsaverage5';
    case 6
        fsaverage_name = 'fsaverage6';
    case 7
        fsaverage_name = 'fsaverage';
end

% Path to annotation files and path to script file
dir_annotations = fullfile(dirFreesurfer, 'subjects', fsaverage_name, 'label');
dir_script      = fullfile(dirFreesurfer, 'matlab');
fname_script    = 'read_annotation.m';

% Left and right hemisphere names
fname_lh = ['lh.', fname, '.annot'];
fname_rh = ['rh.', fname, '.annot'];

% Make sure that annotation files and read_annotation file exist
if ~exist(fullfile(dir_annotations, fname_lh), 'file') || ...
   ~exist(fullfile(dir_annotations, fname_rh), 'file')
    error(['Unable to find ', fname, ' in: ', dir_annotations]);
else
    if ~exist(fullfile(dir_script, fname_script), 'file')
        error(['Unable to find ', fname_script, ' in: ', dir_script]);
    end
end

% Load annotation files
tmp = pwd;
cd(dir_script);
[lh_vertices, lh_label, lh_colortable] = read_annotation(fullfile(dir_annotations, fname_lh)); %#ok<ASGLU>
[rh_vertices, rh_label, rh_colortable] = read_annotation(fullfile(dir_annotations, fname_rh)); %#ok<ASGLU>
cd(tmp);

% What is the "index" of mapping of *_label to entries in
% *_colortable.table(:,5)?
[~, lh_mapping] = ismember(lh_label, lh_colortable.table(:,5));
[~, rh_mapping] = ismember(rh_label, rh_colortable.table(:,5));

% % Extract names
% lh_names = lh_colortable.struct_names(lh_mapping);
% rh_names = rh_colortable.struct_names(rh_mapping);

% Append some integer constant for right hemisphere to maintain unique
% separate mapping: this is not the most robust solution
if offset
    maxVal = max(unique(lh_mapping));
    if maxVal < 100
        toAdd = 100;
    else
        if maxVal < 1000
            toAdd = 1000;
        else
            toAdd = 10000;
        end
    end
else
    toAdd = 0;
end

% Find the location of unknown
lh_loc_unknown = strcmpi(lh_colortable.struct_names, 'unknown');
rh_loc_unknown = strcmpi(rh_colortable.struct_names, 'unknown');

% Get keys
lh_keys = unique(lh_mapping);
rh_keys = unique(rh_mapping);

% Check if medial wall is missing in the keys
if length(lh_keys) ~= length(lh_colortable.struct_names)
    wch = colvec(setdiff(1:length(lh_colortable.struct_names), lh_keys));

    % Append wch
    lh_keys = sort([lh_keys; wch]);
end

% Repeat for the right hemisphere
if length(rh_keys) ~= length(rh_colortable.struct_names)
    wch = colvec(setdiff(1:length(rh_colortable.struct_names), rh_keys));

    % Append wch
    rh_keys = sort([rh_keys; wch]);
end

% Add a constant to rh
rh_keys(~rh_loc_unknown) = rh_keys(~rh_loc_unknown) + toAdd;

% Append lh and rh to ROI names
lh_colortable.struct_names(~lh_loc_unknown) = strcat(lh_colortable.struct_names(~lh_loc_unknown), {'_lh'});
rh_colortable.struct_names(~rh_loc_unknown) = strcat(rh_colortable.struct_names(~rh_loc_unknown), {'_rh'});

% Create a labels table
if splitLR
    left_struct.labels.name  = lh_colortable.struct_names;
    left_struct.labels.key   = lh_keys;
    left_struct.labels.rgba  = lh_colortable.table(:,1:4);
    left_struct.cdata        = lh_mapping;
    
    right_struct.labels.name = rh_colortable.struct_names;
    right_struct.labels.key  = rh_keys;
    right_struct.labels.rgba = rh_colortable.table(:,1:4);
    right_struct.cdata       = rh_mapping;

    out_struct = [left_struct; right_struct];
else
    out_struct.labels.name = [lh_colortable.struct_names; rh_colortable.struct_names(~rh_loc_unknown)];
    out_struct.labels.key  = [lh_keys; rh_keys(~rh_loc_unknown)];
    out_struct.labels.rgba = [lh_colortable.table(:,1:4); rh_colortable.table(~rh_loc_unknown,1:4)];
    out_struct.cdata       = [lh_mapping; rh_mapping];
end