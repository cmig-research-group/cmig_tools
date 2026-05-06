function writeGIfTI_underlay(icoNum, underlay, dirFreesurfer, outName, splitLR)
% Function that writes out GIfTI file(s) given a set of inputs
% Requires the GIfTI toolbox: https://github.com/gllmflndn/gifti to be on
% path such that the function 'gifti' is callable
% 
%% Inputs:
% --------
% icoNum:       scalar      number indicating which icosahedron number
%                           should be used for creating GIfTI images;
%                           defaults to icoNum = 5
%                           should be a number between 5 and 7 (see Notes): 
%                               * icoNum == 5, 10242 vertices, 20480 faces
%                               * icoNum == 6, 40962 vertices, 81920 faces
%                               * icoNum == 7, 163842 vertices, 327680 faces
% 
% underlay:     character   a value indicating which of the following
%                           underlays should be generated as a GIfTI image:
%                               * 'pial'
%                               * 'inflated'
%                               * 'white'
%                               * name of a supported Freesurfer atlas:
%                                   - 'aparc' (refers to aparc)
%                                   - 'aparc.a2009s' (refers to aparc.a2009s)
%                                   - 'yeo17' (refers to Yeo 17 networks)
%                                   - 'yeo7'  (refers to Yeo 7 networks)
% 
% dirFreesurfer: character full path to a location where Freesurfer is
%                           installed (assumes the "base" folder of
%                           Freesurfer, wherein other folders like
%                           "average", "subjects", "matlab", are located;
%                           only required if underlay is one of the
%                           Freesurfer atlas files
% 
% outName:      character   cell string containing the full paths and name
%                           of the output GIfTI files (no extensions); if
%                           empty, uses pwd and generates file names
% 
% splitLR:      logical     yes or no indicating if the GIfTI file should
%                           be written out as separate left and right
%                           hemisphere files or an overall file 
%                           (default: true)
% 
%% Notes:
% Currently, this function sources all the files from the relevant
% fsaverage folders and saves them as GIfTI files
% 
% The supported ico numbers are: 5, 6, and 7
% 
% curvature files are not currently supported
% 
% Following ico numbers are not currently supported:
%   * icoNum == 0, 12 vertices, 20 faces
%   * icoNum == 1, 42 vertices, 80 faces
%   * icoNum == 2, 162 vertices, 320 faces
%   * icoNum == 3, 642 vertices, 1280 faces
%   * icoNum == 4, 2562 vertices, 5120 faces
% 
% Inflated surfaces need to be fixed
% 
%% References:
% https://brainder.org/2016/05/31/downsampling-decimating-a-brain-surface/
% https://www.jiscmail.ac.uk/cgi-bin/wa-jisc.exe?A2=SPM;9b250ca1.1710
% 
%% Requirement(s):
% GIfTI toolbox: https://github.com/gllmflndn/gifti

%% Check inputs
% Check if underlay type is specified
if ~exist('underlay', 'var') || isempty(underlay)
    error('Please specify the underlay type to be saved');
else
    underlay = lower(underlay);

    % Check if valid underlay is requested
    if ~ismember(underlay, {'pial', 'curv', 'curvature', 'inflated', 'white', ...
                            'aparc', 'aparc.a2009s', 'yeo17', 'yeo7'})
        error(['Unknown underlay specified: ', underlay]);
    end
end

% Check ico number
if ~exist('icoNum', 'var') || isempty(icoNum)
    icoNum = 5;
else
    if icoNum < 5 || icoNum > 7
        error('icoNum should be between 5-7 (inclusive)');
    end
end

% Check if Freesurfer path is provided
if ~exist('dirFreesurfer', 'var') || isempty(dirFreesurfer)
    % Path is not provided, try to get from system environment
    [~, dirFreesurfer] = system('$FREESURFER_HOME');

    if isempty(dirFreesurfer)
        error('Please provide full path to where Freesurfer is located');
    else
        % If there is a colon on the path (for example, ":permission denied",
        % strip it out of the path name
        [dirFreesurfer, ~] = strsplit(dirFreesurfer, ':');
        dirFreesurfer      = dirFreesurfer{1};

        % Does this location exist?
        if exist(dirFreesurfer, 'dir')
            warning(['Freesurfer path not provided; queried from system environment: ', dirFreesurfer]);
        else
            error('Please provide full path to where Freesurfer is located');
        end
    end
else
    if ~exist(dirFreesurfer, 'dir')
        error(['Unable to find directory: ', dirFreesurfer]);
    end
end

% Check if output name is provided
if ~exist('outName', 'var') || isempty(outName)
    genNames = true;
else
    genNames = false;

    % Strip extensions out of outName, if provided
    [~, ~, ext] = fileparts(outName);
    outName     = strrep(outName, ext, '');
end

% Do we need to split hemispheres?
if ~exist('splitLR', 'var') || isempty(splitLR)
    splitLR = true;
else
    if ~islogical(splitLR)
        error('splitLR should be either true or false');
    end
end

%% Check for GIfTI toolbox
if ~exist('gifti.m', 'file')
    tmp_gifti_dir = fullfile(fileparts(fileparts(which('writeGIfTI'))), 'gifti_toolbox');
    if exist(tmp_gifti_dir, 'dir')
        warning(['gifti.m not found; adding: ', tmp_gifti_dir, ' to MATLAB path']);
        addpath(tmp_gifti_dir);
    else
        error(['Unable to find gifti.m file on MATLAB path; ', ...
               'please ensure that you have a copy of the GIfTI toolbox on your MATLAB path']);
    end
end

%% Handle atlas files
if ismember(underlay, {'aparc', 'aparc.a2009s', 'yeo17', 'yeo7'})
    switch underlay
        case 'aparc'
            fname = 'aparc';
            Data  = read_annotations(dirFreesurfer, icoNum, fname, splitLR);
            if splitLR
                leftData  = Data(1);
                rightData = Data(2);
            end

        case 'aparc.a2009s'
            fname = 'aparc.a2009s';
            Data  = read_annotations(dirFreesurfer, icoNum, fname, splitLR);
            if splitLR
                leftData  = Data(1);
                rightData = Data(2);
            end

        case 'yeo17'
            fname = 'Yeo2011_17Networks_N1000';
            Data  = read_annotations(dirFreesurfer, icoNum, fname, splitLR);
            if splitLR
                leftData  = Data(1);
                rightData = Data(2);
            end
            
        case 'yeo7'
            fname = 'Yeo2011_7Networks_N1000';
            Data  = read_annotations(dirFreesurfer, icoNum, fname, splitLR);
            if splitLR
                leftData  = Data(1);
                rightData = Data(2);
            end
    end
else
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

    % Working path
    dir_fsFiles = fullfile(dirFreesurfer, 'subjects', fsaverage_name, 'surf');

    switch underlay
        case 'pial'
            left_gii        = gifti(fullfile(dir_fsFiles, 'lh.pial'));
            right_gii       = gifti(fullfile(dir_fsFiles, 'rh.pial'));
            left_vertices   = left_gii.vertices;
            right_vertices  = right_gii.vertices;
            left_faces      = left_gii.faces;
            right_faces     = right_gii.faces; % + size(left_gii.vertices,1);
            if ~splitLR
                allVertices = [left_vertices; right_vertices];
                allFaces    = [left_faces; right_faces + size(left_gii.vertices,1)];
            end

        case 'white'
            left_gii        = gifti(fullfile(dir_fsFiles, 'lh.white'));
            right_gii       = gifti(fullfile(dir_fsFiles, 'rh.white'));
            left_vertices   = left_gii.vertices;
            right_vertices  = right_gii.vertices;
            left_faces      = left_gii.faces;
            right_faces     = right_gii.faces; % + size(left_gii.vertices,1);
            if ~splitLR
                allVertices = [left_vertices; right_vertices];
                allFaces    = [left_faces; right_faces + size(left_gii.vertices,1)];
            end

        case 'inflated'
            left_gii        = gifti(fullfile(dir_fsFiles, 'lh.inflated'));
            right_gii       = gifti(fullfile(dir_fsFiles, 'rh.inflated'));
            left_vertices   = left_gii.vertices;
            right_vertices  = right_gii.vertices;
            left_faces      = left_gii.faces;
            right_faces     = right_gii.faces; % + size(left_gii.vertices,1);
            if ~splitLR
                allVertices = [left_vertices; right_vertices];
                allFaces    = [left_faces; right_faces + size(left_gii.vertices,1)];
            end
    end
end

% Output name
if genNames
    outName = fullfile(pwd, ['Ico', num2str(icoNum), '_', underlay]);
end

% Prepare GIfTI object
if ismember(underlay, {'pial', 'white', 'inflated'})
    if splitLR
        res = gifti(struct('faces', left_faces, 'vertices', left_vertices));
        save(res, [outName, '_lh'], 'GZipBase64Binary'); %#ok<USENS>
            
        res = gifti(struct('faces', right_faces, 'vertices', right_vertices));
        save(res, [outName, '_rh'], 'GZipBase64Binary');
    else
        res = gifti(struct('faces', allFaces, 'vertices', allVertices));
        save(res, outName, 'GZipBase64Binary');
    end
else
    if splitLR
        res = gifti(struct('cdata', leftData.cdata, 'labels', leftData.labels));
        save(res, [outName, '_lh'], 'GZipBase64Binary');

        res = gifti(struct('cdata', rightData.cdata, 'labels', rightData.labels));
        save(res, [outName, '_rh'], 'GZipBase64Binary');
    else
        res = gifti(struct('cdata', Data.cdata, 'labels', Data.labels));
        save(res, outName, 'GZipBase64Binary');
    end
end
end

function out_struct = read_annotations(dirFreesurfer, icoNum, fname, splitLR)
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

% Find the location of unknown
lh_loc_unknown = strcmpi(lh_colortable.struct_names, 'unknown');
rh_loc_unknown = strcmpi(rh_colortable.struct_names, 'unknown');

% Get keys
lh_keys = unique(lh_mapping);
rh_keys = unique(rh_mapping);

% Check if medial wall is missing in the keys
if length(lh_keys) ~= length(lh_colortable.struct_names)
    wch = reshape(setdiff(1:length(lh_colortable.struct_names), lh_keys), [], 1);

    % Append wch
    lh_keys = sort([lh_keys; wch]);
end

% Repeat for the right hemisphere
if length(rh_keys) ~= length(rh_colortable.struct_names)
    wch = reshape(setdiff(1:length(rh_colortable.struct_names), rh_keys), [], 1);

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
end