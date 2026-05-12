function writeGIfTI_underlay(icoNum, underlay, dirFreesurfer, outName, splitLR, offset)
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
% offset:       logical     optional argument that determines if the right
%                           hemisphere labels will be offset by a constant
%                           value (default: false)
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

if ~exist('offset', 'var') || isempty(offset)
    offset = false;
else
    if ~islogical(offset)
        error('offset should be either true or false');
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
            Data = read_fs_annot(dirFreesurfer, icoNum, fname, splitLR, offset);
            if splitLR
                leftData  = Data(1);
                rightData = Data(2);
            end

        case 'aparc.a2009s'
            fname = 'aparc.a2009s';
            Data = read_fs_annot(dirFreesurfer, icoNum, fname, splitLR, offset);
            if splitLR
                leftData  = Data(1);
                rightData = Data(2);
            end

        case 'yeo17'
            fname = 'Yeo2011_17Networks_N1000';
            Data = read_fs_annot(dirFreesurfer, icoNum, fname, splitLR, offset);
            if splitLR
                leftData  = Data(1);
                rightData = Data(2);
            end
            
        case 'yeo7'
            fname = 'Yeo2011_7Networks_N1000';
            Data = read_fs_annot(dirFreesurfer, icoNum, fname, splitLR, offset);
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

            % For inflated surfaces, force a space between left and right
            % hemispheres - say 10mm
            xOffset = max(left_gii.vertices(:,1)) - min(right_gii.vertices(:,1)) + 10;
            right_vertices(:,1) = right_vertices(:,1) + xOffset;

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
        res = gifti(prepAtlasData(leftData, 'CortexLeft'));
        save(res, [outName, '_lh'], 'GZipBase64Binary');

        res = gifti(prepAtlasData(rightData, 'CortexLeft'));
        save(res, [outName, '_rh'], 'GZipBase64Binary');
    else
        % res = gifti(struct('data', Data.cdata, 'label', Data.labels));
        res = gifti(prepAtlasData(Data, 'Cortex'));
        save(res, outName, 'GZipBase64Binary');
    end
end
end

function s = prepAtlasData(data, metaDataValue)

s.metadata   = struct('name',{},'value',{});
s.label.name = data.labels.name';
s.label.key  = int32(data.labels.key');
s.label.rgba = data.labels.rgba;

% Set everything but unknown to opaque
s.label.rgba(~strcmpi(data.labels.name, 'unknown'),4) = 1;

% Update data fields
s.data{1}.metadata                      = struct('name', {'AnatomicalStructurePrimary'}, ...
                                                 'value', {metaDataValue});
s.data{1}.space                         = [];
s.data{1}.attributes.Intent             = 'NIFTI_INTENT_LABEL';
s.data{1}.attributes.DataType           = 'NIFTI_TYPE_INT32';
s.data{1}.attributes.ArrayIndexingOrder = 'ColumnMajorOrder';
s.data{1}.attributes.Dim                = numel(data.cdata);
s.data{1}.attributes.Encoding           = 'GZipBase64Binary';
s.data{1}.attributes.Endian             = 'LittleEndian';
s.data{1}.attributes.ExternalFileName   = '';
s.data{1}.attributes.ExternalFileOffset = '';
s.data{1}.data                          = int32(data.cdata);
end