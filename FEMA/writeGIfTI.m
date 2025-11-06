function writeGIfTI(toWrite, icoNum, outName, underlayType, dirFreesurfer, splitLR)
% Function that writes out GIfTI file(s) given a set of inputs
% Requires the GIfTI toolbox: https://github.com/gllmflndn/gifti to be on
% path such that the function 'gifti' is callable

%% Inputs:
% --------
% toWrite:      [p x v]     vector or matrix of estiamtes for p
%                           coefficients and v voxels
% 
% icoNum:       scalar      number indicating which icosahedron number
%                           should be used for creating GIfTI images;
%                           defaults to icoNum = 5
%                           should be a number between 0 and seven:
%                               * icoNum == 0, 12 vertices, 20 faces
%                               * icoNum == 1, 42 vertices, 80 faces
%                               * icoNum == 2, 162 vertices, 320 faces
%                               * icoNum == 3, 642 vertices, 1280 faces
%                               * icoNum == 4, 2562 vertices, 5120 faces
%                               * icoNum == 5, 10242 vertices, 20480 faces
%                               * icoNum == 6, 40962 vertices, 81920 faces
%                               * icoNum == 7, 163842 vertices, 327680 faces
% 
% outName:      [p x 1]    cell string containing the full paths and names
%                          of the output GIfTI files; if empty, uses pwd
% 
% underlayType  character   a value indicating which of the following
%                           underlays should be generated as a GIfTI image:
%                               * 'pial'
%                               * 'curv' | 'curvature'
%                               * 'inflated'
%                               * name of a supported Freesurfer atlas:
%                                   - 'aparc' (refers to aparc)
%                                   - 'aparc2009' (refers to aparc.a2009s)
%                                   - 'yeo17' (refers to Yeo 17 networks)
%                                   - 'yeo7'  (refers to Yeo 7 networks)
% 
% dirFreesurfer: character full path to a location where Freesurfer is
%                           installed (assumes the "base" folder of
%                           Freesurfer, wherein other folders like
%                           "average", "subjects", "matlab", are located;
%                           only required if underlayType is of the
%                           Freesurfer atlas files
% 
% splitLR:      logical     yes or no indicating if the GIfTI file should
%                           be written out as separate left and right
%                           hemisphere files or an overall file 
%                           (default: true)
% 
%% Notes:
% When generating underlays, toWrite can be set to []; in this case, we
% still need the icoNum and outName to be provided
%
% Currently, the script does not downsample faces but picks the faces from
% the corresponding icsurfs entry
% 
% At the moment, not downsampling pial, inflated, or curvature files
% 
% For atlas files, supported ico numbers are: 3, 4, 5, 6, and 7

%% References:
% https://brainder.org/2016/05/31/downsampling-decimating-a-brain-surface/
% https://www.jiscmail.ac.uk/cgi-bin/wa-jisc.exe?A2=SPM;9b250ca1.1710

%% Check inputs
% Check ico number
if ~exist('icoNum', 'var') || isempty(icoNum)
    icoNum = 5;
else
    if icoNum < 0 || icoNum > 7
        error('icoNum should be between 0-7 (inclusive)');
    end
end

% Check if output names are provided
if ~exist('outName', 'var') || isempty(outName)
    genNames = true;
else
    genNames = false;
end

% Check if underlay type is specified
if exist('underlayType', 'var') && ~isempty(underlayType)
    underlayType = lower(underlayType);

    % First check if valid type of underlay is requested
    if ~ismember(underlayType, {'pial', 'curv', 'curvature', 'inflated', ...
                                'aparc', 'aparc.a2009s', 'yeo17', 'yeo7'})
        error(['Unknown underlayType specified: ', underlayType]);
    else
        genUnderlay = true;
        % If type is valid and an atlas is requested
        if ismember(underlayType, {'aparc', 'aparc.a2009s', 'yeo17', 'yeo7'})
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
        end
    end
else
    genUnderlay = false;
end

% Do we need to split hemispheres?
if ~exist('splitLR', 'var') || isempty(splitLR)
    splitLR = true;
else
    if ~islogical(splitLR)
        error('splitLR should be either true or false');
    end
end

% Finally check toWrite
if genUnderlay
    if exist('toWrite', 'var') && ~isempty(toWrite)
        warning('Since underlayType is specified, first input is ignored');
        toWrite = [];
    end
else
    if ~exist('toWrite', 'var') || isempty(toWrite)
        error('Please provide a vector or matrix of values to write');
    else
        [numRows, numVertices] = size(toWrite);
        if ~genNames
            if ischar(outName)
                outName = cellstr(outName);
            end
            if length(outName) ~= numRows
                error(['Mismatch between number of coefficients: ', num2str(numRows), ...
                       ' and number of output names: ', num2str(length(outName))]);
            end
        end
    end
end

%% Check for GIfTI toolbox
if ~exist('gifti.m', 'file')
    error(['Unable to find gifti.m file on MATLAB path; ', ...
           'please ensure that you have a copy of the GIfTI toolbox on your MATLAB path']);
end

%% Check for SurfView_surfs.mat
toLoad = 'SurfView_surfs.mat';
temp   = fullfile(fileparts(fileparts(which('writeGIfTI'))), 'showSurf');

if ~exist(fullfile(temp, toLoad), 'file')
    error(['Unable to find the file ', toLoad, ' in: ', temp]);
else
    surfaceData = load(fullfile(temp, toLoad));
end

%% Extract all faces from icsurfs (full size)
allFaces = [surfaceData.icsurfs{icoNum+1}.faces; ...
            surfaceData.icsurfs{icoNum+1}.faces + size(surfaceData.icsurfs{icoNum+1}.vertices,1)];

%% Handle undelays
if genUnderlay
    switch underlayType
        case 'pial'
            if splitLR
                left_vertices   = surfaceData.surf_lh_pial.vertices;
                right_vertices  = surfaceData.surf_rh_pial.vertices;
                left_faces      = surfaceData.surf_lh_pial.faces;
                right_faces     = surfaceData.surf_rh_pial.faces + size(surfaceData.surf_lh_pial.vertices,1);
            else
                allVertices = [surfaceData.surf_lh_pial.vertices; surfaceData.surf_rh_pial.vertices];
                allFaces    = [surfaceData.surf_lh_pial.faces; ...
                               surfaceData.surf_rh_pial.faces + size(surfaceData.surf_lh_pial.vertices,1)];
            end

        case {'curv', 'curvature'}
            if splitLR
                left_vertices   = surfaceData.curvvec_lh;
                right_vertices  = surfaceData.curvvec_rh;
                left_faces      = surfaceData.surf_lh_pial.faces;
                right_faces     = surfaceData.surf_rh_pial.faces + size(surfaceData.surf_lh_pial.vertices,1);
            else
                allVertices = [surfaceData.curvvec_lh; surfaceData.curvvec_rh];
                allFaces    = [surfaceData.surf_lh_pial.faces; ...
                               surfaceData.surf_rh_pial.faces + size(surfaceData.surf_lh_pial.vertices,1)];
            end

        case 'inflated'
            if splitLR
                left_vertices   = surfaceData.surf_lh_inflated.vertices;
                right_vertices  = surfaceData.surf_rh_inflated.vertices;
                left_faces      = surfaceData.surf_lh_inflated.faces;
                right_faces     = surfaceData.surf_lh_inflated.faces + size(surfaceData.surf_lh_inflated.vertices,1);
            else
                allVertices = [surfaceData.surf_lh_inflated.vertices; surfaceData.surf_rh_inflated.vertices];
                allFaces    = [surfaceData.surf_lh_inflated.faces; ...
                               surfaceData.surf_lh_inflated.faces + size(surfaceData.surf_lh_inflated.vertices,1)];
            end

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
            Data  = read_annotations(dirFreesurfer, fname, splitLR);
            if splitLR
                leftData  = Data(1);
                rightData = Data(2);
            end
    end

    % How many verticies to keep?
    % endLoc = resampleLocations(icoNum);

    % Downsample
    % if ~isempty(cData)
        % cData = cData(1:endLoc,:);
    % end

    % if ~ismember(underlayType, {'pial', 'inflated', 'curv', 'curvature'})
        % allVertices = allVertices(1:endLoc,:);
    % end

    % Output name
    if genNames
        outName = fullfile(pwd, underlayType);
    end

    % Prepare GIfTI object
    if ismember(underlayType, {'pial', 'inflated', 'curv', 'curvature'})
        if splitLR
            res = gifti(struct('faces', left_faces, 'vertices', left_vertices));
            save(res, [outName, '_lh'], 'GZipBase64Binary');
            
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
else
    % Size in icsurfs
    % numVertices_icsurfs = length(surfaceData.icsurfs{icoNum+1}.vertices);

    % Get underlay information
    allVertices = [surfaceData.icsurfs{icoNum+1}.vertices; surfaceData.icsurfs{icoNum+1}.vertices];
    % allFaces    = [surfaceData.icsurfs{icoNum+1}.faces; ...
    %                surfaceData.icsurfs{icoNum+1}.faces + size(surfaceData.icsurfs{icoNum+1}.vertices,1)];

    % Sanity check
    if length(allVertices) ~= numVertices
        error(['Mismatch between number of vertices in the input: ', num2str(numVertices), ...
               ' and the number of vertices at ico number ', num2str(icoNum), ...
               ' : ', num2str(length(allVertices))]);
    end

    % Generate output names, if required
    if genNames
        outName = fullfile(pwd, strrep(strcat({'FEMA_estimate_'}, num2str((1:numRows)')), ' ', ''));
    end

    % Loop over every coefficient and start saving!
    for rows = 1:numRows
        if splitLR
            res = gifti(struct('cdata', squeeze(toWrite(rows,1:length(surfaceData.icsurfs{icoNum+1}.vertices))')));
            save(res, [outName{rows}, '_lh'], 'GZipBase64Binary');
    
            res = gifti(struct('cdata', squeeze(toWrite(rows,length(surfaceData.icsurfs{icoNum+1}.vertices)+1:end)')));
            save(res, [outName{rows}, '_rh'], 'GZipBase64Binary');
        else
            res = gifti(struct('cdata', squeeze(toWrite(rows,:)')));
            save(res, outName{rows}, 'Base64Binary');
        end
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
    if icoNum < 3
        error('Supported ico numbers for atlas files are: 3, 4, 5, 6, and 7');
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
    error(['Unable to find [l/r]h.aparc.annot file in: ', dir_annotations]);
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
end

function endLoc = resampleLocations(icoNum)
% Trivial function that gives us how many vertices to pick based on
% icosahedron number
switch icoNum
    case 0
        endLoc = 12;

    case 1
        endLoc = 42;
        
    case 2
        endLoc = 162;

    case 3
        endLoc = 642;
        
    case 4
        endLoc = 2562;
        
    case 5
        endLoc = 10242;

    case 6
        endLoc = 40962;

    case 7
        endLoc = 163842;
end
end