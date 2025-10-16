function writeGIfTI(toWrite, icoNum, outDir, basename, colnames, underlayType, dirFreesurfer, splitLR)
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
% outDir:       character   full path to where the GIfTI files should be
%                           written out
% 
% basename:     character   the basename to be used for each of the p GIfTI
%                           images; if empty, defaults to 'FEMA_estimate'
% 
% colnames:     [p x 1]     cell type having the name of each of the
%                           coefficients in toWrite; the names are appended
%                           to the basename to create the name of the
%                           output file; if left empty, '_1', '_2', and so
%                           on are appended to the basename
% 
% underlayType  character   a value indicating which of the following
%                           underlays should be generated as a GIfTI image:
%                               * 'pial'
%                               * 'curv' | 'curvature'
%                               * 'inflated'
%                               * name of a supported Freesurfer atlas:
%                                   - 'aparc' (refers to aparc2009)
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
% splitLR:      logical     yes or no indicating if the underlay should be
%                           written out as separate left and right
%                           hemisphere files or an overall file (only if
%                           underlay is requested; default: false)

%% Notes:
% When generating underlays, toWrite can be set to []; in this case, we
% still need the icoNum, outDir, and basename to be provided; if colnames
% is provided, it is used; otherwise, the 'basename_underlay_icoNum' is
% used as the output name
%
% Currently, the script does not downsample faces but picks the faces from
% the corresponding icsurfs entry
% 
% At the moment, not downsampling pial, inflated, or curvature files
% 
% splitLR needs to be implemented for inflated or atlases

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

% Check output directory and make it if it doesn't exist
if ~exist('outDir', 'var') || isempty(outDir)
    error('Please provide a full path to where the GIfTI images should be saved');
else
    if ~exist(outDir, 'dir')
        mkdir(outDir);
    end
end

% Check output basename
if ~exist('basename', 'var') || isempty(basename)
    basename = 'FEMA_estimate';
end

% Check column names of estimates
if ~exist('colnames', 'var') || isempty(colnames)
    genColNames = true;
else
    genColNames = false;
    if ~iscell(colnames)
        colnames = cellstr(colnames);
    else
        % Ensure that this is a p x 1 vector
        colnames = reshape(colnames, numel(colnames), 1);
    end
end

% Check if underlay type is specified
if exist('underlayType', 'var') && ~isempty(underlayType)
    underlayType = lower(underlayType);

    % First check if valid type of underlay is requested
    if ~ismember(underlayType, {'pial', 'curv', 'curvature', 'inflated', 'aparc', 'yeo17', 'yeo7'})
        error(['Unknown underlayType specified: ', underlayType]);
    else
        genUnderlay = true;
        % If type is valid and an atlas is requested
        if ismember(underlayType, {'aparc', 'yeo17', 'yeo7'})
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
    splitLR = false;
else
    if genUnderlay
        if ~islogical(splitLR)
            error('splitLR should be either true or false');
        end
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
            allVertices = [surfaceData.surf_lh_pial.vertices; surfaceData.surf_rh_pial.vertices];
            allFaces    = [surfaceData.surf_lh_pial.faces; ...
                           surfaceData.surf_rh_pial.faces + size(surfaceData.surf_lh_pial.vertices,1)];
            cData       = [];

        case {'curv', 'curvature'}
            allVertices = [surfaceData.curvvec_lh; surfaceData.curvvec_rh];
            allFaces    = [surfaceData.surf_lh_pial.faces; ...
                          surfaceData.surf_rh_pial.faces + size(surfaceData.surf_lh_pial.vertices,1)];
            cData       = [];

        case 'inflated'
            allVertices = [surfaceData.surf_lh_inflated.vertices; surfaceData.surf_rh_inflated.vertices];
            % allFaces    = [surfaceData.surf_lh_inflated.faces; ...
            %               surfaceData.surf_lh_inflated.faces + size(surfaceData.surf_lh_inflated.vertices,1)];
            cData       = [];

        case 'aparc'
            fname       = 'aparc.a2009s';
            cData       = read_annotations(dirFreesurfer, fname);
            allVertices = [surfaceData.icsurfs{icoNum+1}.vertices; surfaceData.icsurfs{icoNum+1}.vertices];

        case 'yeo17'
            fname       = 'Yeo2011_17Networks_N1000';
            cData       = read_annotations(dirFreesurfer, fname);
            allVertices = [surfaceData.icsurfs{icoNum+1}.vertices; surfaceData.icsurfs{icoNum+1}.vertices];
            
        case 'yeo7'
            fname       = 'Yeo2011_7Networks_N1000';
            cData       = read_annotations(dirFreesurfer, fname);
            allVertices = [surfaceData.icsurfs{icoNum+1}.vertices; surfaceData.icsurfs{icoNum+1}.vertices];
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
    if genColNames
        colnames = cellstr([underlayType, '_', num2str(icoNum)]);
    end    
    outNames = strcat({basename}, {'_'}, colnames, {'.gii'});

    % Prepare GIfTI object
    % if isempty(cData)
    %     res = gifti(struct('faces', allFaces, 'vertices', allVertices));
    % else
    %     res = gifti(struct('faces', allFaces, 'vertices', allVertices, 'cdata', cData));
    % end
    % save(res, fullfile(outDir, outNames{1}), 'Base64Binary');
    if isempty(cData)
        res = gifti(struct('faces', allFaces, 'vertices', allVertices));
        save(res, fullfile(outDir, outNames{1}), 'Base64Binary');
    else
        res = gifti(struct('cdata', cData(1:length(surfaceData.icsurfs{icoNum+1}.vertices))));
        save(res, fullfile(outDir, ['lh.', outNames{1}]), 'Base64Binary');

        res = gifti(struct('cdata', cData(length(surfaceData.icsurfs{icoNum+1}.vertices)+1:end)));
        save(res, fullfile(outDir, ['rh.', outNames{1}]), 'Base64Binary');
    end
    
else
    % Size of inputs
    [numCoeff, numVertices] = size(toWrite);

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
    if genColNames
        colnames = strrep(strcat({'FEMA_estimate_'}, num2str((1:numCoeff)')), ' ', '');
    end
    outNames = strcat(basename, {'_'}, colnames, {'.gii'});

    if length(colnames) ~= numCoeff
        error(['Mismatch between number of coefficients in the input: ', num2str(numCoeff), ...
            ' and the number of entries in colnames: ', num2str(length(colnames))]);
    else
        % Loop over every coefficient and start saving!
        for coeff = 1:numCoeff
            res = gifti(struct('cdata', squeeze(toWrite(coeff,1:length(surfaceData.icsurfs{icoNum+1}.vertices))')));
            save(res, fullfile(outDir, ['lh.', outNames{coeff}]), 'Base64Binary');

            res = gifti(struct('cdata', squeeze(toWrite(coeff,length(surfaceData.icsurfs{icoNum+1}.vertices)+1:end)')));
            save(res, fullfile(outDir, ['rh.', outNames{coeff}]), 'Base64Binary');

            % res = gifti(struct('cdata', squeeze(toWrite(coeff,:)')));
            % save(res, fullfile(outDir, outNames{coeff}), 'Base64Binary');
        end
    end
end
end


function cData = read_annotations(dirFreesurfer, fname)
% Function that reads an annotation file from Freesurfer and returns
% concatenated list of cData

% Path to annotation files and path to script file
dir_annotations = fullfile(dirFreesurfer, 'subjects', 'fsaverage', 'label');
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

% Return concatenated vertices and label
% all_vertices = [lh_vertices; rh_vertices];
cData = [lh_label; rh_label];
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