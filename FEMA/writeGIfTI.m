function writeGIfTI(toWrite, icoNum, outName, splitLR)
% Function that writes out GIfTI file(s) given a set of inputs
% Requires the GIfTI toolbox: https://github.com/gllmflndn/gifti to be on
% path such that the function 'gifti' is callable
% 
%% Inputs:
% --------
% toWrite:      [p x v]     vector or matrix of estimates for p
%                           coefficients and v vertices
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
%                          of the output GIfTI files (no extensions); if
%                          empty, uses pwd and generates file names
% 
% splitLR:      logical     yes or no indicating if the GIfTI file should
%                           be written out as separate left and right
%                           hemisphere files or an overall file 
%                           (default: true)
% 
%% References:
% https://brainder.org/2016/05/31/downsampling-decimating-a-brain-surface/
% https://www.jiscmail.ac.uk/cgi-bin/wa-jisc.exe?A2=SPM;9b250ca1.1710
%
%% Requirement(s):
% GIfTI toolbox: https://github.com/gllmflndn/gifti

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

% Do we need to split hemispheres?
if ~exist('splitLR', 'var') || isempty(splitLR)
    splitLR = true;
else
    if ~islogical(splitLR)
        error('splitLR should be either true or false');
    end
end

% Finally check toWrite
if ~exist('toWrite', 'var') || isempty(toWrite)
    error('Please provide a vector or matrix of values to write');
else
    [numRows, numVertices] = size(toWrite); %#ok<ASGLU>
    if ~genNames
        if ischar(outName)
            outName = cellstr(outName);
        end
        if length(outName) ~= numRows
            error(['Mismatch between number of coefficients: ', num2str(numRows), ...
                   ' and number of output names: ', num2str(length(outName))]);
        else
            % Strip extensions out of outName, if provided
            [~, ~, ext] = fileparts(outName);
            outName     = strrep(outName, ext, '');
        end
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

%% Generate output names, if required
if genNames
    outName = fullfile(pwd, strrep(strcat({'FEMA_estimate_'}, num2str((1:numRows)')), ' ', ''));
end

%% Loop over every coefficient and start saving!
for rows = 1:numRows
    if splitLR
        res = gifti(struct('cdata', squeeze(toWrite(rows,1:numIcoVertices(icoNum))')));
        save(res, [outName{rows}, '_lh'], 'GZipBase64Binary'); %#ok<USENS>
    
        res = gifti(struct('cdata', squeeze(toWrite(rows,numIcoVertices(icoNum)+1:end)')));
        save(res, [outName{rows}, '_rh'], 'GZipBase64Binary');
    else
        res = gifti(struct('cdata', squeeze(toWrite(rows,:)')));
        save(res, outName{rows}, 'GZipBase64Binary');
    end
end
end