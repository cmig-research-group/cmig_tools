function lookup = FEMA_lookupVertices(stats, zValues, seValues, coeff, pThresh, dir_Freesurfer)
% Function to generate vertex-wise lookup tables, given a set of
% FEMA-generated statistics
%% Inputs:
% stats:            [n x v]     beta_hat, sig2mat, or sig2mat .* sig2tvec 
%                               output by FEMA having n number of 
%                               coefficients and v number of vertices
%
% zValues:          [n x v]     Z scores output by FEMA (optional; used for
%                               generating p values); should be empty if 
%                               stats are variance components
%
% seValues:         [n x v]     SE values output by FEMA (optional); should
%                               be empty if stats are variance components
% 
% coeff:            [1 x k]     vector of coefficients of interest (where k 
%                               is a subset of n)
% 
% pThresh:          [1 x 1]     number indicating the p value threshold 
%                               (only used when stats are beta values);
%                               can also be 'none'
% 
% dir_Freesurfer:               full path to a location where the atlas
%                               files i.e., 'lh.aparc.annot' and 
%                               'rh.aparc.annot' are saved
%
%% Output(s):
% lookup:           [k x 1]     cell type having statistics and lookup
%                               information for the k coefficients; the 
%                               contents of the cell are table type
%
%% Notes:
% It is critical that the input stats variable has the correct number of
% vertices as this number is used to subset the vertices in the annotation
% files - v/2 number of vertices are considered
%
% Performs labeling using the DK40 atlas
%
% If dir_Freesurfer is passed, then the Freesurfer function 
% 'read_annotation' should be on the MATLAB path
% 
% p value threshold is applied using a less than operation; specifically, 
% the two sided p-values (if Z values are provided) are calculated as:
% pvals = normcdf(-abs(zValues(coeff,:)))*2;
% and then all p values which are less than the p value threshold are
% considered
% 
% If no threshold operation is applied and/or random effects are being
% summarized, then it is possible that one or more vertices do not have a
% valid lookup label
%
%% Defaults:
% zValues:          []
% coeff:            1:v
% pThresh:          'none'
% dir_Freesurfer:   automatically generated using $FREESURFER_HOME

%% Check inputs and assign defaults
% Check stats
if ~exist('stats', 'var') || isempty(stats)
    error('Please provide an n x v dimensional stats matrix');
else
    [n, v] = size(stats);
end

% Check zValues
if ~exist('zValues', 'var') || isempty(zValues)
    zValues = [];
    doFFX   = false;
else
    [n_z, v_z] = size(zValues);
    if n_z ~= n || v_z ~= v
        error('Mismatch in size of stats and zValues');
    else
        doFFX = true;
    end
end

% Check SE values
if ~exist('seValues', 'var') || isempty(seValues)
    seValues = NaN(size(stats));
else
    if isempty(zValues)
        warning('Disregarding SE values as Z values are not specified');
        seValues = [];
    else
        [n_se, v_se] = size(seValues);
        if n_se ~= n || v_se ~= v
            error('Mismatch in size of stats and SE values');
        end
    end
end

% Check coeff
if ~exist('coeff', 'var') || isempty(coeff)
    coeff    = 1:n;
else
    if length(intersect(1:n, coeff)) ~= length(coeff)
        error('One or more coeff values exceed size of stats');
    end
end
numCoeff = length(coeff);

% Check pThresh
if ~exist('pThresh', 'var') || isempty(pThresh)
    pThresh = 'none';
else
    if ischar(pThresh)
        if ~strcmpi(pThresh, 'none')
            error(['Unknown p value threshold provided: ', pThresh]);
        end
    else
        if isempty(zValues)
            warning('Z scores not specified; disregarding p value threshold');
            pThresh = 'none';
        end
    end
end

% Specify names of annotation files
fname_lh_annot = 'lh.aparc.annot';
fname_rh_annot = 'rh.aparc.annot';
fname_script   = 'read_annotation.m';

% Check dir_Freesurfer
if ~exist('dir_Freesurfer', 'var') || isempty(dir_Freesurfer)
    
    % Call environemnt variable
    [~, dir_Freesurfer] = system('$FREESURFER_HOME');
    
    % If there is a colon on the path (for example, ":permission denied",
    % strip it out of the path name
    [dir_Freesurfer, ~] = strsplit(dir_Freesurfer, ':');
    dir_Freesurfer      = dir_Freesurfer{1};
    
    % Path to annotation files and path to script file
    dir_annotations = fullfile(dir_Freesurfer, 'subjects', 'fsaverage', 'label');
    dir_script      = fullfile(dir_Freesurfer, 'matlab');
    
    % Make sure that annotation files and read_annotation file exist
    if ~exist(fullfile(dir_annotations, fname_lh_annot), 'file') || ...
       ~exist(fullfile(dir_annotations, fname_rh_annot), 'file')
            error(['Unable to find [l/r]h.aparc.annot file in: ', dir_annotations]);
    else
        if ~exist(fullfile(dir_script, fname_script), 'file')
            error(['Unable to find ', fname_script, ' in: ', dir_script]);
        end
    end
else
    % Ensure this directory exists
    if ~exist(dir_Freesurfer, 'dir')
        error(['Unable to find the directory: ', dir_Freesurfer]);
    else
        dir_annotations = dir_Freesurfer;
        
        % Ensure that annotation files exist
        if ~exist(fullfile(dir_annotations, fname_lh_annot), 'file') || ...
           ~exist(fullfile(dir_annotations, fname_rh_annot), 'file')
                error(['Unable to find [l/r]h.aparc.annot file in: ', dir_annotations]);
        else
            % Ensure that read_annotation is accessible
            dir_script = fileparts(which(fname_script));
            if ~exist(dir_script, 'dir')
                error(['Could not find: ', fname_script, '; ensure that the file is on MATLAB path']);
            end
        end
    end
end

%% Load annotation files
tmp = pwd;
cd(dir_script);
[~, lh_label, lh_colortable] = read_annotation(fullfile(dir_annotations, fname_lh_annot));
[~, rh_label, rh_colortable] = read_annotation(fullfile(dir_annotations, fname_rh_annot));
cd(tmp);

%% Label every vertex
% Left hemisphere
lh_labelString = cell(length(lh_label),1);
for clr  = 1:lh_colortable.numEntries
    lh_labelString(lh_label == lh_colortable.table(clr,5)) = strcat({'lh-'}, lh_colortable.struct_names(clr));
end

% Right hemisphere
rh_labelString = cell(length(rh_label),1);
for clr  = 1:rh_colortable.numEntries
    rh_labelString(rh_label == rh_colortable.table(clr,5)) = strcat({'rh-'}, rh_colortable.struct_names(clr));
end

% Truncate left and right hemisphere label string
% Note the half!
trunc_lh_labelString = lh_labelString(1:v/2);
trunc_rh_labelString = rh_labelString(1:v/2);

% Put labels together - left followed by right (this is based on the mask
% information where the 9th vertex is set to zero; this is true for the
% left hemisphere but not for the right hemisphere)
all_labels = [trunc_lh_labelString; trunc_rh_labelString];

%% Create lookup tables
% Initialize
lookup    = cell(numCoeff, 1);

if doFFX
    % Working on fixed effects - results will include stats values, standard
    % error, p values, Z values, and the lookup labels
    for cc = 1:numCoeff
        pvals       = normcdf(-abs(zValues(coeff(cc),:)))*2;
        if ischar(pThresh)
            loc     = true(size(zValues(coeff(cc), :), 2),1);
        else
            loc     = pvals < pThresh;
        end
        lookup{cc}  = cell2table([num2cell(stats(coeff(cc),loc))',      ...
                                  num2cell(seValues(coeff(cc),loc))',   ...
                                  num2cell(pvals(loc))',                ...
                                  num2cell(zValues(coeff(cc),loc))',    ...
                                  all_labels(loc)], ...
                                  'VariableNames', {'Stats', 'SE', 'pValue', 'Zvalue', 'Label'});
    end
else
    % Working on random effects - results will include stats and labels
    for cc = 1:numCoeff
        lookup{cc}  = cell2table([num2cell(stats(coeff(cc),:))', all_labels], ...
                                  'VariableNames', {'Stats', 'Label'});
    end
end