function caller_FEMA_fit(file_X, file_ymat, dirOutput, outPrefix, varargin)
% Caller function for FEMA_fit (compiled version)
% 
%% General specificaton
% For inputs X and y (and associated inputs), rows are assumed to be
% observations and columns are assumed to be variables
% 
% Generally, across all inputs, fid, iid, and eid are assumed to be string
% variables (cell string)
% 
%% To Do:
% Support passing in FamilyStruct
% Regarding contrasts: allow re-analysis on existing models
% Save summary: n obs, ny, Rsquared, singularity, df, model fit?
% Should we separate out permutation results? (problematic)
% Support NIfTI (if variable of interest is passed, then do only those;
%               default: all)
% 
%% Mandatory inputs:
% ------------------
% file_X:           full path to a .mat or a tab-delimited .csv file
%                   containing the design matrix and other information
% 
%                     - if .mat, it should contain the following variables:
%                           * X:        the actual design matrix
%                           * fid:      the family ID for every observation
%                           * iid:      the individual ID for every observation
%                           * eid:      the event ID for every observation
%                           * agevec:   the age values for every observation
%                           * colnames: the names of each column in X
% 
%                     - if .csv, the columns should be named as (in order):
%                           * fid:      the family ID for every observation
%                           * iid:      the individual ID for every observation
%                           * eid:      the event ID for every observation
%                           * agevec:   the age values for every observation
%                           * column 5 onwards should be design matrix variables
%
% file_ymat:        full path to a .mat or a tab-delimited .csv file
%                   containing the outcome variables for which the model
%                   should be fitted (these will be intersected with the
%                   design matrix)
% 
%                     - if .mat, it should contain the following variables:
%                           * iid:      the individual ID for every observation
%                           * eid:      the event ID for every observation
%                           * ymat:     the matrix of outcome variables
% 
%                     - if .csv, the columns should be named as (in order): 
%                           * iid:      the individual ID for every observation
%                           * eid:      the event ID for every observation
%                           * column 3 onwards should be the outcome variables
%
% dir_output:       full path to where the results will be written out
% 
% outPrefix:        prefix to the output file names (without extension;
%                   default: FEMA_fit-yyyyMMMdd-HHmmSS)
%
%% Optional inputs (comma-separated name-value pairs):
% ----------------------------------------------------
% 
% niter:            number of iterations (default: 1)
%
% contrasts:        full path to a tab-delimited csv file that specifies
%                   the contrast (default: [], i.e., no contrasts): 
%                     - first row: name of the variables in the design
%                                  matrix: if file_X was a mat file, these
%                                  names should correspond to colnames; if
%                                  file_X was a csv file, these names
%                                  should match the names of the columns
%                                  from column 5 onwards
%                     - second row onwards should be weight values
%                     - first col: 'Type' which specifies 'u' or 'm'
%                                   followed by a number indicating which
%                                   rows are univariate contrasts and which
%                                   are multivariate contrasts
%                   See, FEMA_parse_contrastFile for specification of the
%                   contrast file
%
% nbins:            bin spacing (default: 20); set to zero to disable
%                   binning; the nbins parameter controls which y variables
%                   will be "binned" when performing the GLS estimation for
%                   the fixed effects; a bin spacing of 20 translates to a 
%                   spacing of 1/20 = 0.05; outcome variables with
%                   estimated random effects in this range will be binned
%                   together for GLS estimation
% 
% GRMfile:          full path to a .mat or .dat file that specifies how
%                   genetically related are the individuals in the study
%                   (default: [], i.e., no GRM file is specified)
%                     - if .mat, it should have the following variables:
%                           * GRM:            matrix of genetic relatedness
%                           * uqObservations: ordering of individuals in the GRM
% 
%                     - if .dat, it should have a a corresponding .csv or
%                       .mat file with the same basename in the same
%                       location; if .csv, this file should be a single
%                       column (no column name) which contains the
%                       iid_list; if .mat, this file should contain a
%                       variable 'uqObservations'
%
% RandomEffects:    list of random effects to estimate (default: {'F', 'S', 'E'}); 
%                   the following random effects are supported:
%                       * F:  family relatedness
%                       * S:  subject effect
%                       * E:  error - always required, added if missing
%                       * A:  additive genetic relatedness - must include
%                             file path to GRM file
%                       * D:  dominant genetic relatedness - square of A
%                       * M:  maternal effect - effect of having same mother
%                       * P:  paternal effect - effect of having same father
%                       * H:  home effect - effect of living at the same address
%                       * T:  twin effect - effect of having the same pregnancy ID
%
%                   If random effect includes 'M', 'P', 'H', or 'T'
%                   effects, an additional tab-delimited two-column csv
%                   file is required for each of these:
%                       * For 'M' effect:
%                           * iid:          individual ID
%                           * MotherID:     mother ID
%                       * For 'P' effect:
%                           * iid:          individual ID
%                           * FatherID:     father ID
%                       * For 'H' effect:
%                           * iid:          individual ID
%                           * HomeID:       address/home ID
%                       * For 'T' effect:
%                           * iid:          individual ID
%                           * PregID:       pregnancy ID
%
%                   The above 'M', 'P', 'H', and 'T' input files should be
%                   passed in as name-value pairs 'MotherID', 'FatherID',
%                   'HomeID', and 'PregID' respectively (which contain the
%                   above columns)
%  
% nperms:           number of permutations to run (default: 0, i.e., no
%                   permutations will be performed)
%
% CovType:          should be one of the following (default: analytic):
%                       * 'analytic':       compound symmetry
%                       * 'unstructured':   unstructured covariance
% 
% FixedEstType:     should be one of the following (default: GLS):
%                       * 'GLS':    generalised least squares
%                       * 'OLS':    ordinary least squares
% 
% RandomEstType:    should be one of the following (default: MoM):
%                       * 'MoM':    method of moments
%                       * 'ML':     maximum likelihood
% 
% GroupByFamType:   should grouping by family type be used (default: true);
%                   only relevant for F, S, and E random effects
% 
% NonnegFlag:       apply non-negativity constraint on random effects
%                   estimation (default: true)
% 
% SingleOrDouble:   control the numerical precision; should be one of the
%                   following (default: double): 
%                       * 'double'
%                       * 'single'
% 
% logLikflag:       flag to computate log-likelihood (default: false)
% 
% PermType:         controls the type of permutation; should be one of:
%                       * 'wildbootstrap':    null bootstrap: creates null
%                                             distribution by randomly
%                                             flipping the sign of each
%                                             observation
%                       * 'wildbootstrap-nn': non-null boostrap: estimates
%                                             distribution around effect of
%                                             interest using sign flipping
%                                             (used for sobel test)
% 
% returnReusable:   controls return of additional variables (default:
%                   true); if true, returns a structure with some
%                   variables that can be reused (primarily by FEMA-GWAS)
% 
% doPar:            should parallel processing be used for GLS solution
%                   (only applicable if CovType is unstructured; default:
%                   false)
% 
% numWorkers:       if doPar is true, how many parallel workers to start
%                   (default: 2) 
% 
% numThreads:       if doPar is true, how many threads per parallel worker
%                   (default: 2)
% 
% saveDesignMatrix: if true, filtered X variable and ID lists are written
%                   out as a separate mat file (default: false)
% 
%% Start
updateString = [char(datetime('now')), ': caller_FEMA_fit: job started'];
disp(updateString);

%% Check mandatory inputs
if ~exist('file_X', 'var') || isempty(file_X)
    error(['Please provide a full path to a file ', ...
           'containing the design matrix and other related information']);
else
    if ~exist(file_X, 'file')
        error(['Unable to find file: ', file_X]);
    end
end

if ~exist('file_ymat', 'var') || isempty(file_ymat)
    error(['Please provide a full path to a file ', ...
           'containining the outcome variables to be analysed']);
else
    if ~exist(file_ymat, 'file')
        error(['Unable to find file: ', file_ymat]);
    end
end

if ~exist('dirOutput', 'var') || isempty(dirOutput)
    error('Please provide a full path to where results should be saved');
else
    if ~exist(dirOutput, 'dir')
        mkdir(dirOutput);
    end
end

if ~exist('outPrefix', 'var') || isempty(outPrefix)
    outPrefix = ['FEMA_fit-', char(datetime('now', 'Format', 'yyyyMMMdd-HHmmSS'))];
end

%% Assign default optional inputs
p = inputParser;
addParameter(p, 'niter', 1);
addParameter(p, 'contrasts', []);
addParameter(p, 'nbins', 20);
addParameter(p, 'GRMFile', [])
addParameter(p, 'RandomEffects', {'F' 'S' 'E'});
addParameter(p, 'nperms', 0);
addParameter(p, 'CovType', 'analytic');
addParameter(p, 'FixedEstType', 'GLS');
addParameter(p, 'RandomEstType', 'MoM');
addParameter(p, 'GroupByFamType', true);
addParameter(p, 'NonnegFlag', true);
addParameter(p, 'SingleOrDouble', 'double');
addParameter(p, 'logLikflag', false);
addParameter(p, 'PermType', 'wildbootstrap');
addParameter(p, 'returnReusable', true);
addParameter(p, 'doPar',false);
addParameter(p, 'numWorkers', 2);
addParameter(p, 'numThreads', 2);
addParameter(p, 'MotherID', []);
addParameter(p, 'FatherID', []);
addParameter(p, 'HomeID', []);
addParameter(p, 'PregID', []);
addParameter(p, 'saveDesignMatrix', false);

%% Parse optional inputs
parse(p, varargin{:})
contrasts        = p.Results.contrasts;
GRMFile          = p.Results.GRMFile;
CovType          = p.Results.CovType;
FixedEstType     = p.Results.FixedEstType;
RandomEstType    = p.Results.RandomEstType;
GroupByFamType   = p.Results.GroupByFamType;
NonnegFlag       = p.Results.NonnegFlag;
SingleOrDouble   = p.Results.SingleOrDouble;
logLikflag       = p.Results.logLikflag;
PermType         = p.Results.PermType;
returnReusable   = p.Results.returnReusable;
doPar            = p.Results.doPar;
MotherID         = p.Results.MotherID;
FatherID         = p.Results.FatherID;
HomeID           = p.Results.HomeID;
PregID           = p.Results.PregID;
saveDesignMatrix = p.Results.saveDesignMatrix;

if ~isnumeric(p.Results.niter)
    niter = str2double(p.Results.niter);
else
    niter = p.Results.niter;
end

if ~isnumeric(p.Results.nbins)
    nbins = str2double(p.Results.nbins);
else
    nbins = p.Results.nbins;
end

if ~isnumeric(p.Results.nperms)
    nperms = str2double(p.Results.nperms);
else
    nperms = p.Results.nperms;
end

if ~isnumeric(p.Results.numWorkers)
    numWorkers = str2double(p.Results.numWorkers);
else
    numWorkers = p.Results.numWorkers;
end

if ~isnumeric(p.Results.numThreads)
    numThreads = str2double(p.Results.numThreads);
else
    numThreads = p.Results.numThreads;
end

if isdeployed
    RandomEffects = strrep(strsplit(p.Results.RandomEffects, ','), ' ', '');
else
    RandomEffects = p.Results.RandomEffects;
end

%% Make some decisions based on optional inputs
if strcmpi(CovType, 'unstructured')
    unstructuredCov = true;
else
    unstructuredCov = false;
end

% If permutation and unstructured covariance, warn the user
if unstructuredCov && nperms > 0
    warning('Permutations not yet implemented for unstructured covariance');
    nperms = 0;
end

% Examine RandomEffects and ensure E is always the last term - relevant for
% unstructured covariance
RandomEffects = rowvec(RandomEffects);
tmp           = strcmpi(RandomEffects, 'E');
if ~any(tmp)
    warning('RandomEffects did not include E term; appending E as the last random effect');
    RandomEffects = [RandomEffects, 'E'];
else
    if find(tmp) ~= length(RandomEffects)
        RandomEffects = [RandomEffects(~tmp), RandomEffects(tmp)];
        logging(['Re-arranging RandomEffects as: ', sprintf('%s ', RandomEffects{:})]);
    end
end

% Grouping by family type is only supported for RandomEffects 'F' 'S' 'E'
% but not supported for unstructured covariance
if ~isempty(setdiff(RandomEffects, {'F' 'S' 'E'})) || unstructuredCov
    GroupByFamType = false;
end

updateString = [char(datetime('now')), ': caller_FEMA_fit: finished input parsing'];
disp(updateString);

%% Load X
updateString = [char(datetime('now')), ': caller_FEMA_fit: loading X'];
disp(updateString);

[~, ~, ext] = fileparts(file_X);
if isempty(ext)
    error('Please provide full path to file_X');
else
    if strcmpi(ext, '.mat')
        temp_X = load(file_X);
        
        % Check to make sure mandatory variables are present
        toCheck_X = {'X', 'fid', 'iid', 'eid', 'agevec', 'colnames'};
        chk_X     = ismember(toCheck_X, fieldnames(temp_X));
        if sum(chk_X) ~= length(toCheck_X)
            missMsg = sprintf('%s, ', toCheck_X{~chk_X});
            error(['The following required variables are missing from ', file_X, ': ', missMsg(1:end-2)]);
        else
            % Assign the right variables
            fid      = colvec(temp_X.fid);
            iid      = colvec(temp_X.iid);
            eid      = colvec(temp_X.eid);
            agevec   = colvec(temp_X.agevec);
            colnames = temp_X.colnames;
            X        = temp_X.X;
        end
    else
        if strcmpi(ext, '.csv')
            temp_X      = readtable(file_X, 'Delimiter', '\t');
            toCheck_X   = {'fid', 'iid', 'eid', 'agevec'};
            temp_names  = temp_X.Properties.VariableNames;
            chk_X       = ismember(toCheck_X, temp_names);
            if sum(chk_X) ~= length(toCheck_X)
                missMsg = sprintf('%s, ', toCheck_X{~chk_X});
                error(['The following required variables are missing from ', file_X, ': ', missMsg(1:end-2)]);
            else
                % Assign the right variables
                fid      = temp_X.fid;
                iid      = temp_X.iid;
                eid      = temp_X.eid;
                agevec   = temp_X.agevec;
                colnames = temp_names(5:end);
                X        = temp_X{:,5:end};            
            end
        else
            error(['Unknown extension for file_X: ', ext]);
        end
    end
end
clear temp_X

%% Load ymat
updateString = [char(datetime('now')), ': caller_FEMA_fit: loading ymat'];
disp(updateString);

[~, ~, ext] = fileparts(file_ymat);
if isempty(ext)
    error('Please provide full path to file_ymat');
else
    if strcmpi(ext, '.mat')
        temp_ymat = load(file_ymat);
    
        % Check to make sure mandatory variables are present
        toCheck_y = {'ymat', 'iid', 'eid'};
        chk_y     = ismember(fieldnames(temp_ymat), toCheck_y);
        if sum(chk_y) ~= length(toCheck_y)
            missMsg = sprintf('%s, ', toCheck_y{~chk_y});
            error(['The following required variables are missing from ', file_ymat, ': ', missMsg(1:end-2)]);
        else
            % Assign the right variables
            iid_y = colvec(temp_ymat.iid);
            eid_y = colvec(temp_ymat.eid);
            ymat  = temp_ymat.ymat;
        end
    else
        if strcmpi(ext, '.csv')
            temp_ymat   = readtable(file_ymat, 'Delimiter', '\t');
            toCheck_y   = {'iid', 'eid'};
            temp_names  = temp_ymat.Properties.VariableNames;
            chk_y       = ismember(toCheck_y, temp_names);
            if sum(chk_y) ~= length(toCheck_y)
                missMsg = sprintf('%s, ', toCheck_y{~chk_y});
                error(['The following required variables are missing from ', file_ymat, ': ', missMsg(1:end-2)]);
            else
                % Assign the right variables
                iid_y       = temp_ymat.iid;
                eid_y       = temp_ymat.eid;
                ymat        = temp_ymat{:,3:end};
                % colnames_y  = temp_names(:,3:end);
            end
        else
            error(['Unknown extension for file_ymat: ', ext]);
        end
    end
end

clear temp_ymat

%% Save original number of observations
ori_nObs_X    = size(X, 1);
ori_nObs_ymat = size(ymat, 1);

%% Intersect data
updateString = [char(datetime('now')), ': caller_FEMA_fit: intersecting data'];
disp(updateString);

% Create an ID list for X and ymat
if isnumeric(eid)
    warning('Converting eid to string type (X)');
    eid = cellstr(num2str(eid));
end
idList_X = strcat(iid, '_', eid);

if isnumeric(eid_y)
    warning('Converting eid to string type (ymat)');
    eid_y = cellstr(num2str(eid_y));
end
idList_ymat = strcat(iid_y, '_', eid_y);

[~, i_X, i_ymat] = intersect(idList_X, idList_ymat, 'stable');

% Only keep common data points
fid         = fid(i_X);
iid         = iid(i_X);
eid         = eid(i_X);
agevec      = agevec(i_X);
X           = X(i_X, :);

ymat        = ymat(i_ymat, :);
% iid_y       = iid_y(i_ymat);
% eid_y       = eid_y(i_ymat);

% Get updated information
[nObs_X, nCols_X] = size(X);
nObs_ymat         = size(ymat,1);

%% Display report before and after intersection
updateString = [char(datetime('now')), ': caller_FEMA_fit: ', num2str(ori_nObs_X), ' observations in X before intersection'];
disp(updateString);

updateString = [char(datetime('now')), ': caller_FEMA_fit: ', num2str(nObs_X), ' observations in X after intersection'];
disp(updateString);

updateString = [char(datetime('now')), ': caller_FEMA_fit: ', num2str(ori_nObs_ymat), ' observations in ymat before intersection'];
disp(updateString);

updateString = [char(datetime('now')), ': caller_FEMA_fit: ', num2str(nObs_ymat), ' observations in ymat after intersection'];
disp(updateString);

%% Perform a check for NaN and Inf
if logical(sum(any(isnan(X)))) || logical(sum(any(isnan(ymat)))) || ...
   logical(sum(any(isinf(X)))) || logical(sum(any(isinf(ymat))))
    error('X and/or ymat have NaN or Inf; please check your data');
end

%% Perform a check for rank deficiency
if rank(X) < nCols_X
    warning('Design matrix X is rank deficient');
end

%% Gather GRM
disp(RandomEffects);
if ismember({'A'}, RandomEffects)
    updateString = [char(datetime('now')), ': caller_FEMA_fit: gathering GRM'];
    disp(updateString);

    % Check if the user has passed in a GRMFile
    if isempty(GRMFile)
        error('RandomEffects included an A effect but GRMFile was not specified');
    else
        if ~exist(GRMFile, 'file')
            error(['Unable to find: ', GRMFile]);
        else
            [~, ~, ext] = fileparts(GRMFile);
            if isempty(ext)
                error('Please provide full path to GRMFile');
            else
                if strcmpi(ext, '.mat')
                    temp_GRM = load(GRMFile);

                    % Check to make sure necessary variables are present
                    toCheck_GRM = {'GRM', 'uqObservations'};
                    chk_GRM     = ismember(fieldnames(temp_GRM), toCheck_GRM);
                    if sum(chk_GRM) ~= length(toCheck_GRM)
                        missMsg = sprintf('%s, ', toCheck_GRM{~chk_GRM});
                        error(['The following required variables are missing from ', GRMFile, ': ', missMsg(1:end-2)]);
                    else
                        % Assign the right variables
                        GRM      = temp_GRM.GRM;
                        iid_list = temp_GRM.uqObservations;
                        clear temp_GRM

                        % Does GRM need re-ordering?
                        [a, ~, ~] = unique(iid, 'stable');
                        if length(a) ~= length(iid_list)
                            reorderGRM = true;
                        else
                            if sum(strcmpi(iid_list, a)) ~= length(a)
                                reorderGRM = true;
                            else
                                reorderGRM = false;
                            end
                        end

                        % Perform reordering, if required
                        if reorderGRM
                            [~, wch] = ismember(a, iid_list);
                            GRM      = GRM(wch, wch);
                        end
                    end
                else
                    if strcmpi(ext, '.dat')
                        basename = strrep(GRMFile, '.dat', '');
                        % Check if there is a corresponding csv file
                        if exist([basename, '.csv'], 'file')
                            iid_list = readtable(GRMFile, 'ReadVariableNames', false);
                        else
                            % Check if a .mat file exists
                            if exist([basename, '.mat'], 'file')
                                temp_GRM = load([basename, '.mat']);

                                % Make sure uqObservations exists
                                if isfield(temp_GRM, 'uqObservations')
                                    iid_list = temp_GRM.uqObservations;
                                    clear temp_GRM
                                else
                                    error(['Unable to find field uqObservations in: ', basename, '.mat']);
                                end
                            else
                                error(['Corresponding csv/mat file missing for: ', GRMFile]);
                            end
                        end

                        % Now read binary file
                        f   = fopen(GRMFile, 'r');
                        GRM = fread(f, [length(iid_list) length(iid_list)], 'double');
                        fclose(f);

                        % Make GRM full matrix
                        GRM = GRM + GRM.';

                        % Make diagonal equal to 1
                        tmp = size(GRM,1);
                        GRM(1:1+tmp:tmp*tmp) = 1;

                        % Convert to single
                        GRM = single(GRM);

                        % Does GRM need re-ordering?
                        [a, ~, ~] = unique(iid, 'stable');
                        if length(a) ~= length(iid_list)
                            reorderGRM = true;
                        else
                            if sum(strcmpi(iid_list, a)) ~= length(a)
                                reorderGRM = true;
                            else
                                reorderGRM = false;
                            end
                        end

                        % Perform reordering, if required
                        if reorderGRM
                            [~, wch] = ismember(a, iid_list);
                            GRM      = GRM(wch, wch);
                        end
                    else
                        error(['Unknown file format specified for GRMFile: ', ext]);
                    end
                end
            end
        end
    end
else
    GRM = [];
end

%% Gather mother ID
if ismember({'M'}, RandomEffects)
    updateString = [char(datetime('now')), ': caller_FEMA_fit: gathering mother ID'];
    disp(updateString);

    % Check if the user has passed in a mother ID fule
    if isempty(MotherID)
        error('RandomEffects included an M effect but MotherID file was not specified');
    else
        if ~exist(MotherID, 'file')
            error(['Unable to find: ', MotherID]);
        else
            [~, ~, ext] = fileparts(MotherID);
            if isempty(ext)
                error('Please provide full path to MotherID file');
            else
                if strcmpi(ext, '.csv')
                    temp_mid    = readtable(MotherID, 'Delimiter', '\t');
                    toCheck_mid = {'iid', 'MotherID'};
                    temp_names  = temp_mid.Properties.VariableNames;
                    chk_mid     = ismember(toCheck_mid, temp_names);
                    if sum(chk_mid) ~= length(toCheck_mid)
                        missMsg = sprintf('%s, ', toCheck_mid{~chk_mid});
                        error(['The following required variables are missing from ', MotherID, ': ', missMsg(1:end-2)]);
                    else
                        % Assign the right variables
                        iid_mid     = temp_mid.iid;
                        MotherID    = temp_mid.MotherID;
                        clear temp_mid

                        % Does MotherID need filtering?
                        [a, ~, ~] = unique(iid, 'stable');
                        if length(a) ~= length(iid_mid)
                            filterMID = true;
                        else
                            if sum(strcmpi(iid_mid, a)) ~= length(a)
                                filterMID = true;
                            else
                                filterMID = false;
                            end
                        end

                        % Perform filtering, if required
                        if filterMID
                            [~, wch] = ismember(a, iid_mid);
                            MotherID = MotherID(wch);
                        end
                    end
                else
                    error(['Unknown file format specified for MotherID: ', ext]);
                end
            end
        end
    end
end

%% Gather father ID
if ismember({'P'}, RandomEffects)
    updateString = [char(datetime('now')), ': caller_FEMA_fit: gathering father ID'];
    disp(updateString);

    % Check if the user has passed in a father ID fule
    if isempty(FatherID)
        error('RandomEffects included an P effect but FatherID file was not specified');
    else
        if ~exist(FatherID, 'file')
            error(['Unable to find: ', FatherID]);
        else
            [~, ~, ext] = fileparts(FatherID);
            if isempty(ext)
                error('Please provide full path to FatherID file');
            else
                if strcmpi(ext, '.csv')
                    temp_fid    = readtable(FatherID, 'Delimiter', '\t');
                    toCheck_fid = {'iid', 'FatherID'};
                    temp_names  = temp_fid.Properties.VariableNames;
                    chk_fid     = ismember(toCheck_fid, temp_names);
                    if sum(chk_fid) ~= length(toCheck_fid)
                        missMsg = sprintf('%s, ', toCheck_fid{~chk_fid});
                        error(['The following required variables are missing from ', FatherID, ': ', missMsg(1:end-2)]);
                    else
                        % Assign the right variables
                        iid_fid     = temp_fid.iid;
                        FatherID    = temp_fid.FatherID;
                        clear temp_fid

                        % Does FatherID need filtering?
                        [a, ~, ~] = unique(iid, 'stable');
                        if length(a) ~= length(iid_fid)
                            filterFID = true;
                        else
                            if sum(strcmpi(iid_fid, a)) ~= length(a)
                                filterFID = true;
                            else
                                filterFID = false;
                            end
                        end

                        % Perform filtering, if required
                        if filterFID
                            [~, wch] = ismember(a, iid_fid);
                            FatherID = FatherID(wch);
                        end
                    end
                else
                    error(['Unknown file format specified for FatherID: ', ext]);
                end
            end
        end
    end
end

%% Gather home ID
if ismember({'H'}, RandomEffects)
    updateString = [char(datetime('now')), ': caller_FEMA_fit: gathering home ID'];
    disp(updateString);

    % Check if the user has passed in a home ID file
    if isempty(HomeID)
        error('RandomEffects included an H effect but HomeID file was not specified');
    else
        if ~exist(HomeID, 'file')
            error(['Unable to find: ', HomeID]);
        else
            [~, ~, ext] = fileparts(HomeID);
            if isempty(ext)
                error('Please provide full path to HomeID file');
            else
                if strcmpi(ext, '.csv')
                    temp_hid    = readtable(HomeID, 'Delimiter', '\t');
                    toCheck_hid = {'iid', 'HomeID'};
                    temp_names  = temp_hid.Properties.VariableNames;
                    chk_hid     = ismember(toCheck_hid, temp_names);
                    if sum(chk_hid) ~= length(toCheck_hid)
                        missMsg = sprintf('%s, ', toCheck_hid{~chk_hid});
                        error(['The following required variables are missing from ', HomeID, ': ', missMsg(1:end-2)]);
                    else
                        % Assign the right variables
                        iid_hid   = temp_hid.iid;
                        HomeID    = temp_hid.HomeID;
                        clear temp_hid

                        % Does FatherID need filtering?
                        [a, ~, ~] = unique(iid, 'stable');
                        if length(a) ~= length(iid_hid)
                            filterHID = true;
                        else
                            if sum(strcmpi(iid_hid, a)) ~= length(a)
                                filterHID = true;
                            else
                                filterHID = false;
                            end
                        end

                        % Perform filtering, if required
                        if filterHID
                            [~, wch] = ismember(a, iid_hid);
                            HomeID   = HomeID(wch);
                        end
                    end
                else
                    error(['Unknown file format specified for HomeID: ', ext]);
                end
            end
        end
    end
end

%% Gather pregnancy ID
if ismember({'T'}, RandomEffects)
    updateString = [char(datetime('now')), ': caller_FEMA_fit: gathering pregnancy ID'];
    disp(updateString);

    % Check if the user has passed in a pregnancy ID file
    if isempty(PregID)
        error('RandomEffects included a T effect but PregID file was not specified');
    else
        if ~exist(PregID, 'file')
            error(['Unable to find: ', PregID]);
        else
            [~, ~, ext] = fileparts(PregID);
            if isempty(ext)
                error('Please provide full path to PregID file');
            else
                if strcmpi(ext, '.csv')
                    temp_pregid    = readtable(PregID, 'Delimiter', '\t');
                    toCheck_pregid = {'iid', 'PregID'};
                    temp_names     = temp_pregid.Properties.VariableNames;
                    chk_pregid     = ismember(toCheck_pregid, temp_names);
                    if sum(chk_pregid) ~= length(toCheck_pregid)
                        missMsg = sprintf('%s, ', toCheck_pregid{~chk_pregid});
                        error(['The following required variables are missing from ', PregID, ': ', missMsg(1:end-2)]);
                    else
                        % Assign the right variables
                        iid_pregID = temp_pregid.iid;
                        PregID     = temp_pregid.PregID;
                        clear temp_pregid

                        % Does FatherID need filtering?
                        [a, ~, ~] = unique(iid, 'stable');
                        if length(a) ~= length(iid_pregID)
                            filterPregID = true;
                        else
                            if sum(strcmpi(iid_pregID, a)) ~= length(a)
                                filterPregID = true;
                            else
                                filterPregID = false;
                            end
                        end

                        % Perform filtering, if required
                        if filterPregID
                            [~, wch] = ismember(a, iid_pregID);
                            PregID   = PregID(wch);
                        end
                    end
                else
                    error(['Unknown file format specified for PregID: ', ext]);
                end
            end
        end
    end
end

%% Read contrast file
if ~isempty(contrasts)
    updateString = [char(datetime('now')), ': caller_FEMA_fit: reading contrast file'];
    disp(updateString);
    [contrasts, hypValues] = FEMA_parse_contrastFile(contrasts, colnames);
else
    hypValues = NaN;
end

%% Run FEMA_fit
updateString = [char(datetime('now')), ': caller_FEMA_fit: calling FEMA_fit'];
disp(updateString);

[beta_hat,      beta_se,        zmat,        logpmat,                       ...
 sig2tvec,      sig2mat,        Hessmat,     logLikvec,                     ...
 beta_hat_perm, beta_se_perm,   zmat_perm,   sig2tvec_perm,                 ...
 sig2mat_perm,  logLikvec_perm, binvec_save, nvec_bins,                     ...
 tvec_bins,     FamilyStruct,   coeffCovar,  reusableVars] =                ...
 FEMA_fit(X, iid, eid, fid, agevec, ymat, niter, contrasts, nbins,          ...
          GRM, 'RandomEffects', RandomEffects, 'nperms', nperms,            ...
          'CovType', CovType, 'FixedEstType', FixedEstType,                 ...
          'RandomEstType', RandomEstType, 'GroupByFamType', GroupByFamType, ...
          'NonnegFlag', NonnegFlag, 'SingleOrDouble', SingleOrDouble,       ...
          'logLikflag', logLikflag, 'PermType', PermType,                   ...
          'returnReusable', returnReusable, 'doPar', doPar,                 ...
          'numWorkers', numWorkers, 'numThreads', numThreads,               ...
          'MotherID', MotherID, 'FatherID', FatherID, 'HomeID', HomeID,     ...
          'PregID', PregID);

%% Determine how to save the results
updateString = [char(datetime('now')), ': caller_FEMA_fit: saving results'];
disp(updateString);

% Get a sense for all variables in the workspace
tmpInfo = whos;

% First: save main (permuted and non-permuted) statistics
saveName = fullfile(dirOutput, [outPrefix, '_estimates.mat']);
toSave   = {'beta_hat', 'beta_se', 'zmat', 'logpmat', 'sig2tvec', 'sig2mat', ...
            'beta_hat_perm', 'beta_se_perm', 'zmat_perm', 'sig2tvec_perm',   ...
            'sig2mat_perm', 'logLikvec', 'logLikvec_perm', 'Hessmat',        ...
            'coeffCovar', 'binvec_save', 'nvec_bins', 'tvec_bins',           ...
            'reusableVars', 'contrasts', 'hypValues'};

if sum([tmpInfo(ismember({tmpInfo(:).name}', toSave)).bytes]) > 2^31
    save(saveName, 'beta_hat', 'beta_se', 'zmat', 'logpmat', 'sig2tvec', 'sig2mat', ...
            'beta_hat_perm', 'beta_se_perm', 'zmat_perm', 'sig2tvec_perm',          ...
            'sig2mat_perm', 'logLikvec', 'logLikvec_perm', 'Hessmat',               ...
            'coeffCovar', 'binvec_save', 'nvec_bins', 'tvec_bins',                  ...
            'reusableVars', 'contrasts', 'hypValues', '-v7.3');
else
    save(saveName, 'beta_hat', 'beta_se', 'zmat', 'logpmat', 'sig2tvec', 'sig2mat', ...
            'beta_hat_perm', 'beta_se_perm', 'zmat_perm', 'sig2tvec_perm',          ...
            'sig2mat_perm', 'logLikvec', 'logLikvec_perm', 'Hessmat',               ...
            'coeffCovar', 'binvec_save', 'nvec_bins', 'tvec_bins',                  ...
            'reusableVars', 'contrasts', 'hypValues');
end

% Second: save filtered design matrix and related variables
if saveDesignMatrix
    saveName = fullfile(dirOutput, [outPrefix, '_designMatrix.mat']);
    toSave   = {'X', 'fid', 'iid', 'eid', 'agevec', 'FamilyStruct', ...
                'MotherID', 'FatherID', 'HomeID', 'PregID'};
    if sum([tmpInfo(ismember({tmpInfo(:).name}', toSave)).bytes]) > 2^31
        save(saveName, 'X', 'fid', 'iid', 'eid', 'agevec', 'FamilyStruct', ...
                       'MotherID', 'FatherID', 'HomeID', 'PregID', '-v7.3');
    else
        save(saveName, 'X', 'fid', 'iid', 'eid', 'agevec', 'FamilyStruct', ...
                       'MotherID', 'FatherID', 'HomeID', 'PregID');
    end
end

% Third: save settings (this does not require v7.3 check)
saveName = fullfile(dirOutput, [outPrefix, '_settings.mat']);
save(saveName, 'niter', 'contrasts', 'nbins', 'RandomEffects', 'nperms',      ...
               'CovType', 'FixedEstType', 'RandomEstType', 'GroupByFamType',  ...
               'NonnegFlag', 'SingleOrDouble', 'logLikflag', 'PermType',      ...
               'returnReusable', 'doPar', 'numWorkers', 'numThreads', 'hypValues');