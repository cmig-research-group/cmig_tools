function [GRM, iid_list] = FEMA_read_GRM(GRMFile, fmt)
% Function to read genetic relatedness matrix (GRM)
%% Inputs:
% GRMFile:  character   full path to the GRM file (with extension)
% 
% fmt:      character   one of the format specifier (optional; see Notes):
%                           * mat
%                           * dat
%                           * gcta
% 
%% Outputs:
% GRM:      single      n x n genetic relatedness matrix
% 
% iid_list: cell        n x 1 list of subjects in the order of the GRM
% 
%% Notes:
% If fmt is not specified:
% - If extension is .mat, assume MATLAB file
% 
% - If extension is .dat, assume binary file and look for supporting
%   <basename>.mat file which has the list of subjects
% 
% - If extension is .bin, assume GCTA-style binary file and look for
%   supporting <basename>.id file which has the list of subjects
% 
% Format specification:
% ---------------------
% .mat file: a single file containing the following variables:
%   - GRM:            n x n full matrix of GRM values
%   - iid_list:       n x 1 cell array of IDs OR
%   - uqObservations: n x 1 cell array of IDs
%   - if both 'iid_list' and 'uqObservations' are present, 'iid_list' takes
%     precedence
% 
% .dat file: binary file with supporting <basename>.mat file:
%   - the GRM matrix is assumed to be lower triangle, saved with double
%     precision
%   - the supporting file contains a cell type variable either named
%     'iid_list' or 'uqObservations' which has the list of unique
%     participants in the order in which they are saved in the binary file
%   - if both 'iid_list' and 'uqObservations' are present, 'iid_list' takes
%     precedence
% 
% .bin file: GCTA-style file with supporting <basename>.id file:
%   - the GRM matrix is assumed to be upper triangle, saved with single
%     precision
%   - the supporting file will have a list of unique participants in the
%     order in which they are saved in the binary file
%   - the supporting .id file has no header
%   - the supporting .id file has two columns: family and participant id;
%     the output iid_list uses participant id

%% Check inputs
if ~exist('GRMFile', 'var') || isempty(GRMFile)
    error('Please input the full path to the GRM file');
end

if ~exist('fmt', 'var') || isempty(fmt)
    [~, ~, c] = fileparts(GRMFile);
    if isempty(c)
        error('Please input the full path to the GRM file');
    else
        if strcmpi(c, '.mat')
            fmt = 'mat';
        else
            if strcmpi(c, '.dat')
                fmt = 'dat';
            else
                if strcmpi(c, '.bin')
                    fmt = 'gcta';
                else
                    error('Unknown format and/or missing extension in the GRM file');
                end
            end
        end
    end
else
    fmt = lower(fmt);
    if ~ismember(fmt, {'mat', 'dat', 'gcta'})
        error('fmt should be one of: mat, dat, or gcta');
    end
end

%% Read GRM
switch fmt
    case 'mat'
        % First, make sure that the variables exist inside the mat file
        variableInfo = who('-file', GRMFile);
        if ismember('GRM', variableInfo)
            if ismember('iid_list', variableInfo)
                results  = load(GRMFile, 'GRM', 'iid_list');
                iid_list = results.iid_list;
            else
                if ismember('uqObservations', variableInfo)
                    results  = load(GRMFile, 'GRM', 'uqObservations');
                    iid_list = results.uqObservations;
                else
                    error(['iid_list and uqObservations are missing in GRM mat file; ' ...
                           'either should be present in the GRM mat file']);
                end
            end
            GRM = single(results.GRM);
        else
            error('GRM variable is missing in the GRM mat file');
        end

    case 'dat'
        % First, make sure that the corresponding .mat file exists
        idFile = strrep(GRMFile, '.dat', '.mat');
        if ~exist(idFile, 'file')
            error('.dat file detected but corresponding .mat file not found');
        else
            % Make sure that iid_list exists
            variableInfo = who('-file', idFile);
            if ismember('iid_list', variableInfo)
                % Load the mat file and get iid_list
                tmp      = load(idFile, 'iid_list');
                iid_list = tmp.iid_list;
            else
                if ismember('uqObservations', variableInfo)
                    tmp      = load(idFile, 'uqObservations');
                    iid_list = tmp.uqObservations;
                else
                    error(['iid_list and uqObservations are missing in the supporting mat file; ' ...
                           'either should be present in the supporting mat file']);
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
        end

    case 'gcta'
        % First, make sure that the corresponding .id file exists
        idFile = strrep(GRMFile, '.bin', '.id');
        if ~exist(idFile, 'file')
            error('.bin file detected but corresponding .id file not found');
        else
            % The .id file has two columns: family ID and participant id;
            % use the second column to get only participant id
            iid_list = readtable(idFile, 'FileType', 'text', 'ReadVariableNames', false);
            iid_list = iid_list{:,2};

            % Determine number of subjects to read
            nSubjs    = height(iid_list);
            sz_toRead = height(iid_list)*(height(iid_list)+1)/2;

            % Get upper triangle index
            idx_triu = find(triu(true(nSubjs, nSubjs)));

            % Read binary file
            fid = fopen(GRMFile, 'r');
            dat = fread(fid, sz_toRead, 'single');
            fclose(fid);

            % Make GRM full
            GRM           = zeros(nSubjs, nSubjs);
            GRM(idx_triu) = dat; %#ok<FNDSB>
            GRM           = triu(GRM) + triu(GRM, 1).';
        end
end