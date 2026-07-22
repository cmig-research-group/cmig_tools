function [ymat, iid_concat, eid_concat, ivec_mask, mask, GRM, preg, address, info] = FEMA_process_data(fstem_imaging,dirname_imaging,datatype,varargin)
%
% ABCD specific function to load and process imaging data 
%
% INPUTS
%	 fstem_imaging <char>		:	name of vertex/voxel-mapped phenotype (e.g., 'thickness-sm16', 'FA') or 
%                                   table name for tabulated ROI data, can be emtpy for datatype 'external'
%	 dirname_imaging <char>		:	path to imaging data directory, if datatype is external than need the 
%                                   full path to the file. 
%	 datatype <char>			:	'voxel','vertex','corrmat', 'roi', 'external'
%									Other then 'external' all code is written to expect ABCD data  
%                                   
% Optional input arguments:
%	 ico <num>					:	ico-number for vertexwise analyses (0-based, default 5)
%	 ranknorm_wholeSample <boolean>			:	rank normalise imaging data (default 0)
%	 standardize_wholeSample <boolean>		:	variance normalise imaging data (default 0)
%	 GRM_file <char>			:	path to genetic relatedness data (GRM) - default [] - only required if A random effect specified
%	 preg_file <char>			:	path to pregnancy data - default [] - only required if T random effect specified
%	 address_file <char>		:	path to address data - default [] - only required if H random effect specified
%   fname_qc <char>			    :	path to QC file - default [] - only supplied if want to filter on some QC measure
%                                   must include participant_id and session_id columns 
%                                   file format can be csv, tsv, parquet 
%   qc_var <char>			    :	name of ONLY ONE variable in QC - default [] - only required if QC random effect specified
% 
% OUTPUTS
%	 ymat						:	matrix of imaging data (n x v)
%	 iid_concat					:	participant_id
%	 eid_concat					:	session_id
%	 ivec_mask					:	vector mask for voxelwise data (155179x1) --> not yet available for vertexwise
%	 mask						:	volume mask mask for voxelwise data (100x100x130) --> not yet available for vertexwise
%	 ymat_names			        :	imaging column labels for external data inputs
%	 GRM						:	intersected genetic relatedness matrix
%

%
% This software is Copyright (c) 2021 The Regents of the University of California. All Rights Reserved.
% See LICENSE.

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% TO DOs: 
% - how can we get study and release for roi data? 

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% allow empty inputs to 
allowEmpty = @(f) @(x) isempty(x) || f(x);
p = inputParser;
addRequired(p, 'fstem_imaging', @(x) (iscell(x) || ischar(x)) && ~isempty(x));
addRequired(p, 'dirname_imaging', @(x) ischar(x) && ~isempty(x));
addRequired(p, 'datatype', @(x) ischar(x) && ~isempty(x));
addParameter(p, 'iid', [], allowEmpty(@(x) iscell(x) || ischar(x)));
addParameter(p, 'eid', [], allowEmpty(@(x) iscell(x) || ischar(x)));
addParameter(p, 'fname_qc', [], allowEmpty(@(x) ischar(x)));
addParameter(p, 'qc_var', [], allowEmpty(@(x) ischar(x)));
addParameter(p, 'ico', 5);
addParameter(p, 'GRM_file', [], allowEmpty(@(x) ischar(x)));
addParameter(p, 'preg_file', [], allowEmpty(@(x) ischar(x)));
addParameter(p, 'address_file', [], allowEmpty(@(x) ischar(x)));
addParameter(p, 'corrvec_thresh', 0.8);
addParameter(p, 'nframes_min', 375);
addParameter(p, 'study', [], allowEmpty(@(x) ischar(x)));
addParameter(p, 'release', [], allowEmpty(@(x) ischar(x)));
addParameter(p, 'wholeSampleTransform', false, @(x) isscalar(x) && (islogical(x) || isnumeric(x)));
addParameter(p, 'transformY', 'none', @(x) ischar(x) && ...
                                           ismember(x, {'none' 'center' 'centre' 'demean' ...
                                           'std' 'standardize' 'normalize' ...
                                           'logn' 'log10' 'inverseranknorm' ...
                                           'ranknorm' 'int'}));
addParameter(p, 'nonDEAP', false, @(x) isscalar(x) && (islogical(x) || isnumeric(x)));

parse(p, fstem_imaging, dirname_imaging, datatype, varargin{:})
iid = p.Results.iid;
eid = p.Results.eid;
fname_qc = p.Results.fname_qc;
qc_var = p.Results.qc_var;
ico = str2num_amd(p.Results.ico);
% icnum = ico + 1;
fname_GRM = p.Results.GRM_file;
fname_preg = p.Results.preg_file;
fname_address = p.Results.address_file;
corrvec_thresh = p.Results.corrvec_thresh;
nframes_min = p.Results.nframes_min;
study = p.Results.study;
release = p.Results.release;
wholeSampleTransform = p.Results.wholeSampleTransform;
transformY = p.Results.transformY;
nonDEAP = p.Results.nonDEAP;

logging(FEMA_info)
% log overall timing 
tOverall = tic; 

% save settings 
info.settings = rmfield(p.Results, {'iid', 'eid'});
info.settings.iid_filter = ~isempty(iid);
info.settings.eid_filter = ~isempty(eid);
% basic info 
info.FEMA_version       = FEMA_info;
info.provenance = 'FEMA_process_data'; 

% check if QC file and variable are both provided
if ~isempty(fname_qc)
    if isempty(qc_var)
        error('QC variable must be specified if QC file is provided');
    end
end

ymat_names = [];
missingness = [];

%%%%% load data %%%%%
switch datatype
    case {'voxel', 'voxelwise'}
        %Load voxelwise imaging data
		logging('Reading voxelwise %s imaging data',fstem_imaging);
		tLoadData = tic; 
		measname = fstem_imaging;
		dirname_volmats = dirname_imaging;
	
		% 6.0 naming convention
		fname_volinfo = sprintf('%s/vol_info.mat',dirname_volmats);
		if ~exist(fname_volinfo,'file') % if older release version 
			fname_volinfo = sprintf('%s/volinfo.mat',dirname_volmats);
		end
		tmp_volinfo = load(fname_volinfo);

		fname_volmat = sprintf('%s/%s.mat',dirname_volmats,lower(measname));
		if ~exist(fname_volmat,'file') % if older release version 
			fname_volmat = sprintf('%s/volmat_%s.mat',dirname_volmats,measname);
		end	
		ymat = getfield(load(fname_volmat),'volmat');
        ymat_names = fstem_imaging;
		mask = tmp_volinfo.vol_mask_sub;
		ivec_mask = find(mask>0.5);

        %% old code 
            % dirlist = tmp_volinfo.dirlist;
            %subjidvec = cell(size(dirlist)); 
		    %sitevec = cell(size(dirlist)); 
		    %datevec = cell(size(dirlist)); 
		    %visitidvec = cell(size(dirlist));
            %if any(contains(dirlist, 'DTIREG_'))
            %    for diri = 1:length(dirlist)
            %    	tmp = regexp(dirlist{diri}, '^DTIREG_(?<site>\w+)_(?<SubjID>\w+)_(?<event>\w+)_(?<date>\d+).', 'names');
            %    	subjidvec{diri} = tmp.SubjID;
            %    	sitevec{diri} = tmp.site;
            %    	eventvec{diri} = tmp.event;
            %    	datevec{diri} = tmp.date;
            %    	visitidvec{diri} = sprintf('%s_%s_%s',tmp.site,tmp.SubjID,tmp.event);
            %    end
            %    iid_concat = colvec(strcat('NDAR_',subjidvec)); % Not sure why imaging subject IDs are missing the NDAR_ part
            %    eid_concat = eventvec';
            %else 
            %    iid_concat = tmp_volinfo.participant_id;
            %    eid_concat = tmp_volinfo.session_id;
            %end 
        iid_concat = tmp_volinfo.participant_id;
        eid_concat = tmp_volinfo.session_id;
        idevent = strcat(iid_concat(:),'_',eid_concat(:));
        corrmat_concat = tmp_volinfo.corrmat;
        
        info.timing.tLoadData = toc(tLoadData);
        
	case {'vertex', 'vertexwise'}
		% Read in vertexwise imaging data
		logging('Reading vertexwise %s imaging data',fstem_imaging);
		tLoadData = tic;
        % No longer reading icsurfs; using utility function numIcoVertices
		%	SurfView_loadsurfs; % Shouldn't be neccessary, if data saved out pre-truncated
		% load('SurfView_surfs.mat'); %this matfile is included in the executable when compiling using -a
        % try
        %     load('SurfView_surfs.mat', 'icsurfs');
        % catch
        %     tmpFile = fullfile(fileparts(fileparts(which('FEMA_process_data'))), 'showSurf', 'SurfView_surfs.mat');
        %     if exist(tmpFile, 'file')
        %         load(tmpFile, 'icsurfs');
        %     else
        %         error(['Could not find SurfView_surfs.mat on the path or in: ', tmpFile]);
        %     end
        % end

		hemistrings = {'lh','rh'};

		measmat = [];
		for hemii = 1:2
			hemi = hemistrings{hemii};
			fname = fullfile(dirname_imaging,[fstem_imaging '_' hemi '.mat']);
			if ~exist(fname,'file') % if older release version
				fname = sprintf('%s/%s-%s.mat',dirname_imaging,fstem_imaging,hemi);  
			end
			% specified icnum
			logging('Reading vertexwise imaging data from %s', fname);
			tmp = load(fname);
            numVertices = numIcoVertices(ico); % Not icnum; icnum is for indexing icsurfs
            measmat = cat(2, measmat,tmp.measmat(:,1:numVertices,1));
			% measmat = cat(2, measmat,tmp.measmat(:,1:size(icsurfs{icnum}.vertices,1)));
		end

		% deal with masking (replace NaNs w. zeros for all subjects, but not for subjects/visits with all missing data)
		ivec_nan = find(mean(isfinite(measmat),2)==0);
		measmat(~isfinite(measmat)) = 0;
		measmat(ivec_nan,:) = NaN;

		ivec_mask=find(nanstd(measmat)~=0);
		if length(ivec_mask)~=size(measmat,2)
			mask=zeros(1,size(measmat,2));
			measmat=measmat(:,ivec_mask);
			mask(ivec_mask)=1;
		else
			mask=[];
		end
		ymat = measmat;
        ymat_names = fstem_imaging;

        % load vol_info 
		fname_volinfo = sprintf('%s/vol_info.mat',dirname_imaging);
		tmp_volinfo = load(fname_volinfo); 

        %% old code 
            %dirlist = tmp_volinfo.dirlist;
		    %subjidvec = cell(size(dirlist)); 
		    %sitevec = cell(size(dirlist)); 
		    %datevec = cell(size(dirlist)); 
		    %eventvec = cell(size(dirlist)); 
		    %timevec = cell(size(dirlist)); 
		    %visitidvec = cell(size(dirlist));
            %
		    %for diri = 1:length(dirlist)
		    %	tmp = regexp(dirlist{diri}, '^[^_]*_(?<site>\w+)_(?<SubjID>\w+)_(?<event>[^_]*)_(?<date>\d+).(?<time>[^_]+)_', 'names'); % AMD -- 3.0 fails for task fMRI data
		    %	sitevec{diri} = tmp.site;
		    %	subjidvec{diri} = tmp.SubjID;
		    %	datevec{diri} = tmp.date;
		    %	if isfield(tmp,'event')
		    %		eventvec{diri} =	tmp.event;
		    %	else
		    %		eventvec{diri} = '';
		    %	end
		    %	visitidvec{diri} = sprintf('%s_%s_%s',tmp.site,tmp.SubjID,tmp.event);
		    %end
            %
            %iid_concat = colvec(strcat('NDAR_',subjidvec));
		    %eid_concat = eventvec;
		    %idevent = strcat(iid_concat(:),'_',eid_concat(:));

        % use new participant_id 
        iid_concat = tmp_volinfo.participant_id;
		eid_concat = tmp_volinfo.session_id;
		idevent = strcat(iid_concat(:),'_',eid_concat(:));
		info.timing.tLoadData = toc(tLoadData);

	case 'corrmat'
		logging('Processing corrmat %s data', fstem_imaging);
		tLoadData = tic;
        % vol info
        fname_volinfo = sprintf('%s/vol_info.mat',dirname_imaging);
        tmp_volinfo = load(fname_volinfo);
		if isstruct(dirname_imaging)
			tmp_corrmat = dirname_imaging;
		else
			fname_corrmat = sprintf('%s/%s.mat',dirname_imaging,fstem_imaging);
			tmp_corrmat = load(fname_corrmat);
		end
		% if isfield(tmp,'nsumvec')	%	Release 4.0
		% 	ivec_tmp = find(tmp.nsumvec>=nframes_min); 
		% 	measmat = tmp.measmat(ivec_tmp,:);
		% elseif isfield(tmp,'dirlist') %	Release 5.1
		% 	ivec_tmp = find(tmp.ntpointvec>=nframes_min); 
		% 	measmat = tmp.corrmat(ivec_tmp,:);
		% else 
		% 	ivec_tmp = find(tmp.ntpointvec>=nframes_min); 
		% 	measmat = tmp.corrmat(ivec_tmp,:);
		% 	tmp.dirlist = tmp_volinfo.dirlist;
		% end
		% dirlist = tmp.dirlist(ivec_tmp);
		% dims = size(measmat); 
        % measmat = reshape(measmat,[dims(1) prod(dims(2:end))]);
		% ymat = measmat;
        ymat = tmp_corrmat.corrmat;
        ymat_names = tmp_corrmat.roinames;
        % ymat_names = [tmp_corrmat.roinames(tmp_corrmat.roi1vec)', tmp_corrmat.roinames(tmp_corrmat.roi2vec)'];

        %% old code 
		    %subjidvec = cell(size(dirlist)); 
		    %sitevec = cell(size(dirlist));`
		    %datevec = cell(size(dirlist)); 
		    %visitidvec = cell(size(dirlist));
		    %for diri = 1:length(dirlist)
		    %	tmp = regexp(dirlist{diri}, '^BOLDPROC_(?<site>\w+)_(?<SubjID>\w+)_(?<event>[^_]*)_(?<date>[^_.]+).(?<time>[^_]+)_', 'names');
		    %	subjidvec{diri} = tmp.SubjID;
		    %	sitevec{diri} = tmp.site;
		    %	eventvec{diri} = tmp.event;
		    %	datevec{diri} = tmp.date;
		    %	visitidvec{diri} = sprintf('%s_%s_%s',tmp.site,tmp.SubjID,tmp.event);
		    %end
            %
		    %iid_concat = colvec(strcat('NDAR_',subjidvec)); % Not sure why imaging subject IDs are missing the NDAR_ part
		    %eid_concat = eventvec;
		
        iid_concat = tmp_volinfo.participant_id;
        eid_concat = tmp_volinfo.session_id;
        idevent = strcat(iid_concat(:),'_',eid_concat(:));

		ivec_mask=[];
		mask=[];

        info.timing.tLoadData = toc(tLoadData);

    case 'roi'
        logging('Reading ROI tabulated data from %s', fstem_imaging);
        tLoadData = tic;
        fname_roi = fullfile(dirname_imaging, fstem_imaging);
        [dirname_roi, fstem_roi, ext_roi] = fileparts(fname_roi);
        if ~exist(fname_roi, 'file')
            if isempty(ext_roi)
                ext_roi = '.parquet';
                fname_roi = [fname_roi, '.parquet'];
                if ~exist(fname_roi, 'file')
                    error('ROI file %s does not exist.', fname_roi);
                end
            else
                error('ROI file %s does not exist.', fname_roi);
            end
        end
        switch ext_roi 
            case '.parquet'
                roidat = parquetread(fname_roi);
                roidat.session_id = cellstr(roidat.session_id);
            case {'.csv' '.tsv' '.txt'}
                opts = detectImportOptions(fname_roi, 'FileType', 'text');
                roidat = readtable(fname_roi, opts);
            otherwise
                error('Unsupported file format for ROI data: %s. Must be .parquet, .csv, .tsv, or .txt', ext_roi);
        end 
        % Check which ID + event columns are present for the different releases
	    %if any(strcmp(roidat.Properties.VariableNames, 'src_subject_id'))
        %	iid_concat = roidat.src_subject_id;
        %	eid_concat = roidat.eventname;
	    %elseif any(strcmp(roidat.Properties.VariableNames, 'participant_id'))
        %	iid_concat = roidat.participant_id;
        %	eid_concat = roidat.session_id;
	    %else
        %	error('No recognized ID/event columns found in the external file.');
	    %end
        iid_concat = roidat.participant_id;
        eid_concat = roidat.session_id;
        roidat = removevars(roidat, {'participant_id', 'session_id'});
	    ymat = table2array(roidat);
	    ymat_names = roidat.Properties.VariableNames;
	    idevent = strcat(iid_concat(:),'_',eid_concat(:));

	    ivec_mask=[];
	    mask=[];

        info.timing.tLoadData = toc(tLoadData);
    case 'external'
	    logging('Reading tabulated imaging data from %s', dirname_imaging);
        tLoadData = tic;
        if nonDEAP
            fname_roidat = fullfile(dirname_imaging, fstem_imaging);
            [~, ~, ext_imaging] = fileparts(fname_roidat);
            switch ext_imaging
                case '.parquet'
                    roidat = parquetread(fname_roidat);
                case {'.csv' '.tsv' '.txt'}
                    opts = detectImportOptions(fname_roidat, 'FileType', 'text');
                    roidat = readtable(fname_roidat, opts);
                otherwise
                    error('Unsupported file format for tabulated data: %s. Must be .parquet, .csv, .tsv, or .txt', ext_imaging);
            end
        else 
            [~, ~, ext_imaging] = fileparts(dirname_imaging);
            if ~isempty(fstem_imaging)
                switch ext_imaging
                    case '.parquet'
                        opts = parquetinfo(dirname_imaging);
                        varInc = [rowvec(fstem_imaging) {'participant_id' 'session_id'}];
                        varPresent = ismember(varInc, opts.VariableNames);
                        if any(~varPresent)
                            missingVars = varInc(~varPresent);
                            error(['Could not find variables in table ', fname_qc, ': ', strjoin(missingVars, ', ')]);
                        else
                            try 
                                roidat = parquetread(dirname_imaging, 'SelectedVariableNames', varInc);
                                roidat.session_id = cellstr(roidat.session_id);
                            catch ME
                                sprintf('Cannot load tabulated data from: %s\n%s', dirname_imaging, rethrow(ME));
                                % error('Cannot load tabulated data from %s', dirname_imaging);
                            end
                        end 
                    case {'.csv' '.tsv' '.txt'}
                        opts = detectImportOptions(dirname_imaging, 'FileType', 'text');
                        varInc = [rowvec(fstem_imaging) {'participant_id' 'session_id'}];
                        varPresent = ismember(varInc, opts.VariableNames);
                        if any(~varPresent)
                            missingVars = varInc(~varPresent);
                            error(['Could not find variables in file ', dirname_imaging, ': ', strjoin(missingVars, ', ')]);
                        else
                            opts.SelectedVariableNames = varInc;  
                            try
                                roidat = readtable(dirname_imaging, opts);
                            catch ME
                                sprintf('Cannot load tabulated data from: %s\n%s', dirname_imaging, rethrow(ME));
                                % error('Cannot load tabulated data from %s', dirname_imaging);
                            end 
                        end
                    otherwise
                        error('Unsupported file format for tabulated data: %s. Must be .parquet, .csv, .tsv, or .txt', ext_imaging); 
                end
            else 
                error('No dependent variable specified');
            end 
        end 
	    % Check which ID + event columns are present for the different releases
	    %if any(strcmp(roidat.Properties.VariableNames, 'src_subject_id'))
        %	iid_concat = roidat.src_subject_id;
        %	eid_concat = roidat.eventname;
	    %elseif any(strcmp(roidat.Properties.VariableNames, 'participant_id'))
        %	iid_concat = roidat.participant_id;
        %	eid_concat = roidat.session_id;
	    %else
        %	error('No recognized ID/event columns found in the external file.');
	    %end
        iid_concat = roidat.participant_id;
        eid_concat = roidat.session_id;

        % remove non-numeric columns from roidat - this should only needed for non-DEAP data
        dataVarNames = setdiff(roidat.Properties.VariableNames, {'participant_id', 'session_id'}, 'stable');
        isNumericData = cellfun(@(x) isnumeric(roidat.(x)), dataVarNames);
        toDropNonNumeric = dataVarNames(~isNumericData);
        if ~isempty(toDropNonNumeric)
            fprintf('%d column(s) [%s] are not numeric. Dropping from analysis.\n', ...
                numel(toDropNonNumeric), strjoin(cellstr(toDropNonNumeric(:)), ', '));
            roidat = removevars(roidat, toDropNonNumeric);
        end

        %define dependent variable, ymat
        ymat = table2array(removevars(roidat, {'participant_id', 'session_id'}));
	    ymat_names = removevars(roidat, {'participant_id', 'session_id'}).Properties.VariableNames;
	    idevent = strcat(iid_concat(:),'_',eid_concat(:));
	    ivec_mask=[];
	    mask=[];

        info.timing.tLoadData = toc(tLoadData);
end

%%%%% back up ymat and lists %%%%%
% ymat_bak = ymat;
% iid_concat_bak = iid_concat;
% eid_concat_bak = eid_concat;
% idevent_bak = idevent;

%%%%% missingness %%%%%
% number of nans in ymat  
tMissingNaN = tic;
defvec_nan = isfinite(sum(ymat, 2)); 
missingness.n_nan = sum(~defvec_nan);
missingness.id_nan = idevent(~defvec_nan);
if strcmp(datatype, 'roi') || strcmp(datatype, 'external')
    % missingness in ymat
    logging('Removing observations based on missingness');
    % Record the number of observations prior to removing
    missingness.n_obs_ymat_preCompCases = height(ymat);

    tmp = ismissing(ymat);

    % Make a report of the missingness pattern
    missingness.variableWise = cell2table([ymat_names',   ...
                                                num2cell(sum(tmp))'], 'VariableNames', ...
                                                {'Variable', 'n_missing'});
end 
info.timing.tMissingNaN = toc(tMissingNaN);
% Generate an error if no data is left
if sum(defvec_nan) == 0
    error('No observations left in the data table; please check variables for missingness');
end

%  how many don't pass threshold on correlation 
if strcmp(datatype, 'voxel')
    tCorrVec = tic;
	corrvec = mean(corrmat_concat,2);
    thresh = corrvec_thresh;
	defvec_corrvec = corrvec>=thresh;
    missingness.n_fail_corrvec = sum(~defvec_corrvec);
    missingness.id_fail_corrvec = idevent(~defvec_corrvec);
    info.timing.tCorrVec = toc(tCorrVec);
else
    defvec_corrvec = true(size(ymat, 1), 1);
end 

% Check for number of time points (corrmat)
if strcmpi(datatype, 'corrmat')
    tnFramesVec = tic;
    defvec_nframes = tmp_corrmat.ntpointvec >= nframes_min;
    missingness.n_fail_nframes  = sum(~defvec_nframes);
    missingness.id_fail_nframes = idevent(defvec_nframes);
    info.timing.tnFramesVec     = toc(tnFramesVec);
else
    defvec_nframes = true(size(ymat, 1), 1);
end
    
% filter on qc if qc file is given 
if ~isempty(fname_qc)
    tQC = tic;
    % check file type of fname_qc
    [dirname_qc, fstem_fname_qc, ext_qc] = fileparts(fname_qc);
    switch ext_qc
        case '.parquet'
            opts = parquetinfo(fname_qc);
            varInc = {'participant_id', 'session_id', qc_var};
            varPresent = ismember(varInc, opts.VariableNames);
            if any(~varPresent)
                missingVars = varInc(~varPresent);
                error(['Could not find variables in table ', fname_qc, ': ', strjoin(missingVars, ', ')]);
            else
                try 
                    qc_tbl = parquetread(fname_qc, 'SelectedVariableNames', varInc);
                catch ME
                    sprintf('Cannot load QC file %s with QC variable %s\n%s', fname_qc, qc_var, rethrow(ME));
                    % error('Cannot load QC file %s with QC variable %s', fname_qc, qc_var);
                    % sprintf('Error: %s', rethrow(ME));
                end
            end 
        case {'.csv' '.tsv' '.txt'}
            opts = detectImportOptions(fname_qc, 'FileType', 'text');
            varInc = {'participant_id', 'session_id', qc_var};
            varPresent = ismember(varInc, opts.VariableNames);
            if any(~varPresent)
                missingVars = varInc(~varPresent);
                error(['Could not find variables in table ', fname_qc, ': ', strjoin(missingVars, ', ')]);
            else
                opts.SelectedVariableNames = {'participant_id', 'session_id', qc_var};  % choose specific columns
                try
                    qc_tbl = readtable(fname_qc, opts);
                catch ME
                    sprintf('Cannot load QC file %s with QC variable %s\n%s', fname_qc, qc_var, rethrow(ME));
                    % sprintf('Error: %s', rethrow(ME));
                end 
            end 
    end
    passQC_idevent = strcat(qc_tbl.participant_id, '_', string(qc_tbl.session_id));
    passQC_idevent = passQC_idevent((qc_tbl.(qc_var) == '1'));
    defvec_qc = ismember(idevent, passQC_idevent);
    missingness.n_fail_qc = sum(~defvec_qc);    
    missingness.id_fail_qc = idevent(~defvec_qc);
    info.timing.tQC = toc(tQC);
else 
    defvec_qc = true(size(ymat, 1), 1);
end 

% Apply transform to whole sample
if wholeSampleTransform
    logging('Transforming whole sample Y');
    tWholeSampleTransform = tic;
    [ymat settingsTransform] = doTransformation(ymat, transformY);
    info.timing.tWholeSampleTransform = toc(tWholeSampleTransform);
end

% filter on iid and eid if given 
% iid 
if ~isempty(iid) 
    tFilterIID = tic;
    if ischar(iid)  
        if isfile(iid) % check if it's a file 
            [~, ~, ext_iid] = fileparts(iid);
            switch ext_iid % load file 
                case '.parquet'
                    iid = parquetread(iid);                        
                case {'.csv', '.tsv', '.txt'}
                    iid = readtable(iid, 'FileType', 'text', 'ReadVariableNames', false);
            end 
            iid = table2cell(iid); % convert to cell 
        else
            iid = cellstr(iid); % if not a file, convert to cell array
        end
    end
    defvec_iid = ismember(iid_concat, iid); 
    missingness.n_rm_iid = sum(~defvec_iid);
    missingness.id_rm_iid = idevent(~defvec_iid);
    info.timing.tFilterIID = toc(tFilterIID);
else 
    defvec_iid = true(size(iid_concat));
end
% eid 
if ~isempty(eid)
    tFilterEID = tic;
    if ischar(eid)  
        if isfile(eid) % check if it's a file 
            [~, ~, ext_eid] = fileparts(eid);
            switch ext_eid % load file 
                case '.parquet'
                    eid = parquetread(eid);
                    
                case {'.csv', '.tsv', '.txt'}
                    eid = readtable(eid, 'FileType', 'text', 'ReadVariableNames', false);
            end 
            eid = table2cell(eid); % convert to cell 
        else
            eid = cellstr(eid); % if not a file, convert to cell array
        end
    end
    defvec_eid = ismember(eid_concat, eid); 
    missingness.n_rm_eid = sum(~defvec_eid);
    missingness.id_rm_eid = idevent(~defvec_eid);
    info.timing.tFilterEID = toc(tFilterEID);
else 
    defvec_eid = true(size(eid_concat));
end
    
%%%%% load GRM %%%%%
if ~isempty(fname_GRM)
	logging('Reading GRM');
    tGRM = tic;
    [GRM.GRM, GRM.iid_list] = FEMA_read_GRM(fname_GRM);
	% GRM = load(fname_GRM);
    iid_grm = GRM.iid_list; 
    [keep, IA, IB] = intersect(iid_concat, iid_grm, 'stable');
    defvec_grm = ismember(iid_concat, keep);
    GRM.GRM = GRM.GRM(IB,IB);
    GRM.iid_list = iid_grm(IB);
    % missingness
    missingness.n_fail_grm = sum(~defvec_grm);
    missingness.id_fail_grm = idevent(~defvec_grm);
    info.timing.tGRM = toc(tGRM);
else 
	defvec_grm = true(size(iid_concat));
    GRM.GRM = [];
    GRM.iid_list = {};
end

%%%%% load pregnancy data %%%%%
if ~isempty(fname_preg)
    logging('Reading pregnancy data');
    tPreg = tic;
    fid = fopen(fname_preg);
    varNames = strsplit(fgetl(fid), {' ' '"'});
    fclose(fid);
    opts=detectImportOptions(fname_preg);
    opts.SelectedVariableNames = [2:5];
    preg = readtable(fname_preg, opts);
    preg = renamevars(preg,1:width(preg),varNames(2:5));
    info.timing.tPreg = toc(tPreg);
else
	defvec_preg = true(size(iid_concat));
    preg = []; 
end

%%%%% load home address data %%%%%
if ~isempty(fname_address)
	logging('Reading home address data');
    tAddress = tic;
	address = readtable(fname_address);
    info.timing.tAddress = toc(tAddress);
else
	defvec_address = true(size(iid_concat));
    address = [];
end

%%%%% filter on all defvecs %%%%%  
tFilterAll = tic;
defvec = defvec_nan & defvec_corrvec & defvec_qc & defvec_nframes & ...
         defvec_iid & defvec_eid & defvec_grm & ...
         defvec_preg & defvec_address;

if sum(defvec) == 0
    error('No observations left in the data table; please check variables for missingness');
end

iid_concat = cellstr(iid_concat(defvec));
eid_concat = cellstr(eid_concat(defvec));
idevent = cellstr(idevent(defvec));
ymat = ymat(defvec, :);
info.timing.tFilterAll = toc(tFilterAll);

% save missingness and final sample size
info.missingness = missingness;
info.n_obs_ymat = size(ymat, 1);
info.n_cols_ymat = size(ymat, 2);
info.ymat_names = ymat_names;

% save timing info 
info.timing.tOverall = toc(tOverall);

logging('Final sample for analysis: %d\n', size(ymat, 1));
logging('***End***');
logging('Elapsed time: %s seconds', num2str(info.timing.tOverall));

end 



