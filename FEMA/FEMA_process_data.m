function [ymat, iid_concat, eid_concat, ivec_mask, mask, colnames_imaging, GRM, preg, address, missingness] = FEMA_process_data(fstem_imaging,dirname_imaging,datatype,varargin)
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
%	 standardize_wholeSample <boolean>			:	variance normalise imaging data (default 0)
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
%	 iid_concat					:	src_subject_id
%	 eid_concat					:	eventname
%	 ivec_mask					:	vector mask for voxelwise data (155179x1) --> not yet available for vertexwise
%	 mask						:	volume mask mask for voxelwise data (100x100x130) --> not yet available for vertexwise
%	 colnames_imaging			:	imaging column labels for external data inputs
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
addRequired(p, 'fstem_imaging', @(x) ischar(x) && ~isempty(x));
addRequired(p, 'dirname_imaging', @(x) ischar(x) && ~isempty(x));
addRequired(p, 'datatype', @(x) ischar(x) && ~isempty(x));
addParameter(p, 'iid', [], allowEmpty(@(x) iscell(x) || ischar(x)));
addParameter(p, 'eid', [], allowEmpty(@(x) iscell(x) || ischar(x)));
addParameter(p, 'fname_qc', [], allowEmpty(@(x) ischar(x)));
addParameter(p, 'qc_var', [], allowEmpty(@(x) ischar(x)));
addParameter(p, 'ranknorm_wholeSample', false);
addParameter(p, 'standardize_wholeSample', false);
addParameter(p, 'wholeSampleTransform', false);
addParameter(p, 'ico', 5);
addParameter(p, 'GRM_file', [], allowEmpty(@(x) ischar(x)));
addParameter(p, 'preg_file', [], allowEmpty(@(x) ischar(x)));
addParameter(p, 'address_file', [], allowEmpty(@(x) ischar(x)));
addParameter(p, 'corrvec_thresh', 0.8);
addParameter(p, 'nframes_min', 375);
addParameter(p, 'study', [], allowEmpty(@(x) ischar(x)));
addParameter(p, 'release', [], allowEmpty(@(x) ischar(x)));

parse(p, fstem_imaging, dirname_imaging, datatype, varargin{:})
iid = p.Results.iid;
eid = p.Results.eid;
fname_qc = p.Results.fname_qc;
qc_var = p.Results.qc_var;
ico = str2num_amd(p.Results.ico);
icnum = ico + 1;
ranknorm_wholeSample = str2num_amd(p.Results.ranknorm_wholeSample);
standardize_wholeSample = str2num_amd(p.Results.standardize_wholeSample);
wholeSampleTransform = p.Results.wholeSampleTransform;
fname_GRM = p.Results.GRM_file;
fname_preg = p.Results.preg_file;
fname_address = p.Results.address_file;
corrvec_thresh = p.Results.corrvec_thresh;
nframes_min = p.Results.nframes_min;
study = p.Results.study;
release = p.Results.release;

logging(FEMA_info)

% check if QC file and variable are both provided
if ~isempty(fname_qc)
    if isempty(qc_var)
        error('QC variable must be specified if QC file is provided');
    end
end

colnames_imaging=[];
missingness = [];

%%%%% load data %%%%%
switch datatype
    case 'voxel'
        %Load voxelwise imaging data
		logging('Reading voxelwise %s imaging data',fstem_imaging);
		tic
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
		mask = tmp_volinfo.vol_mask_sub;
		ivec_mask = find(mask>0.5);
		
		dirlist = tmp_volinfo.dirlist;

        %% old code 
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
        toc 

	case 'vertex'
		% Read in vertexwise imaging data
		logging('Reading vertexwise %s imaging data',fstem_imaging);
		tic
		%	SurfView_loadsurfs; % Shouldn't be neccessary, if data saved out pre-truncated
		% load('SurfView_surfs.mat'); %this matfile is included in the executable when compiling using -a
        load(fullfile(fileparts(fileparts(which('FEMA_process_data'))), 'showSurf', 'SurfView_surfs.mat'), 'icsurfs');

		hemistrings = {'lh','rh'};

		measmat = [];
		for hemii = 1:2
			hemi = hemistrings{hemii};
			fname = sprintf('%s/%s_%s.mat',dirname_imaging,fstem_imaging,hemi); 
			if ~exist(fname,'file') % if older release version
				fname = sprintf('%s/%s-%s.mat',dirname_imaging,fstem_imaging,hemi);  
			end
			% specified icnum
			logging('Reading vertexwise imaging data from %s',fname);
			tmp = load(fname);

			if ~isfield(tmp,'measmat') % Handle task fMRI
				% tmp.measmat = cat(2,measmat,tmp.betamat(:,1:size(icsurfs{icnum}.vertices,1))); % This fails miserably -- must be something off with the beta values
				% tmp.measmat = cat(2,measmat,tmp.betamat(:,1:size(icsurfs{icnum}.vertices,1))./tmp.semat(:,1:size(icsurfs{icnum}.vertices,1))); % Use z-scores per subject / run
				tmp.measmat = tmp.betamat./tmp.semat; % Use z-scores per subject / run
			end
			measmat = cat(2,measmat,tmp.measmat(:,1:size(icsurfs{icnum}.vertices,1)));
		end

		if 1 % Deal with Don's masking (replace NaNs w. zeros for all subjects,but not for subjects/visits with all missing data) - should fix concatenation script
			ivec_nan = find(mean(isfinite(measmat),2)==0);
			measmat(~isfinite(measmat)) = 0;
			measmat(ivec_nan,:) = NaN;
		else % Old version setting all NaNs to zero
			measmat(~isfinite(measmat)) = 0;
		end

		ivec_mask=find(nanstd(measmat)~=0);
		if length(ivec_mask)~=size(measmat,2)
			mask=zeros(1,size(measmat,2));
			measmat=measmat(:,ivec_mask);
			mask(ivec_mask)=1;
		else
			mask=[];
		end
		ymat = measmat;

        % load vol_info 
		fname_volinfo = sprintf('%s/vol_info.mat',dirname_imaging);
		tmp_volinfo = load(fname_volinfo); 
		dirlist = tmp_volinfo.dirlist;

        %% old code 
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
		toc

	case 'corrmat'
		logging('Processing corrmat %s data', fstem_imaging);
		if isstruct(dirname_imaging)
			tmp = dirname_imaging;
		else
			fname_corrmat = sprintf('%s/%s.mat',dirname_imaging,fstem_imaging);
			tmp = load(fname_corrmat);
		end
		if isfield(tmp,'nsumvec')	%	Release 4.0
			ivec_tmp = find(tmp.nsumvec>=nframes_min); 
			measmat = tmp.measmat(ivec_tmp,:);
		elseif isfield(tmp,'dirlist') %	Release 5.1
			ivec_tmp = find(tmp.ntpointvec>=nframes_min); 
			measmat = tmp.corrmat(ivec_tmp,:);
		else 
			ivec_tmp = find(tmp.ntpointvec>=nframes_min); 
			measmat = tmp.corrmat(ivec_tmp,:);
			tmp.dirlist = tmp_volinfo.dirlist;
		end
        % vol info
		fname_volinfo = sprintf('%s/vol_info.mat',dirname_imaging);
		tmp_volinfo = load(fname_volinfo);

		dirlist = tmp.dirlist(ivec_tmp);
		dims = size(measmat); 
        measmat = reshape(measmat,[dims(1) prod(dims(2:end))]);
		ymat = measmat;

        %% old code 
		    %subjidvec = cell(size(dirlist)); 
		    %sitevec = cell(size(dirlist));
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

    case 'roi'
        logging('Reading ROI tabulated data from %s', fstem_imaging);
        fname_roi = fullfile(dirname_imaging, fstem_imaging);
        [dirname_roi, fstem_roi, ext_roi] = fileparts(fname_roi);
        switch ext_roi 
            case '.parquet'
                roidat = parquetread(fname_roi);
            case {'.csv' '.tsv' '.txt'}
                opts = detectImportOptions(fname_roi, 'FileType', 'text');
                roidat = readtable(fname_roi, opts, 'FileType', 'text', 'ReadVariableNames', false);
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
	    ymat = table2array(removevars(roidat, {'participant_id', 'session_id'}));
	    colnames_imaging = removevars(roidat, {'participant_id', 'session_id'}).Properties.VariableNames;
	    idevent = strcat(iid_concat(:),'_',eid_concat(:));

	    ivec_mask=[];
	    mask=[];

    case 'external'
	    logging('Reading tabulated imaging data from %s/%s', dirname_imaging, fstem_imaging);
        [~, ~, ext_imaging] = fileparts(dirname_imaging);
        if ~isempty(fstem_imaging)
            switch ext_imaging
                case '.parquet'
                    opts = parquetinfo(dirname_imaging);
                    varInc = fstem_imaging; % add iid eid to varInc 
                    varPresent = ismember(varInc, opts.VariableNames);
                    if any(~varPresent)
                        missingVars = varInc(~varPresent);
                        error(['Could not find variables in table ', fname_qc, ': ', strjoin(missingVars, ', ')]);
                    else
                        opts.SelectedVariableNames = fstem_imaging;  % choose specific columns
                        try 
                            roidat = parquetread(dirname_imaging);
                        catch ME
                            error(sprintf('Cannot load tabulated data from %s', dirname_imaging));
                            sprintf('Error: %s', rethrow(ME));
                        end
                    end 
                case {'.csv' '.tsv' '.txt'}
                    opts = detectImportOptions(dirname_imaging, 'FileType', 'text');
                    varInc = fstem_imaging; % add iid eid to varInc
                    varPresent = ismember(varInc, opts.VariableNames);
                    if any(~varPresent)
                        missingVars = varInc(~varPresent);
                        error(['Could not find variables in file ', dirname_imaging, ': ', strjoin(missingVars, ', ')]);
                    else
                        opts.SelectedVariableNames = fstem_imaging;  % choose specific columns
                        try
                            roidat = readtable(dirname_imaging, opts);
                        catch ME
                            error(sprintf('Cannot load tabulated data from %s', dirname_imaging));
                            sprintf('Error: %s', rethrow(ME));
                        end 
                    end 
            end
        else 
            switch ext_imaging
                case '.parquet'
                    roidat = readtable(dirname_imaging);
                case {'.csv' '.tsv' '.txt'}
                    opts = detectImportOptions(dirname_imaging, 'FileType', 'text');
                    roidat = readtable(dirname_imaging, opts);
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
        ymat = table2array(removevars(roidat, {'participant_id', 'session_id'}));
	    colnames_imaging = removevars(roidat, {'participant_id', 'session_id'}).Properties.VariableNames;
	    idevent = strcat(iid_concat(:),'_',eid_concat(:));
	    ivec_mask=[];
	    mask=[];

    end

    %%%%% missingness %%%%%
    % number of nans in ymat  
    defvecNaN = isfinite(sum(ymat, 2)); 
    numNaN = length(find(defvecNaN == 0));
    idNaN = strcat(iid_concat(~defvecNaN), '_', eid_concat(~defvecNaN));
    missingness.numNaN = numNaN;
    missingness.idNaN = idNaN; 

    %  how many don't pass threshold on correlation 
    if strcmp(datatype, 'voxel')
    	corrvec = mean(corrmat_concat,2);
        thresh = corrvec_thresh;
    	defvec_corrvec = find(corrvec>=thresh);
        numFailCorr = length(find(defvec_corrvec == 0));
        idFailCorr = strcat(iid_concat(~defvec_corrvec), '_', eid_concat(~defvec_corrvec));
        missingness.numCorr = numFailCorr; % number removed dur to poor registration
        missingness.idCorr = idFailCorr;
    else
        defvec_corrvec = true(size(ymat, 1), 1);
    end 

    % filter on qc if qc file is given 
    if ~isempty(fname_qc)
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
                    opts.SelectedVariableNames = {'participant_id', 'session_id', qc_var};  % choose specific columns
                    try 
                        qc_tbl = parquetread(fname_qc);
                    catch ME
                        error(sprintf('Cannot load QC file %s with QC variable %s', fname_qc, qc_var));
                        sprintf('Error: %s', rethrow(ME));
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
                        error(sprintf('Cannot load QC file %s with QC variable %s', fname_qc, qc_var));
                        sprintf('Error: %s', rethrow(ME));
                    end 
                end 
        end
        passQC_idevent = strcat(qc_tbl.participant_id, '_', string(qc_tbl.session_id));
        passQC_idevent = passQC_idevent((qc_tbl.(qc_var) == '1'));
        defvecQC = ismember(idevent, passQC_idevent);
        idFailQC = idevent(~defvecQC);
        missingness.numFailQC = sum(~defvecQC);    
        missingness.idFailQC = idFailQC;
    else 
        defvecQC = true(size(ymat, 1), 1);
    end 

    % remove nans, fail corrvec, fail qc 
    defvec =  defvecNaN & defvec_corrvec & defvecQC; 
    iid_concat = iid_concat(defvec);
    eid_concat = eid_concat(defvec);
    idevent = idevent(defvec);
    ymat = ymat(defvec, :);
    missingness.numFailAny = sum(~defvec);
    missingness.idFailAny = idevent(~defvec);

    % inverse ranknorm or standardize on qc'd data 
    if ranknorm_wholeSample
    	logging('Rank-norming Y');
    	%ymat=rank_based_INT(ymat);
        [ymat settingsTransform] = doTransformation(ymat, 'ranknorm');
    end

    if standardize_wholeSample
    	logging('Standardizing Y');
    	%ymat=normalize(ymat);
        [ymat settingsTransform] = doTransformation(ymat, 'standardize');
    end

    % filter on iid and eid if given 
    % iid 
    if ~isempty(iid) 
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
    else 
        defvec_iid = true(size(iid_concat));
    end
    % eid 
    if ~isempty(eid)
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
    else 
        defvec_eid = true(size(eid_concat));
    end

    % filter on iid and eid if supplied
    defvec = defvec_iid & defvec_eid;
    iid_concat = iid_concat(defvec);
    eid_concat = eid_concat(defvec);
    idevent = idevent(defvec);
    ymat = ymat(defvec, :);

    % how many nans in iid- and eid- filtered data
    % nans in filtered data
    if ~isempty(eid) || ~isempty(iid)
        defvecNaN_filtered = defvecNaN(defvec) & defvec_iid & defvec_eid;
        missingness.numNaN_filtered = sum(~defvecNaN_filtered);
        missingness.idNaN_filtered = idevent(~defvecNaN_filtered);
        % fail corrvec in filtered data
        defvec_corrvec_filtered = defvec_corrvec(defvec) & defvec_iid & defvec_eid;
        missingness.numFailCorrvec_filtered = sum(~defvec_corrvec_filtered);
        missingness.idFailCorrvec_filtered = idevent(~defvec_corrvec_filtered);
        % fail qc in filtered data
        defvecQC_filtered = defvecQC(defvec) & defvec_iid & defvec_eid;
        missingness.numFailQC_filtered = sum(~defvecQC_filtered);
        missingness.idFailQC_filtered = idevent(~defvecQC_filtered);
        %  fail any in filtered data
        defvecAny_filtered = defvecNaN(defvec) & defvec_corrvec(defvec) & ...
                             defvecQC(defvec) & defvec_iid & defvec_eid;
        missingness.numFailAny_filtered = sum(~defvecAny_filtered);
        missingness.idFailAny_filtered = idevent(~defvecAny_filtered);
    end 

    % for roi and external? user input?  
        %%%%% load GRM %%%%%
    if ~isempty(fname_GRM)
        tic
    	logging('Reading GRM');
    	GRM = load(fname_GRM);
        iid_grm = GRM.iid_list; 
        [keep, IA, IB] = intersect(iid_concat, iid_grm, 'stable');
        defvecGRM = ismember(iid_concat, keep);
        % missingness
        missingness.numGRM = sum(~defvecGRM);
        missingness.idGRM = idevent(~defvecGRM);
        % filter on GRM
        iid_concat = iid_concat(defvecGRM);
        eid_concat = eid_concat(defvecGRM);
        idevent = idevent(defvecGRM);
        ymat = ymat(defvecGRM, :);
        GRM.GRM = GRM.GRM(IB,IB);
        GRM.iid_list = iid_grm(IB);
        toc 
    else 
    	GRM=[];
    end

    %%%%% load pregnancy data %%%%%
    if ~isempty(fname_preg)
	    logging('Reading pregnancy data');
	    fid = fopen(fname_preg);
	    varNames = strsplit(fgetl(fid), {' ' '"'});
	    fclose(fid);
	    opts=detectImportOptions(fname_preg);
	    opts.SelectedVariableNames = [2:5];
	    preg = readtable(fname_preg, opts);
	    preg = renamevars(preg,1:width(preg),varNames(2:5));
    else
    	preg=[];
    end

    %%%%% load home address data %%%%%
    if ~isempty(fname_address)
    	logging('Reading home address data');
    	address = readtable(fname_address);
    else
    	address=[];
    end

    % final sample 
    fprintf('Final sample size: %d\n', size(ymat, 1));
end 



