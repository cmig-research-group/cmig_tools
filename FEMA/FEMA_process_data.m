function [ymat, iid_concat, eid_concat, ivec_mask, mask, colnames_imaging, pihat, preg, address] = FEMA_process_data(fstem_imaging,dirname_imaging,datatype,varargin)
%
% ABCD specific function to load and process imaging data from abcd-sync
%
% INPUTS
%	 fstem_imaging <char>		:	name of vertex/voxel-mapped phenotype (e.g., 'thickness-sm16', 'FA')
%	 dirname_imaging <char>		:	path to imaging data directory
%	 datatype <char>			:	'voxel','vertex','external', 'corrmat'
%									Other then 'external' all code is written to expect ABCD data in same format as abcd-sync
%
% Optional input arguments:
%	 ico <num>					:	ico-number for vertexwise analyses (0-based, default 5)
%	 ranknorm <boolean>			:	rank normalise imaging data (default 0)
%	 varnorm <boolean>			:	variance normalise imaging data (default 0)
%	 pihat_file <char>			:	path to genetic relatedness data (pihat) - default [] - only required if A random effect specified
%	 preg_file <char>			:	path to pregnancy data - default [] - only required if T random effect specified
%	 address_file <char>		:	path to address data - default [] - only required if H random effect specified
%
% OUTPUTS
%	 ymat						:	matrix of imaging data (n x v)
%	 iid_concat					:	src_subject_id
%	 eid_concat					:	eventname
%	 ivec_mask					:	vector mask for voxelwise data (155179x1) --> not yet available for vertexwise
%	 mask						:	volume mask mask for voxelwise data (100x100x130) --> not yet available for vertexwise
%	 colnames_imaging			:	imaging column labels for external data inputs
%	 pihat						:	intersected genetic relatedness matrix
%

%
% This software is Copyright (c) 2021 The Regents of the University of California. All Rights Reserved.
% See LICENSE.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

p = inputParser;
addParamValue(p,'ranknorm',0);
addParamValue(p,'varnorm',0);
addParamValue(p,'ico',5);
addParamValue(p,'pihat_file',[]);
addParamValue(p,'preg_file',[]);
addParamValue(p,'address_file',[]);
addParamValue(p,'corrvec_thresh',0.8);

parse(p,varargin{:})
ico = str2num_amd(p.Results.ico);
icnum = ico + 1;
ranknorm = str2num_amd(p.Results.ranknorm);
varnorm = str2num_amd(p.Results.varnorm);
fname_pihat = p.Results.pihat_file;
fname_preg = p.Results.preg_file;
fname_address = p.Results.address_file;
corrvec_thresh = p.Results.corrvec_thresh;

readvolumesflag = true;
if fstem_imaging(1)=='-'
	readvolumesflag = false;
	fstem_imaging = fstem_imaging(2:end);
end

if ~strcmpi(datatype,'external') %differences between releases not relevant for external data

	colnames_imaging=[];

	if strcmpi(datatype, 'voxel')
		%Load voxelwise imaging data
		logging('Reading VOXELWISE %s imaging data',fstem_imaging);
		tic
		measname = fstem_imaging;
		dirname_volmats = dirname_imaging;
	
		% 6.0 naming convention
		fname_volinfo = sprintf('%s/vol_info.mat',dirname_volmats);
		if ~exist(fname_volinfo,'file') % if older release version 
			fname_volinfo = sprintf('%s/volinfo.mat',dirname_volmats);
		end
		tmp_volinfo = load(fname_volinfo);

		if readvolumesflag
			fname_volmat = sprintf('%s/%s.mat',dirname_volmats,lower(measname));
			if ~exist(fname_volmat,'file') % if older release version 
				fname_volmat = sprintf('%s/volmat_%s.mat',dirname_volmats,measname);
			end	
			ymat = getfield(load(fname_volmat),'volmat');
			mask = tmp_volinfo.vol_mask_sub;
			ivec_mask = find(mask>0.5);
		else
			ymat = NaN([length(tmp_volinfo.dirlist) 0]);
			mask = [];
			ivec_mask =	[];
		end
		ymat_bak = ymat; 

		dirlist = tmp_volinfo.dirlist;
        subjidvec = cell(size(dirlist)); 
		sitevec = cell(size(dirlist)); 
		datevec = cell(size(dirlist)); 
		visitidvec = cell(size(dirlist));
        if any(contains(dirlist, 'DTIREG_'))
            for diri = 1:length(dirlist)
            	tmp = regexp(dirlist{diri}, '^DTIREG_(?<site>\w+)_(?<SubjID>\w+)_(?<event>\w+)_(?<date>\d+).', 'names');
            	subjidvec{diri} = tmp.SubjID;
            	sitevec{diri} = tmp.site;
            	eventvec{diri} = tmp.event;
            	datevec{diri} = tmp.date;
            	visitidvec{diri} = sprintf('%s_%s_%s',tmp.site,tmp.SubjID,tmp.event);
            end
            iid_concat = colvec(strcat('NDAR_',subjidvec)); % Not sure why imaging subject IDs are missing the NDAR_ part
            eid_concat = eventvec';
        else 
            iid_concat = tmp_volinfo.participant_id;
            eid_concat = tmp_volinfo.session_id;
        end 
        idevent = strcat(iid_concat(:),'_',eid_concat(:));
        corrmat_concat = tmp_volinfo.corrmat;

		%QC FOR VOXELWISE DATA
		corrvec=mean(corrmat_concat,2);
        thresh = corrvec_thresh;
		defvec=find(corrvec>=thresh); 
		ymat=ymat(defvec,:);
		idevent=idevent(defvec,:);
		iid_concat=iid_concat(defvec,:);
		eid_concat=eid_concat(defvec,:);
		toc

	elseif strcmpi(datatype, 'vertex')
		% Read in vertexwise imaging data
		logging('Reading VERTEXWISE %s imaging data',fstem_imaging);
		tic
		%	SurfView_loadsurfs; % Shouldn't be neccessary, if data saved out pre-truncated
		% load('SurfView_surfs.mat'); %this matfile is included in the executable when compiling using -a
        load(fullfile(fileparts(fileparts(which('FEMA_process_data'))), 'showSurf', 'SurfView_surfs.mat'), 'icsurfs');

		hemistrings = {'lh','rh'};

		measmat = [];
		for hemii = 1:2
			hemi = hemistrings{hemii};
			fname = sprintf('%s/%s_%s.mat',dirname_imaging,fstem_imaging,hemi); % Should save these pre-truncated to 
			if ~exist(fname,'file') % if older release version
				fname = sprintf('%s/%s-%s.mat',dirname_imaging,fstem_imaging,hemi); % Should save these pre-truncated to 
			end
			% specified icnum
			logging('Reading vertexwise imaging data from %s',fname);
			tmp = load(fname);

			if ~isfield(tmp,'measmat') % Handle task fMRI
				%									 tmp.measmat = cat(2,measmat,tmp.betamat(:,1:size(icsurfs{icnum}.vertices,1))); % This fails miserably -- must be something off with the beta values
				%										tmp.measmat = cat(2,measmat,tmp.betamat(:,1:size(icsurfs{icnum}.vertices,1))./tmp.semat(:,1:size(icsurfs{icnum}.vertices,1))); % Use z-scores per subject / run
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
		fname_volinfo = sprintf('%s/vol_info.mat',dirname_imaging);
		tmp_volinfo = load(fname_volinfo); 
		try 
			dirlist = tmp_volinfo.dirlist;
		catch 
			dirlist = tmp.dirlist;
		end 
		subjidvec = cell(size(dirlist)); 
		sitevec = cell(size(dirlist)); 
		datevec = cell(size(dirlist)); 
		eventvec = cell(size(dirlist)); 
		timevec = cell(size(dirlist)); 
		visitidvec = cell(size(dirlist));

		for diri = 1:length(dirlist)
			tmp = regexp(dirlist{diri}, '^[^_]*_(?<site>\w+)_(?<SubjID>\w+)_(?<event>[^_]*)_(?<date>\d+).(?<time>[^_]+)_', 'names'); % AMD -- 3.0 fails for task fMRI data
			sitevec{diri} = tmp.site;
			subjidvec{diri} = tmp.SubjID;
			datevec{diri} = tmp.date;
			if isfield(tmp,'event')
				eventvec{diri} =	tmp.event;
			else
				eventvec{diri} = '';
			end
			visitidvec{diri} = sprintf('%s_%s_%s',tmp.site,tmp.SubjID,tmp.event);
		end

		iid_concat = colvec(strcat('NDAR_',subjidvec)); % Not sure why imaging subject IDs are missing the NDAR_ part
		eid_concat = eventvec;
		idevent = strcat(iid_concat(:),'_',eid_concat(:));
		toc

	elseif strcmpi(datatype, 'corrmat')
		logging('Processing CORRMAT');
		if isstruct(dirname_imaging)
			tmp = dirname_imaging;
		else
			fname_corrmat = sprintf('%s/%s.mat',dirname_imaging,fstem_imaging);
			tmp = load(fname_corrmat);
		end
		nframes_min = 375; % This should be optional	input param
		if isfield(tmp,'nsumvec')	%	Release 4.0
			ivec_tmp = find(tmp.nsumvec>=nframes_min); 
			measmat = tmp.measmat(ivec_tmp,:);
		elseif isfield(tmp,'dirlist') %	Release 5.1
			ivec_tmp = find(tmp.ntpointvec>=nframes_min); 
			measmat = tmp.corrmat(ivec_tmp,:);
		else 
			ivec_tmp = find(tmp.ntpointvec>=nframes_min); 
			measmat = tmp.corrmat(ivec_tmp,:);
			fname_volinfo = sprintf('%s/vol_info.mat',dirname_imaging);
			tmp2 = load(fname_volinfo);
			tmp.dirlist = tmp2.dirlist;
		end
		dirlist = tmp.dirlist(ivec_tmp);
		dims = size(measmat); measmat = reshape(measmat,[dims(1) prod(dims(2:end))]);
		ymat = measmat;
		subjidvec = cell(size(dirlist)); 
		sitevec = cell(size(dirlist));
		datevec = cell(size(dirlist)); 
		visitidvec = cell(size(dirlist));
		for diri = 1:length(dirlist)
			tmp = regexp(dirlist{diri}, '^BOLDPROC_(?<site>\w+)_(?<SubjID>\w+)_(?<event>[^_]*)_(?<date>[^_.]+).(?<time>[^_]+)_', 'names');
			subjidvec{diri} = tmp.SubjID;
			sitevec{diri} = tmp.site;
			eventvec{diri} = tmp.event;
			datevec{diri} = tmp.date;
			visitidvec{diri} = sprintf('%s_%s_%s',tmp.site,tmp.SubjID,tmp.event);
		end
		iid_concat = colvec(strcat('NDAR_',subjidvec)); % Not sure why imaging subject IDs are missing the NDAR_ part
		eid_concat = eventvec;
		idevent = strcat(iid_concat(:),'_',eid_concat(:));

		ivec_mask=[];
		mask=[];

	end

elseif strcmpi(datatype,'external') %does not need to be intersected with abcd mri info - extract info then go straight to intersecting with design matrix

	logging('Reading EXTERNAL imaging data');
	roidat = readtable(dirname_imaging);
	% Check which ID + event columns are present for the different releases
	if any(strcmp(roidat.Properties.VariableNames, 'src_subject_id'))
    	iid_concat = roidat.src_subject_id;
    	eid_concat = roidat.eventname;
	elseif any(strcmp(roidat.Properties.VariableNames, 'participant_id'))
    	iid_concat = roidat.participant_id;
    	eid_concat = roidat.session_id;
	else
    	error('No recognized ID/event columns found in the external file.');
	end
	ymat=table2array(roidat(:,3:end));
	colnames_imaging=roidat.Properties.VariableNames(3:end);
	idevent = strcat(iid_concat(:),'_',eid_concat(:));

	ivec_mask=[];
	mask=[];

end


if ~isempty(fname_pihat)
	logging('Reading PI HAT');
	pihat = load(fname_pihat);
else
	pihat=[];
end

if ~isempty(fname_preg)
	logging('Reading PREGNANCY ID');
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

if ~isempty(fname_address)
	logging('Reading HOME ID');
	address = readtable(fname_address);
else
	address=[];
end

if ranknorm==1
	logging('Rank-norming Y');
	ymat=rank_based_INT(ymat);
end

if varnorm==1
	logging('Variance-norming Y');
	ymat=normalize(ymat);
end