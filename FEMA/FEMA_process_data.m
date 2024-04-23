function [ymat, iid_concat, eid_concat, ivec_mask, mask, colnames_imaging, pihat, preg, address] = FEMA_process_data(fstem_imaging,dirname_tabulated,dirname_imaging,datatype,varargin)
%
% ABCD specific function to load and process imaging data from abcd-sync
%
% INPUTS
%   fstem_imaging <char>       :  name of vertex/voxel-mapped phenotype (e.g., 'thickness-sm16', 'FA')
%   dirname_tabulated <char>   :  path to tabulated data directory with NDA downloaded txt files
%   dirname_imaging <char>     :  path to imaging data directory
%   datatype <char>            :  'voxel','vertex','external', 'corrmat'
%                                    %Other then 'external' all code is written to expect ABCD data in same format as abcd-sync
%
% Optional input arguments:
%   ico <num>                  :  ico-number for vertexwise analyses (0-based, default 5)
%   ranknorm <boolean>         :  rank normalise imaging data (default 0)
%   pihat_file <char>          :  path to genetic relatedness data (pihat) - default [] - only required if A random effect specified
%   preg_file <char>           :  path to pregnancy data - default [] - only required if T random effect specified
%   address_file <char>        :  path to address data - default [] - only required if H random effect specified
%
% OUTPUTS
%   ymat                       :  matrix of imaging data (n x v)
%   iid_concat                 :  src_subject_id
%   eid_concat                 :  eventname
%   ivec_mask                  :  vector mask for voxelwise data (155179x1) --> not yet available for vertexwise
%   mask                   :  volume mask mask for voxelwise data (100x100x130) --> not yet available for vertexwise
%   colnames_imaging           :  imaging column labels for external data inputs
%   pihat                      :  intersected genetic relatedness matrix
%

%
% This software is Copyright (c) 2021 The Regents of the University of California. All Rights Reserved.
% See LICENSE.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

p = inputParser;
addParamValue(p,'ranknorm',0);
addParamValue(p,'ico',5);
addParamValue(p,'pihat_file',[]);
addParamValue(p,'preg_file',[]);
addParamValue(p,'address_file',[]);

parse(p,varargin{:})
ico = str2num_amd(p.Results.ico);
icnum = ico + 1;
ranknorm = str2num_amd(p.Results.ranknorm);
fname_pihat = p.Results.pihat_file;
fname_preg = p.Results.preg_file;
fname_address = p.Results.address_file;

readvolumesflag = true;
if fstem_imaging(1)=='-'
  readvolumesflag = false;
  fstem_imaging = fstem_imaging(2:end);
end

if ~strcmpi(datatype,'external') %differences between releases not relevant for external data

  colnames_imaging=[];

  %get release version
  reltmp = regexp(dirname_tabulated,'abcd-sync/(?<release>[^/]*)','names');
  if isempty(reltmp)
    error(['The directory structure does not follow the expected scheme.\n',...
      '  The `dirname_tabulated` should contain `.../abcd-sync/#VERSION_NUMBER/...`,\n',...
      '  e.g. `dirname_tabulated=/path/to/abcd-sync/4.0/tabulated/img`.\n',...
      '  yours is: %s'], dirname_tabulated)
    exit
  end
  dataRelease = reltmp.release;

  if strcmpi(datatype, 'voxel')
    %Load voxelwise imaging data
    logging('Reading VOXELWISE %s imaging data',fstem_imaging);
    tic
    measname = fstem_imaging;
    dirname_volmats = dirname_imaging;
    fname_volinfo = sprintf('%s/volinfo.mat',dirname_volmats);
    tmp_volinfo = load(fname_volinfo);

    if readvolumesflag
      fname_volmat = sprintf('%s/volmat_%s.mat',dirname_volmats,measname);
      ymat = getfield(load(fname_volmat),'volmat');
      mask = tmp_volinfo.vol_mask_sub;
      ivec_mask = find(mask>0.5);
    else
      ymat = NaN([length(tmp_volinfo.dirlist) 0]);
      mask = [];
      ivec_mask =  [];
    end

    switch dataRelease
      case '3.0'
        iid_concat = colvec(strcat('NDAR_',tmp_volinfo.subjidvec)); % Not sure why imaging subject IDs are missing the NDAR_ part
        visitid_concat = tmp_volinfo.visitidvec;
        datenumlist_concat = datenum(tmp_volinfo.datevec,'yyyymmdd');
        corrmat_concat = tmp_volinfo.corrmat;
      case {'4.0','5.0', '6.0'}
        dirlist = tmp_volinfo.dirlist;
        subjidvec = cell(size(dirlist)); sitevec = cell(size(dirlist)); datevec = cell(size(dirlist)); visitidvec = cell(size(dirlist));
        for diri = 1:length(dirlist)
          tmp = regexp(dirlist{diri}, '^DTIREG_(?<site>\w+)_(?<SubjID>\w+)_(?<event>\w+)_(?<date>\d+).', 'names');
          subjidvec{diri} = tmp.SubjID;
          sitevec{diri} = tmp.site;
          eventvec{diri} = tmp.event;
          datevec{diri} = tmp.date;
          visitidvec{diri} = sprintf('%s_%s_%s',tmp.site,tmp.SubjID,tmp.event);
        end
        iid_concat = colvec(strcat('NDAR_',subjidvec)); % Not sure why imaging subject IDs are missing the NDAR_ part
        visitid_concat = visitidvec;
        corrmat_concat = tmp_volinfo.corrmat;
      otherwise
        warning('Data release %s may not work properly. Check the code that parses voxelwise data directory names')
    end
    toc

  elseif strcmpi(datatype, 'vertex')
    % Read in vertexwise imaging data
    logging('Reading VERTEXWISE %s imaging data',fstem_imaging);
    tic
    %  SurfView_loadsurfs; % Shouldn't be neccessary, if data saved out pre-truncated
    load('SurfView_surfs.mat'); %this matfile is included in the executable when compiling using -a

    hemistrings = {'lh','rh'};

    measmat = [];
    for hemii = 1:2
      hemi = hemistrings{hemii};
      fname = sprintf('%s/%s-%s.mat',dirname_imaging,fstem_imaging,hemi); % Should save these pre-truncated to specified icnum
      logging('Reading vertexwise imaging data from %s',fname);
      tmp = load(fname);

      if ~isfield(tmp,'measmat') % Handle task fMRI
        %                   tmp.measmat = cat(2,measmat,tmp.betamat(:,1:size(icsurfs{icnum}.vertices,1))); % This fails miserably -- must be something off with the beta values
        %                    tmp.measmat = cat(2,measmat,tmp.betamat(:,1:size(icsurfs{icnum}.vertices,1))./tmp.semat(:,1:size(icsurfs{icnum}.vertices,1))); % Use z-scores per subject / run
        tmp.measmat = tmp.betamat./tmp.semat; % Use z-scores per subject / run
      end

      if 0
        ic = NaN;
        for ic = 1:length(icsurfs)
          if size(icsurfs{ic}.vertices,1) == size(tmp.measmat,2)
            icnum_true = ic;
          end
        end
        if icnum~=icnum_true
          icnum=icnum_true;
          warning('Wrong ico entered for vertexwise data loaded. This has been edited to the correct ico.')
        end
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
    dirlist = tmp.dirlist;
    subjidvec = cell(size(dirlist)); sitevec = cell(size(dirlist)); datevec = cell(size(dirlist)); eventvec = cell(size(dirlist)); timevec = cell(size(dirlist)); visitidvec = cell(size(dirlist));

    for diri = 1:length(dirlist)
      switch dataRelease
        case '3.0'
          tmp = regexp(dirlist{diri}, '^FSURF_(?<site>\w+)_(?<SubjID>\w+)_(?<date>\d+)_', 'names'); % 3.0 directory naming convention
        case {'4.0', '5.0', '6.0'}
          tmp = regexp(dirlist{diri}, '^[^_]*_(?<site>\w+)_(?<SubjID>\w+)_(?<event>[^_]*)_(?<date>\d+).(?<time>[^_]+)_', 'names'); % 4.0 directory naming convention (with 3.0 data, tmp.event has duplicate of date
      end

      tmp = regexp(dirlist{diri}, '^[^_]*_(?<site>\w+)_(?<SubjID>\w+)_(?<event>[^_]*)_(?<date>\d+).(?<time>[^_]+)_', 'names'); % AMD -- 3.0 fails for task fMRI data
      sitevec{diri} = tmp.site;
      subjidvec{diri} = tmp.SubjID;
      datevec{diri} = tmp.date;
      if isfield(tmp,'event')
        eventvec{diri} =  tmp.event;
      else
        eventvec{diri} = '';
      end
      switch dataRelease
        case '3.0'
          visitidvec{diri} = sprintf('%s_%s_%s',tmp.site,tmp.SubjID,tmp.date);
        case {'4.0', '5.0', '6.0'}
          visitidvec{diri} = sprintf('%s_%s_%s',tmp.site,tmp.SubjID,tmp.event);
      end
    end

    iid_concat = colvec(strcat('NDAR_',subjidvec)); % Not sure why imaging subject IDs are missing the NDAR_ part
    visitid_concat = visitidvec;
    toc


  elseif strcmpi(datatype, 'corrmat')
    logging('Processing CORRMAT');
    if isstruct(dirname_imaging)
      tmp = dirname_imaging;
    else
      fname_corrmat = sprintf('%s/%s.mat',dirname_imaging,fstem_imaging);
      tmp = load(fname_corrmat);
    end
    nframes_min = 375; % This should be optional  input param
    if isfield(tmp,'nsumvec')  %  Release 4.0
      ivec_tmp = find(tmp.nsumvec>=nframes_min); 
      measmat = tmp.measmat(ivec_tmp,:);
    else %  Release 5.1
      ivec_tmp = find(tmp.ntpointvec>=nframes_min); 
      measmat = tmp.corrmat(ivec_tmp,:);
    end
    dirlist = tmp.dirlist(ivec_tmp);
    dims = size(measmat); measmat = reshape(measmat,[dims(1) prod(dims(2:end))]);
    ymat = measmat;
    subjidvec = cell(size(dirlist)); sitevec = cell(size(dirlist)); datevec = cell(size(dirlist)); visitidvec = cell(size(dirlist));
    for diri = 1:length(dirlist)
      %    tmp = regexp(dirlist{diri}, '^BOLDPROC_(?<site>\w+)_(?<SubjID>\w+)_(?<date>\d+)_', 'names');
      tmp = regexp(dirlist{diri}, '^BOLDPROC_(?<site>\w+)_(?<SubjID>\w+)_(?<event>[^_]*)_(?<date>[^_.]+).(?<time>[^_]+)_', 'names');
      %    tmp = regexp(dirlist{diri}, '^[^_]*_(?<site>\w+)_(?<SubjID>\w+)_(?<event>[^_]*)_(?<date>[^_.]+).(?<time>[^_]+)_', 'names'); % 4.0 directory naming convention (with 3.0 data, tmp.event has duplicate of date
      subjidvec{diri} = tmp.SubjID;
      sitevec{diri} = tmp.site;
      eventvec{diri} = tmp.event;
      datevec{diri} = tmp.date;
      visitidvec{diri} = sprintf('%s_%s_%s',tmp.site,tmp.SubjID,tmp.event);
    end
    iid_concat = colvec(strcat('NDAR_',subjidvec)); % Not sure why imaging subject IDs are missing the NDAR_ part
    visitid_concat = visitidvec;

    ivec_mask=[];
    mask=[];

  end

  switch dataRelease
    case '3.0'
      tic
      files=dir([dirname_tabulated '/abcd_mri01*']);
      fname_tabulated=[dirname_tabulated '/' files.name];
      logging('Reading tabulated imaging data from %s',fname_tabulated);
      imgtable = readtable(fname_tabulated);
      iid_imgtable = imgtable.subjectkey;
      eid_imgtable = imgtable.eventname;
      visitid_imgtable = imgtable.mri_info_visitid;
      toc

      % Merge vertex/voxelwise and tabulated imaging data
      tmplist_concat = visitid_concat;
      tmplist_imgtable = visitid_imgtable;
      [dummy IA IB] = intersect(tmplist_concat,tmplist_imgtable,'stable');
      eid_concat = eid_imgtable(IB);
      iid_concat = iid_concat(IA);
      idevent = strcat(iid_concat(:),'_',eid_concat(:));
      imgtable=imgtable(IB,:);
      ymat=ymat(IA,:);

    case '4.0'
      tic
      files=dir([dirname_tabulated '/abcd_mri01*']);
      fname_tabulated=[dirname_tabulated '/' files.name];
      logging('Reading tabulated imaging data from %s',fname_tabulated);
      imgtable = readtable(fname_tabulated);
      imgtable.idevent=strcat(imgtable.src_subject_id,'_',imgtable.eventname);

      files=dir([dirname_tabulated '/abcd_imgincl01*']);
      fname_incflag=[dirname_tabulated '/' files.name];
      inctable = readtable(fname_incflag);
      inctable.idevent=strcat(inctable.src_subject_id,'_',inctable.eventname);

      if ismember('dataset_id',imgtable.Properties.VariableNames)
        imgtable=removevars(imgtable,'dataset_id');
        inctable=removevars(inctable,'dataset_id');
      end

      [dummy IA IB] = intersect(imgtable.idevent,inctable.idevent,'stable');
      imgtable=join(imgtable(IA,:),inctable(IB,:));

      iid_imgtable = imgtable.subjectkey;
      eid_imgtable = imgtable.eventname;
      date_imgtable = imgtable.mri_info_studydate;
      visitid_imgtable = imgtable.mri_info_visitid;

      %for dati=1:length(date_imgtable)
      %tmp = regexp(imgtable.mri_info_studydate{2}, '^(?<month>\d+)/(?<day>\d+)/(?<year>\d+)', 'names');
      %datevec{dati}=sprintf('%s%s%s',tmp.year,tmp.month,tmp.day);
      %      visitid_imgtable{dati}=sprintf('%s_%d',visitid_imgtable{dati},date_imgtable(dati));
      %end
      toc

      % Merge vertex/voxelwise and tabulated imaging data
      tmplist_concat = visitid_concat;
      tmplist_imgtable = visitid_imgtable;
      [dummy IA IB] = intersect(tmplist_concat,tmplist_imgtable,'stable');
      eid_concat = eid_imgtable(IB);
      iid_concat = iid_concat(IA);
      idevent = strcat(iid_concat(:),'_',eid_concat(:));
      imgtable=imgtable(IB,:);
      ymat=ymat(IA,:);

    case '5.0'
      tic
      % if using "released"
      if contains(dirname_tabulated, 'released')
        fname_mri_info = 'mri_y_adm_info';
        fname_imgincl = 'mri_y_qc_incl';
      elseif contains(dirname_tabulated, 'img')
        fname_mri_info = 'abcd_mri01';
        fname_imgincl = 'abcd_imgincl01';
      end
      % mri_info file
      files=dir([dirname_tabulated '/', fname_mri_info, '*']);
      fname_tabulated=[dirname_tabulated '/' files.name];
      logging('Reading tabulated imaging data from %s',fname_tabulated);
      imgtable = readtable(fname_tabulated);
      imgtable.idevent=strcat(imgtable.src_subject_id,'_',imgtable.eventname);
      % imgincl file
      files=dir([dirname_tabulated '/', fname_imgincl, '*']);
      fname_incflag=[dirname_tabulated '/' files.name];
      inctable = readtable(fname_incflag);
      inctable.idevent=strcat(inctable.src_subject_id,'_',inctable.eventname);

      if ismember('dataset_id',imgtable.Properties.VariableNames)
        imgtable=removevars(imgtable,'dataset_id');
        inctable=removevars(inctable,'dataset_id');
      end

      [dummy IA IB] = intersect(imgtable.idevent,inctable.idevent,'stable');
      imgtable=join(imgtable(IA,:),inctable(IB,:));

      iid_imgtable = imgtable.src_subject_id;
      eid_imgtable = imgtable.eventname;
      date_imgtable = imgtable.mri_info_studydate;
      visitid_imgtable = imgtable.mri_info_visitid;

      %for dati=1:length(date_imgtable)
      %tmp = regexp(imgtable.mri_info_studydate{2}, '^(?<month>\d+)/(?<day>\d+)/(?<year>\d+)', 'names');
      %datevec{dati}=sprintf('%s%s%s',tmp.year,tmp.month,tmp.day);
      %      visitid_imgtable{dati}=sprintf('%s_%d',visitid_imgtable{dati},date_imgtable(dati));
      %end
      toc

      % Merge vertex/voxelwise and tabulated imaging data
      tmplist_concat = visitid_concat;
      tmplist_imgtable = visitid_imgtable;
      [dummy IA IB] = intersect(tmplist_concat,tmplist_imgtable,'stable');
      eid_concat = eid_imgtable(IB);
      iid_concat = iid_concat(IA);
      idevent = strcat(iid_concat(:),'_',eid_concat(:));
      imgtable=imgtable(IB,:);
      ymat=ymat(IA,:);
    
    case '6.0'
      tic
      % if using "released"
      if contains(dirname_tabulated, 'released')
        fname_mri_info = 'mri_y_adm_info';
        fname_imgincl = 'mri_y_qc_incl';
      elseif contains(dirname_tabulated, 'img')
        fname_mri_info = 'abcd_mri01';
        fname_imgincl = 'abcd_imgincl01';
      end
      % mri_info file
      files=dir([dirname_tabulated '/', fname_mri_info, '*']);
      fname_tabulated=[dirname_tabulated '/' files.name];
      logging('Reading tabulated imaging data from %s',fname_tabulated);
      imgtable = readtable(fname_tabulated);
      % imgtable.idevent=strcat(imgtable.src_subject_id,'_',imgtable.eventname);
      % imgincl file
      files=dir([dirname_tabulated '/', fname_imgincl, '*']);
      fname_incflag=[dirname_tabulated '/' files.name];
      inctable = readtable(fname_incflag);
      % inctable.idevent=strcat(inctable.src_subject_id,'_',inctable.eventname);

      [dummy IA IB] = intersect(imgtable.VisitID,inctable.VisitID,'stable');
      imgtable=join(imgtable(IA,:),inctable(IB,:));

      iid_imgtable = imgtable.src_subject_id;
      eid_imgtable = imgtable.eventname;
      date_imgtable = imgtable.mri_info_studydate;
      visitid_imgtable = imgtable.mri_info_visitid;

      %for dati=1:length(date_imgtable)
      %tmp = regexp(imgtable.mri_info_studydate{2}, '^(?<month>\d+)/(?<day>\d+)/(?<year>\d+)', 'names');
      %datevec{dati}=sprintf('%s%s%s',tmp.year,tmp.month,tmp.day);
      %      visitid_imgtable{dati}=sprintf('%s_%d',visitid_imgtable{dati},date_imgtable(dati));
      %end
      toc

      % Merge vertex/voxelwise and tabulated imaging data
      tmplist_concat = visitid_concat;
      tmplist_imgtable = visitid_imgtable;
      [dummy IA IB] = intersect(tmplist_concat,tmplist_imgtable,'stable');
      eid_concat = eid_imgtable(IB);
      iid_concat = iid_concat(IA);
      idevent = strcat(iid_concat(:),'_',eid_concat(:));
      imgtable=imgtable(IB,:);
      ymat=ymat(IA,:);

  end

  %QC FOR VOXELWISE DATA
  if strcmpi(datatype,'voxel')
    corrvec=mean(corrmat_concat(IA,:),2);
    thresh = 0.8;
    switch dataRelease
      case '3.0'
        defvec=find(corrvec>=thresh);
      case '4.0'
        defvec=find(corrvec>=thresh & str2double(imgtable.imgincl_dmri_include)==1 & str2double(imgtable.imgincl_t1w_include)==1);
      case {'5.0', '6.0'}
        defvec=find(corrvec>=thresh & imgtable.imgincl_dmri_include==1 & imgtable.imgincl_t1w_include==1); 
    end
    ymat=ymat(defvec,:);
    idevent=idevent(defvec,:);
    iid_concat=iid_concat(defvec,:);
    eid_concat=eid_concat(defvec,:);
  end


elseif strcmpi(datatype,'external') %does not need to be intersected with abcd mri info - extract info then go straight to intersecting with design matrix

  logging('Reading EXTERNAL imaging data');
  roidat = readtable(dirname_imaging);
  iid_concat = roidat.src_subject_id; % Not sure why imaging subject IDs are missing the NDAR_ part
  eid_concat = roidat.eventname;
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



