function FEMA_DEAP_job_assembly(fstem_imaging,dirname_out,dirname_imaging,dirname_cache,dirname_jobs,nfrac,X,colnames_model,varargin)

  if isdeployed
    nfrac = str2num(nfrac);
    if isnan(nfrac)
      error('One or more input arguments could not be converted to a number.');
    end
  end

  if nargin < 9
    logging('Usage: FEMA_DEAP_wrapper(fstem_imaging,dirname_out,dirname_imaging,dirname_cache,dirname_jobs,nfrac,X,colnames_model)');
    error('Incorrect number of input arguments')
  end

  % code only support voxel datatype for now
  datatype = 'voxel';

  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% Look into how include the following code as "macro"

  inputs = inputParser;
  % used variables
  addParamValue(inputs,'RandomEffects',{'F' 'S' 'E'}); % Default to Family, Subject, and eps
  addParamValue(inputs,'synth',0); % AMD - put back synth option
  addParamValue(inputs,'colsinterest',1);
  addParamValue(inputs,'contrasts',[]);
  addParamValue(inputs,'output','nifti'); % but we only support this format...
  addParamValue(inputs,'ivnames','');

  parse(inputs,varargin{:})
  contrasts = str2num_amd(inputs.Results.contrasts);
  if ~isfinite(contrasts)
    fname_contrasts = inputs.Results.contrasts;
    logging('Reading contrast matrix from %s',fname_contrasts);
    contrasts = readtable(fname_contrasts);
  end
  outputFormat = inputs.Results.output;
  if ~(contains(outputFormat,'mat') || contains(outputFormat,'nifti'))
    error('Incorrect output format (%s)',outputFormat)
  end

  if ~isempty(inputs.Results.ivnames)
    ivnames = split(inputs.Results.ivnames,',');
  else     
    ivnames = {}; 
  end

  RandomEffects = rowvec(split(strrep(strrep(inputs.Results.RandomEffects,'{',''),'}','')));
  synth = str2num_amd(inputs.Results.synth);
  colsinterest = str2num_amd(inputs.Results.colsinterest);

  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

  logging('***Start');

  starttime = now();

  dirname_results = sprintf('%s/results/%s_%s_%s',dirname_jobs,designid,datatype,fstem_imaging);

  if exist(dirname_results,'dir')
    cmd = sprintf('rm -rf %s',dirname_results);
    [s r e] = jsystem(cmd);
  end

  dirname_volmats = dirname_imaging;
  fname_volinfo = sprintf('%s/volinfo.mat',dirname_volmats);
  mask = getfield(load(fname_volinfo,'vol_mask_sub'),'vol_mask_sub'); 
  ivec_mask = find(mask>0.5); % Should save and read ivec_mask instead
  zmat = zeros(size(X,2),length(ivec_mask));
  logpmat = zeros(size(X,2),length(ivec_mask));
  beta_hat = zeros(size(X,2),length(ivec_mask));
  beta_se = zeros(size(X,2),length(ivec_mask));
  sig2mat = zeros(length(RandomEffects),length(ivec_mask));
  sig2tvec = zeros(1,length(ivec_mask));

  % load results from each worker
  for fraci = 1:nfrac
    if exist(fname_results,'file')
      existvec(fraci) = true;    
      fname_results = sprintf('%s/job_%02d_%d.mat',dirname_results,fraci,nfrac);
      try
        tmp = load(fname_results);
      catch
        pause(2);
        tmp = load(fname_results);
      end

      fname_out = sprintf('%s/%s_%s_%02d_%02d.mat',dirname_cache,datatype,fstem_imaging,fraci,nfrac);
      try
        load(fname_out,'ivec_frac'); % Should perhaps check to touch file insetad, to make sure that worker is done writing?
      catch
      end

      zmat(:,ivec_frac) = tmp.zmat;
      logpmat(:,ivec_frac) = tmp.logpmat;
      beta_hat(:,ivec_frac) = tmp.beta_hat;
      beta_se(:,ivec_frac) = tmp.beta_se;
      sig2mat(:,ivec_frac) = tmp.sig2mat;
      sig2tvec(:,ivec_frac) = tmp.sig2tvec;
    end
  end

  % SAVE OUTPUT
  save_params = struct('fstem_imaging',fstem_imaging,'datatype',datatype,'outdir',dirname_out,'synth',synth);
  base_variables_to_save = {'X','iid','eid','colnames_model','contrasts','datatype','inputs','zmat','logpmat','beta_hat','beta_se','sig2mat','sig2tvec','save_params','mask'};

  if ~exist(dirname_out,'dir')
    mkdir(dirname_out); 
  end

  fpath_out = sprintf('%s/FEMA_wrapper_output_%s_%s.mat',dirname_out,datatype,fstem_imaging);

  %write column names to json for DEAP
  fname_col = sprintf('%s/FEMA_results_colnames.json',dirname_out);
  out = struct('colnames_model',{colnames_model},'RandomEffects',{RandomEffects});
  jsonStr = jsonencode(out);
  fid = fopen(fname_col,'w');
  fprintf(fid,'%s\n',jsonStr);
  fclose(fid);

  % =========================================================================
  % Write VOXEL results (mat, nifti, or deap)
  % =========================================================================
  if strcmpi(datatype,'voxel')
    vol_z = zeros([size(mask) size(zmat,1)]);
    vol_logp = zeros([size(mask) size(zmat,1)]);
    vol_beta_hat = zeros([size(mask) size(zmat,1)]);
    vol_beta_se = zeros([size(mask) size(zmat,1)]);
    for j = 1:size(zmat,1)
      vol_z(:,:,:,j) = single(fullvol(zmat(j,:),mask));
      vol_logp(:,:,:,j) = single(fullvol(logpmat(j,:),mask));
      vol_beta_hat(:,:,:,j) = single(fullvol(beta_hat(j,:),mask));
      vol_beta_se(:,:,:,j) = single(fullvol(beta_se(j,:),mask));
    end

    vol_sig2t = zeros([size(mask) 1]);
    vol_sig2t(ivec_mask) = single(sig2tvec);
    vol_sig2 = zeros([size(mask) size(sig2mat,1)]);
    for j = 1:size(sig2mat,1)
      vol_sig2(:,:,:,j) = single(fullvol(sig2mat(j,:),mask));
    end

    if contains(lower(outputFormat), 'nifti')
      results = struct('beta_hat',vol_beta_hat,'beta_se',vol_beta_se,'zmat',vol_z,'logpmat',vol_logp,'sig2tvec',vol_sig2t,'sig2mat',vol_sig2);
      writeNIFTI(results, dirname_out, fstem_imaging, ivnames, colnames_model(colsinterest)); % Should have option to write only subset of colnames_model annd / or stats volumes
    end
  else
    fprintf(1,'Error: datatype=%s not implemented\n',datatype);
  end


  % Need to add code to gather outputs
  % FEMA_wrapper_script_bottom;

  endtime = now();
  logging('***Done*** (%0.2f seconds)',(endtime-starttime)*3600*24);


  % ToDos
  %   Make sure colsinterest is passed from DEAP

