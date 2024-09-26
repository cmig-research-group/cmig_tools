function FEMA_DEAP_wrapper(fstem_imaging,fname_design,dirname_out,dirname_tabulated,dirname_imaging,datatype,dirname_cache,dirname_jobs,designid,nfrac,varargin)

  if isdeployed
    nfrac = str2num(nfrac);
    if isnan(nfrac)
      error('One or more input arguments could not be converted to a number.');
    end
  end

  if nargin < 9
    logging('Usage: FEMA_DEAP_wrapper(fstem_imaging,fname_design,dirname_out,dirname_tabulated,dirname_imaging,datatype,dirname_jobs,designid,nfrac)');
    error('Incorrect number of input arguments')
  end


  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% Look into how include the following code as "macro"

  inputs = inputParser;
  % addParamValue(inputs,'ranknorm',0);
  % addParamValue(inputs,'ico',5);
  % addParamValue(inputs,'pihat_file',[]);
  % addParamValue(inputs,'preg_file',[]);
  % addParamValue(inputs,'address_file',[]);
  % addParamValue(inputs,'nperms',0);
  % addParamValue(inputs,'mediation',0);
  % addParamValue(inputs,'tfce',0);

  %FEMA_fit variable inputs
  % addParamValue(inputs,'niter',1);
  % addParamValue(inputs,'nbins',20);
  % addParamValue(inputs,'CovType','analytic');
  % addParamValue(inputs,'FixedEstType','GLS');
  % addParamValue(inputs,'RandomEstType','MoM');
  % addParamValue(inputs,'GroupByFamType',true);
  % addParamValue(inputs,'NonnegFlag',true); % Perform lsqnonneg on random effects estimation
  % addParamValue(inputs,'SingleOrDouble','double');
  % addParamValue(inputs,'logLikflag',0);
  % addParamValue(inputs,'Hessflag',false);
  % addParamValue(inputs,'ciflag',false);
  % addParamValue(inputs,'permtype','wildbootstrap');

  % used variables
  addParamValue(inputs,'RandomEffects',{'F' 'S' 'E'}); % Default to Family, Subject, and eps
  addParamValue(inputs,'synth',0); % AMD - put back synth option
  addParamValue(inputs,'colsinterest',1);
  addParamValue(inputs,'contrasts',[]);
  addParamValue(inputs,'output','nifti'); % but we only support this format...
  addParamValue(inputs,'ivnames','');

  parse(inputs,varargin{:})
  % Display input arguments for log
  % niter = str2num_amd(inputs.Results.niter);
  % ico = str2num_amd(inputs.Results.ico);
  % ranknorm = str2num_amd(inputs.Results.ranknorm);
  % nbins = str2num_amd(inputs.Results.nbins);
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

  %fname_pihat = inputs.Results.pihat_file;
  %fname_address = inputs.Results.address_file;
  %fname_pregnancy = inputs.Results.preg_file;
  %CovType = inputs.Results.CovType;     
  %FixedEstType = inputs.Results.FixedEstType;
  %RandomEstType = inputs.Results.RandomEstType;
  %GroupByFamType = inputs.Results.GroupByFamType;
  %NonnegFlag = str2num_amd(inputs.Results.NonnegFlag);
  %SingleOrDouble = inputs.Results.SingleOrDouble;
  %OLSflag = ismember(lower(FixedEstType),{'ols'});
  %GLSflag = ismember(lower(FixedEstType),{'gls'});
  %logLikflag = str2num_amd(inputs.Results.logLikflag);
  %Hessflag = str2num_amd(inputs.Results.Hessflag);
  %ciflag = str2num_amd(inputs.Results.ciflag);
  %nperms = str2num_amd(inputs.Results.nperms);
  %permtype = inputs.Results.permtype;     
  %mediation = str2num_amd(inputs.Results.mediation);
  %tfce = str2num_amd(inputs.Results.tfce);

  RandomEffects = rowvec(split(strrep(strrep(inputs.Results.RandomEffects,'{',''),'}','')))
  synth = str2num_amd(inputs.Results.synth);
  colsinterest = str2num_amd(inputs.Results.colsinterest);

  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


  if ~ismember(lower(datatype),{'voxel' 'vertex' 'external' 'corrmat'})
    error('Input error: invalid datatype')
  end

  logging('***Start');

  starttime = now();

  rng shuffle % Set random number generator so different every time -- should allow for option to control this

  dirname_results = sprintf('%s/results/%s_%s_%s',dirname_jobs,designid,datatype,fstem_imaging);

  if exist(dirname_results,'dir')
    cmd = sprintf('rm -rf %s',dirname_results);
    [s r e] = jsystem(cmd);
  end

  % Create analysis job
  [colnames_model X] = FEMA_DEAP_wrapper_submit(fstem_imaging,fname_design,dirname_tabulated,dirname_imaging,'voxel',designid,dirname_jobs);
  % dispatch work to workers now that job files exists

  % Load design
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

  % zmat = NaN(size(X,2),);

  % function to process results as each worker completes
  function process_results(fraci)
    MAX_RETRIES = 10;
    fname_results = sprintf('%s/job_%02d_%02d.mat',dirname_results,fraci,nfrac);
    retry = 0;
    tmp = struct;
    while retry < MAX_RETRIES
      try
        tmp = load(fname_results, 'zmat','logpmat','beta_hat','beta_se','sig2mat','sig2tvec');
        retry = MAX_RETRIES + 1;
        % logging(fieldnames(tmp));
        fname_out = sprintf('%s/%s_%s_%02d_%02d.mat',dirname_cache,datatype,fstem_imaging,fraci,nfrac);
        cache = load(fname_out,'ivec_frac'); % Should perhaps check to touch file insetad, to make sure that worker is done writing?
        % logging(fieldnames(tmp));
        zmat(:,cache.ivec_frac) = tmp.zmat;
        logpmat(:,cache.ivec_frac) = tmp.logpmat;
        beta_hat(:,cache.ivec_frac) = tmp.beta_hat;
        beta_se(:,cache.ivec_frac) = tmp.beta_se;
        sig2mat(:,cache.ivec_frac) = tmp.sig2mat;
        sig2tvec(:,cache.ivec_frac) = tmp.sig2tvec;
      catch
        retry = retry + 1;
        pause(5);
        logging('Error: Could not load %s, retry %d',fname_results, retry);
      end
    end
    if retry == MAX_RETRIES
      logging('Error: Could not load %s, giving up',fname_results);
    end
  end

  % Create client to to dispatch work to workers
  client = FEMA_DEAP_client(designid, fstem_imaging, @process_results, nfrac);

  % block until all workers are done, will call the process_results function as each worker completes
  client.wait_for_completion();

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
end
