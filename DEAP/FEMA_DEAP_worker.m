function FEMA_DEAP_worker(fstem_imaging,datatype,fraci,nfrac,dirname_cache,dirname_jobs)

  if isdeployed
    fraci = str2num(fraci);
    nfrac = str2num(nfrac);
    if isnan(fraci) || isnan(nfrac)
      error('One or more input arguments could not be converted to a number.');
    end
  end

  % nested function makes this code static and therefore I need to explicetly access the variables from the loaded mat file
  fname_cache = sprintf('%s/%s_%s.mat',dirname_cache,datatype,fstem_imaging);
  fprintf("Loading cache %s\n", fname_cache);
  tic 
  cache = load(fname_cache);
  iid = cache.iid;
  eid = cache.eid;
  fid = cache.fid;
  agevec = cache.agevec;
  pihatmat = cache.pihatmat;
  RandomEffects = cache.RandomEffects;
  PregID = cache.PregID;
  HomeID = cache.HomeID;
  FamilyStruct = cache.FamilyStruct;
  toc

  fname_frac_cache = sprintf('%s/%s_%s_%02d_%02d.mat',dirname_cache,datatype,fstem_imaging,fraci,nfrac);
  fprintf("Loading frac cache %s\n", fname_frac_cache);
  tic % Takes ~1.4sec for nfrac = 20
  frac_cache = load(fname_frac_cache);
  ymat_frac = frac_cache.ymat_frac;
  ivec_frac = frac_cache.ivec_frac;
  toc

  colsinterest = [1]; % Should allow this to be passed in

  % Check for pending analyses

  ymat_bak = ymat_frac;

  function run_job(job_name, job_fstem_imaging)
    if job_fstem_imaging ~= fstem_imaging
      error('job_fstem_imaging does not match fstem_imaging');
    end
    starttime = now();
    fprintf(1,'\n***Start*** %s %s\n',mfilename,datestr(starttime));
    job_details = sprintf('%s_%s_%s',job_name,datatype,fstem_imaging);
    part_name = sprintf('job_%02d_%02d',fraci,nfrac);
    fname_job_design = sprintf('%s/design/%s.mat',dirname_jobs,job_name); 
    dirname_results = sprintf('%s/results/%s',dirname_jobs,job_details); 
    tic
    taskinfo = load(fname_job_design); % This takes most of the time -- should use fread instead?
    toc

    idvec = strcat(iid,'_',eid);
    idvec2 = strcat(taskinfo.iid,'_',taskinfo.eid);
    
    ymat_frac = ymat_bak;
    if ~isequal(idvec,idvec2)
      [dummy IA IB] = intersect(idvec,idvec2,'stable');
      X = taskinfo.X(IB,:);
      ymat = ymat_frac(IA,:);
      % Should re-compute FamilyStruct, etc.
    else
      X = taskinfo.X;
      ymat = ymat_frac;
      IA = [1:length(idvec)];
      IB = [1:length(idvec)];
    end
    
    jvec_covars = setdiff([1:size(X,2)],colsinterest);
    X_tmp = X - X(:,jvec_covars)*pinv(X(:,jvec_covars))*X;
    X_reasid = cat(2,X_tmp(:,jvec_covars),X(:,end));
    ymat = ymat - X(:,jvec_covars)*pinv(X(:,jvec_covars))*ymat;
    % BP-MOD: [beta_hat beta_se zmat logpmat sig2tvec sig2mat] = FEMA_DEAP_fit(X,iid,eid(IA),fid(IA),agevec(IA),ymat,[],[],[],[],'FamilyStruct',FamilyStruct); % Work with Pravesh to update to use integreated FEMA_DEAP_fit -- need to re-compute FamilyStruct if iid and eid change
    [beta_hat beta_se zmat logpmat sig2tvec sig2mat] = FEMA_fit(X,iid,eid(IA),fid(IA),agevec(IA),ymat,[],[],[],[],'FamilyStruct',FamilyStruct); % Work with Pravesh to update to use integreated FEMA_DEAP_fit -- need to re-compute FamilyStruct if iid and eid change

    if ~exist(dirname_results,'dir')
      cmd = sprintf('mkdir -p %s',dirname_results);
      [s r e] = jsystem(cmd);
    end

    fname_results = sprintf('%s/%s.mat',dirname_results,part_name);
    tic
    save(fname_results,'beta_hat','beta_se','zmat','logpmat','sig2tvec','sig2mat');
    toc
    endtime = now();
    fprintf(1,'***Done*** (%0.2f seconds)\n',(endtime-starttime)*3600*24);
  end

  FEMA_DEAP_server(fraci, @run_job)
end
% ToDo
%   Should check if analysis already run -- generate warning
%   Limit to colsofinterest -- should save this in taskinfo and/or design
%   Perform residualization of X in FEMA_DEAP_wrapper.m, avoid re-doing residualization if isequal(idvec,idvec2)
