function [colnames_model X] = FEMA_DEAP_wrapper_submit(fstem_imaging,fname_design,dirname_tabulated,dirname_imaging,datatype,designid,dirname_jobs)
  if nargin < 5
      fprintf(1,'Usage: FEMA_DEAP_wrapper(fstem_imaging,fname_design,dirname_tabulated,dirname_imaging,datatype)\n',mfilename);
      error('Incorrect number of imput arguments')
  end

  starttime = now();
  logging('\n***Start*** %s %s\n',mfilename,datestr(starttime));

  if ~ismember(lower(datatype),{'voxel'})
    error('Input error: invalid datatype')
  end

  % LOAD AND PROCESS IMAGING DATA FOR ANALYSIS - ABCD specific function unless datatype='external'
  [ymat, iid_concat, eid_concat, ivec_mask, mask, colnames_imaging, pihat, preg, address] = FEMA_process_data(sprintf('-%s',fstem_imaging),dirname_tabulated,dirname_imaging,datatype); % Look for way to speed this up (currently takes ~2.5 sec)

  logging('Checkpoint 1 (%0.2f seconds)',(now-starttime)*3600*24); % This typically takes ~2.76 seconds

  [X,iid,eid,fid,agevec,ymat,contrasts,colnames_model,pihatmat,PregID,HomeID] = FEMA_intersect_design(fname_design, ymat, iid_concat, eid_concat); % Make version of this that accepts empty ymat S

  logging('Checkpoint 2 (%0.2f seconds)',(now-starttime)*3600*24);

  dirname_design = sprintf('%s/design',dirname_jobs);
  if ~exist(dirname_design,'dir')
    jsystem(sprintf("mkdir -p %s", dirname_design));
  end
  fname_job_design = sprintf('%s/design/%s.mat',dirname_jobs,designid); 

  save(fname_job_design,'X','iid','eid','fid','agevec');

  logging('Checkpoint 3 (%0.2f seconds)',(now-starttime)*3600*24);

  %   incorporate computation and saving of iid, eid, X_resid from SSE_DEAP_worker.m

