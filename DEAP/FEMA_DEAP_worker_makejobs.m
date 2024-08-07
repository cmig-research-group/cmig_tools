function FEMA_DEAP_worker(fstemlist,datatype,nfrac,dirname_cache,dirname_jobs,batchdir)

if ~exist(batchdir)
  mkdir(batchdir)
end
jsystem(sprintf('rm -rf %s/job*.m %s/pbsout',batchdir,batchdir));

fname_list = sprintf('%s/scriptlist.txt',batchdir);
fid_list = fopen(fname_list,'wb');
jobnum = 0;
for fsi = 1:length(fstemlist)
  fstem_imaging = fstemlist{fsi};
  for fraci = 1:nfrac
    jobnum = jobnum+1; jobname = sprintf('job_%05d',jobnum);
    fname = sprintf('%s/%s.m',batchdir,jobname);
    fid = fopen(fname,'wb');
    fprintf(fid,'FEMA_DEAP_worker(''%s'',''%s'',%d,%d,''%s'',''%s'');\n',fstem_imaging,datatype,fraci,nfrac,dirname_cache,dirname_jobs);
    fprintf(fid,'exit\n');
    fclose(fid);
    fprintf(fid_list,'%s\n',jobname);
  end
end
fclose(fid_list);

% FEMA_DEAP_worker(fstem_imaging,'voxel',fraci,nfrac,dirname_cache,dirname_jobs);

