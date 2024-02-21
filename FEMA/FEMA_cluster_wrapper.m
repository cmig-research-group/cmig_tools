function FEMA_cluster_wrapper(fstem_imaging,fname_design,dirname_out,dirname_tabulated,dirname_imaging,datatype,batchdir,njobs,varargin)

      if ~exist(batchdir,'dir')
            mkdir(batchdir)
      end
      jsystem(sprintf('rm -rf %s/job*.m %s/pbsout',batchdir,batchdir));
      fname = sprintf('%s/scriptlist.txt',batchdir);
      fid_list = fopen(fname,'wb');
      nperms_per_job=varargin{find(strcmpi(varargin,'nperms'))+1};
      dirname_out = fullfile(dirname_out, sprintf('dt-%s_img-%s_njobs-%d_nperms-%d',datatype,fstem_imaging,njobs,nperms_per_job));
      for outnum=1:length(dirname_out)
            fprintf('\nCluster output for design matrix %d will be written to: %s ',outnum,dirname_out{outnum});
      end
      jobnum = 0;
      for diri = 1:njobs
            jobnum = jobnum+1; jobname = sprintf('job_%05d',jobnum);
                  if iscell(dirname_out)
                        dirname_out_job = cell(1,length(dirname_out));
                        for fi = 1:length(dirname_out)
                              dirname_out_job{fi} = sprintf('%s/%s',dirname_out{fi},jobname);
                        end
                  else
                        dirname_out_job = sprintf('%s/%s',dirname_out,jobname);
                  end
            fname = sprintf('%s/%s.m',batchdir,jobname);
            fid = fopen(fname,'wb');
            fprintf(fid,'FEMA_wrapper(%s,%s,%s,%s,%s,%s',inputForm(fstem_imaging),inputForm(fname_design),inputForm(dirname_out_job),inputForm(dirname_tabulated),inputForm(dirname_imaging),inputForm(datatype));,
            
            for j = 1:length(varargin);
                  tmp = varargin{j};
                  %    disp(tmp);
                  tmp = inputForm(tmp);
                  %    disp(tmp);
                  fprintf(fid,',%s',tmp);
            end

            fprintf(fid,');\n');
            fprintf(fid,'exit\n');
            fclose(fid);
            fprintf(fid_list,'%s\n',jobname);
      end
      
      fclose(fid_list);

end


