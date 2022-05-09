function [zmat_perm, beta_hat_perm, tfce_perm, colnames_interest, save_params, mask] = FEMA_cluster_gather(outdir,clusterjobs_outdir,nperms,varargin)
    %
    % Gathers outputs from FEMA_wrapper if it was run as multiple jobs (to maximise nperms)
    %           1) Gathers outputs and saves them as one file in `outdir`
    %           2) If `tfce` and `mostest_reg` args are passed, computes TFCE and MOSTest statistics
    %           NB IMPORTANT: relies on having only one .mat file in each of the job directories
    %
    % USAGE: ...
    %
    % INPUTS
    %   outdir <char>                 :  path to the output directory to save the gather data file
    %   clusterjobs_outdir  <char>    :  directory (string) containing all job ouput directories to be concatenated
    %   nperms                        :  total number of permutations
    %
    % Optional input arguments:
    %   mostest_reg <int> or 'auto'        :  regularization parameter for the MOSTest, see `FEMA_MOSTest.m` for more info
    %                                               If not specified, MOSTest will not run.
    %                                               If set to 'auto', the best parameter will be run (not yet implemented)
    %                                                     Otherwise set to the chosen custom integer
    %   calc_perm_pvals <bool>, default 1  :  whether to compute rank and extrapolated p-values from the permuted empirical distribution
    %                                          uses `FEMA_perm_significance.m`
    %
    % OUTPUTS
    %   zmat_perm                     :  permutation matrix of z-scores
    %   beta_hat_perm                 :  permutation matrix of beta coefficients
    %   tfce_perm                     :  if TFCE has been run, will concatenate its values
    %   save_params <struct>          :  structure output from FEMA_wrapper with parameters to save file with FEMA naming and directory conventions
    %                                           N.B. - outputs will only be saved if save_params provided as an input
    %
    
    %
    % This software is Copyright (c) 2021 The Regents of the University of California. All Rights Reserved.
    % See LICENSE
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
          significance_args = inputParser;
          addParamValue(significance_args,'mostest_reg',[]);
          addParamValue(significance_args,'calc_perm_pvals',1);
          addParamValue(significance_args,'fstem_pheno','');
    
          parse(significance_args,varargin{:})
          significance_args = significance_args.Results;
          fstem_pheno = significance_args.Results.fstem_pheno;
          tic
          dlist = dir(sprintf('%s/*',clusterjobs_outdir));
          dlist = dlist([dlist.isdir]); % only dirs, not files
          dlist = dlist(3:end); % skip . and ..
          jobdirs = strcat(clusterjobs_outdir,'/',{dlist.name});
    
          njobs = length(jobdirs);
          starttime = now();
          logging('Gathering FEMA outputs from %d directories...',njobs);
    
          digits_njobs = max(ceil(log10(njobs+1)),1);
          loop_timer_start = now();
          for job=1:njobs
                estimated_time_remaining = (now()-loop_timer_start)*3600*24/job * (njobs - job);
                logging('Reading job %0*d of %d (%.2fs - remaining %.0fs)',digits_njobs,job,njobs,toc,estimated_time_remaining);
                tic
                % Find .mat file in job subdirectory, exits if more than one
                dirlist = dir(sprintf('%s/*%dperms',jobdirs{job},nperms));
                dirlist = {dirlist.name};
                if length(dirlist) ~= 1
                      error(sprintf('Found no matching directory with nperms %d in the job subdirectories.',nperms));
                end 
                dirname = sprintf('%s/%s',jobdirs{job},dirlist{1}); 
                mat_files = dir(sprintf('%s/*%s*.mat',dirname,fstem_pheno));
                job_to_load = sprintf('%s/%s', mat_files(1).folder, mat_files(1).name);
                if length(mat_files) ~= 1
                      error(sprintf(['Found several `.mat` files in the job subdirectories. ', ...
                            'Please erase/move the ones that are not relevant.\n', ...
                            '  Job directory: %s\n', ...
                            '  Found %d results: %s\n', ...
                            'Exiting now.'], clusterjobs_outdir, length(mat_files), job_to_load))
                      return
                end
                % Load job results
                if job==1
                      load(job_to_load); % This seems really slow (8GB/66sec ~= 120MB/sec) -- also use lower iconum
                      assert(size(zmat_perm,1)==length(colsinterest),'First dimension of `zmat_perm` must correspond to `length(colsinterest)`')
                      files_to_load = {'beta_hat_perm', 'zmat_perm'}; %speed up loading time for job>1
                      % Create output file name
                      [fdir, fname] = fileparts(job_to_load);
                      [~, perm_info] = fileparts(fdir);
                      outdir_perms = fullfile(outdir, perm_info);
                      fpath_out = fullfile(outdir_perms, strcat(fname, ".mat"));
                      
                      % Initialize matrices and add unpermuted element (1st column in 3rd dimension)
                      zmat_perm_cat=NaN(length(colsinterest),size(zmat_perm,2),nperms+1,'single');
                      zmat_perm_cat(:,:,1) = zmat_perm(:,:,1); 
                      beta_hat_perm_cat=NaN(length(colsinterest),size(zmat_perm,2),nperms+1,'single');
                      beta_hat_perm_cat(:,:,1) = beta_hat_perm(:,:,1); 
                      if exist('tfce_perm','var')
                            tfce_perm_cat=NaN(length(colsinterest),size(zmat_perm,2),nperms+1,'single');
                            tfce_perm_cat(:,:,1) = tfce_perm(:,:,1);
                            files_to_load = [files_to_load, 'tfce_perm'];
                      end
                else
                      load(job_to_load, files_to_load{:});
                end
                % Copy over the permutations to gathering variables
                nperms_per_job = size(zmat_perm,3)-1;
                zmat_perm_cat(:,:,1+(job-1)*nperms_per_job+[1:nperms_per_job]) = zmat_perm(:,:,2:end);
                beta_hat_perm_cat(:,:,1+(job-1)*nperms_per_job+[1:nperms_per_job]) = beta_hat_perm(:,:,2:end);
                if exist('tfce_perm','var')
                      tfce_perm_cat(:,:,1+(job-1)*nperms_per_job+[1:nperms_per_job]) = tfce_perm(:,:,2:end);
                end
                
                % PrintMemoryUsage
          end
          logging('Gathering completed in %.0fs...',(now()-starttime)*3600*24);
          zmat_perm = zmat_perm_cat;
          beta_hat_perm = beta_hat_perm_cat;
          if exist('tfce_perm_cat','var')
                tfce_perm = tfce_perm_cat;
          end
    
          % get `fstem_imaging` and `datatype` from file name
          tmp = split(fname, 'output_');
          tmp = split(tmp(2), '_');
          datatype = tmp{1};
          tmp = join(tmp(2:end), '_');
          tmp = split(tmp, '_synth');
          fstem_imaging = tmp{1};
          synth = contains(fname, 'synth');
          save_params = struct('fstem_imaging',fstem_imaging,'datatype',datatype,'outdir',outdir_perms,'synth',synth);
    
          tic
          variables_to_save = {'X','contrasts','inputs','datatype','colnames_model','colnames_interest',...
                'colsinterest','zmat','beta_hat','beta_se','sig2mat','sig2tvec','zmat_perm','beta_hat_perm','save_params','mask'};
          if ~exist(outdir_perms), mkdir(outdir_perms); end
          if inputs.Results.tfce==0
                if strcmpi(datatype,'voxel')
                      save(fpath_out,variables_to_save{:},'vol_z','vol_beta_hat','-v7.3')
                elseif strcmpi(datatype, 'vertex') | strcmpi(datatype, 'external') | strcmpi(datatype, 'corrmat')
                      save(fpath_out,variables_to_save{:},'logpmat','-v7.3')
                end
          elseif inputs.Results.tfce==1
                if strcmpi(datatype,'voxel')
                      save(fpath_out,variables_to_save{:},'vol_z','vol_beta_hat','tfce_perm','-v7.3')
                elseif strcmpi(datatype, 'vertex') | strcmpi(datatype, 'external') | strcmpi(datatype, 'corrmat')
                      save(fpath_out,variables_to_save{:},'logpmat','tfce_perm','-v7.3')
                end
          end
          logging('Concatenated file saved to %s (%.1fs)',fpath_out,toc)
          
          if significance_args.calc_perm_pvals==1
                      stattype='z';
                      logging('Start significance computation for %s', stattype)
                      [fpath_out,rank_log10pval_uncorr, rank_log10pval_fwecorr]=FEMA_perm_significance(zmat_perm, stattype, colnames_interest,'save_params',save_params,'extrapolate',0,'mask',mask);
                      if inputs.Results.tfce==1
                            stattype='tfce';
                            logging('Start significance computation for %s', stattype)
                            [fpath_out,rank_log10pval_uncorr, rank_log10pval_fwecorr]=FEMA_perm_significance(tfce_perm, stattype, colnames_interest,'save_params',save_params,'extrapolate',0,'mask',mask);
                      end
          end
          if ~isempty(significance_args.mostest_reg)
                [fpath_out,rank_log10pval_mostest, extrap_log10pval_mostest] = FEMA_MOSTest(zmat_perm, colnames_interest,'save_params',save_params,'k', significance_args.mostest_reg);
          end
    
          logging('Completed');
    
    end
    