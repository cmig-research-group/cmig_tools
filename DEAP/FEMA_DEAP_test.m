% ToDos
%   Write code to break ymat into multiple segments -- insert into FEMA_wrapper
%   Write code to save model / analysis meta data to file -- insert into FEMA_wrapper
%   Write worker process code to
%     startup: read in specific segment of data (data type & segment # of total segments)
%     loop: wait for next analyis, read in saved model meta data, save results for segment
%   Write worker process to assemble processing results

addpath(genpath('/home/dale/matlab/FEMA_local'));

dirname_tabulated = '/space/amdale/1/tmp/ABCD_cache/abcd-sync/4.0/tabulated/released';
dirname_voxelwise = '/space/amdale/1/tmp/ABCD_cache/abcd-sync/4.0/imaging_concat/voxelwise/ABCD2_cor10';
dirname_vertexwise = '/space/amdale/1/tmp/ABCD_cache/abcd-sync/4.0/imaging_concat/vertexwise';
dirname_corrmat = '/space/amdale/1/tmp/ABCD_cache/abcd-sync/4.0/imaging_concat/corrmat';


%  Cache voxelwise phenotypes

nfrac = 20;
dirname_cache = '/home/dale/FEMA_cache';
dirname_jobs = '/space/amdale/1/tmp/jobdir';

measlist = {{'dmri' 'RNI'} {'dmri' 'RND'} {'smri' 'JA'}}; % List of phenotypes to cache
for measi = 1:length(measlist)
  subdir  = measlist{measi}{2};
  voxelwise_pheno = measlist{measi}{2};
  fstem_pheno_voxel = sprintf('%s',voxelwise_pheno);
  dirname_out = '/home/dale/FEMA_output/test'; jsystem(sprintf('rm -rf %s',dirname_out)); mkdir(dirname_out); % This shouldn't be needed by FEMA_DEAP_init
  fname_design = '/space/amdale/1/tmp/ABCD_cache/abcd-sync/4.0/support_files/design_matrices/ABCD_rel4.0_cross_desmat_PCs_SES_motion_interview_age.txt'; colsinterest = [1]; % Age design from Clare -- should use only standard covariates
  FEMA_DEAP_init(fstem_pheno_voxel,fname_design,dirname_out,dirname_tabulated,dirname_voxelwise_tmp,'voxel'); % Should add dirname_cache and nfrac as arguments, get rid of dirname_out, dirname_tabulated? -- should write out imaging phenotypes only
end


% Queue up specific analysis
fstem_imaging = 'JA'; 
subdir = 'smri';
fname_design = '/space/amdale/1/tmp/ABCD_cache/abcd-sync/4.0/support_files/design_matrices/ABCD_rel4.0_cross_desmat_PCs_SES_motion_interview_age.txt'; 
colsinterest = [1]; 
designid = 'des00001'; % Age design from Clare -- should be uip[dated to Release 5.0
dirname_out = '/home/dale/FEMA_output/test'; 
jsystem(sprintf('rm -rf %s',dirname_out)); 
mkdir(dirname_out); % This shouldn't be needed by FEMA_DEAP_init

dirname_imaging = sprintf('%s/%s',dirname_voxelwise,subdir);
FEMA_DEAP_wrapper(fstem_imaging,fname_design,dirname_out,dirname_tabulated,dirname_imaging,'voxel',dirname_cache,dirname_jobs,designid,nfrac);




% Test individual worker (these should be running as concurrent parallel jobs)
batchdir = '/home/dale/batchdirs/DEAP_workers';
fstemlist = {'JA'};
FEMA_DEAP_worker_makejobs(fstemlist,'voxel',nfrac,dirname_cache,dirname_jobs,batchdir);

fstem_imaging = 'JA';

fraci = 5; 
FEMA_DEAP_worker(fstem_imaging,'voxel',designid,fraci,nfrac,dirname_cache,dirname_jobs);

% Write code to gather results from workers


%tmp_volinfo = load('/space/syn65/1/data/abcd-sync/5.0/imaging_concat/voxelwise/ABCD3_cor10/smri/volinfo.mat');
tbl_tabulated = readtable('/space/syn65/1/data/abcd-sync/5.0/tabulated/released/abcd_mri01.txt');

% ToDo
%   Update FEMA_worker to read from dirname_jobs





% Vertexwise
%icnum = 6; % 10242 vertices per hemisphere -- should match with what DEAP is using
%icnum = 5; % 2562 vertices per hemisphere -- should match with what DEAP is using
icnum = 4; % 642 vertices per hemisphere -- should match with what DEAP is using -- this is not currently supported by FEMA_wrapper for other than tfmri
sm = 16; % Should be user settable

vertexwise_pheno = 'thickness'; subdir = 'smri';
%vertexwise_pheno = 'area'; subdir = 'smri';
%vertexwise_pheno = 'sulc'; subdir = 'smri';
%vertexwise_pheno = 'RNI-gm'; subdir = 'dmri';
%vertexwise_pheno = 'RND-gm'; subdir = 'dmri';
%vertexwise_pheno = 'FA-gm'; subdir = 'dmri';
%vertexwise_pheno = 'nBack_1_2_back_vs_0_back'; subdir = 'tfmri';

if strcmp(subdir,'tfmri')
  fstem_pheno_vertex = sprintf('%s',vertexwise_pheno);
else
%  fstem_pheno_vertex = sprintf('%s_ic%d_sm%d',vertexwise_pheno,icnum,sm);
  fstem_pheno_vertex = sprintf('%s_ic%d_sm%d',vertexwise_pheno,5,sm);
end
dirname_vertexwise_tmp = sprintf('%s/%s',dirname_vertexwise,subdir);
[fnames_out_vertex] = FEMA_wrapper(fstem_pheno_vertex,fname_design,dirname_out,dirname_tabulated,dirname_vertexwise_tmp,'vertex','ico',icnum-1,'output','mat','niter',niter,'pihat_file',fname_pihat,'preg_file',fname_preg,'address_file',fname_address,'nperms',nperms,'RandomEffects',RandomEffects,'RandomEstType',RandomEstType,'colsinterest',colsinterest);





% Old junk

addpath(genpath('/home/dale/matlab/SSE_local'));

% Test script for calling SSE_wrapper from matlab and command line

fname_design = '/space/amdale/1/tmp/ABCD_cache/SSE_examples/ABCD_3.0_design_example1.txt'; design_id = 'des00001';
fname_contrasts = '/space/amdale/1/tmp/ABCD_cache/SSE_examples/contrast_example1.txt';
dirname_tabulated = '/space/amdale/1/tmp/ABCD_cache/abcd-sync/3.0/tabulated/released';
dirname_jobs = '/space/amdale/1/tmp/jobdir';
dirname_out = '/space/amdale/1/tmp/SSE_results/design_id';
dirname_vertexwise = '/space/amdale/1/tmp/ABCD_cache/abcd-sync/3.0/vertexwise';
dirname_voxelwise = '/space/amdale/1/tmp/ABCD_cache/abcd-sync/3.0/voxelwise/ABCD1_cor10/volmats_subsamp';
winsor_std = 0;
preresid = 0;
ranknorm = 0;


% Test vertexwise
fstem_vertex_pheno = 'thickness-sm256';
%fstem_vertex_pheno = 'area-sm256';
%fstem_vertex_pheno = 'sulc-sm256';
%fstem_vertex_pheno = 'N0-gwc_ic5_sm256';
%fstem_vertex_pheno = 'ND-gwc_ic5_sm256';
icnum = 6; % 10242 vertices per hemisphere -- should match with what DEAP is using

SSE_DEAP_worker(fstem_vertex_pheno,'vertex',dirname_vertexwise,dirname_tabulated,dirname_jobs,'ico',icnum-1); % Should be launched as a persistent process / demon, and restarted as needed

SSE_DEAP_wrapper(fstem_vertex_pheno,fname_design,dirname_out,dirname_tabulated,dirname_vertexwise,'vertex',dirname_jobs,design_id,'winsor_std',winsor_std,'preresid',preresid,'ranknorm',ranknorm,'output','matnifti'); % Return X, ymat, iid, fid; Get rid of histogram plotting


% Test voxelwise 
fstem_voxel_pheno = 'JA';
%fstem_voxel_pheno = 'ND';
%fstem_voxel_pheno = 'N0';

SSE_DEAP_worker(fstem_voxel_pheno,'voxel',dirname_voxelwise,dirname_tabulated,dirname_jobs);

SSE_DEAP_wrapper(fstem_voxel_pheno,fname_design,dirname_out,dirname_tabulated,dirname_voxelwise,'voxel',dirname_jobs,design_id,'winsor_std',winsor_std,'preresid',preresid,'ranknorm',ranknorm,'output','mat','niter',1,'nbins',10);

% Test corrmat

