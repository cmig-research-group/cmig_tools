%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% FEMA RESAMPLING DEMO USING FEMA_wrapper.m
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% This script is designed to provide an end-to-end demo of using the resampling methods and
% associated functions within the FEMA package primarily using `FEMA_wrapper.m`.

% This demo will only work for ABCD INVESTIGATORS with access to abcd-sync.
% You need to have abcd-sync mirrored onto your server in the exact same directory structure in order
% for FEMA_wrapper.m to find the correct files for the demo analysis.

% We are in the process of creating demo scripts with dummy data for non-ABCD Investigators.
% FEMA_fit.m the core function of FEMA is NOT ABCD specific and can be run with any data in the 
% FEMA specified format.

% When using FEMA_wrapper.m, vertexwise, voxelwise and corrmat analyses using FEMA_wrapper.m require abcd-sync data,
% therefore this only works for ABCD Investigators. However, anyone can run FEMA_wrapper.m using
% external data of any kind by specifying 'datatype' as external.

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% Resampling procedure %%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% FEMA uses a residual, wild bootstrap (WB) resampling procedure that can be used to generate a null
% distribution in order to calculate standard errors and confidence intervals that are not dependent on
% the shape of the empirical distribution (Wu, 1986). This is useful if the dependent variable is not
% normally distributed or if the model is heteroskedastic. FEMA generates a synthetic y (dependent variable) 
% by multiplying residuals by a random ‘wild’ weight (v) and computing y like so: Yi = Yhat_i + EiVi

% The type of resampling method can be set using `PermType`:
%     PermType = 'wildbootstrap'  (null)
%     PermType = 'wildbootstrap-nn' (non-null)

% Resampling statistics is needed to compute standard errors and p-values when the imaging data is
% not normally distributed, to run threshold free cluster enhancement (TFCE) and to run Sobel's test
% for mediation analysis.

% TFCE and Sobel's test are only supported for vertexwise and voxelwise analysis.

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% Single computer vs cluster? %%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% For demo purposes this script will only run 10 permutations on a single computer.
% For valid statistics at the appropriate alpha level for your analyses >1000 permutations should
% be run. For computational efficiency we advise that high numbers of permutations are run on a
% cluster. The FEMA package is designed such that `FEMA_wrapper.m` can be run on several instances
% in parallel. The outputs can then be gathered into a single file using `FEMA_cluster_gather.m`.
% The associated functions that need to be run using the full set of permutations can then be applied
% either to the output from a single iteration of `FEMA_wrapper.m` or to the output from
% `FEMA_cluster_gather.m`

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% This script uses example design matrices within the `support_files`
% directory on `abcd-sync`, and concatenated imaging data within the
% `imaging_concat` directory on `abcd-sync`.

% After running this demo and saving the output from FEMA, the results can
% be visualised using our visualisation tools: `showSurf` and `showVol`.
% We have included demos for these tools, which will run off the results from
% this demo. There is one for vertexwise data (`FEMA_showSurf_demo.m`) and
% one for voxelwise data (`FEMA_showVol_demo.m`).

% Code written by Clare E Palmer, Pierre Nedelec, John Iversen, Diliana Pecheva & Anders Dale (2022)

% This software is Copyright (c) 2021 The Regents of the University of California. All Rights Reserved.
% See LICENSE.

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% DEMO SCRIPT WITH TUTORIAL

% Run `FEMA_resampling_demo.m` in Matlab as a script for a full end-to-end demo.
% The demo will prompt you to give user specific inputs in the command window.

% If you would like to edit the script manually to change some settings then SAVE A LOCAL
% COPY of this script before running.

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% REQUIREMENTS TO RUN FEMA_wrapper.m
%
% CLONE GitHub directory for cmig_tools
% ADD to MATLAB path:

% e.g. if cloned into ~/github: 
% addpath(genpath('~/github/cmig_tools'))

% If wanting to run threshold free cluster enhancement (TFCE) you need to additionally download PALM and add directory to path: https://fsl.fmrib.ox.ac.uk/fsl/fslwiki/PALM
% addpath(genpath('~/PALM'))

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% 1) Add cmig_tools and PALM directories to your path

git_dir=input('Please specify the path to your cmig_tools directory. If this is already in your path leave blank and press ENTER: ','s');

if ~isempty(git_dir)
fprintf('\nAdding cmig_tools directory to MATLAB path... \n')
addpath(genpath(git_dir))
end

fprintf('If wanting to run TFCE, you will need to download the PALM toolbox (https://fsl.fmrib.ox.ac.uk/fsl/fslwiki/PALM) \n')
palm_dir=input('Please specify the path to PALM toolbox. If you are not planning to run TFCE leave blank and press ENTER: ','s');
if ~isempty(palm_dir)
      fprintf('\nAdding PALM directory to MATLAB path... \n')
      addpath(genpath(palm_dir))
end

% 1) Get path to local copy of abcd-sync.  This can be done by specifying
% the path manually or by using abcdConfig to specify paths.

% `abcdConfig` is a function that users can run once from the MATLAB command
% line during set up which will prompt the user to select directory paths for
% those required to run FEMA.  This will then store the paths in a JSON
% file in your home directory for use everytime you want to run FEMA to save you
% inputting the path to abcd-sync everytime.
%%
abcd_config=input('\nDo you want to use abcdConfig to provide the path to abcd-sync? (enter numeric value: 1 or 0) '); %Set to 1 to use abcdConfig to specify path to abcd-sync

if abcd_config==1
      cfg = abcdConfig('FEMA');
      abcd_sync_path=cfg.data.abcd_sync;
elseif abcd_config==0
      abcd_sync_path=input('Please input path to local abcd-sync mirror: ', 's');
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% 2) Specify data release

dataRelease = '4.0'; %'3.0' or '4.0'

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% 3) Specify where to store results

outdir_root = fullfile('~/FEMA_demo_results',dataRelease);
fprintf('\n\nStoring results here: %s \n', outdir_root)
outcheck=input('Would you like to specify a DIFFERENT output directory? (enter numeric value: 1 or 0): ');
if outcheck==1
      outdir_root=input('Please specify output directory: ','s');
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% 4) Specify which imaging analyses to demo

% FEMA can be run using four datatypes, which determine the data used for Y
% in the linear mixed effects analysis:

%   datatype = 'vertex' --> loads concatenated .mat file of vertexwise
%   imaging data (nobs x nvertices) in specified modality (e.g. cortical thickness, surface area)
%   and resolution i.e. number of vertices (e.g. ico = 5)

%   datatype = 'voxel' --> loads concatenated .mat file of voxelwise
%   imaging data (nobs x nvoxels) in specified modality (e.g. FA, MD, RNI)

%   datatype = 'corrmat' --> loads concatenated .mat file of a correlation
%   matrix derived from the imaging data (XxXxX???) in specified modality
%   (e.g. task fMRI, resting state)

%   datatype = 'external' --> loads .txt file including any
%   data, imaging (e.g. ROI data) or behavioral data. These data will be
%   the dependent variables (Y) for each model. FEMA will run mass
%   univariate mixed effects associations across every column of this .txt
%   file

% FEMA also has functionality to generate null synthesized data in the same size
% as the datatype specified by inputting synth=1 as input to FEMA_wrapper.m

% For this demo, set the below variables to 1 for all datatypes you'd
% like to demo.  The output from FEMA for each datatype is saved with a
% unique output filename within the output directory specified above.
%%
fprintf('\n\nFor demo purposes nperms=10. For valid statistics run with much higher nperms and on the cluster to increase efficiency.\n')

fprintf('\nPlease specify which imaging data to demo: \n')

doVertexwise = input('Run vertexwise demo? (enter numeric value: 1 or 0) '); % run vertexwise analysis (datatype = 'vertex')
doTFCEvertex = input('Run TFCE on the FEMA output? (enter numeric value: 1 or 0) ');
doSobelvertex = input('Run mediation analysis vertexwise? (enter numeric value: 1 or 0) ');

doVoxelwise = input('\nRun voxelwise demo? (enter numeric value: 1 or 0) '); % run voxelwise analysis (datatype = 'voxel')
doTFCEvoxel = input('Run TFCE on the FEMA output? (enter numeric value: 1 or 0) ');
doSobelvoxel = input('Run mediation analysis voxelwise? (enter numeric value: 1 or 0) ');

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Check if PALM toolbox in path for TFCE analysis
if doTFCEvertex==1 | doTFCEvoxel==1
      if ~exist('palm_calcarea')
      error('User has selected to demo TFCE, but PALM toolbox is not in MATLAB path. Re-run specifying PALM directory, or select 0 to demo TFCE analyses.')
      end
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% 5) If editting FEMA_resampling_demo.m SAVE LOCAL COPY

% 6) RUN SCRIPT FOR DEMO

% For the demo, the code below this line does not need to be edited.

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% See below for detailed annotations about the inputs to FEMA_wrapper.m for
% each datatype used for this demo and examples of how to run FEMA_wrapper.m

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% INPUTS TO FEMA_wrapper.m

dirname_tabulated = fullfile(abcd_sync_path,dataRelease,'tabulated/released'); % directory to tabulated imaging data on abcd-sync 

switch dataRelease
  case '4.0'
    atlasVersion = 'ABCD2_cor10';
    dirname_tabulated = fullfile(abcd_sync_path,'4.0','tabulated/released'); 
    fname_design = fullfile(abcd_sync_path,'4.0','support_files/design_matrices','FEMA_demo_rel4.0_designmatrix_age.txt');
    fname_pihat = fullfile(abcd_sync_path, '4.0','genomics','ABCD_rel4.0_pihat.mat'); 
end

[root,filename]=fileparts(fname_design);
dirname_out=strcat(outdir_root,'/',filename);

% NOTE: `fname_design` can also be a cell array of multiple design matrices
% (txt file paths) that will be looped over in `FEMA_wrapper.m`
% Use makeDesign.R from ~/github/cmig_utils/r to make design matrix compatible with `FEMA_wrapper.m`
% See makeDesign_demo.R for tutorial on how to make design matrices
% compatible with `FEMA_wrapper.m`

% Optional inputs for `FEMA_wrapper.m` depending on analysis
contrasts=[]; % Contrasts relate to columns in design matrix e.g. [1 -1] will take the difference between cols 1 and 2 in your design matrix (X).  This needs to be padded with zeros at the beginning but not the end.
ranknorm = 0; % Rank normalizes dependent variables (Y) (default = 0)
nperms = 20; % Number of permutations - if wanting to use resampling methods nperms>0
RandomEffects = {'F','S','E'}; % Random effects to include: family, subject, error
mediation = 0; % If wanting to use outputs for a mediation analysis set mediation=1 - ensures same resampling scheme used for each model in fname_design
PermType = 'wildbootstrap'; %Default resampling method is null wild-bootstrap - to run mediation analysis need to use non-null wild-bootstrap ('wildboostrap-nn')
tfce = 0; % If wanting to run threshold free cluster enhancement (TFCE) set tfce=1 (default = 0)
colsinterest=[1]; % Only used if nperms>0. Indicates which IVs (columns of X) the permuted null distribution and TFCE statistics will be saved for (default 1, i.e. column 1)


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% DEMO VERTEXWISE ANALYSIS

% This demo produces vertexwise age associations with cortical thickness controlling for 
% sociodemographic information and genetic ancestry.  The results can be visualised using
% `showSurf_demo.m`

if doVertexwise

      datatype='vertex'; % imaging modality selected
      modality='smri'; % concatenated imaging data stored in directories based on modality (smri, dmri, tfmri, roi)

      % Uses path structure in abcd-sync to automatically find data
      dirname_imaging = fullfile(abcd_sync_path, dataRelease, 'imaging_concat/vertexwise/', modality); % filepath to imaging data
  
      switch dataRelease
            case '4.0'
                  fstem_imaging = 'thickness_ic5_sm256'; % name of imaging phenotype - data already saved as ico=5
                  %a few of the many other choices in 4.0, including now RSI: 
                  %   'thickness_ic5_sm1000', 'area_ic5_sm1000',
                  %   'N0-gm_ic5_sm256','N0-gwc_ic5_sm256','N0-wm_ic5_sm256',...
                  %   'ND-gm_ic5_sm256','ND-gwc_ic5_sm256','ND-wm_ic5_sm256'
      end
  
      ico = 5; % icosahedral number
  
      % Once all filepaths and inputs have been specified FEMA_wrapper.m can be run in one line

      if doTFCEvertex==1

            tfce=1;
            PermType='wildbootstrap';

      elseif doTFCEvertex==0

            tfce=0;

      end

      % RUN FEMA
      [fpaths_out beta_hat beta_se zmat logpmat sig2tvec sig2mat beta_hat_perm beta_se_perm zmat_perm sig2tvec_perm sig2mat_perm inputs mask tfce_perm colnames_interest save_params] = FEMA_wrapper(fstem_imaging, fname_design, dirname_out, dirname_tabulated, dirname_imaging, datatype,...
      'ico', ico, 'ranknorm', ranknorm, 'contrasts', contrasts, 'RandomEffects', RandomEffects, 'pihat_file', fname_pihat, 'nperms', nperms, 'mediation',mediation,'PermType',PermType,'tfce',tfce,'colsinterest',colsinterest);

      % Estimate p-values based on the empirical null distribution for the z statistics
      stattype='z';
      [fpath_out,rank_log10pval_uncorr, rank_log10pval_fwecorr]=FEMA_perm_significance(zmat_perm, stattype, colnames_interest, 'save_params',save_params,'extrapolate',0,'mask',mask);

      if ~isempty(tfce_perm)

            % If TFCE run, estimate p-values for TFCE statistics based on the empirical null distribution
            stattype='tfce';
            [fpath_out,rank_log10pval_uncorr, rank_log10pval_fwecorr]=FEMA_perm_significance(tfce_perm, stattype, colnames_interest, 'save_params',save_params,'extrapolate',0,'mask',mask);

      end

      if doSobelvertex==1

            tfce=0; % Cannot run TFCE and mediation at the same time
            mediation=1; % Set mediation to 1
            PermType='wildbootstrap-nn'; % Requires non-null WB

            % For mediation analysis, FEMA requires two nested models: reduced and full. These can be specified as a cell-array in fname_design as below.

            switch dataRelease
            case '4.0'
                  atlasVersion = 'ABCD2_cor10';
                  fname_design = {fullfile(abcd_sync_path,'4.0','support_files/design_matrices','FEMA_demo_rel4.0_designmatrix_mediation_reduced.txt');
                                    fullfile(abcd_sync_path,'4.0','support_files/design_matrices','FEMA_demo_rel4.0_designmatrix_mediation_full.txt')};
            end

            % As two design matrices are provided. Two separate output directories are required.
            [~,filename]=fileparts(fname_design{1});
            dirname_out_med{1}=strcat(outdir_root,'/',filename);
            [~,filename]=fileparts(fname_design{2});
            dirname_out_med{2}=strcat(outdir_root,'/',filename);

            % RUN FEMA
            [fpaths_out beta_hat beta_se zmat logpmat sig2tvec sig2mat beta_hat_perm beta_se_perm zmat_perm sig2tvec_perm sig2mat_perm inputs mask tfce_perm colnames_interest save_params] = FEMA_wrapper(fstem_imaging, fname_design, dirname_out_med, dirname_tabulated, dirname_imaging, datatype,...
            'ico', ico, 'ranknorm', ranknorm, 'contrasts', contrasts, 'RandomEffects', RandomEffects, 'pihat_file', fname_pihat, 'nperms', nperms, 'mediation',mediation,'PermType',PermType,'tfce',tfce,'colsinterest',colsinterest);

            FEMA_outfile_reduced=fpaths_out{1};
            FEMA_outfile_full=fpaths_out{2};

            % Uses save_params from full model output (last model in fname_design)

            % For mediation analysis, the user must specify the independent variable of interest that is present in both design matrices as `IVname`.
            % Sobel's test is testing for a mediation effect on the association between the specified IV of interest (`IVname`) and the
            % imaging phenotype (`fstem_imaging`). The mediator is the variable that is present in the full model, but absent from the reduced model.
            IVname='interview_age';
            
            alpha=0.05; %generates 95% confidence intervals

            [fpath_out, mediation_tstat, mediation_pvals]=FEMA_sobel_test(FEMA_outfile_reduced,FEMA_outfile_full,IVname,alpha, 'save_params',save_params);


      end

end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% DEMO VOXELWISE ANALYSES

% This demo produces voxelwise age associations with restricted normalised isotropic (RNI) diffusio controlling for 
% sociodemographic information and genetic ancestry.  The results can be visualised using
% `showVol_demo.m`

if doVoxelwise

      datatype = 'voxel'; % imaging modality selected
      modality='dmri'; % concatenated imaging data stored in directories based on modality (smri, dmri, tfmri, roi)

      % uses path structure in abcd-sync to automatically find data
      dirname_imaging = fullfile(abcd_sync_path, dataRelease, '/imaging_concat/voxelwise/', atlasVersion, modality); % filepath to imaging data
      fname_design = fullfile(abcd_sync_path,'4.0','support_files/design_matrices','FEMA_demo_rel4.0_designmatrix_age.txt');
      
      switch dataRelease
            case '4.0'
                  fstem_imaging = 'RNI'; %Name of imaging phenotype - possible values: DTI measures (FA, MD) or RSI measures (FNI, RDF, RD, RIF, RI, RND, RNI, RNT + same for 'H' or Jacobian (JA)
      end
  
      % Once all filepaths and inputs have been specified FEMA_wrapper.m can be run in one line

      if doTFCEvoxel==1

            tfce=1;
            mediation=0;
            PermType='wildbootstrap';

      elseif doTFCEvoxel==0

            tfce=0;
            mediation=0;

      end

      % RUN FEMA
      [fpaths_out beta_hat beta_se zmat logpmat sig2tvec sig2mat beta_hat_perm beta_se_perm zmat_perm sig2tvec_perm sig2mat_perm inputs mask tfce_perm colnames_interest save_params] = FEMA_wrapper(fstem_imaging, fname_design, dirname_out, dirname_tabulated, dirname_imaging, datatype,...
      'ico', ico, 'ranknorm', ranknorm, 'contrasts', contrasts, 'RandomEffects', RandomEffects, 'pihat_file', fname_pihat, 'nperms', nperms, 'mediation',mediation,'PermType',PermType,'tfce',tfce,'colsinterest',colsinterest);

      % Estimate p-values based on the empirical null distribution for the z statistics
      stattype='z';
      [fpath_out,rank_log10pval_uncorr, rank_log10pval_fwecorr]=FEMA_perm_significance(zmat_perm, stattype, colnames_interest, 'save_params',save_params,'extrapolate',0, 'mask',mask);

      if ~isempty(tfce_perm)

            % If TFCE run, estimate p-values for TFCE statistics based on the empirical null distribution
            stattype='tfce';
            [fpath_out,rank_log10pval_uncorr, rank_log10pval_fwecorr]=FEMA_perm_significance(tfce_perm, stattype, colnames_interest, 'save_params',save_params,'extrapolate',0, 'mask',mask);

      end

      if doSobelvoxel==1

            tfce=0; % Cannot run TFCE and mediation at the same time
            mediation=1; % Set mediation to 1
            PermType='wildbootstrap-nn'; % Requires non-null WB

            % For mediation analysis, FEMA requires two nested models: reduced and full. These can be specified as a cell-array in fname_design as below.

            switch dataRelease
            case '4.0'
                  atlasVersion = 'ABCD2_cor10';
                  fname_design = {fullfile(abcd_sync_path,'4.0','support_files/design_matrices','FEMA_demo_rel4.0_designmatrix_mediation_reduced.txt');
                                    fullfile(abcd_sync_path,'4.0','support_files/design_matrices','FEMA_demo_rel4.0_designmatrix_mediation_full.txt')};
            end

            % As two design matrices are provided. Two separate output directories are required.
            [~,filename]=fileparts(fname_design{1});
            dirname_out_med{1}=strcat(outdir_root,'/',filename);
            [~,filename]=fileparts(fname_design{2});
            dirname_out_med{2}=strcat(outdir_root,'/',filename);

            % RUN FEMA
            [fpaths_out beta_hat beta_se zmat logpmat sig2tvec sig2mat beta_hat_perm beta_se_perm zmat_perm sig2tvec_perm sig2mat_perm inputs mask tfce_perm colnames_interest save_params] = FEMA_wrapper(fstem_imaging, fname_design, dirname_out_med, dirname_tabulated, dirname_imaging, datatype,...
            'ico', ico, 'ranknorm', ranknorm, 'contrasts', contrasts, 'RandomEffects', RandomEffects, 'pihat_file', fname_pihat, 'nperms', nperms, 'mediation',mediation,'PermType',PermType,'tfce',tfce,'colsinterest',colsinterest);

            FEMA_outfile_reduced=fpaths_out{1};
            FEMA_outfile_full=fpaths_out{2};

            % Uses save_params from full model output (last model in fname_design)

            % For mediation analysis, the user must specify the independent variable of interest that is present in both design matrices as `IVname`.
            % Sobel's test is testing for a mediation effect on the association between the specified IV of interest (`IVname`) and the
            % imaging phenotype (`fstem_imaging`). The mediator is the variable that is present in the full model, but absent from the reduced model.
            IVname='interview_age';
            
            alpha=0.05; %generates 95% confidence intervals

            [fpath_out, mediation_tstat, mediation_pvals]=FEMA_sobel_test(FEMA_outfile_reduced,FEMA_outfile_full,IVname,alpha, 'save_params',save_params,'mask',mask);


      end

end
