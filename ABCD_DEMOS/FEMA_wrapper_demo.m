%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% FEMA DEMO USING FEMA_wrapper.m
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% This script is designed to provide an end-to-end demo of the
% FEMA code and can be used as an example script for how to develop your
% own scripts to analyse ABCD data using FEMA via `FEMA_wrapper.m`

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
% This script uses example design matrices within the `support_files`
% directory on `abcd-sync`, and concatenated imaging data within the
% `imaging_concat` directory on `abcd-sync`.  For the voxelwise analysis, this
% demo script will replicate the findings from Palmer et al., (2022),
% Developmental Cognitive Neuroscience
% (ncbi.nlm.nih.gov/pmc/articles/PMC8671104/)

% In this paper, Palmer et al. showed voxelwise associations between
% restricted diffusion metrics derived from the diffusion MRI data and age
% from 9-14 years in the ABCD sample.  This demo will conduct these same
% developmental associations using the same design matrix and analysis inputs
% as in this paper.  This demo also includes age associations with other imaging modalities: cortical thickness, resting-state functional connectivity.

% After running this demo and saving the output from FEMA, the results can
% be visualised using our visualisation tools: `showSurf` and `showVol`.
% We have included demos for these tools, which will run off the results from
% this demo. There is one for vertexwise data (`FEMA_showSurf_demo.m`) and
% one for voxelwise data (`FEMA_showVol_demo.m`).

% This software is Copyright (c) 2021 The Regents of the University of California. All Rights Reserved.
% See LICENSE.

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% ORGANISATION OF THE FEMA PACKAGE

% FEMA_wrapper.m is a wrapper function designed to simplify analysing data
% from the ABCD Study.  FEMA_wrapper.m calls the three main functions
% required to run FEMA with ABCD data:

%   1) FEMA_process_data.m = loads imaging/behavioral data (Y) and
%   intersects with imaging info relevant for ABCD analysis

%   2) FEMA_intersect_design.m = loads design matrix (X) and intersects
%   with imaging/behavioral data (Y)

%   3) FEMA_fit.m = this is the core FEMA function that runs the linear
%   mixed effects analysis.  If wanting to use FEMA to analyse other data
%   not related to ABCD you can directly input data into FEMA_fit.m and
%   write your own wrapper script to save the output.

% FEMA_wrapper.m saves the output from FEMA_fit.m into a specified
% output directory

% There are several inputs that can be given to FEMA_wrapper.m that will
% specify the analyses you would like to run.  Examples of these will be 
% shown below.  There is more detail on these in the documentation for FEMA_wrapper.m

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% DEMO SCRIPT WITH TUTORIAL

% Run `FEMA_wrapper_demo.m` in Matlab as a script for a full end-to-end demo.
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

dirname_out = fullfile('~/FEMA_demo_results',dataRelease);
fprintf('\n\nStoring results here: %s \n', dirname_out)
outcheck=input('Would you like to specify a DIFFERENT output directory? (enter numeric value: 1 or 0): ');
if outcheck==1
      dirname_out=input('Please specify output directory: ','s');
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
fprintf('\n\nFor this demo, please specify which imaging data to demo: \n')
doVertexwise = input('Run vertexwise demo? (enter numeric value: 1 or 0) '); % run vertexwise analysis (datatype = 'vertex')
doVoxelwise = input('Run voxelwise demo? (enter numeric value: 1 or 0) '); % run voxelwise analysis (datatype = 'voxel')
doCorrmat  = input('Run corrmat demo? (enter numeric value: 1 or 0) '); % run corrmat analysis (datatype = 'corrmat')
doExternal = input('Run external demo? (enter numeric value: 1 or 0) '); % run external analysis (datatype = 'external') 
doSynthTest = input('Run synthesized vertexwise demo? (enter numeric value: 1 or 0) '); % run example synthesized analysis (datatype = 'vertex'; synth = 1)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% 5) If editting FEMA_wrapper_demo.m SAVE LOCAL COPY

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
  case '3.0'
    atlasVersion = 'ABCD1_cor10';
    dirname_tabulated = fullfile(abcd_sync_path,'3.0','tabulated/released');
    fname_design = fullfile(abcd_sync_path,'3.0','support_files/design_matrices','FEMA_demo_rel3.0_designmatrix_age.txt');
    fname_pihat = fullfile(abcd_sync_path, '3.0','genomics','ABCD_rel3.0_pihat.mat');
  case '4.0'
    atlasVersion = 'ABCD2_cor10';
    dirname_tabulated = fullfile(abcd_sync_path,'4.0','tabulated/released'); 
    fname_design = fullfile(abcd_sync_path,'4.0','support_files/design_matrices','FEMA_demo_rel4.0_designmatrix_age.txt');
    fname_pihat = fullfile(abcd_sync_path, '4.0','genomics','ABCD_rel4.0_pihat.mat'); 
end

[root,filename]=fileparts(fname_design);
dirname_out=strcat(dirname_out,'/',filename);

% NOTE: `fname_design` can also be a cell array of multiple design matrices
% (txt file paths) that will be looped over in `FEMA_wrapper.m`
% Use makeDesign.R from ~/github/cmig_utils/r to make design matrix compatible with `FEMA_wrapper.m`
% See makeDesign_demo.R for tutorial on how to make design matrices
% compatible with `FEMA_wrapper.m`

% Optional inputs for `FEMA_wrapper.m` depending on analysis
contrasts=[]; % Contrasts relate to columns in design matrix e.g. [1 -1] will take the difference between cols 1 and 2 in your design matrix (X).  This needs to be padded with zeros at the beginning but not the end.
ranknorm = 1; % Rank normalizes dependent variables (Y) (default = 0)
nperms = 0; % Number of permutations - if wanting to use resampling methods nperms>0
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
            case '3.0'
                  %fstem_imaging = 'area-sm256';
                  fstem_imaging = 'thickness-sm256'; % name of imaging phenotype
                  %fstem_imaging = 'sulc-sm256';
            case '4.0'
                  fstem_imaging = 'thickness_ic5_sm256'; % name of imaging phenotype - data already saved as ico=5
                  %a few of the many other choices in 4.0, including now RSI: 
                  %   'thickness_ic5_sm1000', 'area_ic5_sm1000',
                  %   'N0-gm_ic5_sm256','N0-gwc_ic5_sm256','N0-wm_ic5_sm256',...
                  %   'ND-gm_ic5_sm256','ND-gwc_ic5_sm256','ND-wm_ic5_sm256'
      end
  
      ico = 5; % icosahedral number
  
      % Once all filepaths and inputs have been specified FEMA_wrapper.m can be run in one line

      % RUN FEMA
      [fpaths_out beta_hat beta_se zmat logpmat sig2tvec sig2mat beta_hat_perm beta_se_perm zmat_perm sig2tvec_perm sig2mat_perm inputs mask tfce_perm analysis_params] = FEMA_wrapper(fstem_imaging, fname_design, dirname_out, dirname_tabulated, dirname_imaging, datatype,...
      'ico', ico, 'ranknorm', ranknorm, 'contrasts', contrasts, 'RandomEffects', RandomEffects, 'pihat_file', fname_pihat, 'nperms', nperms, 'mediation',mediation,'PermType',PermType,'tfce',tfce,'colsinterest',colsinterest);

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

      switch dataRelease
            case '3.0'
                  fstem_imaging = 'N0'; % Name of imaging phenotype - possible values: DTI measures (FA, MD) or RSI measures (N0, ND) or Jacobian (JA)
            case '4.0'
                  fstem_imaging = 'RNI'; %Name of imaging phenotype - possible values: DTI measures (FA, MD) or RSI measures (FNI, RDF, RD, RIF, RI, RND, RNI, RNT + same for 'H' or Jacobian (JA)
      end
  
      % Once all filepaths and inputs have been specified FEMA_wrapper.m can be run in one line

      % RUN FEMA
      [fpaths_out beta_hat beta_se zmat logpmat sig2tvec sig2mat beta_hat_perm beta_se_perm zmat_perm sig2tvec_perm sig2mat_perm inputs mask tfce_perm analysis_params] = FEMA_wrapper(fstem_imaging, fname_design, dirname_out, dirname_tabulated, dirname_imaging, datatype,...
      'ranknorm', ranknorm, 'contrasts', contrasts, 'RandomEffects', RandomEffects, 'pihat_file', fname_pihat, 'nperms', nperms, 'mediation',mediation,'PermType',PermType,'tfce',tfce,'colsinterest',colsinterest);

end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% DEMO CORRMAT ANALYSIS

if doCorrmat

      datatype = 'corrmat'; % imaging modality selected
      modality = 'restingstate'; % concatenated corrmat data stored in directories based on modality (fmritask, restingstate)

      % uses path structure in abcd-sync to automatically find data
      dirname_imaging = fullfile(abcd_sync_path, '4.0', '/imaging_concat/corrmat',modality); % filepath to imaging data

      fstem_imaging = 'rsfmri_fd0.20'; % name of imaging phenotype
      
      % Once all filepaths and inputs have been specified FEMA_wrapper.m can be run in one line

      % RUN FEMA
      [fpaths_out beta_hat beta_se zmat logpmat sig2tvec sig2mat beta_hat_perm beta_se_perm zmat_perm sig2tvec_perm sig2mat_perm inputs mask tfce_perm analysis_params] = FEMA_wrapper(fstem_imaging, fname_design, dirname_out, dirname_tabulated, dirname_imaging, datatype,...
      'ranknorm', ranknorm, 'contrasts', contrasts, 'RandomEffects', RandomEffects, 'pihat_file', fname_pihat, 'nperms', nperms, 'mediation',mediation,'PermType',PermType,'tfce',tfce,'colsinterest',colsinterest);

end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% DEMO EXTERNAL ANALYSIS

if doExternal

      datatype='external'; % can use txt with columns of ANY data type (e.g. ROIs, behavior) - runs mass univaraite LME across every column
  
      %External txt file (Y) must be created in the following way:
      % Y(:,1) = subject ID --> colname: 'src_subject_id'
      % Y(:,2) = event --> colname: 'eventname'
      % Y(:,3:end) = imaging data
      fstem_imaging='RNI';
      dirname_imaging = fullfile(abcd_sync_path, '4.0', '/imaging_concat/external/',sprintf('ALL_%s_ROIs_bysubj_mean.txt',fstem_imaging)); % example external data file

      % Once all filepaths and inputs have been specified FEMA_wrapper.m can be run in one line
      
      % RUN FEMA
      [fpaths_out beta_hat beta_se zmat logpmat sig2tvec sig2mat beta_hat_perm beta_se_perm zmat_perm sig2tvec_perm sig2mat_perm inputs mask tfce_perm analysis_params] = FEMA_wrapper(fstem_imaging, fname_design, dirname_out, dirname_tabulated, dirname_imaging, datatype,...
      'ranknorm', ranknorm, 'contrasts', contrasts, 'RandomEffects', RandomEffects, 'pihat_file', fname_pihat, 'nperms', nperms, 'mediation',mediation,'PermType',PermType,'tfce',tfce,'colsinterest',colsinterest);

end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% DEMO NULL SYNTHESIZED VOXELWISE/VERTEXWISE DATA

if doSynthTest

      datatype = 'vertex'; % depends on what type of data you want to synthesize
      modality='smri'; % concatenated imaging data stored in directories based on modality (smri, dmri, tfmri, roi)

      % Uses path structure in abcd-sync to automatically find data
      dirname_imaging = fullfile(abcd_sync_path, dataRelease, 'imaging_concat/vertexwise/', modality); % filepath to imaging data
  
      switch dataRelease
            case '3.0'
                  %fstem_imaging = 'area-sm256';
                  fstem_imaging = 'thickness-sm256'; % name of imaging phenotype
                  %fstem_imaging = 'sulc-sm256';
            case '4.0'
                  fstem_imaging = 'thickness_ic5_sm256'; % name of imaging phenotype
                  %a few of the many other choices in 4.0, including now RSI: 
                  %   'thickness_ic5_sm1000', 'area_ic5_sm1000',
                  %   'N0-gm_ic5_sm256','N0-gwc_ic5_sm256','N0-wm_ic5_sm256',...
                  %   'ND-gm_ic5_sm256','ND-gwc_ic5_sm256','ND-wm_ic5_sm256'
      end
  
      ico = 5; % ico number of vertexwise data

      synth = 1; % will intiate FEMA_synthesize.m inside FEMA_wrapper.m to generate null synthesized data based on Y
  
      % Once all filepaths and inputs have been specified FEMA_wrapper.m can be run in one line
      
      % RUN FEMA
      [fpaths_out beta_hat beta_se zmat logpmat sig2tvec sig2mat beta_hat_perm beta_se_perm zmat_perm sig2tvec_perm sig2mat_perm inputs mask tfce_perm analysis_params] = FEMA_wrapper(fstem_imaging, fname_design, dirname_out, dirname_tabulated, dirname_imaging, datatype,...
       'ico', ico, 'ranknorm', ranknorm, 'contrasts', contrasts, 'RandomEffects', RandomEffects, 'pihat_file', fname_pihat, 'nperms', nperms, 'mediation',mediation,'PermType',PermType,'tfce',tfce,'colsinterest',colsinterest,'synth',synth);

end




