%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% DEMO USING showVol.m WITH FEMA OUTPUT
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% This script is designed to provide an end-to-end demo of how to
% visualize VOXELWISE results generated from running FEMA.

% `showVol.m` is a tool to plot whole brain voxelwise statistics. It uses
% an interactive GUI that allows the user to move around the brain, explore
% effect from different planes and overlay regions of interest (ROIs).
% For the ABCD Study, subcortical parcellations from freesurfer and
% other external atlases have been registered to the ABCD atlas space. For
% these to provide accurate spatial labelling, statistical volumes also need
% to be within the ABCD atlas space. The function `convertFEMAVols` ensures
% volume data are aligned with these ROIs. This function also allows the user to
% set color overlays of statistical effects on top of a background image
% (default: T1w image).  The below demo shows the user how to
% visualize ACBD analysis results that used data registered the ABCD atlas space.

% Code written by Anders Dale, John Iversen, Clare E Palmer and Diliana Pecheva, 2021
%
% This software is Copyright (c) 2021 The Regents of the University of California. All Rights Reserved.
% See LICENSE.

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% DEMO SCRIPT

% Run `FEMA_showVol_demo.m` in MATLAB as a script for a full end-to-end demo.
% The demo will prompt you to give user specific inputs in the command window.

% If you would like to edit the script manually to change some settings then SAVE A LOCAL
% COPY of this script before running.

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% REQUIREMENTS TO RUN showVol.m
%
% CLONE GitHub directory for cmig_tools
% ADD to MATLAB path:

% e.g. if cloned into ~/github: 
% addpath(genpath('~/github/cmig_tools'))

git_dir=input('Please specify the path to your cmig_tools directory. If this is already in your path leave blank and press ENTER: ','s');

if ~isempty(git_dir)
fprintf('\nAdding cmig_tools directory to MATLAB path... \n')
addpath(genpath(git_dir))
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% VISULATION OF VOXELWISE FEMA OUTPUT USING 'showVol'


% 1) Specify where FEMA_wrapper output is saved and load into MATLAB workspace

dirname_out=input('Specify directory of where FEMA output saved: ','s'); % directory of where FEMA output saved
fstem_imaging='RNI'; %imaging phenotype used for analysis
fname_results = sprintf('%s/FEMA_wrapper_output_voxel_%s.mat',dirname_out,fstem_imaging);
load(fname_results,'vol_beta_hat','vol_z','colnames_model'); % load FEMA output - only need some variables

% 2) Specify ABCD release version as the atlas used for voxelwise registration is different from 3.0 to 4.0

dataRelease='4.0';

switch dataRelease
case '3.0'
  atlasVersion = 'ABCD1_cor10';
case '4.0'
  atlasVersion = 'ABCD2_cor10';
end
  
loadPrerenderedIfNeeded(atlasVersion) % loads atlas specific data
global PRI % saves atlas data as a global variable so does not need to be loaded every time you run showVol

vol_stat = vol_z; %specify output to plot e.g. vol_z, vol_beta_hat etc
stat_name = 'z stat'; %name to label output

nvox=156662;
thresh_bonfcorr=abs(norminv(0.05/nvox));

% 3) Use a helper function to convert FEMA output into volumes for showVol. 

%  for each IV, creates a raw grayscale map and colorized map overlaid on T1 
limits = [-40 40 thresh_bonfcorr];  % Limits for colormap: (1) min (2) max (3) threshold (leave third value empty to have unthresholded maps)
index = [1 2 3];                    % Indicates which IVs to plot: colnames_model(index)
showPval = false;                   % If showPval==true, creates a map of log10 p-values
cmap = redblackblue_alpha;          % Colormap
cmap(:,4) = 1;
bg_vol = [];                        % Set background image to overlay statistics. Default: atlas T1
CSFprob = [];                       % Probabilistic threshold for masking CSF e.g. if CSFprob=0.8, only voxels in which 80% of participants labelled that voxel as CSF will be masked (Default: no masking)
interpMethod = [];                  % By default uses linear interpolation to go from 2mm to 1mm voxels. Specify 'nearest' to show actual resolution

vols = convertFEMAVols(vol_stat, fstem_imaging, stat_name, colnames_model, index, limits, showPval, cmap, bg_vol, CSFprob, interpMethod, atlasVersion);

% 4) Visualize using showVol
switch atlasVersion
      case {'ABCD1', 'ABCD1_cor10'}
            %for backwards compatibility, PRI naming is different for ABCD1
            showVol(vols, PRI.vol_T1,PRI.vol_T2,PRI.vol_aseg_rgb,PRI.vol_fiber_rgb,PRI.vol_CO,PRI.vol_FOD_1mm)
      
      case {'ABCD2', 'ABCD2_cor10'}
            coords = [-7 -13 -4]; % Opens showVol to a specific location e.g. R NAcc
            showVol(vols, PRI.ABCD2.T1, PRI.ABCD2.CO, PRI.ABCD2.FOD, struct('roiatlas','ABCD2'), coords)
end
