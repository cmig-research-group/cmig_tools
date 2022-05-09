%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% DEMO USING showSurf.m WITH FEMA OUTPUT
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% This script is designed to provide an end-to-end demo of how to
% visualize VERTEXWISE results generated from running `FEMA_wrapper_demo.m`

% `showSurf.m` is a tool to plot surface based statistics onto a pial or
% inflated cortical surface based on FreeSurfer segmentation.  This tool
% produces MATLAB figures across different views of the brain (lateral and
% medial; left and right hemispheres) that can be used for publications.

% For publication high resolution quality images we recommend running
% vertexwise analyses at the highest icosahedral order (ico) available 7.

% Code written by Anders Dale, John Iversen, Clare E Palmer and Diliana Pecheva, 2021
%
% This software is Copyright (c) 2021 The Regents of the University of California. All Rights Reserved.
% See LICENSE.

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% DEMO SCRIPT

% Run `FEMA_showSurf_demo.m` in MATLAB as a script for a full end-to-end demo.
% The demo will prompt you to give user specific inputs in the command window.

% If you would like to edit the script manually to change some settings then SAVE A LOCAL
% COPY of this script before running.

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% REQUIREMENTS TO RUN showSurf.m
%
% CLONE GitHub directories for cmig_utils, showSurf
% ADD ALL ABCD CMIG tools directories to MATLAB path:

% e.g. if cloned into ~/github: 
% addpath(genpath('~/github/cmig_utils'))
% addpath(genpath('~/github/showSurf'))

git_dir=input('Please specify the path to your cmig_tools directory. If this is already in your path leave blank and press ENTER: ','s');

if ~isempty(git_dir)
fprintf('\nAdding cmig_tools directory to MATLAB path... \n')
addpath(genpath(git_dir))
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% VISUALIZATION OF VERTEXWISE FEMA OUTPUT USING 'showSurf.m'

% 1) Specify where FEMA_wrapper output is saved and load into MATLAB workspace

dirname_out=input('Specify directory of where FEMA output saved: ','s'); % directory of where FEMA output saved
fstem_imaging='thickness_ic5_sm256'; % imaging phenotype used for analysis
fname_results = sprintf('%s/FEMA_wrapper_output_vertex_%s.mat',dirname_out,fstem_imaging); % FEMA output filename
load(fname_results) % load FEMA output

dataRelease='4.0'; % ABCD data release

% 2) Load surface templates for plotting

% In order to plot statistics on a cortical surface, we require a template
% of the cortical surface to project the results onto.  Below we load these
% surface templates.  `SurfView_surfs.mat` contains surfaces for all icos.

load SurfView_surfs.mat % load surface templates
ico=5; % ico number
icnum=ico+1; % index for ico number (icnum = ico + 1)
icnvert = size(icsurfs{icnum}.vertices,1); % indices of vertices for specified icosahedral order

% 3) Plot surface-wise statistics for different IVs of interest

% The demo below will produce figures for the IVs (columns of X) from
% the FEMA analysis specfiied by `ncoeff`

ncoeff=1;

for coeffnum = ncoeff

      statname = 'z stat'; % specify name of statistic plotting for figure label
      vertvals = zmat(coeffnum,:); % specify statistics to plot
      vertvals_lh = vertvals(1:icnvert); % divide statistics by hemisphere for plotting
      vertvals_rh = vertvals(icnvert+[1:icnvert]);
    
      % specify limits for plot based on vertvals
      fmax = min(300,max(abs(vertvals))); % max limit for colorbar
      fmin = 0.0; % min limit of colorbar
      fmid = fmax/2; % middle of colorbar
      fvals = [fmin fmid fmax];
      clim = [-fmax fmax]; % set colorbar limits
      
      cm = blueblackred(); % set colormap

      curvcontrast = [0.2 0.2]; % contrast of gyri/sulci

      bgcol = [0 0 0]; % change to [1 1 1] for white background
      
      fh = figure(coeffnum + 100*(str2num(dataRelease(1))-3)); clf; % number matlab figure window
      set(fh,'Color',bgcol); fh.InvertHardcopy = 'off';
      subplot(2,2,1); SurfView_show_new(surf_lh_pial,surf_rh_pial,vertvals_lh,vertvals_rh,fvals,cm,'left', [1 0],curvvec_lh,curvvec_rh,icsurfs{icnum},[],curvcontrast,bgcol); set(gca,'visible','off'); axis tight;
      subplot(2,2,2); SurfView_show_new(surf_lh_pial,surf_rh_pial,vertvals_lh,vertvals_rh,fvals,cm,'right',[0 1],curvvec_lh,curvvec_rh,icsurfs{icnum},[],curvcontrast,bgcol); set(gca,'visible','off'); axis tight;
      subplot(2,2,3); SurfView_show_new(surf_lh_pial,surf_rh_pial,vertvals_lh,vertvals_rh,fvals,cm,'right',[1 0],curvvec_lh,curvvec_rh,icsurfs{icnum},[],curvcontrast,bgcol); set(gca,'visible','off'); axis tight;
      subplot(2,2,4); SurfView_show_new(surf_lh_pial,surf_rh_pial,vertvals_lh,vertvals_rh,fvals,cm,'left', [0 1],curvvec_lh,curvvec_rh,icsurfs{icnum},[],curvcontrast,bgcol); set(gca,'visible','off'); axis tight;
      titleAx = axes;
      set(titleAx,'position',[0 0 1 1],'units','normalized');axis off;
      text(titleAx, 0.5,1,sprintf('%s ~ %s [%s]', fstem_imaging, colnames_model{coeffnum}, statname),'color','w','fontweight','bold','interpreter','none','verticalalignment','top','horizontalalignment','center','fontsize',14)
      
      % Set colorbar
      colormap(cm);
      cb = colorbar('color', 'w');
      % cb.Label.Interpreter = 'latex';
      cb.Label.String = strcat('z-score');
      cb.Label.FontSize = 10;
      cb.Box = 'off';
      cb.Position = [.92 .08 .02 .8150];
      caxis(clim);
end
