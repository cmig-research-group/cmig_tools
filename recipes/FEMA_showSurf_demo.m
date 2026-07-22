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

% 2) Specify visual preferences for plotting

legendPosition = 'South'; % 'South' or 'East'
title = 1; % whether to include title at top of plot
cm = blueblackred(); % set colormap (preferred: blueblackred or fire)
curvcontrast = [0.2 0.2]; % contrast of gyri/sulci
bgcol = [0 0 0]; % change to [1 1 1] for white background
polarity = 2; % set to 1 for unipolar values e.g. values that range 0 to 1

% 3) Load surface templates for plotting

% In order to plot statistics on a cortical surface, we require a template
% of the cortical surface to project the results onto.  Below we load these
% surface templates.  `SurfView_surfs.mat` contains surfaces for all icos.

load SurfView_surfs.mat % load surface templates
ico=5; % ico number
icnum=ico+1; % index for ico number (icnum = ico + 1)
icnvert = size(icsurfs{icnum}.vertices,1); % indices of vertices for specified icosahedral order

% 4) Plot surface-wise statistics for different IVs of interest

% The demo below will produce figures for the IVs (columns of X) from
% the FEMA analysis specfiied by `ncoeff`

ncoeff=1;

for coeffnum = ncoeff

      statname = 'z stat'; % specify name of statistic plotting for figure label
      vertvals = zmat(coeffnum,:); % specify statistics to plot
      vertvals_lh = vertvals(1:icnvert); % divide statistics by hemisphere for plotting
      vertvals_rh = vertvals(icnvert+[1:icnvert]);
      
      % specify limits for plot based on vertvals
      fmax = min(300,max(abs(vertvals))); % max limit for plotting purposes
      fmin = 0.0; % min limit for plotting purposes
      fmid = fmax/2; % middle value for plotting purposes
      fvals = [fmin fmid fmax]; % this will be passed to the SurfView_show_new function

      % set colorbar limits - usually [fmin fmax] or [-fmax fmax]
      clim = [-fmax fmax]; 
      
      % Create figure
      % The first two position arguments specify where the figure will
      % appear on the screen; the next two arguments (16 and 10) specify
      % the size of the figure in centimeters - 16 cm wide and 10 cm long
      % is a reasonable choice for an A4 sized paper
      fh = figure('Units', 'centimeters', 'Position', [10 10 16 10], 'Color', bgcol, 'InvertHardcopy', 'off');
      
      % If user wants to number the figure, these two lines should suffice
      % fh = figure(coeffnum + 100*(str2double(dataRelease(1))-3));
      % set(fh, 'Units', 'centimeters', 'Position', [10 10 16 10], 'Color', bgcol, 'InvertHardcopy', 'off');
      
      % Define spacing for axes
      % hvgap controls the horizontal and vertical spaces between the axes
      % btgap controls the space from the bottom of the figure and the top
      % of the figure respectively
      % lrgap controls the space from the left of the figure and the right
      % of the figure respectively
      hvgap = [0.02 0.02];
      if strcmpi(legendPosition, 'south')
          lrgap = [0.02 0.02];
          if title
              btgap = [0.12 0.08];
          else
              btgap = [0.12 0.01];
          end
      else
          if strcmpi(legendPosition, 'east')
              lrgap = [0.02 0.138];
              if title
                  btgap = [0.018 0.08];
              else
                  btgap = [0.018 0.018];
              end
          end
      end
      
      % Create axes
      allH = tight_subplot(2, 2, hvgap, btgap, lrgap);
      %hold(allH(:), 'on');

      axes(allH(1)); SurfView_show_new(surf_lh_pial,surf_rh_pial,vertvals_lh,vertvals_rh,fvals,cm,'left', [1 0],curvvec_lh,curvvec_rh,icsurfs{icnum},polarity,curvcontrast,bgcol); set(gca,'visible','off'); axis tight;
      axes(allH(2)); SurfView_show_new(surf_lh_pial,surf_rh_pial,vertvals_lh,vertvals_rh,fvals,cm,'right',[0 1],curvvec_lh,curvvec_rh,icsurfs{icnum},polarity,curvcontrast,bgcol); set(gca,'visible','off'); axis tight;
      axes(allH(3)); SurfView_show_new(surf_lh_pial,surf_rh_pial,vertvals_lh,vertvals_rh,fvals,cm,'right',[1 0],curvvec_lh,curvvec_rh,icsurfs{icnum},polarity,curvcontrast,bgcol); set(gca,'visible','off'); axis tight;
      axes(allH(4)); SurfView_show_new(surf_lh_pial,surf_rh_pial,vertvals_lh,vertvals_rh,fvals,cm,'left', [0 1],curvvec_lh,curvvec_rh,icsurfs{icnum},polarity,curvcontrast,bgcol); set(gca,'visible','off'); axis tight;
      
      if(title)
          titleAx = axes;
          set(titleAx,'position',[0 0 1 1],'units','normalized');axis off;
          text(titleAx, 0.5,1,sprintf('%s ~ %s [%s]', fstem_imaging, colnames_model{coeffnum}, statname),'color','w','fontweight','bold','interpreter','none','verticalalignment','top','horizontalalignment','center','fontsize',14)
      end

      % Set colorbar
      colormap(cm);
      cb                    = colorbar('color', 'w');
      cb.FontSize           = 10;
      cb.Label.String       = strcat('z-score');
      cb.Label.FontSize     = 12;
      cb.Label.FontWeight   = 'bold';   
      cb.Box                = 'off';

      if strcmpi(legendPosition, 'south')
          cb.Location = 'south';
          if title
              cb.Position(1)      = allH(1).Position(1);
              cb.Position(2)      = cb.Position(2) - hvgap(1);
              cb.Position(3)      = allH(1).Position(3)*2 + hvgap(1);
          else
              cb.Position(1)      = allH(1).Position(1);
              cb.Position(2)      = cb.Position(2) - btgap(1);
              cb.Position(3)      = allH(1).Position(3)*2 + hvgap(1);
          end
      else
          if strcmpi(legendPosition, 'east')
              cb.Location = 'eastoutside';
              if title
                  cb.Position(1)      = allH(4).Position(1) + allH(4).Position(3) + 0.01;
                  cb.Position(2)      = allH(3).Position(2);
                  cb.Position(4)      = allH(1).Position(4)*2 + hvgap(1);
              else
                  cb.Position(1)      = allH(4).Position(1) + allH(4).Position(3) + 0.16;
                  cb.Position(2)      = allH(3).Position(2);
                  cb.Position(4)      = allH(1).Position(4)*2 + hvgap(1);
              end
          end
      end
      caxis(clim);   
end