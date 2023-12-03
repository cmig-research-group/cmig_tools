% ABCD_anatomy_visualization_atlas_template
% catch-all example of visualizing ABCD prerendered volumes for anatomy work

%% select atlas version and load prerendered images and FOD. 
%  If not at UCSD, showVolData files can be downloaded here: 
%  https://drive.google.com/drive/folders/1Uq04VEKRkgQ5SlWDfV9T4KXbwNuCY0nR?usp=sharing
%  Also, see https://github.com/cmig-research-group/showVol

% Select version of ABCD atlas

atlasVersion = 'ABCD2_cor10';
fodtype = {'v','p80'}; %load several for now, to compare!

% New in Dec 2023 -- not totally complete (no FOD, need some double-checking)
atlasVersion = 'ABCD3_cor10';
fodtype = [];

% older
% atlasVersion = 'ABCD1_cor10';
% fodtype={''};

% Load prerendered images and FODs
cfg = abcdConfig('showVol'); %will ask you to locate the showVolData directory on first run
loadPrerenderedIfNeeded(atlasVersion, fodtype) %loads prerendered images into cache only if necessary to save time
global PRI %this global variable contains loaded prerendered images

%% Example of loading in other volumes e.g. a sharper T1 from thalamus atlas
fname = fullfile(cfg.data.showVolData, ['Atlas/T1_thalamus_' atlasVersion '.mat']);
if exist(fname,'file')
  tmp = load(fname); %var has same name as file, so T1_thalamus_ABCD2
  try
    thalamus_T1 = tmp.vol_T1; %FIXME--naming conventions across atlas versions need standardizing
  catch
    thalamus_T1 = tmp.(['T1_thalamus_' atlasVersion(1:5)]);
  end
end


%% examples of using config structure to customize showVol
cfg2=[];
cfg2.roiatlas = atlasVersion(1:5);          %which atlas's ROIs to use. Defaults to ABCD1
% cfg2.RCS       = [100 116 132];      %open with crosshairs at specific coordinate, e.g. for Left GPi
% cfg2.winpos    = [0.5 0.5 200 100];  %customize this to your liking (do "get(gca,'Position')" to learn the position
% cfg2.linewidth = 3;                  %thicker crosshair lines
% cfg2.link      = false;              %don't link this figures cursor position to other windows
% cfg2.roifile   = atlas_file;         %ROI probability volumes
% cfg2.roiset    = filename            %open a saved ROI display settings file (saved from an earlier run of showVol)

%% running showVol
% shows two examples for different atlas versions -- the PRI naming hasn't been quite standardized yet. For now, just inspect the 
%   contents of the PRI variable to see what is available and what it is called...
switch cfg2.roiatlas
  case 'ABCD1'
    showVol(atlas_T1('ABCD1'), thalamus_T1, PRI.ABCD1.T2, PRI.ABCD1.aseg_rgb, PRI.ABCD1.fiber_rgb, ...
      PRI.ABCD1.CO, PRI.ABCD1.FOD, cfg2);
    
  case 'ABCD2'
    showVol(atlas_T1('ABCD2'), atlas_T1('ABCD1'),thalamus_T1, PRI.ABCD2.aparcaseg_rgb, PRI.ABCD2.fiber_rgb, PRI.ABCD2.CO, ...
       PRI.ABCD2.FOD, PRI.ABCD2.FOD_p80, cfg2);
     
  case 'ABCD3'
    showVol(atlas_T1('ABCD3'), atlas_T1('ABCD3', true), atlas_T1('ABCD2'), atlas_T1('ABCD1'), PRI.ABCD3.aseg_rgb, PRI.ABCD3.aparcaseg_rgb,...
      PRI.ABCD3.fiber_rgb, PRI.ABCD3.CO, cfg2);
end

% TO DO: standardize variable and volume names
% fix T1/C0 mixup for ABCD2 in prepareAtlases



