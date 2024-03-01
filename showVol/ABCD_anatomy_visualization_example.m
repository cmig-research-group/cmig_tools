% ABCD_anatomy_visualization_atlas_template
% catch-all example of visualizing ABCD prerendered volumes for anatomy work

%% select atlas version and load prerendered images and FOD. 
%  If not at UCSD, showVolData files can be downloaded here: 
%  https://drive.google.com/drive/folders/1Uq04VEKRkgQ5SlWDfV9T4KXbwNuCY0nR?usp=sharing
%  Also, see https://github.com/cmig-research-group/showVol
%
% TODO: will move showVolData to OSF. Link to cmig-tools on github

% Select version of ABCD atlas

atlasVersion = '5.0_ABCD3_cor10';
fodtype = [];

% older Atlases
% atlasVersion = 'ABCD1_cor10';
% fodtype={''};

% atlasVersion = 'ABCD2_cor10';
% fodtype = {'v','p80'}; %load several for now, to compare!

% Load prerendered images and FODs
cfg = abcdConfig('showVol'); %ensure ~/abcdConfig.json is configured to point to your showValData directory
loadPrerenderedIfNeeded(atlasVersion, fodtype) %loads prerendered images into cache on first call. This was most needed with FODs, but saves a little time regardles
global PRI %this global variable contains loaded prerendered images

%% Example of loading in other volumes e.g. a sharper T1 from thalamus atlas
fname = fullfile(cfg.data.showVolData, 'Atlas', atlasVersion, ['thalamus_T1_' atlasVersion '.mat']);
%fname = fullfile(cfg.data.showVolData, ['Atlas/T1_thalamus_' atlasVersion '.mat']); %older atlas versions pre-5.0_ABCD3

if exist(fname,'file')
  tmp = load(fname);
  try
    thalamus_T1 = tmp.vol_T1; %FIXME--naming conventions across atlas versions need standardizing
  catch
    try
      thalamus_T1 = tmp.(['T1_thalamus_' atlasVersion(1:5)]);
    catch
      thalamus_T1 = tmp.thalamus_T1;
    end
  end
end


%% examples of using config structure to customize showVol
cfg2=[];
cfg2.roiatlas = atlasVersion;          %which atlas's ROIs to use. MANDATORY
% cfg2.RCS       = [100 116 132];      %open with crosshairs at specific coordinate, e.g. for Left GPi
% cfg2.winpos    = [0.5 0.5 200 100];  %customize this to your liking (do "get(gca,'Position')" to learn the position
% cfg2.linewidth = 3;                  %thicker crosshair lines
% cfg2.link      = false;              %don't link this figures cursor position to other windows
% cfg2.roifile   = atlas_file;         %ROI probability volumes
% cfg2.roiset    = filename            %open a saved ROI display settings file (saved from an earlier run of showVol)

%% running showVol
% shows examples for different atlas versions -- the PRI naming hasn't been quite standardized yet. For now, just inspect the 
%   contents of the PRI variable to see what is available and what it is called...
if contains(atlasVersion, 'ABCD1')
    showVol(atlas_T1('ABCD1'), thalamus_T1, PRI.ABCD1.T2, PRI.ABCD1.aseg_rgb, PRI.ABCD1.fiber_rgb, ...
      PRI.ABCD1.CO, PRI.ABCD1.FOD, cfg2);
end
    
if contains(atlasVersion,'ABCD2')
    showVol(atlas_T1('ABCD2'), atlas_T1('ABCD1'),thalamus_T1, PRI.ABCD2.aparcaseg_rgb, PRI.ABCD2.fiber_rgb, PRI.ABCD2.CO, ...
       PRI.ABCD2.FOD, PRI.ABCD2.FOD_p80, cfg2);
end

if contains(atlasVersion, '5.0_ABCD3')
    showVol(atlas_T1('5.0_ABCD3'), atlas_T1('5.0', true), atlas_T1('ABCD2',true), atlas_T1('ABCD1',true),thalamus_T1,...
      PRI.ABCD3_50.aseg_rgb, PRI.ABCD3_50.aparcaseg_rgb, PRI.ABCD3_50.fiber_rgb, PRI.ABCD3_50.CO, cfg2);
end

% TO DO: standardize variable and volume names across atlas versions
% fix T1/C0 mixup for ABCD2 in prepareAtlases



