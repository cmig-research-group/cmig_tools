function prepareAtlases_ABCD3_cor10
% prepareAtlases_ABCD3_cor10 prepare atlases for use with showVol
%
% Dec 2023, relevant ROI files are now in 
%   /space/syn65/1/data/abcd-sync/5.0/imaging_concat/voxelwise/roi (atlases) and /smri (T1, etc)
%
% ShowVolData still in
%   /space/amdale/1/tmp/ABCD_cache/ABCD-Atlas/showVolData
%
% atlas_dspace link in
%   /space/amdale/1/tmp/ABCD_cache/ABCD-Atlas/showVolData/atlas_dspace
%  to
%   /home/dale/atlas_dspace/atlas_dspace_ABCD3_cor10.mat 
%   (also seems to be present at /space/amdale/1/data/atlas_dspace/ -- identical?)
%
% no atlas_dspace_ABCD3_cor10_FOD/
%
% OUTPUTS
%   outputs to cfg.data.showVolData (/space/amdale/1/tmp/ABCD_cache/ABCD-Atlas/Atlas/showVolData)
%   ROIs: showVolAtlases_ABCD2_cor10.mat
%   rgb_<atlasName>_atlas,for atlasName in {aseg, aparc, fiber, thalamus, pauli} as of May 2021
%   T1_,T2_, C0 as available"
%

% optimization:remove all-zero volumes and compress volumes
% this is for ABCD2, building on atlas_dspace_ABCD2_cor10 and using
%   code from TractQuant/TractQuant_ABCD_makejobs.m to generate rgb and labelprobs data structures
%   relies on fs_colorlut, now in /utils

% 9/28/2021: return to practice of normalized ROIs for Pauli/Thalamus atlases (which are at base binary ROIs)
% 12/2/2023: update for ABCD3

persistent AD

timestamp = datestr(now);
atlasVersion = 'ABCD3_cor10';
useAtlasDspace = true; % previously for T1, T2 etc we used 'means' calculated by Diliana--for now take directly from atlas_dspace

% KEY: these are the new preprocessed data from Diliana's pipeline, versioned by date
processDate = '2021-09-05';
source = ['/space/amdale/1/tmp/ABCD_cache/ABCD-Atlas/showVolDataSource_' processDate '/'];
meandir = 'thresh_mean0.80_s0.50/';

processDate = '2023/04/02'; %most recent imaging_concat/voxelwise/roi modification date. (atlas_dspace is 2022-01-02)
source = '/space/syn65/1/data/abcd-sync/5.0/imaging_concat/voxelwise/roi';
sourceAnat = '/space/syn65/1/data/abcd-sync/5.0/imaging_concat/voxelwise/smri';
meandir = '';

cfg = abcdConfig('showVol');
outdir = fullfile(cfg.data.showVolData, 'Atlas','');
%outdir = fullfile('~/Images', 'Atlas','') %for TESTING


% FIXME: thalamaus and pauli atlases for ABCD3?
atlases = {'aseg','aparcaseg','fiber'} %, 'thalamus','pauli'}; %aseg+aparc in a single volume, but we'll split for conveniences
atlas_dspace = '/space/amdale/1/tmp/ABCD_cache/ABCD-Atlas/showVolData/atlas_dspace';

files = {fullfile(source,meandir,'volmat_aseg_prob200x200x260.mat'), ...
    fullfile(source,meandir,'volmat_aparcaseg_prob200x200x260.mat'), ...
    fullfile(source,meandir,'volmat_fiber_prob200x200x260.mat'), ...
    fullfile(source,'MNI_Atlases','Thalamus/ABCD2_cor10_snap.mat'), ...
    fullfile(source,'MNI_Atlases','Pauli_SubcortNuclei/ABCD2_cor10_snap.mat')};
  
notes = {['ABCD3_cor10 ASEG ROIs, ' processDate], ...
  ['ABCD3_cor10 APARCASEG ROIs, ' processDate], ...
  ['ABCD3_cor10 FIBER ROIs, ' processDate], ...
  ['ABCD3 Thalamus atlas (add citation), ' processDate], ...
  ['ABCD3 Pauli atlas (add citation), ' processDate]};

doNormalize = [0 0 0 1 1];

overlayColors = {[1 0 1], [1 0 .5], [0 1 .5],[],[]};
defaultThreshold = [0.7 0.7 0.25 0.5 0.5];

atlasFname = fullfile(outdir, ['showVolAtlases_' atlasVersion '.mat']);

%process aseg and fibers (in same file)
if isempty(AD)
  disp('Loading atlas_dspace aseg, aparc, wmparc & fiber ROIs')
  load(fullfile(atlas_dspace,'atlas_dspace_ABCD3_cor10.mat'));
  AD = atlas_dspace; %for convenience
else
  disp('Using Cached atlas_dspace')
end

% volume anatomical parameters
Mvxl2lph = M_RAS_TO_LPH * AD.M_atl; %M_atl is in RAS

%shift to actual anatomical center of new 1mm atlas
center_1mm = [100 100 130];
center_1mm = center_1mm - 1; % ** Emprically, the atlas_dspace T1 has its center at 99 99 129. This does not match M_atl--check w/ A&D
Mvxl2lph(1,4) = -center_1mm(2);
Mvxl2lph(2,4) = center_1mm(3);
Mvxl2lph(3,4) = center_1mm(1);
voxSize = 1;

%% ROI Probabilities %%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% loop over ABCD ROIs
for iA = 1:3
  atlas = atlases{iA};
  disp(atlas)
  
  A = load(files{iA});
  
  %remove empty volumes (normalization no longer needed)
   try
    probMax = max(A.volmat_prob,[],1:3);
  catch
    probMax = max(max(max(A.volmat_prob,[],1),[],2),[],3);
  end
  goodvols = squeeze(probMax)>0;
  %A.volmat_prob = A.volmat_prob ./ probMax;
  A.roiabbreviations = [];
  
  %remove derivative & redundant ROIs
  
  if strcmp(atlas,'fiber')
    goodvols(38:end)=false;
    disp('Removing allFibers ROIs from fiber atlas')
    %remove derivative ROIs:
    % ROIs 36-37 are R/L frontal fornix--keep now
    % 38-42 are all fibers, R/L all fibers (no cc), and R/L all fibers

    %FIXME: note the fiber has 115 rgb but only 35 valid ROIs--looks like aparacaseg colors by mistake?
    %[roicodes_fiber,roinames_fiber,rgbv_fiber] = fs_colorlut(fullfile(fundir, 'MMIL_FSColorLUT.txt'));
    %[~,~,ind_fiber] = intersect(AD.indlist_fiber+10000, roicodes_fiber); %index of filist into roicodes
    %roicolors       = rgbv_fiber(ind_fiber,1:3);
    %A.roirgb = roicolors(:,1:35);
    % NB roirgb and roicolors don't actually match
  end
  
  %remove all of the aseg rois from this. In principle, does fitting jointly change them from the aseg-only data?
  %  checking the probMax for aseg ROIs, they are identical in the aseg and aparcaseg atlases, from
  %  which we assume the probs are the same (except for CC_All and Ventricles_All...)
  if strcmp(atlas,'aparcaseg')
    ctxRois = contains(A.roinames,'ctx-');
    goodvols(~ctxRois) = false;
  end
  
  %remove empty/unneded volumes
  if sum(~goodvols)
    fprintf(' Removing %d empty/unneeded vols:\n',sum(~goodvols))
    disp(A.roinames(~goodvols(:)))
    A.volmat_prob = A.volmat_prob(:,:,:,goodvols);
    A.roicodes    = A.roicodes(goodvols);
    A.roinames    = A.roinames(goodvols);
    A.roirgb      = A.roirgb(goodvols,:);
  end  
  
  %re-sort ROIs, add abbreviations
  if strcmp(atlas,'aparcaseg')
    disp('re-sorting aparc cortical ROIs')
    aparcAbbrev = {'bSTS','cACC','cMFC','Cun','Ent','Fus','IPC','ITC','isthCing','lOcc','lOFC',...
      'Ling','mOFC','mTem','paraHipp','paraCen','POperc','POrb','PTri','periCalc','postCen','PCC',...
      'preCen','preCun','rACC','rMFC','SFC','SPC','STC','SMar','FP','TP','tTem','Ins'};
    roiabbreviations = [aparcAbbrev aparcAbbrev];
    order = [1:34;35:68];
    order = order(:);
    A.volmat_prob = A.volmat_prob(:,:,:,order);
    A.roicodes    = A.roicodes(order);
    A.roinames    = A.roinames(order);
    A.roirgb      = A.roirgb(order,:);
    A.roiabbreviations = roiabbreviations(order);
  end
    

  %generate roinames and roicodes for aseg, aparcaseg and wmparc
  % fundir  = fileparts(which(mfilename));
  % [roicodes_fs,roinames_fs,rgbv_fs] = fs_colorlut(fullfile(fundir,'FreeSurferColorLUT.txt'));
  % [~,~,ind_aseg]      = intersect(AD.indlist_aseg, roicodes_fs); %index of filist into roicodes (goodvols only)
  % [~,~,ind_aparcaseg] = intersect(AD.indlist_aparcaseg, roicodes_fs); %index of filist into roicodes
  % [~,~,ind_wmparc]    = intersect(AD.indlist_wmparc, roicodes_fs); %index of filist into roicodes
  %
  
  %set up anat struct
  % TODO
  % Using a constant color for each atlas, rather than individually coloring as we do for Thalamus and Pauli.
  %   This can be made an option later?
  % make abbreviations for standard ROIs
 
  atl=[];
  atl.overlayColor = overlayColors{iA};
  atl.prob = single(A.volmat_prob);
  atl.roicodes = A.roicodes;
  atl.roinames = A.roinames;
  atl.roiabbreviations = A.roiabbreviations;
  atl.roicolors = A.roirgb;
  atl.subjlist = A.subjlist;
  atl.size = size(atl.prob(:,:,:,1));
  atl.roi_text = [];
  atl.file = files{iA};
  atl.probabilityThreshold = defaultThreshold(iA);
  atl.showNames = false;
  atl.outlineColor = 'w';
  
  %make a list of ROI names for use UI for highlighting ROIs
  atl.uiNames = {'None'};
  atl.uiRoiIdx = 0;
  for iV = 1:size(atl.prob,4)
    atl.uiNames{end+1} = atl.roinames{iV};
    atl.uiRoiIdx(end+1) = iV;
  end
  %for UI state
  atl.uiAtlasName = [atlas ' [' atlasVersion ']'];
  atl.uiRoiOverlayIdx = []; %indices into roi volume of selected ROIs to overlay
  atl.uiRoiOverlayImg = [];
  atl.uiRoiOverlaySelected = []; %indices into uiNames
  
  %add additional shape info
  atl.Mvxl2lph = Mvxl2lph;
  atl.voxelSize = voxSize;
  atl.dimr=size(atl.prob,1); atl.dimc = size(atl.prob,2); atl.dimd = size(atl.prob,3);
  atl.vx=atl.voxelSize; atl.vy=atl.voxelSize; atl.vz=atl.voxelSize;
  atl.rcscent = center_1mm(:);
  lphcent = atl.Mvxl2lph * [atl.rcscent(:); 1];
  atl.lphcent = lphcent(1:3);
  
  atl.datePrepared = timestamp;
  atl.notes = notes{iA};
  atl.atlasVersion = atlasVersion;
  
  anat.(atlas) = atl;
  
  % write out RGB image
  rgbFname = [files{iA}(1:end-4) '_rgbsum.mat'];
  tmp = load(rgbFname);
  rgb_vol = ctx_mgh2ctx(single(zeros(atl.size)), M_LPH_TO_RAS*Mvxl2lph); %starter volume
  rgb_vol.imgs = tmp.volmat_prob_rgbsum;
  
  %FIXME: fiber rgbsum is 'overexposed', but this method doesn't work--look into
  if 0 && strcmp(atlas,'fiber')
    prob = anat.(atlas).prob;
    c = anat.(atlas).roicolors / 255;
    rgb_sum = 0 * rgb_vol.imgs; %zero it out
    %add colored ROIs to image
    for iR = 1:size(prob,4)
      roiprob = prob(:,:,:,iR);
      roi_rgb = cat(4,c(iR,1)*roiprob, c(iR,2)*roiprob, c(iR,3)*roiprob);
      rgb_sum = rgb_sum + roi_rgb;
    end
    rgb_sum(rgb_sum>1) = 1;
    rgb_vol.imgs = rgb_sum;
  end
  
  %add additional shape info
  rgb_vol.Mvxl2lph = Mvxl2lph;
  rgb_vol.voxelSize = voxSize;
  rgb_vol.dimr=size(atl.prob,1); atl.dimc = size(atl.prob,2); atl.dimd = size(atl.prob,3);
  rgb_vol.vx=atl.voxelSize; atl.vy=atl.voxelSize; atl.vz=atl.voxelSize;
  rgb_vol.rcscent = center_1mm(:);
  lphcent = atl.Mvxl2lph * [atl.rcscent(:); 1];
  rgb_vol.lphcent = lphcent(1:3);
  
  rgb_vol.name = [atlas ' rgb [' atlasVersion ']'];
  rgb_vol.datePrepared = timestamp;
  rgb_vol.notes = [notes{iA} ' RGB'];
  rgb_vol.atlasVersion = atlasVersion;
  varname = [atlas '_rgb_' atlasVersion];
  eval([varname ' = rgb_vol;']); %assign to atlas-specific name
  save(fullfile(outdir,[atlas '_rgb_' atlasVersion '.mat']), varname)
  
  
end %loop over ABCD ROI atlases

%% ABCD PRERENDERED IMAGES - T1, T2, C0 %%
images = {'T1','T2','FA', 'CO'};

if useAtlasDspace
  
  disp(' save T1 (atlas_dspace)')
  clear vol
  vol.imgs = AD.muvols_T1;
  vol.Mvxl2lph = Mvxl2lph;
  vol.dimc = 200; vol.dimr=200; vol.dimd=260;
  vol.vx=1;vol.vy=1;vol.vz=1;
  vol.lphcent=[0 0 0]';
  vol.atlasVersion = atlasVersion;
  vol.name = 'T1 [ABCD3]';
  T1_ABCD3_cor10 = vol;
  save(fullfile(outdir,['T1_' atlasVersion '.mat']), 'T1_ABCD3_cor10')
  
  disp(' save T2 (atlas_dspace)')
  clear vol
  vol.imgs = AD.muvols_T2;
  vol.Mvxl2lph = Mvxl2lph;
  vol.dimc = 200; vol.dimr=200; vol.dimd=260;
  vol.vx=1;vol.vy=1;vol.vz=1;
  vol.lphcent=[0 0 0]';
  vol.atlasVersion = atlasVersion;
  vol.name = 'T2 [ABCD3]';
  T2_ABCD3_cor10 = vol;
  save(fullfile(outdir,['T2_' atlasVersion '.mat']), 'T2_ABCD3_cor10')
  
else
  
  %T1
  disp(' save T1')
  fname = fullfile(sourceAnat,'volmat_nu_mean200x200x260.mat');
  tmp = load(fname);
  vol = ctx_mgh2ctx(single(tmp.vol_mean), M_LPH_TO_RAS*Mvxl2lph);
  vol.lphcent = [0 0 0]'; %have to do this, because we are using 0-based M, while the function assumes 1-based
  vol.atlasVersion = atlasVersion;
  vol.name = 'T1 [ABCD3]';
  T1_ABCD3_cor10 = vol;
  save(fullfile(outdir,['T1_' atlasVersion '.mat']), 'T1_ABCD3_cor10')
  
  %T2
  disp(' save T2')
  fname = fullfile(sourceAnat,'volmat_T2_mean200x200x260.mat');
  tmp = load(fname);
  vol = ctx_mgh2ctx(single(tmp.vol_mean), M_LPH_TO_RAS*Mvxl2lph);
  vol.lphcent = [0 0 0]'; %have to do this, because we are using 0-based M, while the function assumes 1-based
  vol.atlasVersion = atlasVersion;
  vol.name = 'T2 [ABCD3]';
  T2_ABCD3_cor10  = vol;
  save(fullfile(outdir,['T2_' atlasVersion '.mat']), 'T2_ABCD3_cor10')
  
  %FA
  disp(' save FA')
  fname = fullfile(sourceAnat,'volmat_FA_mean200x200x260.mat');
  tmp = load(fname);
  vol = ctx_mgh2ctx(single(tmp.vol_mean), M_LPH_TO_RAS*Mvxl2lph);
  vol.lphcent = [0 0 0]'; %have to do this, because we are using 0-based M, while the function assumes 1-based
  vol.atlasVersion = atlasVersion;
  vol.name = 'FA [ABCD3]';
  FA_ABCD3_cor10  = vol;
  save(fullfile(outdir,['FA_' atlasVersion '.mat']), 'FA_ABCD3_cor10')
  
end

%CO
disp(' save CO (atlas_dspace)')
vol = AD.vol_CO_rgb;
vol.Mvxl2lph = Mvxl2lph;
vol.atlasVersion = atlasVersion;
vol.name = 'diffusion direction [ABCD3]';
CO_ABCD3_cor10 = vol;
save(fullfile(outdir,['CO_' atlasVersion '.mat']), 'CO_ABCD3_cor10')


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Thalamus & Pauli (Atlases imported from other studies)
%hand-code abbreviations for Thalamus
thalamusAbbrev = {'PUL','A','MD','VLD','CllpmPUL','VA','VLV','PUL','A','MD','VLD','CllpmPUL','VA','VLV'}';

for iA = 4:5
  if ~exist(files{iA},'file'), continue, end

    atlas = atlases{iA};
    disp(atlas)
    A = load(files{iA});
    nroi = size(A.vol_labels_reg,4);
    
    prob = single(A.vol_labels_reg);
    if doNormalize(iA)
      disp('Normalizing Probabilities')
      prob = prob ./ max(prob,[],1:3);
    end
    
    anat.(atlas).prob = prob;
    try
        name = A.T.Name;
    catch
        name = A.T.LabelName;
    end
    anat.(atlas).roinames = name;
    if strcmp(atlas,'thalamus')
        anat.(atlas).roiabbreviations = thalamusAbbrev;
    else
        anat.(atlas).roiabbreviations = A.T.Abbreviation;
    end
    anat.(atlas).roicolors = [A.T.R A.T.G A.T.B];
    m = max(anat.(atlas).roicolors(:));
    if m > 1
        anat.(atlas).roicolors = anat.(atlas).roicolors / 255; %rescale to matlab [0 1]
    end
    anat.(atlas).size = size(anat.(atlas).prob(:,:,:,1));
    anat.(atlas).roi_text = []; %for text display of current ROI names in figure
    anat.(atlas).file = files{iA};
    anat.(atlas).overlayColor = []; %unused
    anat.(atlas).probabilityThreshold = 0.5;
    anat.(atlas).showNames = true;

    anat.(atlas).uiNames = {'None'}; anat.(atlas).uiRoiIdx = 0;
    for iV = 1:nroi
        anat.(atlas).uiNames{end+1} = sprintf('%s (%s)',anat.(atlas).roinames{iV}, anat.(atlas).roiabbreviations{iV}) ;
        anat.(atlas).uiRoiIdx(end+1) = iV;
    end
    anat.(atlas).uiRoiOverlayIdx = []; %indices into roi volume of selected ROIs to overlay
    anat.(atlas).uiRoiOverlayImg = [];
    anat.(atlas).uiRoiOverlaySelected = []; %indices into uiNames
    
    %add additional shape info
    anat.(atlas).Mvxl2lph = Mvxl2lph;
    anat.(atlas).voxelSize = voxSize;
    anat.(atlas).dimr=size(anat.(atlas).prob,1); anat.(atlas).dimc = size(anat.(atlas).prob,2); anat.(atlas).dimd = size(anat.(atlas).prob,3);
    anat.(atlas).vx=anat.(atlas).voxelSize; anat.(atlas).vy=anat.(atlas).voxelSize; anat.(atlas).vz=anat.(atlas).voxelSize;
    anat.(atlas).rcscent = center_1mm(:);
    lphcent = anat.(atlas).Mvxl2lph * [anat.(atlas).rcscent(:); 1];
    anat.(atlas).lphcent = lphcent(1:3);
    
    anat.(atlas).notes = notes{iA};
    anat.(atlas).dataPrepared = timestamp;
    anat.(atlas).atlasVersion = atlasVersion;
    
    %% also save some high-res T1 and T2 images and rgb atlas if present
    if isfield(A,'vol_T1_reg')
        disp(' save T1')
        vol = ctx_mgh2ctx(single(A.vol_T1_reg), M_LPH_TO_RAS*Mvxl2lph);
        vol.lphcent = [0 0 0]'; %have to do this, because we are using 0-based M, while the function assumes 1-based
        vol.atlasVersion = atlasVersion;
        vol.name = [atlas ' atlas T1 [ABCD3]'];
        varname = [atlas '_T1'];
        eval([varname ' = vol;']); %assign to atlas-specific name
        save(fullfile(outdir,[atlas '_T1_' atlasVersion '.mat']), varname);
    end
    if isfield(A,'vol_T2_reg')
        disp(' save T2')
        vol = ctx_mgh2ctx(single(A.vol_T2_reg), M_LPH_TO_RAS*Mvxl2lph);
        vol.lphcent = [0 0 0]';
        vol.atlasVersion = atlasVersion;
        vol.name = [atlas ' atlas T2 [ABCD3]'];
        varname = [atlas '_T2'];
        eval([varname ' = vol;']); %assign to atlas-specific name
        save(fullfile(outdir,[atlas '_T2_' atlasVersion '.mat']), varname);
        
    end
    if isfield(A,'vol_labels_rgb')
        disp(' save rgb')
        vol = A.vol_labels_rgb;
        vol.Mvxl2lph = Mvxl2lph;
        vol.lphcent = [0 0 0]';
        vol.atlasVersion = atlasVersion;
        vol.name = [atlas ' atlas rgb [ABCD3]'];
        varname = [atlas '_rgb'];
        eval([varname ' = vol;']); %assign to atlas-specific name
        save(fullfile(outdir,[atlas '_rgb_' atlasVersion '.mat']),varname);
    end
     if isfield(A,'vol_labels_T1_rgb')
        disp(' save T1 + rgb')
        vol = A.vol_labels_T1_rgb;
        vol.Mvxl2lph = Mvxl2lph;
        vol.lphcent = [0 0 0]';
        vol_labels_rgb.atlasVersion = atlasVersion;
        vol_labels_rgb.name = [atlas ' atlas rgb+T1 [ABCD3]'];
        varname = [atlas '_rgb_T1'];
        eval([varname ' = vol;']); %assign to atlas-specific name
        save(fullfile(outdir,[atlas '_rgb_T1_' atlasVersion '.mat']),varname);
     end
     if isfield(A,'vol_labels_T2_rgb')
        disp(' save T2 + rgb')
        vol = A.vol_labels_T2_rgb;
        vol.Mvxl2lph = Mvxl2lph;
        vol.lphcent = [0 0 0]';
        vol_labels_rgb.atlasVersion = atlasVersion;
        vol_labels_rgb.name = [atlas ' atlas rgb+T2 [ABCD2]'];
        varname = [atlas '_rgb_T2'];
        eval([varname ' = vol;']); %assign to atlas-specific name
        save(fullfile(outdir,[atlas '_rgb_T2_' atlasVersion '.mat']),varname);
    end
    
end

%% compress probabilities and unpack atlases into individual variables (so we can load as needed)
% FIXME: could also save with '-struct' argument in newer matlab
for iA = 1:length(atlases)
  fprintf('%s: ', atlases{iA})
  anat.(atlases{iA}).prob = compressVol(anat.(atlases{iA}).prob);
  anat.(atlases{iA}).atlasVersion = atlasVersion;
  eval([atlases{iA} ' = anat.' atlases{iA} ';'])
end

%% set up ROI overlay defaults
params.showOverlay = true; %toggle for overlay
params.binarizeOverlay = false; %entire ROI > threshold rendered same color, false: show probability
params.showOutline = true; %outline ROIs
params.outlineColor = [1 1 1];
params.outlineWidth = 1;
params.overlayAlpha = 0.4;
params.imageFade = 1; %start with no fade
params.showRoiLabel = true;
params.currentAtlas = '';
params.generatedDate = timestamp;
params.notes = ['ABCD3_cor10, '  processDate];

%% SAVEt
disp('save atlases')
save(atlasFname,'atlases','params',atlases{:},'-v7.3')

%% create symmetrized atlases
% TODO

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Helper Functions

%% compression/expansion
%simples compression method: only store non-zero values. since ROIs are very sparse
%for the most part, the space savings (thus load time improvement) is
%large. E.g. aseg only 1.4% of voxels are non-zero
function c = compressVol(vol)
s = size(vol);
if prod(s) >= 2^32, error('Volume too large to compress with int32 indexes'), end
nz = cast(find(vol(:)~=0),'uint32');

c = {s nz(:) single(vol(nz))};

fprintf('(%s%% compression)\n',num2str(100*2*length(nz)/prod(s),2)); %factor of 2 since we save index and values

% companion expansion function when loading atlas
%  e.g. 
% if iscell(aseg.prob)
%     aseg.prob = expandVol(aseg.prob);
% end
function vol = expandVol(c)
vol = zeros(c{1},'single');
vol(c{2}) = c{3};
