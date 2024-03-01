function img = maskBrain(img, atlas)
% maskBrain mask skull and background
%
%   img = maskBrain(img, atlas)
%
%  atlas is one of 3.0_ABCD1_cor10, 4.0_ABCD2_cor10 or 5.0_ABCD3_cor10 (or unambiguous abbreviations, e.g. 5.0 or ABCD2: see validateAtlasVersion.m)
%
%   if input img is [] just returns the mask
%
% mask is computed by summing all aseg rois and taking any voxel with summed probability > 0.01%
%
% demo, compare brain masks for atlas versions:
%   showVol(maskBrain([],'ABCD1'),maskBrain([],'ABCD2'),maskBrain([],'5.0_ABCD3'))
%
% TODO: check that this method matches vol_mask_aseg
% this caches the loaded brainMask

persistent brainMask
persistent lastAtlas

if ~exist('atlas','var') || isempty(atlas)
    error('maskBrain: must specify an atlas version including release e.g. 5.0_ABCD3_cor10. See validateAtlasVersion.m')
end

atlas = validateAtlasVersion(atlas); %expand to canonical form

sameAtlas = strncmp(atlas, lastAtlas, 5);

%cached load
if isempty(brainMask) || ~sameAtlas
  cfg = abcdConfig('showVol');
  
  lastAtlas = atlas;
  
  switch atlas
    case '3.0_ABCD1_cor10'
      
      load(fullfile(cfg.data.showVolData,'atlas_dspace','rgbsnap_ABCD1_cor10.mat'),'vol_mask_aseg');
      brainMask = vol_mask_aseg > 0.01;
      
    case '4.0_ABCD2_cor10'
      
      load(fullfile(cfg.data.showVolData,'Atlas',['showVolAtlases_ABCD2_cor10.mat']),'aseg')
      c = aseg.prob;
      if iscell(c) %decompress
        aseg.prob = zeros(c{1}, 'single');
        aseg.prob(c{2}) = c{3};
      end
      brainMask = sum(aseg.prob, 4);
      brainMask = brainMask > 0.01;
      
    case '5.0_ABCD3_cor10'
      
      disp('Loading ABCD3 aseg to cache brainMask')
      load(fullfile(cfg.data.showVolData,'Atlas',atlas,['showVolAtlases_5.0_ABCD3_cor10.mat']),'aseg')
      c = aseg.prob;
      if iscell(c) %decompress
        aseg.prob = zeros(c{1}, 'single');
        aseg.prob(c{2}) = c{3};
      end
      brainMask = sum(aseg.prob, 4);
      brainMask = brainMask > 0.01;
      
  end
end

if ~exist('img','var') || isempty(img)
  img = brainMask;
else
  img(~brainMask) = 0;
end