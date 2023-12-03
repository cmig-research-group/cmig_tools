function img = maskBrain(img, atlas)
% maskBrain mask skull and background
%
%   img = maskBrain(img, atlas)
%
% ` atlas is ABCD1_cor10 or ABCD2_cor10 (or just ABCD1 or ABCD2) - ABCD1 by default if not specified
%
%   if input img is [] just returns the mask
%
% demo, comparing brain masks for atlas versions:
%   showVol(atlas_T1,atlas_T1('ABCD1',true),atlas_T1('ABCD2'),atlas_T1('ABCD2',true),maskBrain([],'ABCD1'),maskBrain([],'ABCD2'))
%
% TODO: check that this method of finding aseg_mask matches
% this caches the loaded brainMask

persistent brainMask
persistent lastAtlas

if ~exist('atlas','var') || isempty(atlas)
  atlas = 'ABCD1_cor10';
  disp('maskBrain; Using ABCD1 atlas by default. Specify atlas if you want something else.')
end

sameAtlas = strncmp(atlas, lastAtlas, 5);
if sameAtlas && ~isempty(brainMask), return, end %mask is cached

%cached load
if isempty(brainMask) || ~sameAtlas
  cfg = abcdConfig('showVol');
  
  lastAtlas = atlas;
  
  switch atlas
    case {'ABCD1_cor10','ABCD1'}
      
      load(fullfile(cfg.data.showVolData,'atlas_dspace','rgbsnap_ABCD1_cor10.mat'),'vol_mask_aseg');
      brainMask = vol_mask_aseg > 0.01;
      
    case {'ABCD2_cor10','ABCD2'}
      
      load(fullfile(cfg.data.showVolData,'Atlas',['showVolAtlases_ABCD2_cor10.mat']),'aseg')
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