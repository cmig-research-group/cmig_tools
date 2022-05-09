function vol_T1 = atlas_T1(atlas, mask)
% atlas_T1 return the T1 in ABCD 1mm atlas space
%
%  Convenience function. Info is also available in PRI
%
%  vol_T1 = atlas_T1(atlas, mask)
%
% optional input:
%   mask:     if true, return skull stripped brain-only T1
%
% demo:
%     showVol(atlas_T1,atlas_T1('ABCD1',true),atlas_T1('ABCD2'),atlas_T1('ABCD2',true))

if ~exist('mask','var')
  mask = false;
end

if ~exist('atlas','var') || isempty(atlas)
  atlas = 'ABCD1_cor10';
  disp('atlas_T1; Using ABCD1 atlas by default. Specify atlas if you want something else.')
end

cfg = abcdConfig('showVol');

switch atlas
  case {'ABCD1','ABCD1_cor10'}
    pre_rendered_path = fullfile(cfg.data.showVolData,'atlas_dspace','rgbsnap_ABCD1_cor10.mat');
    
    load(pre_rendered_path,'vol_T1_mu','M_atl');
    
    %convert LPH center to 0-based centering in default ABCD atlas space
    Mvxl2lph = M_RAS_TO_LPH * M_atl; %M_atl is in RAS, apparently
    Mvxl2lph(1,4) = -100;
    Mvxl2lph(2,4) = 130;
    Mvxl2lph(3,4) = 100;
    
    vol_T1 = ctx_mgh2ctx(vol_T1_mu, M_LPH_TO_RAS*Mvxl2lph);
    vol_T1.name = 'Atlas T1';
    
  case {'ABCD2_cor10','ABCD2'}
    fname = fullfile(cfg.data.showVolData, 'Atlas','T1_ABCD2_cor10.mat');
    tmp = load(fname);
    vol_T1 = tmp.T1_ABCD2_cor10;
    
  otherwise
    error('Invalid atlas value (%s)',atlas)
end

%mask it
if mask
  vol_T1.imgs = maskBrain(vol_T1.imgs, atlas);
end
