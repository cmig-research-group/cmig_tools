function vol_T1 = atlas_T1(atlas, mask)
% atlas_T1 return the T1 in ABCD 1mm atlas space
%
%  Convenience function. T1 is also available in PRI
%
%  vol_T1 = atlas_T1(atlas, mask)
%
%  atlas is one of 3.0_ABCD1_cor10, 4.0_ABCD2_cor10 or 5.0_ABCD3_cor10 (or unambiguous abbreviations, e.g. 5.0 or ABCD2: see validateAtlasVersion.m)
%
% optional input (default true):
%   mask:     if false, return full T1, else returns skull-stripped T1 by default
%
% demo:
%     showVol(atlas_T1('ABCD1',true),atlas_T1('ABCD2',true),atlas_T1('5.0_ABCD3',true))

if ~exist('mask','var')
  mask = true;
end

if ~exist('atlas','var') || isempty(atlas)
  error('atlas_T1: must specify an atlas version including release e.g. 5.0_ABCD3_cor10. See validateAtlasVersion.m')
end

atlas = validateAtlasVersion(atlas); %expand to canonical form

cfg = abcdConfig('showVol');

switch atlas
  case 'ABCD1_cor10' %we did things differently back then...
    pre_rendered_path = fullfile(cfg.data.showVolData,'atlas_dspace','rgbsnap_ABCD1_cor10.mat');
    
    load(pre_rendered_path,'vol_T1_mu','M_atl');
    
    %convert LPH center to 0-based centering in default ABCD atlas space
    Mvxl2lph = M_RAS_TO_LPH * M_atl; %M_atl is in RAS, apparently
    Mvxl2lph(1,4) = -100;
    Mvxl2lph(2,4) = 130;
    Mvxl2lph(3,4) = 100;
    
    vol_T1 = ctx_mgh2ctx(vol_T1_mu, M_LPH_TO_RAS*Mvxl2lph);
    vol_T1.name = 'Atlas T1';
    
  %case {'ABCD2_cor10','5.0_ABCD3_cor10', '6.0_ABCD3_cor10'} %this pattern should work for anything from ABCD2_cor10 onwards
  otherwise
    fname = fullfile(cfg.data.showVolData, 'Atlas',atlas,['T1_' atlas '.mat']);
    tmp = load(fname);
    atlasFieldname = atlas(atlas~='.');  %convert to valid field name
    vol_T1 = tmp.(['T1_' atlasFieldname ]);
end

%mask it
if mask
  vol_T1.imgs = maskBrain(vol_T1.imgs, atlas);
end
