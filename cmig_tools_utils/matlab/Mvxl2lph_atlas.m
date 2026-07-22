function M = Mvxl2lph_atlas(atlas)
% Mvxl2lph_atlas return the canonical Mvxl2lph transform for current atlas
%
%   M = Mvxl2lph_atlas(atlas)
%
% (Note, the M_atl field in many prerendered atlases is different: 
%   a) converts from voxel to RAS
%   b) has a center that is off by 1 mm
%
%   to convert this back to that format: M_LPH_TO_RAS*Mvxl2lph (e.g. to use if calling ctx_mgh2ctx)
%

if nargin < 1 || isempty(atlas)
  error('You must specify an atlas version. See validateAtlasVersion.m for valid choices')
end

atlas = validateAtlasVersion(atlas);

switch atlas
  case {'3.0_ABCD1_cor10', '4.0_ABCD2_cor10'}
    
    rcs_center1mm = [100 100 130];
    M =  ... %new 1mm voxel space
      [0  1  0 -rcs_center1mm(2); ...
       0  0 -1  rcs_center1mm(3);...
      -1  0  0  rcs_center1mm(1); ...
       0  0  0  1];
     
%   case {'5.0_ABCD3_cor10'} %Before Nov 2024, I adjusted the center. 
%     
%     rcs_center1mm = [99 99 129];
%     M =  ... %new 1mm voxel space
%       [0  1  0 -rcs_center1mm(2); ...
%        0  0 -1  rcs_center1mm(3);...
%       -1  0  0  rcs_center1mm(1); ...
%        0  0  0  1];
     
  case {'5.0_ABCD3_cor10', '6.0_ABCD3_cor10'} %FIXME: empirically the center is 100 99 130! 
    
    rcs_center1mm = [101 101 131];
    M =  ... %new 1mm voxel space
      [0  1  0 -rcs_center1mm(2); ...
       0  0 -1  rcs_center1mm(3);...
      -1  0  0  rcs_center1mm(1); ...
       0  0  0  1];
    
  otherwise
      error('Atlas ''%s'' not recognized', atlas)
end
