function M = Mvxl2lph_atlas(atlas)
% Mvxl2lph_atlas return the canonical Mvxl2lph transform for current atlas
%
% (Note, the M_atl field in many prerendered atlases is different: 
%   a) converts from voxel to RAS
%   b) has a center that is off by 1 mm
%
% by default, atlas is ABCD1_cor10; will add others over time

if nargin < 1
  atlas = 'ABCD1_cor10';
end

switch atlas
  case {'ABCD1_cor10' 'ABCD1' 'ABCD2_cor10' 'ABCD2'}
    
    center1mm = [100 100 130];
    M =  ... %new 1mm voxel space
      [0  1  0 -center1mm(2); ...
       0  0 -1  center1mm(3);...
      -1  0  0  center1mm(1); ...
       0  0  0  1];
    
  otherwise
      error('Atlas ''%s'' not recognized', atlas)
end
