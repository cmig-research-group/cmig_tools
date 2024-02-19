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
% As of Dec 2023, default atlas is now ABCD3_cor10, but ABCD2 and ABCD1 are options

if nargin < 1 || isempty(atlas)
  atlas = 'ABCD3_cor10';
  fprintf(2, '%s using default %s atlas\n',mfilename,atlas)
end

switch atlas
  case {'ABCD1_cor10' 'ABCD1' 'ABCD2_cor10' 'ABCD2'}
    
    center1mm = [100 100 130];
    M =  ... %new 1mm voxel space
      [0  1  0 -center1mm(2); ...
       0  0 -1  center1mm(3);...
      -1  0  0  center1mm(1); ...
       0  0  0  1];
     
  case {'ABCD3_cor10' 'ABCD3'  } %FIXME: Dec 23, this is tentative, but empirically the center has changed by 1mm
    
    center1mm = [99 99 129];
    M =  ... %new 1mm voxel space
      [0  1  0 -center1mm(2); ...
       0  0 -1  center1mm(3);...
      -1  0  0  center1mm(1); ...
       0  0  0  1];
    
  otherwise
      error('Atlas ''%s'' not recognized', atlas)
end
