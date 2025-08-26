function niftiwrite_amd(vol,fname_nii,M,symmetric_cLim)

if ~exist('symmetric_cLim', 'var') || isempty(symmetric_cLim)
    symmetric_cLim = false;
end

if strcmp(fname_nii(max(1,(length(fname_nii)-2)):end),'.gz')
  compressedflag = true;
  fname_nii = fname_nii(1:end-3);
else
  compressedflag = false;
end

T = M'; T(4,:) =  M*[1 1 1 1]'; Transform = struct('Dimensionality',3,'T',T); 

switch ndims(vol) % check if input vol is 3D or 4D
  case 3
    PixelDimensions = sum(M(:,1:3).^2,1).^.5;
  case 4
    PixelDimensions = [sum(M(:,1:3).^2,1).^.5, 1];
end

% Calculate min and max
cal_min = min(vol, [], 'all', 'omitmissing');
cal_max = max(vol, [], 'all', 'omitmissing');

% Handle scenario where the min or max value is Inf
if isinf(cal_max)
    cal_max = max(vol(~isinf(vol)));
end

if isinf(cal_min)
    cal_min = min(vol(~isinf(vol)));
end

% If user wants symmetric colour limits, set to the largest value
if symmetric_cLim
    max_all = max(abs(cal_min), abs(cal_max));
    cal_min = -max_all;
    cal_max = max_all;
end

% Setting the scale and intercept: is this required?
% This is already set as defaults, so no need
% scl_slope     = 1;
% scl_intercept = 0;

% niftiwrite simplifies cal_min and cal_max to DisplayIntensityRange
info = struct('Transform', Transform, 'Datatype', class(vol), 'ImageSize', size(vol), ...
              'Description', '', 'Version', 'NIfTI1', 'Qfactor', 1,                   ...
              'PixelDimensions', PixelDimensions, 'SpaceUnits', 'Millimeter',         ...
              'TimeUnits', 'Second', 'SliceCode', 'Unknown', 'AdditiveOffset', 0,     ...
              'MultiplicativeScaling', 1, 'TimeOffset', 0, 'FrequencyDimension', 0,   ...
              'PhaseDimension', 0, 'SpatialDimension', 0,                             ...
              'DisplayIntensityRange', [cal_min cal_max], 'TransformName','Sform');
niftiwrite(vol,fname_nii,info,'Compressed',compressedflag);