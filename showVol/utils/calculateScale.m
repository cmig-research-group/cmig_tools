function sf = calculateScale(handles, vol)

%calculate scaling factors, using correct voxel sizes from vol struct
% FIXME as currently written, image size is taken from first image, which is how the r,c,s axes are defined, so
% make scales relative to that.
vol1 = handles.vols{1};
if nargin < 2
  vol = handles.vols{handles.currentVol};
end
try
    sf.r = vol1.vx/vol.vx;
    sf.c = vol1.vy/vol.vy;
    sf.s = vol1.vz/vol.vz;
catch
    sf.r = 1;
    sf.c = 1;
    sf.s = 1;
end
sf.scale = mean([sf.r sf.c sf.s]);

return

% FIXME: what sets rrMax? rrMax gets unreliable with multi-scaled voxels
if isfield(handles.vols{handles.currentVol},'imgs')
  ima = handles.vols{handles.currentVol}.imgs;
  sf.r = size(ima,1)/handles.rrMax;
  sf.c = size(ima,2)/handles.ccMax;
  sf.s = size(ima,3)/handles.ssMax;
elseif isfield(handles.vols{handles.currentVol},'imgs1')
  ima = handles.vols{handles.currentVol}.imgs3{handles.ss};
  sf.r = size(ima,1)/handles.rrMax;
  sf.c = size(ima,2)/handles.ccMax;
  ima = handles.vols{handles.currentVol}.imgs2{handles.cc};
  sf.s = size(ima,2)/handles.ssMax;
else %simply an image volume
  ima = handles;
  sf.r = size(ima,1)/handles.rrMax;
  sf.c = size(ima,2)/handles.ccMax;
  sf.s = size(ima,3)/handles.ssMax;
end
sf.scale = mean([sf.r sf.c sf.s]);
