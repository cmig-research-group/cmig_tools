function vol_out = mmil_dilate_mask(vol_in,varargin)
%function vol_out = mmil_dilate_mask(vol_in,[options])
%
% Purpose: fill and dilate a mask volume through smoothing and thresholding
%
% Usage: vol_out = mmil_dilate_mask(vol_in,'key1', value1,...);
%
% Required Input:
%   vol_in: input mask volume (ctx format)
%     Can be [] if 'fname_in' is specified
%
% Optional Input:
%  'fname_in': full path of input mask volume file (e.g. brain.mgz, aseg.mgz)
%    If supplied, vol_in is ignored and fname_in is loaded instead
%    {default: []}
%  'fname_out': full path of output mask volume file
%    If supplied, vol_out is saved as fname_out
%    {default: []}
%  'thresh0': threshold applied to mask before first smoothing step
%    {default: 0.01}
%  'smooth1': smoothing sigma (mm) for initial fill step
%    (set to 0 for no fill)
%    {default: 10}
%  'thresh1': threshold applied to mask after first smoothing step
%    (if smooth1=0, this is ignored)
%    {default: 0.5}
%  'smooth2': smoothing sigma (mm) for second dilation step
%    (set to 0 for no dilation)
%    {default: 30}
%  'thresh2': threshold applied to mask after second smoothing step
%    (if smooth2=0, this is ignored)
%    {default: 0.5}
%  'smooth3': smoothing sigma (mm) for third dilation step
%    (set to 0 for no dilation)
%    {default: 5}
%  'thresh3': threshold applied to mask after third smoothing step
%    sets small values to zero, but does not binarize mask
%    (if smooth3=0, this is ignored)
%    {default: 0.01}
%  'erode_flag': allow mask erosion by not resetting original values
%     after smoothing and thresholding
%    {default: 0}
%  'forceflag': [0|1] whether to overwrite existing fname_out
%    Ignored if fname_out is empty
%    {default: 0}
%
% Output:
%   vol_out: output mask volume (ctx format)
%
% Created:  06/08/09 by Don Hagler
% Last Mod: 02/15/11 by Don Hagler
%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% parse input parameters

if ~mmil_check_nargs(nargin,1), return; end;
parms = mmil_args2parms(varargin, { ...
  'fname_in',[],[],...
  'fname_out',[],[],...
  'thresh0',0.01,[0,1000],...
  'smooth1',10,[0,100],...
  'thresh1',0.5,[0,1],...
  'smooth2',30,[0,100],...
  'thresh2',0.5,[0,1],...
  'smooth3',5,[0,100],...
  'thresh3',0.01,[0,1000],...
  'erode_flag',false,[false true],...
  'forceflag',false,[false true],...
});

if ~isempty(parms.fname_out) &&...
   exist(parms.fname_out,'file') && ~parms.forceflag
  if nargout>0
    vol_out = ctx_load_mgh(parms.fname_out);
  end;
  return;
end;

if ~isempty(parms.fname_in)
  if ~exist(parms.fname_in,'file')
    error('file %s not found',parms.fname_in);
  end;
  [vol_in,mr_parms] = ctx_load_mgh(parms.fname_in);
else
  mr_parms = [];
end;

vol_out = [];

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% select non-zero voxels
volmask = vol_in;
volmask.imgs = 1.0*(vol_in.imgs>parms.thresh0);
volmask.maxI = 1.0;
vol_in = volmask;

% convert smoothing kernels from mm to voxels
ps = sqrt(sum(volmask.Mvxl2lph(1:3,1:3).^2,1)); % voxel dimensions in mm
parms.smoothvec1 = parms.smooth1./ps;
parms.smoothvec2 = parms.smooth2./ps;
parms.smoothvec3 = parms.smooth3./ps;

% smooth
if parms.smooth1>0
  volmask.imgs = ...
    mmil_smooth3d(volmask.imgs,...
      parms.smoothvec1(1),parms.smoothvec1(2),parms.smoothvec1(3));
  % truncate
  volmask.imgs = 1.0*(volmask.imgs>=parms.thresh1);
  % refresh from input
  if ~parms.erode_flag
    volmask.imgs(vol_in.imgs>0)=1;
  end;
end

% smooth
if parms.smooth2>0
  volmask.imgs = ...
    mmil_smooth3d(volmask.imgs,...
      parms.smoothvec2(1),parms.smoothvec2(2),parms.smoothvec2(3));
  % truncate
  volmask.imgs = 1.0*(volmask.imgs>=parms.thresh2);
  % refresh from input
  if ~parms.erode_flag
    volmask.imgs(vol_in.imgs>0)=1;
  end;
end

% smooth
if parms.smooth3>0
  volmask.imgs = max(0,...
    mmil_smooth3d(volmask.imgs,...
      parms.smoothvec3(1),parms.smoothvec3(2),parms.smoothvec3(3)));
  % set very small values to 0
  volmask.imgs(volmask.imgs<parms.thresh3) = 0;
  % refresh from input
  if ~parms.erode_flag
    volmask.imgs(vol_in.imgs>0)=1;
  end;
end;

if ~isempty(parms.fname_out)
  ctx_save_mgh(volmask,parms.fname_out,mr_parms);
end;

if nargout>0
  vol_out = volmask;
end;

