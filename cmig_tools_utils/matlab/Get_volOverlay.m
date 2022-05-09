function [vols_logp_rgb, cmap] = Get_volOverlay(vol_tstat, tmin, tmax, vol_T1, tthresh, stat)

% This is a wrapper to get the blended volume overlays
% Require input:
%   vol_tstat: the statistical map generated as the voxelwise analysis
%              assuming it is already masked and has consistent orientation as vol_T1 (100 x 120 x 120)
%   tmax: the maximun value of tstats for the coloring map
%   tmin: the minimun value of tstats for the coloring map
%         NOTE -- recommend to do symmetric tmax and tmin if the stat is signed. 
%   vol_T1: the mean T1 in atlas space. It is oriented as 100 x 120 x 120 now
%           can load from the pre-saved volumn .mat
%   tthresh: optional, if specified suppresses -stats with |value| < tthresh
%     note, this needs a new version of vol_color_overlay from the showVol
%     repository.
%    stat: optional, default is false (will give lop pvals), if stat=true,
%       will output color overlay of original statistic
%
%   output: vols_logp_rgb, opyionally the colormap cmap

if nargin < 4
  vol_T1 = atlas_T1; %load default ABCD T1
  stat=false;
end

fprintf('%s -- %s.m: Generating the T1 overlays \r\n', datestr(now), mfilename);
vol_T1_rgb = vol_color_overlay(vol_T1,gray,[0 100]);

fprintf('%s -- %s.m: Overlaying the statistic map, using color range %.2f to %.2f \r\n', datestr(now), mfilename, tmin, tmax);
tmp = ctx_mgh2ctx(vol_tstat, eye(4)); 
tmp.imgs(~isfinite(tmp.imgs)) = 0; 
if stat==false
tmp.imgs = -sign(tmp.imgs).*log10(2*normcdf(-abs(tmp.imgs)));  
end
range = [tmin tmax];
if nargin > 4 && ~isempty(tthresh)
  range = [tmin tmax tthresh];
end
[vols_logp_rgb, cm] = vol_color_overlay(tmp, mmil_cmap_blueblackred(), range, vol_T1_rgb);


if nargout > 1
  cmap = cm;
end