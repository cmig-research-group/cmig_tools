function [vol_rgb, cmap] = vol_color_overlay(vol_ovl,cm,crange,vol_bg_rgb)
%
%   vol_rgb = vol_color_overlay(vol_ovl, cm, crange, vol_bg_rgb)
%
% 10/20: JRI: if crange has 3rd value [tmin tmax tthresh], use it as a threshold for displaying stat results only > abs(crange(3))
% colormap can have 4th, alpha column
%
%   [vol_rgb, cmap] = vol_color_overlay(vol_ovl, cm, crange, vol_bg_rgb)
%     optionally return the colormap

if size(cm,2) == 4
  alpha = cm(:,4);
  cm = cm(:,1:3);
else
  alpha = ones(size(cm,1),1);
end

vol_ovl.imgs(isnan(vol_ovl.imgs(:))) = 0; %replace any masked regions with 0 (this presumes that 0 will be rendered as transparent)

vol_ind = max(1,min(size(cm,1),interp1(crange(1:2),[1 size(cm,1)],vol_ovl.imgs,'linear','extrap')));

%use threshold to render stats around zero as totally transparent
if length(crange)>2
  zero_ind = floor(interp1(crange(1:2),[1 size(cm,1)],abs(crange(3))*[-1 1],'linear','extrap'));
  alpha(zero_ind(1):zero_ind(2)) = 0;
end

vol_ovl_rgb = vol_ovl;
vol_ovl_rgb.imgs = reshape(cm(round(vol_ind(:)),:),[size(vol_ind) 3]);

if exist('vol_bg_rgb','var') && ~isempty(vol_bg_rgb)
  if ~isempty(alpha)
    vol_alpha = reshape(alpha(round(vol_ind(:))), size(vol_ind) );
    vol_alpha = repmat(vol_alpha, [1 1 1 3]);
  else
    vol_alpha = repmat(max(vol_ovl_rgb.imgs,[],4),[1 1 1 3]);
  end
  vol_rgb = vol_bg_rgb;
  vol_rgb.imgs = (1-vol_alpha).*vol_bg_rgb.imgs + vol_alpha.*vol_ovl_rgb.imgs;
else
  vol_rgb = vol_ovl_rgb;
end

if nargout > 1
  cmap = cat(2,cm,alpha);
end

%showVol(vol_rgb)

