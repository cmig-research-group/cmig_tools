function [vol_rgb_ctx, colormap] = convertVolColor(vol, colormap, varargin)
%
%	 vol_rgb = vol_color_overlay(vol, colormap, crange, vol_bg_rgb)
%
% 10/20: JRI: if crange has 3rd value [tmin tmax tthresh], use it as a threshold for displaying stat results only > abs(crange(3))
% colormap can have 4th, alpha column
%
%	 [vol_rgb, colormap] = vol_color_overlay(vol, colormap, crange, vol_bg_rgb)
%		 optionally return the colormap

	p = inputParser;
	addParamValue(p, 'crange', [min(vol, [], 'all') max(vol, [], 'all')]);
	addParamValue(p, 'cmid', []);
	addParamValue(p, 'thresh', []);
	addParamValue(p, 'M', eye(4));
	addParamValue(p, 'name', []); 

	parse(p, varargin{:})
	crange = p.Results.crange;
	cmid = p.Results.cmid;
	thresh = p.Results.thresh;
	M = p.Results.M; 
	name = p.Results.name; 

	if size(colormap,2) == 4
		alpha = colormap(:,4);
		colormap = colormap(:,1:3);
	else
		alpha = ones(size(colormap,1),1);
	end

	% !! DO I WANT THIS? !!
	%%replace any masked regions with 0 (this presumes that 0 will be rendered as transparent)
	%vol(isnan(vol(:))) = 0;

	if isempty(cmid)
		vol_ind = max(1, min(size(colormap,1), interp1(crange(1:2),[1 size(colormap,1)], vol,'linear','extrap')));
	else 
		idx_pos =  find(vol>=cmid);
		idx_neg = find(vol<cmid); 

		vol_ind = nan(size(vol));
	 
		vol_ind(idx_neg) = max(1, min(size(colormap,1), interp1([crange(1), cmid], [1 size(colormap,1)/2], vol(idx_neg),'linear','extrap')));
		vol_ind(idx_pos) = max(1, min(size(colormap,1), interp1([cmid, crange(2)], [size(colormap,1)/2 size(colormap,1)], vol(idx_pos),'linear','extrap')));
	end 
	vol_rgb = reshape(colormap(round(vol_ind(:)),:),[size(vol_ind) 3]);	

	vol_rgb_ctx = ctx_mgh2ctx(vol_rgb, M); 
	vol_rgb_ctx.name = name; 
	vol_rgb_ctx.limits = crange;
	vol_rgb_ctx.colormap = colormap; 
	vol_rgb_ctx.statname = '';
	if ~isempty(cmid)
		vol_rgb_ctx.cmid = cmid; 
	end 

	%use threshold to render stats around zero as totally transparent
	if ~isempty(thresh)
		zero_ind = floor(interp1(crange(1:2), [1 size(colormap,1)], abs(thresh)*[-1 1],'linear','extrap'));
		alpha(zero_ind(1):zero_ind(2)) = 0;
	end

	if exist('vol_bg_rgb','var') && ~isempty(vol_bg_rgb)
		if ~isempty(alpha)
			vol_alpha = reshape(alpha(round(vol_ind(:))), size(vol_ind) );
			vol_alpha = repmat(vol_alpha, [1 1 1 3]);
		else
			vol_alpha = repmat(max(vol_rgb_ctx.imgs,[],4),[1 1 1 3]);
		end
		vol_rgb = vol_bg_rgb;
		vol_rgb_ctx.imgs = (1-vol_alpha).*vol_bg_rgb.imgs + vol_alpha.*vol_rgb_ctx.imgs;
	else
		vol_rgb_ctx = vol_rgb_ctx;
	end

	if nargout > 1
		colormap = cat(2,colormap,alpha);
	end

	%showVol(vol_rgb)
end 
