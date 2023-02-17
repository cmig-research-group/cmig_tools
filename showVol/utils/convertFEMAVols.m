function vols = convertFEMAVols(vol_stat, measure, stat_name, colnames_model, index, limits, pval, cmap, bg_vol, CSFprob, method, atlasVersion, rgb)
%convertFEMAVols Take output from FEMA pipeline and convert to 1mm volumes usable with new showVol
%
%    vols = convertFEMAVols(vol_stat, measure, statname, colnames_model, index, limits, pval?, cmap, bg_vol, CSFprob, method, atlasVersion)
%
% example: z_vols = convertFEMAVols(vol_z, 'N0', 'tstat', colnames_model, 1:5, true, [-10 10 5]); 
%
%  vol_stat & 
%  colnames_model   output from FEMA_wrapper (stacked 2mm stats volumes for each covariate and covariate names)
%  measure          dependent measure name (e.g. 'N0')
%  stat_name        name of statistic (e.g. 'tstat' or 'beta', etc)
%
% optional:
%  index            indices of factor maps to display (not specified or empty -> all)
%  limits           use this uniform stat limits for all plots. (not given or empty: use max abs for each covariate.)
%                    [min max] or to add a 'transparency threshold': [min max threshold] (see Get_volOverlay2)
%                   further, if a scalar < 1, make this proportion of max transparent
%  pval             default is false (will give raw stat), if pval=true, will overlay signed log of pvalue
%  cmap             use your own RGB or RGBA colormap, which will override the default red/black/blue map. 
%                    RGB map will use threshold to determine transparency. RGBA map has
%                    additional column for transparency (1=opaque, 0=transparent, showing the background volume)
%  bg_vol           background volume for stats overlay. If unspecified or empty, defaults 
%                    to standard skull-stripped ABCD T1 as returned by
%                    atlas_T1(true). Can be volume structure or filename
%  CSFprob          mask CSF regsions with prob > CSFprob; if not specified
%                    or [], no masking
%  method           by default uses linear interpoltion, 'nearest' will use the nearest neighbor
%  atlasVersion     ABCD1 or ABCD2
%  rgb              optional, if not empty then will return statmap with RGB
%                    values corresponding to the three dimensions in vol_stat (size(vol_stat,
%                    4) must be 3)
%
% returns a cell array of vols, which can be passed to showVol. In the template atlas geometry (Fall, 2020):
%
%   showVol(vols,...)

if ~exist('atlasVersion','var')
  atlasVersion ='ABCD1';
  disp('Defaulting to ABCD1 atlas. set atlasVersion argument to ''ABCD2'' to use it instead.')
end

if ~exist('index','var') || isempty(index)
  index = 1:size(vol_stat, 4);
end
colorize = true;
if ~exist('limits','var')
  limits = [];
end
if ~exist('pval','var')
  pval = false;
  cmap = redblackblue;
end
if ~exist('cmap','var') || isempty(cmap)
  cmap = mmil_cmap_blueblackred;
end
if ~exist('bg_vol','var') || isempty(bg_vol)
  bg_vol = atlas_T1(atlasVersion,true); %default background image (brain-only T1)
else
  if ischar(bg_vol)
    bg_vol = load(bg_vol);
  end   
end
if ~exist('CSFprob','var') || isempty(CSFprob)
  CSFprob = []; %no masking
end
if ~exist('method','var') || isempty(method)
  method = 'linear';
end

template_vol = atlas_T1(atlasVersion); %used only to get the geometry metadata right
keyboard;
if exist('rgb')
    assert(length(rgb)==3,'Argument rgb must have length 3.');
    vols = template_vol;

    switch method
        case 'linear'
          statimgs_r = upsample_volume(vol_stat(:,:,:,1));
          statimgs_g = upsample_volume(vol_stat(:,:,:,2));
          statimgs_b = upsample_volume(vol_stat(:,:,:,3));
        case 'nearest'
          statimgs_r = upsample_volume_nearest(vol_stat(:,:,:,1));
          statimgs_g = upsample_volume_nearest(vol_stat(:,:,:,2));
          statimgs_b = upsample_volume_nearest(vol_stat(:,:,:,3));
        otherwise
          error('unknown method')
    end
    %mask it
    if ~isempty(CSFprob)
        statimgs_r = maskCSF(statimgs_r, CSFprob, atlasVersion);
        statimgs_g = maskCSF(statimgs_g, CSFprob, atlasVersion);
        statimgs_b = maskCSF(statimgs_b, CSFprob, atlasVersion);
    end
    
    imgs_r = maskBrain(statimgs_r, atlasVersion);
    imgs_g = maskBrain(statimgs_g, atlasVersion);
    imgs_b = maskBrain(statimgs_b, atlasVersion);


    vols.imgs = cat(4, imgs_r, imgs_g, imgs_b);

    formula = [measure ' ~ ' strjoin(rgb, '_')];
    vols.name = sprintf('%s [RGB]',formula);
    vols.statname = strjoin(rgb, '_');

else

    i=1;
    for c = 1:length(index)
      iC = index(c);
      vols{i} = template_vol;
      %statimgs = imresize3(vol_stat(:,:,:,iC), 2);%,'nearest'); %make nearest an option
      switch method
        case 'linear'
          statimgs = upsample_volume(vol_stat(:,:,:,iC));
        case 'nearest'
          statimgs = upsample_volume_nearest(vol_stat(:,:,:,iC));
        otherwise
          error('unknown method')
      end
      %mask it
      if ~isempty(CSFprob)
        statimgs = maskCSF(statimgs, CSFprob, atlasVersion);
      end
      vols{i}.imgs = maskBrain(statimgs, atlasVersion);
      
      
      formula = [measure ' ~ ' colnames_model{iC}];
      vols{i}.name = sprintf('%s [%s]',formula,stat_name);
      vols{i}.statname = stat_name;
      
      thresh = [];
      if isempty(limits) || numel(limits)==1
        statmax = nanmax(abs(statimgs(:)));
        statmin = -statmax;
        thislimits = [statmin statmax];
        if numel(limits)==1 && limits>0 && limits<1 %single element gives transparency _proportion_ relative to max
          thislimits(3) = thislimits(2)*limits;
        end
      else
        thislimits = limits;
      end
      
      if numel(thislimits) == 3
        thresh = thislimits(3);
      end
      
      vols{i}.limits = thislimits;
      vols{i}.colormap = gray(128);
      
      %add a colorized stat map as well?
      if nargin > 5 && colorize==true
        i = i+1;
        [vol_stat_rgb, cmap] = Get_volOverlay2(statimgs, thislimits, bg_vol, pval, cmap);
        vols{i} = template_vol;
        vols{i}.imgs = vol_stat_rgb.imgs;
        if ~isempty(thresh)
          threshStr = sprintf(', > %s', num2str(thresh,2));
        else
          threshStr = '';
        end
        vols{i}.limits = thislimits;
        if pval % pvals, so convert color limits. wrinkle is some Inf values
          vols{i}.limits = -sign(vols{i}.limits).*log10(2*normcdf(-abs(vols{i}.limits)));
          if any(isinf(vols{i}.limits))
            vols{i}.limits = [-322  322]; %beyond z=38.4 normcdf returns 0
          end
          statStr = 'log p';
        else
          statStr = stat_name;
        end
        vols{i}.colormap = cmap;
        vols{i}.statname = statStr;
        vols{i}.name = sprintf('%s [%s +/- %s%s]',formula, statStr, num2str(vols{i}.limits(2),3), threshStr);
    
      end
      
      i = i+1;
    
    end
end





