
%% ------------------------------------------------------------
%% --- add ROI overlays to image and/or compute outlines
%% ------------------------------------------------------------
function [ima, outline] = anatomyAddRoiOverlay(handles, ima, rr, cc, ss, mode)
% overlay probability for selected ROIs & calculate contour outlines
%mode is 'fill' (default) or 'outline' or 'filloutline' for both
% default is fill for lowres images and outline for FOD.

%return if nothing to do
if ~handles.hasABCDBrain || ~handles.anat.params.showOverlay
  outline = [];
  return
end

atlases = handles.anat.atlases;
shown = false;
for iA = 1:length(atlases)
  atlas = atlases{iA};
  if isfield(handles.anat,atlas) && ~isempty(handles.anat.(atlas).uiRoiOverlayImg)
    shown = true;
    break
  end
end
if ~shown
  outline = [];
  return
end

%add outline based on checkbox
if strcmp(mode,'fill') && handles.anat.params.showOutline
  mode = 'filloutline';
end

%for scalar volumes, convert to RGB if filling ROI, but preserve ability to adjust brightness and contrast
% FIXME: one inconvenience is that the contrast/darkness sliders don't change the image look when
% ROIs are displayed. Need to switch ROI off then back on to update it
if size(ima,3) == 1 && contains(mode,'fill')
  handles = updateDisplay_clims(handles);
  CLim = get(gca,'CLim');
  ima = (ima - CLim(1))/(CLim(2) - CLim(1));
  ima(isnan(ima)) = 0;
  ima = repmat(ima, 1, 1, 3);
end

isCoronal = length(ss)==1;
isSagittal = length(cc)==1;
isAxial = length(rr)==1; %x/y are transposed in axial view

outlineCoord = cell(0);
outlineColor = [];
outlineLabel = cell(0,3);
mask = [];
roiImg = [];

M_img = handles.vols{handles.currentVol}.Mvxl2lph;

%loop over atlases
for iA = 1:length(atlases)
  atlas = atlases{iA};
  if ~isfield(handles.anat, atlas), continue; end
  
  A = handles.anat.(atlas);
  if isempty(A.uiRoiOverlayImg); continue; end
  
  aimg = A.uiRoiOverlayImg;
  M_atlas = A.Mvxl2lph;
  
  %resample atlas into current image coordinates
  [ra, ca, sa] = rcs2rcs(M_img, M_atlas, A.size, rr, cc, ss); %image rcs to atlas rcs
  ra = floor(ra); ca = floor(ca); sa = floor(sa); % convert to matrix indices
  
  %collapse rois into single rgb image & mask
  if contains(mode, 'fill')
    
    overlay = squeeze(aimg(ra, ca, sa, :));
    if isAxial, overlay = permute(overlay, [2 1 3]); end
        
    if isempty(mask)
      mask = ( sum(overlay,3) > 0 );
      roiImg = zeros(size(overlay,1), size(overlay,2));
    else
      mask = (mask | sum(overlay,3) > 0);
    end
    for iR = 1:size(overlay,3)
      o = overlay(:,:,iR);
      if any(o(:))
        if handles.anat.params.binarizeOverlay
          o(o>A.probabilityThreshold) = 1;
        end
        c = A.roicolors(A.uiRoiOverlayIdx(iR),:);
        %c = A.roicolors(iR,:); %somehow this recreates colors in Anders 5/23/21 email ??
        m = max(c(:)); %FIXME: temporary fix, will be fixed permanentyly in prepareAtlases.m
        if m > 1
          c = c/255;
        end
        roiImg = roiImg + repmat(o,[1 1 3]) .* permute(c, [3 1 2]);
      end
    end
  end
  
  %make outlines
  if contains(mode, 'outline')
    
    %keep full atlas slice for outline calculation, at original scales
    if isCoronal
      overlayUnthresholded = squeeze(aimg(:, :, sa, :));
    elseif isSagittal
      overlayUnthresholded = squeeze(aimg(:, ca, :, :));
    elseif isAxial
      overlayUnthresholded = squeeze(aimg(ra, :, :, :));
      overlayUnthresholded = permute(overlayUnthresholded, [2 1 3]);
    end
    
    level = A.probabilityThreshold;
    for iR = 1:size(overlayUnthresholded,3)
      o = overlayUnthresholded(:,:,iR);
      if ~any(o(:)), continue, end
      cmat = contourc(double(o),[level level]);
      if isempty(cmat), continue, end
      cmat(:,1) = [];
      breaks = find(cmat(1,:)==level);
      cmat(:,breaks) = nan;
      
      % NOTE: We are using atlas coords for axis now, even for FODs so no
      % need to scale up outline coords any longer
      %convert cmat to coords of current image, from atlas coords
%       if isCoronal
%         [r2, c2, ~] = rcs2rcs(M_atlas, M_img, [size(ima) ss], cmat(1,:), cmat(2,:),ss);
%         cmat = [r2; c2];
%       elseif isSagittal
%         [r2, ~, s2] = rcs2rcs(M_atlas, M_img, [size(ima,1) cc size(ima,2)], cmat(1,:), cc, cmat(2,:));
%         cmat = [r2; s2];
%       elseif isAxial
%         [~, c2, s2] = rcs2rcs(M_atlas, M_img, [rr size(ima)], rr, cmat(1,:), cmat(2,:));
%         cmat = [c2; s2];
%       end
%       cmat(:,breaks) = nan;
      
      if isempty(breaks), breaks = size(cmat,2);end
      center = nanmean(cmat(:,1:breaks(1)-1),2);
      outlineCoord{end+1} = cmat;
      if contains(mode,'fill') && handles.anat.params.overlayAlpha>0
        c = [1 1 1]; %when showing colored roiprob, use white outline FIXME: make this UI-selectable
      else
        c = A.roicolors(A.uiRoiOverlayIdx(iR),:);
        m = max(c(:)); %FIXME: temporary fix, will be fixed permanentyly in prepareAtlases.m
        if m > 1
          c = c/255;
        end
        c = min(1, [.75 .75 .75] + c); %desaturate towards white
      end
      outlineColor = cat(1, outlineColor, c);
      
      %if we have an abbreviation, add it as a label if showNames is true
      if A.showNames && ~isempty(A.roiabbreviations)
        label = A.roiabbreviations(A.uiRoiOverlayIdx(iR));
      else
        label = '';
      end
      outlineLabel(end+1,:) = {center(1) center(2) label};
      
    end %loop over ROIs
  end %if outline
  
end %loop over atlas

%this was nice, but not sure needed any more
%         %outline _only_ add a second contour at 50%, useful in FOD display,
%         outline50 = [];
%         if ~contains(mode, 'fill') && handles.anat.fiberProbabilityThreshold < 0.4
%             level = 0.5;
%             cmat = [];
%             for iR = 1:size(overlay,3)
%                 cmat = cat(2, cmat, contourc(double(overlayUnthresholded(:,:,iR)),[level level]) );
%             end
%             cmat(:,cmat(1,:)==level) = nan;
%             outline50 = cat(2, outline50, cmat*a2i);
%         end
%         if ~isempty(outline50), outline = {outline outline50};end


%combine original image (possibly faded) with alpha blended roi image
if contains(mode,'fill')
  fade = handles.anat.params.imageFade;
  if isinteger(ima), ima = cast(ima,'single')/255; end
  % experimental: want low probs to be more transparent than higher ones, so that fringes turn more transparent
  fadeProb = 0.15;
  lum = mean(roiImg,3);
  inFringe = (lum<fadeProb);
  shadedMask = cast(mask,'single');
  shadedMask(inFringe) = lum(inFringe)/fadeProb;
  
  %ima = fade * ima.*~mask + (1-handles.anat.params.overlayAlpha) * ima.*mask + handles.anat.params.overlayAlpha * roiImg.*shadedMask;
  ima = ima.*(1-handles.anat.params.overlayAlpha*shadedMask/4) + handles.anat.params.overlayAlpha * roiImg.*shadedMask; %try new method, get rid of 'fade' param
  %ima = ima + handles.anat.params.overlayAlpha * roiImg.*shadedMask; %try new method, get rid of 'fade' param
  ima = min(ima,1); %clamp to 1
end
if contains(mode,'outline')
  outline = {outlineCoord outlineColor outlineLabel};
else
  outline = [];
end

%% ------------------------------------------------------------
