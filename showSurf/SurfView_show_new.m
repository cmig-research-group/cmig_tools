function SurfView_show_new(surf_lh,surf_rh,vertvals_lh,vertvals_rh,fvals,cmap,viewname,includevec,curvvec_lh,curvvec_rh,icsurf,polarity,curvcontrast,bgcol)
%  SurfView_show_new(surf_lh,surf_rh,vertvals_lh,vertvals_rh,fvals,cmap,viewname,includevec,curvvec_lh,curvvec_rh,icsurf,polarity,curvcontrast,bgcol)
%
%
% This software is Copyright (c) 2022 The Regents of the University of California. All Rights Reserved.
% See LICENSE.

if ~exist('curvvec_lh','var'), curvvec_lh = []; end
if ~exist('curvvec_rh','var'), curvvec_rh = []; end
if ~exist('polarity','var') | isempty(polarity), polarity = 2; end
if ~exist('curvcontrast','var'), curvcontrast = []; end
if ~exist('bgcol','var') | isempty(bgcol), bgcol = [0 0 0]; end % Default to black background

if exist('icsurf','var')
  nverts = size(icsurf.vertices,1);
  nfacess = size(icsurf.faces,1);
  surf_lh.vertices = surf_lh.vertices(1:nverts,:);
  surf_rh.vertices = surf_rh.vertices(1:nverts,:);
  surf_lh.faces = icsurf.faces;
  surf_rh.faces = icsurf.faces;
  if size(vertvals_lh,1)==1 | size(vertvals_lh,1)==3, vertvals_lh = vertvals_lh'; end
  if size(vertvals_rh,1)==1 | size(vertvals_rh,1)==3, vertvals_rh = vertvals_rh'; end
  if length(vertvals_lh)>nverts
    vertvals_lh = vertvals_lh(1:nverts,:);
    vertvals_rh = vertvals_rh(1:nverts,:);
  end
  if ~isempty(curvvec_lh), curvvec_lh = curvvec_lh(1:nverts); end
  if ~isempty(curvvec_rh), curvvec_rh = curvvec_rh(1:nverts); end
end
tmp1 = surf_lh; tmp2 = surf_rh; 
if includevec(1) & ~includevec(2)
  surfstruct = tmp1;
  curvvec = curvvec_lh;
  vertvals = vertvals_lh;
elseif ~includevec(1) & includevec(2)
  surfstruct = tmp2;
  curvvec = curvvec_rh;
  vertvals = vertvals_rh;
elseif includevec(1) & includevec(2)
  surfstruct = SurfView_concat(tmp1,tmp2); 
  curvvec = cat(1,curvvec_lh,curvvec_rh);
  vertvals = cat(1,vertvals_lh,vertvals_rh);
end

FaceVertexCData = SurfView_CData_new(vertvals,fvals,cmap,curvvec,polarity,curvcontrast);

p = patch(surfstruct);
set(p, 'FaceColor', 'interp', 'EdgeColor', 'none', 'LineStyle', 'none', 'SpecularStrength', 0, 'SpecularStrength', 0, 'AmbientStrength', 0.4, 'FaceVertexCData', FaceVertexCData);
set(gca,'Color',bgcol,'XTickLabel',[],'YTickLabel',[],'ZTickLabel',[]);

switch lower(viewname)
  case 'left',             camorbit(-90,0); camorbit(0,-90);
  case 'right',            camorbit(90,0); camorbit(0,-90);
  case 'front',            camorbit(180,0); camorbit(0,-90);
  case {'back','behind'},  camorbit(0,-90);
  case {'top','above'},    camorbit(0,0); 
  case {'bottom','below'}, camorbit(180,0); camorbit(0,180);
end
camlight headlight;
daspect([1 1 1]);
material dull;
lighting phong;



% ToDo
%   - Look for surface mask (get rid of non-cortical regions)
%   - Merge w. ic?.tri files for subsampled data (copy coords from high-res surfaces, keep topology from lo-res) 
