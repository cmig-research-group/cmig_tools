function FaceVertexCData = SurfView_CData_new(vertvals,fvals,cmap,curvvals,polarity,curvcontrast)
% This software is Copyright (c) 2022 The Regents of the University of California. All Rights Reserved.
% See LICENSE.
if ~exist('polarity','var'), polarity = 2; end
if isempty(fvals)
  fmax = max(abs(vertvals)); fmin = fmax/10; fmid = (fmin+fmax)/2;
else
  fmin = fvals(1); fmid = fvals(2); fmax = fvals(3);
end
if ~exist('curvcontrast','var') | isempty(curvcontrast), curvcontrast = [0.2 0.2]; end
if min(size(vertvals))==2
  if size(vertvals,1)==2
    vertvals = vertvals';
  end
  wvec = vertvals(:,2);
  vertvals = vertvals(:,1);
elseif size(vertvals,2)~=3
  vertvals = colvec(vertvals);
end
if ~exist('curvvals','var') | isempty(curvvals), curvvals = -ones(size(vertvals)); end
tmp1 = curvcontrast(1)*(1-curvcontrast(2)*colvec(sign(curvvals)))*ones(1,3);
fvec = zeros(size(vertvals));

if size(vertvals,2)==3
    FaceVertexCData = vertvals;
else    
    if polarity == 2
        ivec = find(abs(vertvals)>fmin&abs(vertvals)<=fmid); fvec(ivec) = sign(vertvals(ivec)).*(0.5*(abs(vertvals(ivec))-fmin)/(fmid-fmin));
        ivec = find(abs(vertvals)>fmid); fvec(ivec) = sign(vertvals(ivec)).*(0.5+0.5*(abs(vertvals(ivec))-fmid)/(fmax-fmid));
        fvec = max(-1,min(1,fvec));
        ci = round((1+size(cmap,1))/2);
        indvec = ci+0.5*fvec*(size(cmap,1)-1);
    else
      ivec = find(vertvals<=fmid); fvec(ivec) = 0.5*(vertvals(ivec)-fmin)/(fmid-fmin);
      ivec = find(vertvals>fmid); fvec(ivec) = 0.5+0.5*(vertvals(ivec)-fmid)/(fmax-fmid);
      fvec = max(0,min(1,fvec));
      indvec = fvec*(size(cmap,1)-1);
    end

    %indvec = max(1,min(size(cmap,1),indvec));
    tmp2 = interp1(cmap,indvec,'linear','extrap');
    if ~exist('wvec','var'), wvec = 2*min(abs(fvec),0.5); end
    wmat = repmat(wvec,[1 3]);
    FaceVertexCData = (1-wmat).*tmp1 + wmat.*tmp2;
end
