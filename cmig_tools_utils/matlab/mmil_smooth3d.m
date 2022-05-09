function vol_out = mmil_smooth3d(vol_in,sig1,sig2,sig3)
%function vol_out = mmil_smooth3d(vol_in,sig1,sig2,sig3)
%
% Created:   03/09/10 by Don Hagler based on code from Anders Dale
% Rcnt Mod:  05/21/12 by Don Hagler
% Last Mod:  09/14/12 by Don Hagler
%

if ~mmil_check_nargs(nargin,2), return; end;
if ~exist('sig2','var') || isempty(sig2), sig2 = sig1; end;
if ~exist('sig3','var') || isempty(sig3), sig3 = sig1; end;


orig_type = class(vol_in);
if ~strcmp(orig_type,'single')
  vol_in = single(vol_in);
end;

kdim = size(vol_in,1);
idim = size(vol_in,2);
jdim = size(vol_in,3);
filtvol = 0;
tmp = zeros(kdim,idim,jdim);
for k=1:kdim
  tmp(k,:,:) = k-(kdim/2+1);
end
filtvol = filtvol+(tmp/kdim*sig1).^2;
for i=1:idim
  tmp(:,i,:) = i-(idim/2+1);
end
filtvol = filtvol+(tmp/idim*sig2).^2;
for j=1:jdim
  tmp(:,:,j) = j-(jdim/2+1);
end
filtvol = filtvol+(tmp/jdim*sig3).^2;
filtvol = exp(-1/2*filtvol);
clear('tmp');
vol_out = fftshift(ifftn(vol_in));
vol_out = fftn(ifftshift(vol_out.*filtvol));

vol_out = real(vol_out);
if strcmp(orig_type,'double')
  vol_out = double(vol_out);
end;


