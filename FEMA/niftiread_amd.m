function [vol M] = niftiread_amd(fname_nii)

vol = niftiread(fname_nii);

if nargout>1
  info = niftiinfo(fname_nii);
  M = info.Transform.T';
  M(:,4) = M*[-1 -1 -1 1]';
end

