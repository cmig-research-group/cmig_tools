function vol_out = subsample_volume(vol_in,ds1,ds2,ds3)

if ~exist('ds1','var'), ds1 = 2; end
if ~exist('ds2','var'), ds2 = 2; end
if ~exist('ds3','var'), ds3 = 2; end

ssfacts = [ds1 ds2 ds3];

if isstruct(vol_in)
  if ~all(ssfacts==1)
    [tmp_in M_in] = ctx_ctx2mgh(vol_in);
    tmp_out = tmp_in(1:ssfacts(1):end,1:ssfacts(2):end,1:ssfacts(3):end,:);
    M_out = M_in; M_out(1:3,1:3) = M_in(1:3,1:3)*diag(ssfacts); 
    r0_in = M_in*[1 1 1 1]'; r0_out = M_out*[1 1 1 1]';
    M_out(1:3,4) =  M_out(1:3,4) - (r0_out(1:3)-r0_in(1:3)); % Adjust M_out to make coordinate of first voxel consistent
    vol_out = ctx_mgh2ctx(tmp_out,M_out);
  else % No subsampling needed
    vol_out = vol_in;
  end
else
  vol_out = vol_in(1:ssfacts(1):end,1:ssfacts(2):end,1:ssfacts(3):end,:);
end

