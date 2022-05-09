function vol = ABCD_fullvol(volvec,vol_mask)

ivec_mask = find(vol_mask>=0.5);
nframes = length(volvec)/length(ivec_mask);

vol = zeros([size(vol_mask) nframes],class(volvec));
for fi = 1:nframes
  vol_tmp = zeros([size(vol_mask) 1],class(volvec));
  vol_tmp(ivec_mask) = volvec((fi-1)*length(ivec_mask)+[1:length(ivec_mask)]);
  vol(:,:,:,fi) = vol_tmp;
end

