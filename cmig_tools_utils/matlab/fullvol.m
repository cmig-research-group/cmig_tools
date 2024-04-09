function vol = ABCD_fullvol(volvec,vol_mask)

ivec_mask = find(vol_mask>=0.5);

if size(volvec,1)>1 
    nframes = size(volvec,1);
else
    nframes = length(volvec)/length(ivec_mask);
end

vol = zeros([size(vol_mask) nframes],class(volvec));
for fi = 1:nframes
  vol_tmp = zeros([size(vol_mask) 1],class(volvec));
  if size(volvec,1)>1
      vol_tmp(ivec_mask) = volvec(fi,:);
  else
      vol_tmp(ivec_mask) = volvec((fi-1)*length(ivec_mask)+[1:length(ivec_mask)]);
  end
  vol(:,:,:,fi) = vol_tmp;
end

