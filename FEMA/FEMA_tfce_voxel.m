function tfcestat = FEMA_tfce_voxel(D,data)

% Define default "delta h". Can also be specified by the user
dh = max(D(:))/100;
E = 0.5;
H = 2.0;
conn = 26;

tfcestat = zeros(size(D));
for h = dh:dh:max(D(:));
  CC = bwconncomp(D>=h,conn);
  integ = cellfun(@numel,CC.PixelIdxList).^E * h^H;
  for c = 1:CC.NumObjects,
    tfcestat(CC.PixelIdxList{c}) = tfcestat(CC.PixelIdxList{c}) + integ(c);
  end
end

%Correct for delta h
tfcestat = tfcestat * dh;

