% convert rcs coords of one volume to another. Will need to round to be used as image indices
function [r2, c2, s2] = rcs2rcs(M1, M2, M2_size, r, c, s)

if isequal(M1, M2)
  r2 = r; c2 = c; s2 = s;
  return
end

%FIXME: is this really the simplest way?
vox1 = {r c s};
vox2 = {};
for i = 1:3
  vox = zeros(4,length(vox1{i}));
  vox(i,:) = vox1{i};
  vox(4,:) = 1;
  
  lph = M1 * vox;
  lph(4,:) = 1;
  
  vox = inv(M2)*lph;
  vox2{i} = vox(i,:);
end
[r2, c2, s2] = deal(vox2{:});

%clamp to edges of vol2
r2 = min(M2_size(1), max(1, r2));
c2 = min(M2_size(2), max(1, c2));
s2 = min(M2_size(3), max(1, s2));

