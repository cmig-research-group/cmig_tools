function vol_out = upsample_volume_nearest(vol_in,ds1,ds2,ds3)

if ~exist('ds1','var'), ds1 = 2; end
if ~exist('ds2','var'), ds2 = 2; end
if ~exist('ds3','var'), ds3 = 2; end

ssfacts = [ds1 ds2 ds3];

if isstruct(vol_in)
  [vol M_in] = ctx_ctx2mgh(vol_in);
else
  vol = vol_in;
end

nvox = prod(size(vol));
dim1 = size(vol,1);
dim2 = size(vol,2);
dim3 = size(vol,3);
indvec = [1:nvox]';
[i1,i2,i3] = ind2sub(size(vol),indvec);
%vol_us = zeros(dim1*ds1,dim2*ds2,dim3*ds3,'single');
vol_us = zeros([dim1*ds1 dim2*ds2 dim3*ds3 size(vol,4)],class(vol));
vol_us(1+([1:dim1]-1)*ds1,1+([1:dim2]-1)*ds2,1+([1:dim3]-1)*ds3,:) = vol;
for i1 = 1:(dim1)
  for d1 = 1:(ds1-1)
    ii1 = 1+(i1-1)*ds1;
    if i1<dim1
      f = d1/ds1;
      %fill with nearest original voxel
      if f <= 0.5
        vol_us(ii1+d1,:,:,:) = vol_us(ii1,:,:,:);
      else
        vol_us(ii1+d1,:,:,:) = vol_us(ii1+ds1,:,:,:);
      end
    else
      vol_us(ii1+d1,:,:,:) = vol_us(ii1,:,:,:);
    end
  end
end
for i2 = 1:(dim2)
  for d2 = 1:(ds2-1)
    ii2 = 1+(i2-1)*ds2;
    if i2<dim2
      f = d2/ds2;
      %fill with nearest original voxel
      if f <= 0.5
        vol_us(:,ii2+d2,:,:) = vol_us(:,ii2,:,:);
      else
        vol_us(:,ii2+d2,:,:) = vol_us(:,ii2+ds2,:,:);
      end
    else
      vol_us(:,ii2+d2,:,:) = vol_us(:,ii2,:,:);
    end
  end
end
for i3 = 1:(dim3)
  for d3 = 1:(ds3-1)
    ii3 = 1+(i3-1)*ds3;
    if i3<dim3
      f = d3/ds3;
      %fill with nearest original voxel
      if f <= 0.5
        vol_us(:,:,ii3+d3,:) = vol_us(:,:,ii3,:);
      else
        vol_us(:,:,ii3+d3,:) = vol_us(:,:,ii3+ds3,:);
      end
    else
      vol_us(:,:,ii3+d3,:) = vol_us(:,:,ii3,:);
    end
  end
end

if isstruct(vol_in)
  M_out = M_in; M_out(1:3,1:3) = M_in(1:3,1:3)*diag(1./ssfacts);
  r0_in = M_in*[1 1 1 1]'; r0_out = M_out*[1 1 1 1]';
  M_out(1:3,4) =  M_out(1:3,4) - (r0_out(1:3)-r0_in(1:3)); % Adjust M_out to make coordinate of first voxel consistent
  vol_out = ctx_mgh2ctx(vol_us,M_out);
else
  vol_out = vol_us;
end

