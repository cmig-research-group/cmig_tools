function vols_ctx = ctx_mgh2ctx(vols,M)

% Updated to handle n-d volumes

if ~exist('M','var'), M = eye(4); end

dims = size(vols);
width = dims(1); height=dims(2); depth=dims(3); nvol=prod(dims(4:end)); 

[xsize,ysize,zsize,x_r,x_a,x_s,y_r,y_a,y_s,z_r,z_a,z_s,c_r,c_a,c_s] = mat2mgh(M,width,height,depth);

vols_ctx = struct;
vols_ctx.imgs = double(vols);
vols_ctx.Mvxl2lph = M_RAS_TO_LPH * M;
vols_ctx.dimc = height; % size(vol,2);
vols_ctx.dimr = width; % size(vol,1);
vols_ctx.dimd = depth; % size(vol,3);
vols_ctx.vx = xsize;
vols_ctx.vy = ysize;
vols_ctx.vz = zsize;
vols_ctx.lphcent=[-c_r;-c_a;c_s];
vols_ctx.DirCol=[-x_r;-x_a;x_s];
vols_ctx.DirRow=[-y_r;-y_a;y_s];
vols_ctx.DirDep=[-z_r;-z_a;z_s];

return

% Old code

[width,height,depth,nvol] = size(vols) ;

[xsize,ysize,zsize,x_r,x_a,x_s,y_r,y_a,y_s,z_r,z_a,z_s,c_r,c_a,c_s] = mat2mgh(M,width,height,depth);

vols_ctx = cell(1,nvol);
for vi = 1:nvol
  vol = vols(:,:,:,vi);
  vol_ctx.imgs = double(vol);
  vol_ctx.Mvxl2lph = M_RAS_TO_LPH * M;
  vol_ctx.dimc = size(vol,2);
  vol_ctx.dimr = size(vol,1);
  vol_ctx.dimd = size(vol,3);
  vol_ctx.vx = xsize;
  vol_ctx.vy = ysize;
  vol_ctx.vz = zsize;
  vol_ctx.lphcent=[-c_r;-c_a;c_s];
  vol_ctx.DirCol=[-x_r;-x_a;x_s];
  vol_ctx.DirRow=[-y_r;-y_a;y_s];
  vol_ctx.DirDep=[-z_r;-z_a;z_s];
  if nvol>1
    vols_ctx{vi} = vol_ctx;
  else
    vols_ctx = vol_ctx;
  end
end

% Should modify this sso it can take a cell array of 3-D matrices in vols

