function [xsize,ysize,zsize,x_r,x_a,x_s,y_r,y_a,y_s,z_r,z_a,z_s,c_r,c_a,c_s] = mat2mgh(M,width,height,depth)

% function [xsize,ysize,zsize,x_r,x_a,x_s,y_r,y_a,y_s,z_r,z_a,z_s,c_r,c_a,c_s] = mat2mgh(M,width,height,depth)

M=matlab2c(M);

xsize = norm(M(:,1));
ysize = norm(M(:,2));
zsize = norm(M(:,3));
x_r = M(1,1)/xsize;
x_a = M(2,1)/xsize;
x_s = M(3,1)/xsize;
y_r = M(1,2)/ysize;
y_a = M(2,2)/ysize;
y_s = M(3,2)/ysize;
z_r = M(1,3)/zsize;
z_a = M(2,3)/zsize;
z_s = M(3,3)/zsize;
ci = (width)/2 ; cj = (height)/2 ; ck = (depth)/2 ;
c_r = M(1,4) + (M(1,1)*ci + M(1,2)*cj + M(1,3)*ck) ;
c_a = M(2,4) + (M(2,1)*ci + M(2,2)*cj + M(2,3)*ck);
c_s = M(3,4) + (M(3,1)*ci + M(3,2)*cj + M(3,3)*ck);

return
