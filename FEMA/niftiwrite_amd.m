function niftiwrite_amd(vol,fname_nii,M)

if strcmp(fname_nii(max(1,(length(fname_nii)-2)):end),'.gz')
  compressedflag = true;
  fname_nii = fname_nii(1:end-3);
end

T = M'; T(4,:) =  M*[1 1 1 1]'; Transform = struct('Dimensionality',3,'T',T); 
info = struct('Transform',Transform,'Datatype',class(vol),'ImageSize',size(vol),'Description','','Version','NIfTI1','Qfactor',1,'PixelDimensions',sum(M(:,1:3).^2,1),'SpaceUnits','Millimeter',...
              'TimeUnits','Second','SliceCode','Unknown','AdditiveOffset',0,'MultiplicativeScaling',1,'TimeOffset',0,'FrequencyDimension',0,'PhaseDimension',0,'SpatialDimension',0,'DisplayIntensityRange',[0 0],'TransformName','Sform');
niftiwrite(vol,fname_nii,info,'Compressed',compressedflag);

