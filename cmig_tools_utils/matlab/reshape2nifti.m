function vol_out = reshaoe2nifti(v)

dims = size(v); 
f = factor(dims(1));
if length(f)==3
  dims_tmp = f;
elseif length(f)>3
  dims_tmp = [prod(f(1:(length(f)-2))) f((length(f)-1):end)];
else
  dims_tmp = [1 1 1];
  dims_tmp((end-length(f)+1):end) = f;
%  keyboard % Not sure if this results in valid NIFTI file -- may need to have all non-singleton diemnsions
end

vol_out  = reshape(v,[dims_tmp dims(2)]);

