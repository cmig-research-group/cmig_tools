function writeNifti(results, dirname_out, fstem_imaging, ivnames, colnames_model)
% writeNifti write voxelwise data to nifti, called from SSE_wrapper
%
%   writeNifti(results, dirname_out, fstem_imaging, ivnames, colnames_model)
%
% writes a volume for each variable type in results and each IV in ivnames
   
%new behavior 6/2021: write out only beta_hat and a single file per column, but only for IVs
% FIXME: until we get a list of actual IVs, write out all columns, with exception of "mri_info_*" for scanner info
% we could add to this list, but other covariates are potentially things people would want to model
% 7/2021 added 'ivnames' input--comma separated list of IVs. Also will write an individual .mat

randomFields = {'sig2tvec', 'sig2mat'};

% parse IVs
if isempty(ivnames)
  excludeCol = strmatch('mri_info_',colnames_model);
  nCol = length(colnames_model);
  ivCol = setdiff(1:nCol, excludeCol);
else
  [~,ivCol,~] = intersect(colnames_model,ivnames);
end
if length(ivCol) < 1, error('No IVs found! Not writing nifti.'), end

% =========================================================================
% write out main effects (for IVs of interest)
fieldnamelist = setdiff(fieldnames(results),randomFields);

for fi = 1:length(fieldnamelist)
  fieldname = fieldnamelist{fi};
  vol_nifti = results.(fieldname);
  
  %write one file per IV
  for iv = ivCol(:)'
    colname = colnames_model{iv};
    vol_nifti_col = vol_nifti(:,:,:,iv);
    %vol_nifti_col = upsample_volume(vol_nifti_col); % make this optional?
    
    fname_nii = sprintf('%s/SSE_results_voxelwise_%s_%s__%s.nii',dirname_out,fstem_imaging,fieldname,colname);
    niftiwrite(vol_nifti_col,fname_nii);
    fprintf(1,'File %s written (dims = [%s])\n',fname_nii,num2str(size(vol_nifti_col),'%d '));
  end
end

% =========================================================================
% write out the random effects
fieldnamelist = randomFields;
for fi = 1:length(fieldnamelist)
  fieldname = fieldnamelist{fi};
  vol_nifti = results.(fieldname);
  
  %write one file per IV
  for iv = 1:size(vol_nifti,4)
    colname = sprintf('%02d',iv);
    vol_nifti_col = vol_nifti(:,:,:,iv);
    %vol_nifti_col = upsample_volume(vol_nifti_col); % make this optional?
    
    fname_nii = sprintf('%s/SSE_results_voxelwise_%s_%s__%s.nii',dirname_out,fstem_imaging,fieldname,colname);
    niftiwrite(vol_nifti_col,fname_nii);
    fprintf(1,'File %s written (dims = [%s])\n',fname_nii,num2str(size(vol_nifti_col),'%d '));
  end
end



% =========================================================================
%write column names to json
colnames_model = colnames_model(ivCol); %subset to the IVs actually written
fname_col = sprintf('%s/SSE_results_colnames.json',dirname_out);
out = struct('colnames_model',{colnames_model});
jsonStr = jsonencode(out);
fid = fopen(fname_col,'w');
fprintf(fid,'%s\n',jsonStr);
fclose(fid);