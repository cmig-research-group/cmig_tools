function writeNIFTI(results, dirname_out, fstem_imaging, ivnames, colnames_model)
% writeNifti write voxelwise data to nifti, called from SSE_wrapper
%
%   writeNifti(results, dirname_out, fstem_imaging, ivnames, colnames_model)
%
% writes a volume for each variable type in results and each IV in ivnames
   
%new behavior 6/2021: write out only beta_hat and a single file per column, but only for IVs
% FIXME: until we get a list of actual IVs, write out all columns, with exception of "mri_info_*" for scanner info
% we could add to this list, but other covariates are potentially things people would want to model
% 7/2021 added 'ivnames' input--comma separated list of IVs. Also will write an individual .mat

% These should be saved along with volinfo, and passed along to this function; hardcode for now
M_atl = [0    -1     0   101; 0     0     1  -131; -1     0     0   101; 0     0     0     1];
M_atl_sub = [0    -2     0   102; 0     0     2  -132; -2     0     0   102; 0     0     0     1];

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
%    colname = colnames_model{iv};
    colname = sprintf('%02d',iv);
    vol_nifti_col = vol_nifti(:,:,:,iv); M_tmp = M_atl_sub;
%    vol_nifti_col = upsample_volume(vol_nifti_col); M_tmp = M_atl; % Optionally, upsample volume to full resolution
% Should optionally upsample volume?
    if 0
      fname_nii = sprintf('%s/SSE_results_voxelwise_%s_%s_%s.nii',dirname_out,fstem_imaging,fieldname,colname);
      niftiwrite(vol_nifti_col,fname_nii,'Compressed',true); % Should find a way to save geometry info with niftiwrite -- use fs_save_mgh / mri_convert hack for now
    else
      fname_mgh = sprintf('%s/SSE_results_voxelwise_%s_%s_%s.mgh',dirname_out,fstem_imaging,fieldname,colname); fname_nii = strrep(fname_mgh,'.mgh','.nii.gz');
      fs_save_mgh(vol_nifti_col,fname_mgh,M_tmp);
      cmd = sprintf('mri_convert %s %s',fname_mgh,fname_nii);
      [s r e] = jsystem(cmd); fname_vol_z_nii = fname_nii;
      delete(fname_mgh);
    end
    fprintf(1,'File %s written (dims = [%s])\n',fname_nii,num2str(size(vol_nifti_col),'%d '));
  end
end

% =========================================================================
% write out the random effects
fieldnamelist = randomFields;
for fi = 1:length(fieldnamelist)
  fieldname = fieldnamelist{fi};
  vol_nifti = results.(fieldname);
  for iv = 1:size(vol_nifti,4)
    colname = sprintf('%02d',iv);
    vol_nifti_col = vol_nifti(:,:,:,iv); M_tmp = M_atl_sub;
%    vol_nifti_col = upsample_volume(vol_nifti_col); M_tmp = M_atl; % Optionally, upsample volume to full resolution
    if 0
      fname_nii = sprintf('%s/SSE_results_voxelwise_%s_%s_%s.nii',dirname_out,fstem_imaging,fieldname,colname);
      niftiwrite(vol_nifti_col,fname_nii);
    else
      fname_mgh = sprintf('%s/SSE_results_voxelwise_%s_%s_%s.mgh',dirname_out,fstem_imaging,fieldname,colname); fname_nii = strrep(fname_mgh,'.mgh','.nii.gz');
      fs_save_mgh(vol_nifti_col,fname_mgh,M_tmp);
      cmd = sprintf('mri_convert %s %s',fname_mgh,fname_nii);
      [s r e] = jsystem(cmd); fname_vol_z_nii = fname_nii;
      delete(fname_mgh);
    end
    fprintf(1,'File %s written (dims = [%s])\n',fname_nii,num2str(size(vol_nifti_col),'%d '));
  end
end

% ToDos
%   Write NIFTI files as .nii.gz
%   Write out z-scores and/or logp
%   Make sure coordinate infor is correct (in atlas space) -- does niftiwrite allow for specification of transformation matrix?
%   Write random effects with name of random eefect in filename 
