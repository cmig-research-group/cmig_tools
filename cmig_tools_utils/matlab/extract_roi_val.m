function [roi_mean, varargout] = extract_roi_val(vol, roinames, thresh, fname_parc)

	% Input Variables:
	% vol: 3D or 4D matrix of voxelwise data from which you wish to extract ROI values
	% 
	% roinames: cell array of strings containing the names of the ROIs you wish to extract
	%           - can be laterilised, eg. 'R_CST' for right cortico-spinal tract, and the function will output 
	% 		 	the mean voxel value for the specified hemisphere
	% 			- can be  bilateral, eg. 'CST', and the function will output the mean voxel value across both 			% 			hemispheres
	%      
	% thresh: threshold for binarizing the ROI masks (eg 0.8 for 80% probability)
	% 
	% fname_parc: path to .mat files of parcellations in ABCD3 space (e.g., aseg.mat, aparcaseg.mat, fiber.mat, 
	% 			  thalamus_abcd3_cor10.mat, pauli_subcortnuclei_abcd3_cor10.mat)
	%
	% Output variables
	% roi_mean: mean value of the voxels within the ROI
	% 			if vol is 4D, roi_mean will be a matrix of size length(roinames) x size(vol,4)
	% 			if vol is 3D, roi_mean will be a vector of length(roinames)
	%
	% Optional output variables 
	% vox_inc: indices of the voxels included for each ROI, a cell array of size length(roinames) 

	%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
	
	%$ load parcelltion 
	parc = load(fname_parc);
	if contains(fname_parc, 'thalamus_abcd3_cor10')
		parc.volmat_prob = parc.vol_labels_reg; 
		parc.roinames = parc.thalamus_lut.LabelName;
	end
	if contains(fname_parc, 'pauli_subcortnuclei_abcd3_cor10')
		parc.volmat_prob = parc.vol_labels_reg; 
		parc.roinames = parc.pauli_lut.Abbreviation;
	end
	% check dimnesions of the input volume and parc  (should we downsample parc or upsample input volume?)
	% stats output from FEMA are 100x100x130 whereas the parc is 200x200x260)
	dim_vol = size(vol);
	dim_parc = size(parc.volmat_prob);
	if ~isequal(dim_vol(1:3), dim_parc(1:3))
		% upsample the volume to match the parc dimensions
		if all(dim_vol(1:3)*2 == dim_parc(1:3))
			vol= upsample_volume(vol);
		else 
			error('Volume and parcellation dimensions do not match')
		end
	end

	% find the indices of the ROIs in the parc structure
	roi_mean = NaN(length(roinames), 1);
	vox_inc = cell(length(roinames), 1);
	for ri=1:length(roinames)
		% find all rois that match the input roiname 
		idx_roi = find(contains(parc.roinames, roinames{ri})); 
		% get the voxels inside rois
		ivec_roi = []; 
		for ii = 1:length(idx_roi)
			ivec_roi_tmp = find(parc.volmat_prob(:,:,:,idx_roi(ii)) > thresh);
			ivec_roi = [ivec_roi; ivec_roi_tmp];
		end 
		vox_inc{ri} = ivec_roi;
		if length(dim_vol)>3
			for di=1:dim_vol(4)
				vol_tmp = vol(:,:,:,di);
				vol_roi= vol_tmp(ivec_roi);
				roi_mean(ri,di) = mean(vol_roi);
			end
		else 
			vol_tmp = vol;
			vol_roi= vol_tmp(ivec_roi);
			roi_mean(ri) = mean(vol_roi);
		end 
	end 
	varargout{1} = ivec_roi;
end 
