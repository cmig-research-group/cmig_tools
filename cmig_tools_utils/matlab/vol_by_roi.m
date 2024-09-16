function [roi_sum, roi_vox] = vol_by_roi(vol, varargin)
	% Function to group voxelwise statistics into their ROIs for a given ABCD
	% atlas.
	%% Inputs:
	% vol                   voxelwise statmap (e.g., from FEMA output)
	%% Optional inputs:
	% atlas                 string indicating which ABCD atlas to use (default:
	%                       'ABCD3')
	% incl_aseg_rois        list of aseg ROIs to include; default is to
	%                       include all ROIs in 'annot.aseg.roinames'
	% incl_pauli_rois       list of pauli ROIs to include; default is to
	%                       include all ROIs in 'annot.pauli.roinames'
	% incl_thalamus_rois    list of thalamus ROIs to include; default is to
	%                       include all ROIs in 'annot.thalamus.roinames'

	p = inputParser;
	addParamValue(p,'atlas', '');
	addParamValue(p,'incl_aseg_rois', []);
	addParamValue(p,'incl_aparc_rois', []);
	addParamValue(p,'incl_pauli_rois', []);
	addParamValue(p,'incl_thalamus_rois', []);
	addParamValue(p,'incl_fiber_rois', []);
	addParamValue(p,'thresh', 0.8);
	addParamValue(p,'dirname_abcdsync', []);
	addParamValue(p,'release', '5.0');

	parse(p,varargin{:})
	atlas = p.Results.atlas;
	incl_aseg_rois = p.Results.incl_aseg_rois;
	incl_aparc_rois = p.Results.incl_aparc_rois;
	incl_pauli_rois = p.Results.incl_pauli_rois;
	incl_thalamus_rois = p.Results.incl_thalamus_rois;
	incl_fiber_rois = p.Results.incl_fiber_rois;
	thresh = p.Results.thresh;
	dirname_abcdsync = p.Results.dirname_abcdsync;
	release = p.Results.release;

	% Load atlases
    if ~isempty(dirname_abcdsync)
		fprintf('Loading atlases from %s/%s\n', dirname_abcdsync, release);
		dirname_roi = sprintf('%s/%s/imaging_concat/voxelwise/roi', dirname_abcdsync, release);
		if ~isempty('incl_aseg_rois') 
			fname_roi = fullfile(dirname_roi, 'aseg.mat');
			aseg = load(fname_roi, 'volmat_prob', 'roinames', 'roicodes');
			aseg.volmat_prob = subsample_volume(aseg.volmat_prob);
		end
		if ~isempty(incl_aparc_rois) 
			fname_roi = fullfile(dirname_roi, 'aparcaseg.mat');
			aparcaseg = load(fname_roi, 'volmat_prob', 'roinames', 'roicodes');
			aparcaseg.volmat_prob = subsample_volume(aparcaseg.volmat_prob);
		end
		if ~isempty(incl_pauli_rois) 
			fname_roi = fullfile(dirname_roi, 'pauli_subcortnuclei_abcd3_cor10.mat');
			pauli = load(fname_roi, 'vol_labels_reg', 'T');
			pauli.roinames = pauli.T.Abbreviation;
			pauli.volmat_prob = subsample_volume(pauli.vol_labels_reg);
		end
		if ~isempty(incl_thalamus_rois) 
			fname_roi = fullfile(dirname_roi, 'thalamus_abcd3_cor10.mat');
			thalamus = load(fname_roi, 'vol_labels_reg', 'T');
			thalamus.roinames = thalamus.T.LabelName;
			thalamus.volmat_prob = subsample_volume(thalamus.vol_labels_reg);
		end
		if ~isempty(incl_fiber_rois) 
			fname_roi = fullfile(dirname_roi, 'fiber.mat');
			fiber = load(fname_roi, 'volmat_prob', 'roinames', 'roicodes');
			fiber.volmat_prob = subsample_volume(fiber.volmat_prob);
		end
	else 
		if isempty(atlas)
	        disp("No atlas specified; defaulting to ABCD3");
	        atlas = 'ABCD3';
	    end
		try 
	    	cfg = abcdConfig;                
	    	annot = load(fullfile(cfg.data.showVolData,'Atlas',sprintf('showVolAtlases_%s_cor10.mat',atlas)));
		catch
	    	error('Could not load atlas');
		end
		if ~isempty(incl_aseg_rois)
	    	c = annot.aseg.prob;
	    	prob = zeros(c{1}, 'single');
	    	prob(c{2}) = c{3};
	    	annot.aseg.prob = prob;
	    	prob = subsample_volume(prob);
			aseg.volmat_prob = prob;
			aseg.roinames = annot.aseg.roinames;
	    end 
		if ~isempty(incl_aparc_rois)
	    	c = annot.aparcaseg.prob;
	    	prob = zeros(c{1}, 'single');
	    	prob(c{2}) = c{3};
	    	annot.aparcaseg.prob = prob;
			prob = subsample_volume(prob);
			aparcaseg.volmat_prob = prob;
			aparcaseg.roinames = annot.aparcaseg.roinames;
	    end
		if ~isempty(incl_pauli_rois)
	    	c = annot.pauli.prob;
	    	prob = zeros(c{1}, 'single');
	    	prob(c{2}) = c{3};
	    	annot.pauli.prob = prob;
	    	prob = subsample_volume(prob);
			pauli.volmat_prob = prob;
			pauli.roinames = annot.pauli.roinames;
	    end
		if ~isempty(incl_thalamus_rois)
	    	c = annot.thalamus.prob;
	    	prob = zeros(c{1}, 'single');
	    	prob(c{2}) = c{3};
	    	annot.thalamus.prob = prob;
	    	prob = subsample_volume(prob);
			thalamus.volmat_prob = prob;
			thalamus.roinames = annot.thalamus.roinames;
	    end
		if ~isempty(incl_fiber_rois)
	    	c = annot.fiber.prob;
	    	prob = zeros(c{1}, 'single');
	    	prob(c{2}) = c{3};
	    	annot.fiber.prob = prob;
	    	prob = subsample_volume(prob);
			fiber.volmat_prob = prob;
			fiber.roinames = annot.fiber.roinames;
	    end
	end 
	% Create output table
    roi_sum = struct('sum', [], 'cats', [], 'region_size', []);
    roi_vox = struct('vox', [], 'cats', []);

	if ~isempty(incl_aseg_rois)
		[roi_sum, roi_vox] = extract_roi_val(vol, aseg.roinames, thresh, incl_aseg_rois, aseg.volmat_prob, roi_sum, roi_vox);
	end    
	if ~isempty(incl_aparc_rois) 
		[roi_sum, roi_vox] = extract_roi_val(vol, aparcaseg.roinames, thresh,  incl_aparc_rois, aparcaseg.volmat_prob, roi_sum, roi_vox); 
	end
	if ~isempty(incl_pauli_rois) 
		[roi_sum, roi_vox] = extract_roi_val(vol, pauli.roinames, thresh,  incl_pauli_rois, pauli.volmat_prob, roi_sum, roi_vox); 
	end
	% N.B. threhsold for pauli and thalamus is hard coded at 0.5 because it's not a probability 
    if ~isempty(incl_thalamus_rois) 
		[roi_sum, roi_vox] = extract_roi_val(vol, thalamus.roinames, 0.5,  incl_thalamus_rois, thalamus.volmat_prob, roi_sum, roi_vox); 
	end 
    if ~isempty(incl_fiber_rois) 
		[roi_sum, roi_vox] = extract_roi_val(vol, fiber.roinames, 0.5,  incl_fiber_rois, fiber.volmat_prob, roi_sum, roi_vox); 
	end
end