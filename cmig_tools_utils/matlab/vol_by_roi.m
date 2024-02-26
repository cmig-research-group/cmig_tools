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

    if ~exist('atlas')
        disp("No atlas specified; defaulting to ABCD3");
        atlas = 'ABCD3';
    end
    cfg = abcdConfig;                
    annot = load(fullfile(cfg.data.showVolData,'Atlas',sprintf('showVolAtlases_%s_cor10.mat',atlas)));
    c = annot.aseg.prob;
    prob = zeros(c{1}, 'single');
    prob(c{2}) = c{3};
    annot.aseg.prob = prob;
    prob = prob(1:2:end, 1:2:end, 1:2:end, :);
    if ~exist('incl_aseg_rois') incl_aseg_rois = annot.aseg.roinames; end
    % roinames = incl_aseg_rois;
    
    if exist('annot.pauli')
        pauli_prob = zeros(annot.pauli.prob{1}, 'single');
        pauli_prob(annot.pauli.prob{2}) = annot.pauli.prob{3};
        pauli_prob = pauli_prob(1:2:end, 1:2:end, 1:2:end, :);
        if ~exist('incl_pauli_rois') incl_pauli_rois = annot.pauli.roinames; end
        % roinames = [roinames; incl_pauli_rois];
    end
    
    if exist('annot.thalamus')
        thalamus_prob = zeros(annot.thalamus.prob{1}, 'single');
        thalamus_prob(annot.thalamus.prob{2}) = annot.thalamus.prob{3};
        thalamus_prob = thalamus_prob(1:2:end, 1:2:end, 1:2:end, :);
        if ~exist('incl_thalamus_rois') incl_thalamus_rois = annot.thalamus.roinames; end
    end

    if exist('annot.fiber')
        fiber_prob = zeros(annot.fiber.prob{1}, 'single');
        fiber_prob(annot.thalamus.prob{2}) = annot.fiber.prob{3};
        fiber_prob = fiber_prob(1:2:end, 1:2:end, 1:2:end, :);
        if ~exist('incl_fiber_rois') incl_fiber_rois = annot.fiber.roinames; end
    end

    % Create output table
    roi_sum = struct('sum', [], 'cats', [], 'region_size', []);
    roi_vox = struct('vox', [], 'cats', []);

    % NB: probability hard coded, consider letting user specify
    [roi_sum, roi_vox] = extract_roi_val(vol, annot.aseg.roinames, 0.8, incl_aseg_rois, prob, roi_sum, roi_vox);
    if exist('pauli_names') [roi_sum, roi_vox] = extract_roi_val(vol, annot.pauli.roinames, 0.8,  incl_pauli_rois, pauli_prob, roi_sum, roi_vox); end
    if exist('thalamus_names') [roi_sum, roi_vox] = extract_roi_val(vol, annot.thalamus.roinames, 0.8,  incl_thalamus_rois, thalamus_prob, roi_sum, roi_vox); end
    if exist('fiber_names') [roi_sum, roi_vox] = extract_roi_val(vol, annot.fiber.roinames, 0.8,  incl_fiber_rois, fiber_prob, roi_sum, roi_vox); end

end