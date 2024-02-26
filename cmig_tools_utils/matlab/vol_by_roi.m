function [roi_sum, roi_vox] = vol_by_roi(vol, ...
                                            atlas)
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
    aseg_names = annot.aseg.roinames;
    incl_aseg_rois = {'Cerebellum-White-Matter', 'Cerebellum-Cortex', 'Thalamus-Proper', 'Pallidum'};
    incl_aseg_rois = aseg_names;
    roinames = aseg_names;
    
    if exist('annot.pauli')
        pauli_prob = zeros(annot.pauli.prob{1}, 'single');
        pauli_prob(annot.pauli.prob{2}) = annot.pauli.prob{3};
        pauli_prob = pauli_prob(1:2:end, 1:2:end, 1:2:end, :);
        pauli_names = annot.pauli.roinames;
        roinames = [roinames; pauli_names];
        incl_pauli_rois = {'Caudate', 'Putamen', 'Subthalamic Nucleus', 'Red Nucleus', 'Substantia Nigra pars compacta ', 'Substantia Nigra pars reticulata'};
        incl_pauli_rois = pauli_names;
    end
    
    if exist('annot.thalamus')
        thalamus_prob = zeros(annot.thalamus.prob{1}, 'single');
        thalamus_prob(annot.thalamus.prob{2}) = annot.thalamus.prob{3};
        thalamus_prob = thalamus_prob(1:2:end, 1:2:end, 1:2:end, :);
        thalamus_names = annot.thalamus.roinames;
        incl_thalamus_rois = {'Pulvinar', 'Anterior', 'Medio_Dorsal', 'Ventral_Latero_Dorsal','Central_Lateral-Lateral_Posterior-Medial_Pulvinar','Ventral_Anterior','Ventral_Latero_Ventral'};
        incl_thalamus_rois = thalamus_names;
    end

        if exist('annot.aparcaseg')
        aparcaseg_prob = zeros(annot.aparcaseg.prob{1}, 'single');
        aparcaseg_prob(annot.aparcaseg.prob{2}) = annot.aparcaseg.prob{3};
        aparcaseg_prob = aparcaseg_prob(1:2:end, 1:2:end, 1:2:end, :);
        aparcaseg_names = annot.aparcaseg.roinames;
    end

    % Create output table
    roi_sum = struct('sum', [], 'cats', [], 'region_size', []);
    roi_vox = struct('vox', [], 'cats', []);

    % NB: probability hard coded, consider letting user specify
    [roi_sum, roi_vox] = extract_roi_val(vol, aseg_names, 0.8, incl_aseg_rois, prob, roi_sum, roi_vox);
    if exist('pauli_names') [roi_sum, roi_vox] = extract_roi_val(vol, pauli_names, 0.8,  incl_pauli_rois, pauli_prob, roi_sum, roi_vox); end
    if exist('thalamus_names') [roi_sum, roi_vox] = extract_roi_val(vol, thalamus_names, 0.8,  incl_thalamus_rois, thalamus_prob, roi_sum, roi_vox); end

end