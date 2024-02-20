function [roi_sum, roi_vox] = extract_roi_val(vox_mat, roi_names, thresh, incl_rois, prob, roi_sum, roi_vox)

    for i = 1:length(incl_rois)
        incl_roi = incl_rois{i};
        indices = find(cellfun(@(x) contains(x, incl_roi), roi_names));
        roi_ind = false(size(prob(:,:,:,1)));
        for ind = 1:length(indices)
            % Extract roi indicies
            roi_ind = roi_ind | prob(:, :, :, indices(ind))>thresh;
        end
        roi_vox.cats = [roi_vox.cats; repmat({incl_roi}, sum(roi_ind(:)), 1)];
        roi_sum.cats = [roi_sum.cats, {incl_roi}];
        % Extract effects from ROI
        tmp = vox_mat .* roi_ind;
        roi_sum.sum = [roi_sum.sum, sum(tmp(:))];
        roi_sum.region_size = [roi_sum.region_size, sum(roi_ind(:))];
        roi_vox.vox = [roi_vox.vox; vox_mat(roi_ind)];
    end
end
