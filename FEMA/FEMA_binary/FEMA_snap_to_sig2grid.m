function [binvec, nvec_bins] = FEMA_snap_to_sig2grid(sig2mat, nbins, FamilyStruct, W_1)

nsig2bins = FamilyStruct.nsig2bins;
num_y     = size(sig2mat, 2);
num_RFX   = size(sig2mat, 1);

if nbins == 0
    % No binning: each outcome has its own precision matrix.
    binvec    = 1:num_y;
    nvec_bins = ones(num_y, 1);
else
    % The first num_RFX-1 coordinates are non-residual variance proportions.
    % The last coordinate summarizes the inverse IRLS working weights.
    sig2vals = [double(sig2mat(1:end-1, :)); 1 - 4 ./ mean(double(W_1), 1)];
    sig2vals = min(1, max(0, sig2vals));

    % Convert each [0,1] coordinate to a grid index in 1:nbins.
    sig2gridi = min(nbins, floor(nbins * sig2vals) + 1).';
    binvec    = sub2ind_amd(nbins * ones(1, num_RFX), sig2gridi).';

    % Only occupied bins are stored; the vector still has nsig2bins rows.
    nvec_bins = sparse(binvec(:), 1, 1, nsig2bins, 1);
end

end
