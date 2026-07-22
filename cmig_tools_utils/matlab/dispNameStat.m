function dispName = dispNameStat(statName)
% Return display name for FEMA statistic fields.
switch statName
    case 'beta_hat'
        dispName = 'Beta';
    case 'beta_se'
        dispName = 'Standard Error';
    case 'zmat'
        dispName = 'T statistics';
    case 'logpmat'
        dispName = 'Signed -log10(p)';
    case 'Wald'
        dispName = 'Wald statistics';
    case 'p_Wald'
        dispName = 'Wald p value';
    case 'logp_Wald'
        dispName = 'Wald -log10(p)';
    case 'sig2tvec'
        dispName = 'Total variance';
    case 'sig2mat'
        dispName = 'Variance';
    case 'sig2mat_normalized'
        dispName = 'Normalized variance';
    case 'sig2tvec_perm'
        dispName = 'Total variance (permuted)';
    case 'sig2mat_perm'
        dispName = 'Variance (permuted)';
    case 'sig2mat_normalized_perm'
        dispName = 'Normalized variance (permuted)';
    otherwise
        error('Unknown statName: %s', statName);
end
end
