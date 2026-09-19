function [u_update, prob, deviance] = compute_BLUP(u, X, beta_hat, sig2mat, allWsTerms, ...
                                                RandomEffects, RandomVar, ymat)

    num_obs = size(X, 1);
    num_y   = size(ymat, 2);
    num_RE = sum(~strcmpi(RandomEffects, 'E'));
    
    u_res      = u - X * beta_hat;
    Wu = zeros(num_obs, num_y, class(u));
    for yy = 1:num_y
        Wu(:,yy) = allWsTerms{yy} * u_res(:,yy);
    end

    u_update_r = zeros(num_obs, num_y, class(u));
    for ri = 1:num_RE
        fieldName = sprintf('V_%s', RandomEffects{ri});
        Z = RandomVar.(fieldName);
        Randombeta = (Z' * Wu) .* sig2mat(ri,:);
        u_update_r = u_update_r + Z * Randombeta;
    end

    u_update  = X * beta_hat + u_update_r;
    prob = 1 ./ (1 + exp(-u_update));
    prob = max(min(prob, 1 - 1e-6), 1e-6);
    deviance = -2 * sum(ymat .* log(prob) + (1 - ymat) .* log(1 - prob), 1);

end