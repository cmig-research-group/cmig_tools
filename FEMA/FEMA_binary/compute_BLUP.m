function [u_update, p, deviance] = compute_BLUP(u, X, beta_hat, sig2mat, allWsTerms, ...
                                                RandomEffects, RandomVar, ymat)
    
    u_res      = u - X * beta_hat;
    nobs       = size(X,1);
    Wu         = allWsTerms * u_res;

    % The third dimension indexes the columns of y
    Randombeta = zeros(nobs, length(RandomEffects)-1, size(ymat,2));

    for r = 1:(length(RandomEffects)-1) % excluding the error part
        RandomVar_tmp     = RandomVar.(sprintf('V_%s', RandomEffects{r}));
        Randombeta_tmp    = sig2mat(r) * RandomVar_tmp' * Wu; 
        Randombeta(:,r,:) = RandomVar_tmp * Randombeta_tmp; % random coefficients for all observations
    end

    % Summing across the random effects; squeeze drops the singleton random
    % effects dimension since that is what is being collapsed
    u_update_r = squeeze(sum(Randombeta, 2));

    % X*beta_hat will nobs * nyvar; to each, add u_update_r but cognizant
    % of the third dimeension; result should be nobs * nyvars
    u_update  = X * beta_hat + u_update_r;

    p = 1 ./ (1 + exp(-u_update));

    p = max(min(p, 1 - 1e-6), 1e-6);

    deviance = -2 * sum(ymat .* log(p) + (1 - ymat) .* log(1 - p), 1);

end
