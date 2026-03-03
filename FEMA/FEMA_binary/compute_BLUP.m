function [u_update, p, deviance] = compute_BLUP(u, X, beta_hat, sig2mat, allWsTerms, ...
                                                RandomEffects, RandomVar, ymat)
    
    u_res = u - X * beta_hat;
    nobs = size(X,1);
    Wu                 = allWsTerms * u_res;
    Randombeta   = zeros(nobs, length(RandomEffects)-1);

    for r = 1:(length(RandomEffects)-1) % excluding the error part

        RandomVar_tmp         = RandomVar.(sprintf('V_%s', RandomEffects{r}));
        Randombeta_tmp        = sig2mat(r) * RandomVar_tmp' * Wu; 
        Randombeta(:,r) = RandomVar_tmp * Randombeta_tmp; % random coefficients for all observations

    end


    u_update_r = sum(Randombeta, 2);

    u_update  = X * beta_hat + u_update_r;

    p = 1 ./ (1 + exp(-u_update));

    p = max(min(p, 1 - 1e-6), 1e-6);

    deviance = -2 * sum(ymat .* log(p) + (1 - ymat) .* log(1 - p), 1);

end
