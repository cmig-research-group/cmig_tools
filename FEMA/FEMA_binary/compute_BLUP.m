function [u_update, p, deviance] = compute_BLUP(u, X, beta_hat, sig2mat, allWsTerms, ...
                                                RandomEffects, RandomVar, ymat)
    
    u_res = u - X * beta_hat;
    nobs = size(X,1);
    Randombeta                              = cell(1, size(ymat,2));
    
    for coli = 1:size(ymat,2)

        Wu                 = allWsTerms{coli} * u_res(:,coli);
        Randombeta{coli}   = zeros(nobs, length(RandomEffects)-1);

        for r = 1:(length(RandomEffects)-1) % excluding the error part
    
            RandomVar_tmp         = RandomVar.(sprintf('V_%s', RandomEffects{r}));
            Randombeta_tmp        = sig2mat(r,coli) * RandomVar_tmp' * Wu; 
            Randombeta{coli}(:,r) = RandomVar_tmp * Randombeta_tmp; % random coefficients for all observations

        end

    end

    u_update_r = cell2mat(cellfun(@(x) sum(x,2), Randombeta, 'UniformOutput', false));

    u_update  = X * beta_hat + u_update_r;

    p = 1 ./ (1 + exp(-u_update));

    % p = 0.95 * p + (1 - 0.95) * mean(p);

    p = max(min(p, 1 - 1e-6), 1e-6);

    deviance = -2 * sum(ymat .* log(p) + (1 - ymat) .* log(1 - p), 1);

end