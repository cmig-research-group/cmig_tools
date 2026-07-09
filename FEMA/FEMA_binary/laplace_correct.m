function [tau2_lap, intercept_lap] = laplace_correct(y, X, beta_hat, RandomEffects, RandomVar, tau2_init, IC_fam, IC_subj)


    % user parameters
    mode_max_iter = 10;
    max_iter = 10;
    mode_tol      = 1e-4;
    opt_tol       = 1e-3;
    scale_range   = 5;
    max_laplace_obs = 20000;

    N = size(y,1);

    opt = optimset('Display', 'off', 'TolX', opt_tol,'MaxFunEvals', max_iter);
    
    isE = strcmpi(RandomEffects, 'E');
    if any(isE) && numel(tau2_init) == numel(RandomEffects)
        RandomEffects = RandomEffects(~isE);
        tau2_init = tau2_init(~isE);
    end

    tau2_init = max(tau2_init, 1e-10);

    R = numel(RandomEffects);

    % Large-N - mini-batch
    if size(y,1) > max_laplace_obs
        if any(strcmpi(RandomEffects, 'F'))
            sample_block = IC_fam;
        elseif any(strcmpi(RandomEffects, 'S'))
            sample_block = IC_subj;
        end
    
        J_block = max(sample_block);
        block_size = accumarray(sample_block, 1, [J_block, 1], @sum, 0);
    
        rng_state = rng;
        rng(1, 'twister');
        block_order = randperm(J_block);
        rng(rng_state);
    
        cum_size = cumsum(block_size(block_order));
        n_block_keep = find(cum_size >= max_laplace_obs, 1, 'first');
    
        chosen_block = block_order(1:n_block_keep);
        keep_obs = ismember(sample_block, chosen_block);
    else
        keep_obs = true(N, 1);
    end

    y = y(keep_obs);
    X = X(keep_obs, :);
    
    IC_fam  = IC_fam(keep_obs);
    IC_subj = IC_subj(keep_obs);

    eta_fixed = X * beta_hat;
    delta0 = 0;
    p = 1 ./ (1+exp(-eta_fixed));
    logit_gap = abs(log(mean(y)/(1-mean(y))) - log(mean(p)/(1-mean(p))));
    imbalance = min(mean(y), 1-mean(y));
    do_joint = (imbalance < 0.25) || (logit_gap > 0.1);

    % search interval
    lo = log(max(tau2_init / scale_range, 1e-8));
    hi = log(max(tau2_init * scale_range, 1e-6));
    logtau = log(tau2_init);
    

    Zparts  = cell(1, R);
    ridpart = cell(1, R);

    for r = 1:R
        field_name = sprintf('V_%s', RandomEffects{r});
        Zr = sparse(double(RandomVar.(field_name)(keep_obs, :)));
        used_cols = find(any(Zr, 1));
        Zr = Zr(:, used_cols);
        Zparts{r} = Zr;
        ridpart{r} = r * ones(size(Zr, 2), 1);
    end

    Z = [Zparts{:}];
    rid = vertcat(ridpart{:});  % rid(k) tells which tau2 controls b_k
    q = size(Z, 2);
    b_cache = zeros(q, 1);      % warm start across objective evaluations


    if do_joint

        coord_passes = 2;
        % joint intercept + sigma correction
        lo_delta = -abs(beta_hat(1));
        hi_delta =  abs(beta_hat(1));

        for pass = 1:coord_passes
    
            par_old = [delta0; logtau(:)];
    
            obj_delta = @(d) obj_sparse_global(logtau, d);
            delta0 = fminbnd(obj_delta, lo_delta, hi_delta, opt);
    
            for r = 1:R
                obj_r = @(ar) obj_sparse_global(set_one(logtau, r, ar), delta0);
                logtau(r) = fminbnd(obj_r, lo(r), hi(r), opt);
            end
    
            par_new = [delta0; logtau(:)];
    
            if norm(par_new - par_old) < opt_tol * (1 + norm(par_old))
                break;
            end
        end
    
    else
    
        coord_passes = 1;

        % sigma-only correction
        if R == 1
    
            obj = @(a) obj_sparse_global(a, delta0);
            logtau = fminbnd(obj, lo, hi, opt);
    
        else
    
            for pass = 1:coord_passes
    
                logtau_old = logtau;
    
                for r = 1:R
                    obj_r = @(ar) obj_sparse_global(set_one(logtau, r, ar), delta0);
                    logtau(r) = fminbnd(obj_r, lo(r), hi(r), opt);
                end
    
                if norm(logtau - logtau_old) < opt_tol * (1 + norm(logtau_old))
                    break;
                end
            end
        end
    end

    tau2_lap = exp(logtau);
    delta0_lap = delta0;
    intercept_lap = beta_hat(1) + delta0_lap;


    function f = obj_sparse_global(logtau2,delta0_current)
        % logtau2 can be scalar or vector.
        tau2 = exp(logtau2(:));
        eta_current = eta_fixed + delta0_current * X(:,1);
        [ll, b_cache] = laplace_loglik( ...
            tau2, b_cache, y, eta_current, Z, rid, mode_max_iter);
        f = -ll;
    end


end

function x = set_one(x, r, val)
    x(r) = val;
end



