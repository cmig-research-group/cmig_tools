function [tau2_lap, intercept_lap] = laplace_correct(y, X, beta_hat, eta_fixed, RandomEffects, RandomVar, tau2_init, IC_fam, IC_subj)


    % user parameters
    mode_max_iter   = 10;
    max_iter        = 10;
    opt_tol         = 1e-3;
    scale_range     = 5;
    max_laplace_obs = 20000;

    N = size(y,1);
    num_y = size(y,2);

    opt = optimset('Display', 'off', 'TolX', opt_tol,'MaxFunEvals', max_iter);
    
    isE = strcmpi(RandomEffects, 'E');
    if any(isE) && size(tau2_init,1) == numel(RandomEffects)
        RandomEffects = RandomEffects(~isE);
        tau2_init = tau2_init(~isE,:);
    end

    tau2_init = max(tau2_init, 1e-10);

    R = numel(RandomEffects);

    % Large-N - mini-batch
    keep_obs = [];
    if N > max_laplace_obs
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

        y         = y(keep_obs,:);
        X         = X(keep_obs,:);
        eta_fixed = eta_fixed(keep_obs,:);
    end

    mean_y      = mean(y, 1);

    % search interval
    lo_all     = log(max(tau2_init / scale_range, 1e-8));
    hi_all     = log(max(tau2_init * scale_range, 1e-6));
    logtau_all = log(tau2_init);
    
    % Random-effect loading matrices and their dmperm decomposition are shared
    field_names = strcat('V_', RandomEffects);
    Zparts  = cell(1, R);
    ridpart = cell(1, R);

    for r = 1:R
        Zr = RandomVar.(field_names{r});
        if isempty(keep_obs)
            Zr = sparse(double(Zr));
        else
            Zr = sparse(double(Zr(keep_obs, :)));
        end
        Zr = Zr(:, any(Zr, 1));
        Zparts{r} = Zr;
        ridpart{r} = r * ones(size(Zr, 2), 1);
    end

    Z = [Zparts{:}];
    rid = vertcat(ridpart{:});  % rid(k) tells which tau2 controls b_k

    Z_pattern = spones(Z);
    K_pattern = spones(Z_pattern' * Z_pattern) | speye(size(Z,2));
    [solve_info.perm, ~, solve_info.borders, ~] = dmperm(K_pattern);
    solve_info.component_size = diff(solve_info.borders);
    solve_info.all_scalar = all(solve_info.component_size == 1);
    solve_info.one_general = length(solve_info.component_size) == 1;
    solve_info.Zt = Z';
    solve_info.Z2t = (Z .* Z)';
    solve_info.scalar_idx = [];
    solve_info.block_sizes = [];
    solve_info.block_indices = {};
    if ~solve_info.all_scalar && ~solve_info.one_general
        scalar_components = find(solve_info.component_size == 1);
        solve_info.scalar_idx = solve_info.perm( ...
            solve_info.borders(scalar_components));
        solve_info.block_sizes = unique( ...
            solve_info.component_size(solve_info.component_size > 1))';
        solve_info.block_indices = cell(length(solve_info.block_sizes), 1);
        for ss = 1:length(solve_info.block_sizes)
            block_size = solve_info.block_sizes(ss);
            components = find(solve_info.component_size == block_size);
            indices = zeros(block_size, length(components));
            for gg = 1:length(components)
                cc = components(gg);
                indices(:,gg) = solve_info.perm( ...
                    solve_info.borders(cc):solve_info.borders(cc+1)-1);
            end
            solve_info.block_indices{ss} = indices;
        end
    end

    q             = size(Z, 2);
    tau2_lap      = zeros(R, num_y);
    intercept_lap = zeros(1, num_y);

    for yy = 1:num_y
        y_current         = y(:,yy);
        eta_fixed_current = eta_fixed(:,yy);
        p_current         = 1 ./ (1 + exp(-eta_fixed_current));
        mean_p            = mean(p_current);
        logit_gap = abs(log(mean_y(yy)/(1-mean_y(yy))) - ...
                        log(mean_p/(1-mean_p)));
        imbalance = min(mean_y(yy), 1-mean_y(yy));
        do_joint  = (imbalance < 0.25) || (logit_gap > 0.1);

        delta0 = 0;
        lo      = lo_all(:,yy);
        hi      = hi_all(:,yy);
        logtau  = logtau_all(:,yy);
        b_cache = zeros(q, 1); % warm start across objective evaluations

        if do_joint
    
            % joint intercept + sigma correction
            lo_delta = -abs(beta_hat(1,yy));
            hi_delta =  abs(beta_hat(1,yy));
    
            for pass = 1:2
        
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
    
            % sigma-only correction
            if R == 1
        
                obj = @(a) obj_sparse_global(a, delta0);
                logtau = fminbnd(obj, lo, hi, opt);
        
            else
        
                for r = 1:R
                    obj_r = @(ar) obj_sparse_global(set_one(logtau, r, ar), delta0);
                    logtau(r) = fminbnd(obj_r, lo(r), hi(r), opt);
                end
            end
        end

        tau2_lap(:,yy)    = exp(logtau);
        intercept_lap(yy) = beta_hat(1,yy) + delta0;

    end

    function f = obj_sparse_global(logtau2,delta0_current)
        % logtau2 can be scalar or vector.
        tau2 = exp(logtau2(:));
        eta_current = eta_fixed_current + delta0_current * X(:,1);
        [ll, b_cache] = laplace_loglik( ...
            tau2, b_cache, y_current, eta_current, Z, rid, mode_max_iter, solve_info);
        f = -ll;
    end


end

function x = set_one(x, r, val)
    x(r) = val;
end



