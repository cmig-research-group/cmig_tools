function [ll, b] = laplace_loglik(tau2, b, y, eta0, Z, rid, maxit, solve_info)

    gdiag = tau2(rid);
    Ginv_diag = 1 ./ gdiag;
    q = numel(gdiag);

    if nargin < 8 || isempty(solve_info)
        Z_pattern = spones(Z);
        K_pattern = spones(Z_pattern' * Z_pattern) | speye(q);
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
            n_components = length(solve_info.component_size);
            component_indices = cell(n_components, 1);
            for cc = 1:n_components
                component_indices{cc} = solve_info.perm( ...
                    solve_info.borders(cc):solve_info.borders(cc+1)-1);
            end
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
                    indices(:,gg) = component_indices{components(gg)};
                end
                solve_info.block_indices{ss} = indices;
            end
        end
    end

    if ~solve_info.all_scalar
        Ginv = spdiags(Ginv_diag, 0, q, q);
    end

    for it = 1:maxit
        eta = eta0 + Z * b;
        p = 1 ./ (1 + exp(-eta));
        W = max(p .* (1 - p), 1e-10);

        grad = solve_info.Zt * (y - p) - Ginv_diag .* b;
        if solve_info.all_scalar
            Kdiag = full(solve_info.Z2t * W) + Ginv_diag;
            step = grad ./ Kdiag;
        else
            K = solve_info.Zt * (Z .* W) + Ginv;
            K = (K + K') / 2;
            if solve_info.one_general
                step = K \ grad;
            else
                step = zeros(q, 1);
                scalar_idx = solve_info.scalar_idx;
                step(scalar_idx) = grad(scalar_idx) ./ full(diag(K(scalar_idx,scalar_idx)));

                for ss = 1:length(solve_info.block_sizes)
                    block_size = solve_info.block_sizes(ss);
                    indices = solve_info.block_indices{ss};
                    n_group = size(indices, 2);
                    K_page = zeros(block_size, block_size, n_group);
                    g_page = reshape(grad(indices), block_size, 1, n_group);
                    for gg = 1:n_group
                        idx = indices(:,gg);
                        K_page(:,:,gg) = full(K(idx,idx));
                    end
                    step_page = pagemldivide(K_page, g_page);
                    step(indices) = reshape(step_page, block_size, n_group);
                end
            end
        end
        b = b + step;

        if norm(step) < 1e-4 * (1 + norm(b))
            break;
        end
    end

    eta = eta0 + Z * b;
    p = 1 ./ (1 + exp(-eta));
    W = max(p .* (1 - p), 1e-10);

    if solve_info.all_scalar
        Kdiag = full(solve_info.Z2t * W) + Ginv_diag;
        logdetK = sum(log(Kdiag));
    else
        K = solve_info.Zt * (Z .* W) + Ginv;
        K = (K + K') / 2;
        if solve_info.one_general
            logdetK = logdet_spd(K);
        else
            scalar_idx = solve_info.scalar_idx;
            logdetK = sum(log(full(diag(K(scalar_idx,scalar_idx)))));
            for ss = 1:length(solve_info.block_sizes)
                indices = solve_info.block_indices{ss};
                for gg = 1:size(indices,2)
                    idx = indices(:,gg);
                    logdetK = logdetK + logdet_spd(full(K(idx,idx)));
                end
            end
        end
    end

    eta_log1p = max(eta, 0) + log1p(exp(-abs(eta)));
    cond_ll = sum(y .* eta - eta_log1p);
    quad_penalty = sum((b.^2) .* Ginv_diag);
    logdetG = sum(log(gdiag));

    ll = cond_ll - 0.5 * quad_penalty - 0.5 * logdetG - 0.5 * logdetK;
end
