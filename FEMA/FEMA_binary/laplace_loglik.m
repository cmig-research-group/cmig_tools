function [ll, b] = laplace_loglik(tau2, b, y, eta0, Z, rid, maxit)

    gdiag = tau2(rid);
    Ginv_diag = 1 ./ gdiag;
    q = numel(gdiag);

    for it = 1:maxit
        eta = eta0 + Z * b;
        p = 1 ./ (1 + exp(-eta));
        W = max(p .* (1 - p), 1e-10);

        grad = Z' * (y - p) - Ginv_diag .* b;
        K = Z' * (Z .* W) + spdiags(Ginv_diag, 0, q, q);
        K = (K + K') / 2;

        step = K \ grad;
        b = b + step;

        if norm(step) < 1e-4 * (1 + norm(b))
            break;
        end
    end

    eta = eta0 + Z * b;
    p = 1 ./ (1 + exp(-eta));
    W = max(p .* (1 - p), 1e-10);

    K = Z' * (Z .* W) + spdiags(Ginv_diag, 0, q, q);
    K = (K + K') / 2;

    eta_log1p = max(eta, 0) + log1p(exp(-abs(eta)));
    cond_ll = sum(y .* eta - eta_log1p);
    quad_penalty = sum((b.^2) .* Ginv_diag);
    logdetG = sum(log(gdiag));
    logdetK = logdet_spd(K);

    ll = cond_ll - 0.5 * quad_penalty - 0.5 * logdetG - 0.5 * logdetK;
end