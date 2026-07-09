function ld = logdet_spd(A)
% Numerically stable log determinant for a small SPD matrix.
    A = (A + A') / 2;
    jitter = 0;

    for k = 1:5
        [R, p] = chol(A + jitter * eye(size(A)));
        if p == 0
            ld = 2 * sum(log(diag(R)));
            return;
        end
        jitter = max(1e-10, 10 * max(jitter, eps));
    end

    ev = eig(A);
    ev = max(real(ev), 1e-12);
    ld = sum(log(ev));
end