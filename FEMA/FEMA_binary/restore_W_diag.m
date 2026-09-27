function allV = restore_W_diag(Vis_base, Vs_base, W_1_base, W_1_target, page_limit, num_blocks)
% restore outcome-specific diagonal working variances from
% a shared bin-level covariance inverse. Starting from Vis_base = inv(Vs_base),
% it applies sequential Sherman-Morrison rank-one updates for the diagonal
% differences W_1_target - W_1_base. Each column of W_1_target represents
% one outcome. If an update is numerically unstable, the corresponding
% covariance is inverted directly. The output contains the vectorized
% precision matrix for each outcome, with size m^2-by-num_y.

if nargin < 5 || isempty(page_limit)
    page_limit = 100000;
end

m        = size(W_1_target, 1);
if nargin >= 6 && ~isempty(num_blocks)
    target_2d  = false;
    num_y      = numel(W_1_target) / (m * num_blocks);
elseif ismatrix(W_1_target) % no binning
    target_2d  = true;
    num_y      = size(W_1_target, 2);
    num_blocks = 1;
    W_1_target = reshape(W_1_target, m, 1, num_y);
else
    target_2d  = false;
    num_blocks = size(W_1_target, 2);
    num_y      = size(W_1_target, 3);
end

W_1_base = double(W_1_base(:));
W_1_target = double(W_1_target);
Vs_base  = double(Vs_base);
Vis_base = double(Vis_base);
R_base   = Vs_base - diag(W_1_base);

% dmperm identifies components in one family template.
pattern = spones(sparse(R_base)) | speye(m);
[perm, ~, borders, ~] = dmperm(pattern);
component_size = diff(borders);
n_components   = num_blocks * length(component_size);

if n_components == 1
    % General connected covariance retain a sparse direct-solve fallback.
    target = reshape(W_1_target, m, num_y);
    Vis    = zeros(m, m, num_y);
    R_sparse = sparse(R_base);
    I_sparse = speye(m);
    for jj = 1:num_y
        if all(target(:,jj) == W_1_base)
            Vis(:,:,jj) = Vis_base;
            continue;
        end
        V_target = R_sparse + spdiags(target(:,jj), 0, m, m);
        V_inverse = V_target \ I_sparse;
        if any(~isfinite(nonzeros(V_inverse)))
            V_inverse = pinv(full(V_target));
        end
        Vis(:,:,jj) = full((V_inverse + V_inverse') / 2);
    end
    allV = reshape(Vis, m^2, num_y);
    return;
end

Vis = zeros(m, m, num_blocks, num_y);

% Scalar components require no matrix construction or factorization.
scalar_components = find(component_size == 1);
for cc = scalar_components
    local_idx  = perm(borders(cc));
    value      = R_base(local_idx,local_idx) + ...
                 reshape(W_1_target(local_idx,:,:), num_blocks, num_y);
    inverse_value = 1 ./ value;
    inverse_value(value == 0) = 0;
    Vis(local_idx,local_idx,:,:) = reshape(inverse_value, 1, 1, num_blocks, num_y);
end

% Equal-sized multidimensional components are solved together
block_sizes = unique(component_size(component_size > 1))';
for block_size = block_sizes
    components = find(component_size == block_size);
    n_template  = length(components);
    local_indices = cell(n_template, 1);
    R_template    = zeros(block_size, block_size, n_template);
    for gg = 1:n_template
        local_indices{gg} = perm(borders(components(gg)):borders(components(gg)+1)-1);
        R_template(:,:,gg) = R_base(local_indices{gg},local_indices{gg});
    end
    I_block = eye(block_size);
    y_chunk     = max(1, floor(page_limit / (n_template * num_blocks)));
    for first_y = 1:y_chunk:num_y
        last_y = min(first_y + y_chunk - 1, num_y);
        y_idx  = first_y:last_y;
        num_y_chunk = length(y_idx);
        n_instance = num_blocks * num_y_chunk;
        n_page = n_template * n_instance;
        A_page = zeros(block_size, block_size, n_page);
        for gg = 1:n_template
            local_idx = local_indices{gg};
            page_idx  = gg:n_template:n_page;
            A_page(:,:,page_idx) = repmat(R_template(:,:,gg), ...
                                          1, 1, n_instance);
            for kk = 1:block_size
                A_page(kk,kk,page_idx) = A_page(kk,kk,page_idx) + ...
                    reshape(W_1_target(local_idx(kk),:,y_idx), ...
                            1, 1, n_instance);
            end
        end
        inverse_page = pagemldivide(A_page, I_block);
        bad_page = reshape(any(any(~isfinite(inverse_page), 1), 2), 1, []);
        for page_idx = find(bad_page)
            inverse_page(:,:,page_idx) = pinv(A_page(:,:,page_idx));
        end
        for gg = 1:n_template
            local_idx = local_indices{gg};
            page_idx  = gg:n_template:n_page;
            Vis(local_idx,local_idx,:,y_idx) = ...
                reshape(inverse_page(:,:,page_idx), block_size, block_size, ...
                        num_blocks, num_y_chunk);
        end
    end
end

% Preserve the exact base inverse when no diagonal changed.
unchanged = reshape(all(W_1_target == reshape(W_1_base, m, 1, 1), 1), ...
                    num_blocks, num_y);
for jj = 1:num_y
    for bb = find(unchanged(:,jj))'
        Vis(:,:,bb,jj) = Vis_base;
    end
end

Vis = (Vis + permute(Vis, [2 1 3 4])) / 2;
if target_2d
    allV = reshape(Vis, m^2, num_y);
else
    allV = reshape(Vis, m^2, num_blocks, num_y);
end

end
