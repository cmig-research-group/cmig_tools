function ymat_bootstrap = generate_null_sample(X, iid, fid, agevec, ...
                           beta_hat_null, sig2mat_null, clusterinfo, RandomEffects)
nobs = size(X,1);
ymat_bootstrap = zeros(nobs, 1);

y_RFX_sim = nan(nobs, 1);
nfam = length(clusterinfo);
for fi = 1:nfam
    tmp = 0;
    for ri = 1:(length(RandomEffects)-1) % 排除error term
        tmp = tmp + sqrt(sig2mat_null(ri)) .* ...
                    mvnrnd(zeros(length(clusterinfo{fi}.jvec_fam), 1), ...
                    double(getfield(clusterinfo{fi}, ...
                    sprintf('V_%s',RandomEffects{ri}))),1)';
    end
    y_RFX_sim(clusterinfo{fi}.jvec_fam,:) = tmp;
end

eta_sim = X * beta_hat_null + y_RFX_sim;

p_sim = 1 ./ (1 + exp(-eta_sim));

ymat_bootstrap = rand(nobs, 1) < p_sim;
end