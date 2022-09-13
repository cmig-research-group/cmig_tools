function loglike = FEMA_logLik(sig2vec,X,yvec_res,clusterinfo,Ss)

% Calculate log-likelihood with respect to random random effects.
% Higher values correspont to better fit of the model to the data.
% Note that some previous versions of this function returned negative log
% likelihood (aka cost function - lower is better).
%
% sig2vec         - vector of random effects (variance parameters) for a given phenotype;
%                   typically, sig2vec = sig2mat(:,i) * sig2tvec(i), i - phenotype index
% yvec_res        - vector of residuals
% clusterinfo, Ss - output of FEMA_parse_family
%
%
% FEMA_logLik doesn't have RandomEffects argument
% instead, it works on assumption that the order of random effects in 'x' is the same as in 'Ss'

assert(length(sig2vec) == length(Ss));

Sigma = sig2vec(1) * Ss{1};
for ri = 2:length(Ss)
  Sigma = Sigma + sig2vec(ri) * Ss{ri};
end

loglike = 0;
for fi = 1:length(clusterinfo)
  jvec_fam = clusterinfo{fi}.jvec_fam;
  Sigma_fam = full(Sigma(jvec_fam, jvec_fam));
  Sigma_fam = (Sigma_fam + Sigma_fam')/2;
  loglike = loglike + log(mvnpdf(double(yvec_res(jvec_fam)),0,Sigma_fam));
  
  % ReML adjustment would be like this:
  % Newton-Raphson and EM Algorithms for Linear Mixed-Effects Models for Repeated-Measures Data 
  % MARY J. LINDSTROM and DOUGLAS M. BATES
  % https://www.tandfonline.com/doi/pdf/10.1080/01621459.1988.10478693
  % 
  % Why does REML adjustment scale with X ?
  % X_fam = X(jvec_fam, :);
  % loglike = loglike - 0.5 * log(det(X_fam' * pinv(Sigma_fam) * X_fam));  
end
