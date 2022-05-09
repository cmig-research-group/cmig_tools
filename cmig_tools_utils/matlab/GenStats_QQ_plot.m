function [logqmat hv_logp logqmat_ci] = GenStats_QQ_plot(logpmat,hv_z,ci_alpha)

if exist('hv_z','var') & ~isempty(hv_z)
  hv_logp = -log10(2*normcdf(-abs(hv_z)));
else
  hv_logp = linspace(0,200,10000);
end

if exist('ci_alpha','var') & isempty(ci_alpha), ci_alpha = 0.05; end

for j = 1:size(logpmat,2)
  logpvec = logpmat(:,j);
  hc = hist(logpvec,hv_logp);
  chc = cumsum(hc)/sum(hc);
  if j==1, logqmat = NaN(length(chc),size(logpmat,2)); logqmat_ci = NaN(length(chc),size(logpmat,2),2); end
  logqmat(:,j) = -log10(1-chc);
  if exist('ci_alpha','var')
    [phat pci] = binofit(cumsum(hc),sum(hc)*ones(size(hc)),ci_alpha); % Compute confidence interval, assuming independence (no LD)
    logqmat_ci(:,j,1) = -log10(1-pci(:,1));
    logqmat_ci(:,j,2) = -log10(1-pci(:,2));
  end
end
if nargout==0
  plot([0 7],[0 7],'k:',logqmat,hv_logp,'LineWidth',2); xlim([0 7]); ylim([0 10]); h=xlabel('Empirical -log_1_0(q)'); set(h,'FontSize',20); h=ylabel('Nominal -log_1_0(p)'); set(h,'FontSize',20);
  if exist('ci_alpha','var')
    hold on;
    plot(logqmat_ci(:,1),hv_logp,'k:',logqmat_ci(:,2),hv_logp,'k:','LineWidth',2); xlim([0 7]); ylim([0 10]); h=xlabel('Empirical -log_1_0(q)'); set(h,'FontSize',20); h=ylabel('Nominal -log_1_0(p)'); set(h,'FontSize',20);
  end
end

