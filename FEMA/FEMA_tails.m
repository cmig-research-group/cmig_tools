function [nlog10cdfvec_ecdf nlog10cdfvec_dale nlog10cdfvec_pareto pd_dale pd_pareto xvec] =  FEMA_tails(statvec,pl,pu,xvec,distrib,fh)

statvec = statvec(isfinite(statvec)); % Shouldn't be needed

if ~exist('pl','var') || isempty(pl)
  pl = 0.99;
end

if ~exist('pu','var') || isempty(pu)
  pu = 1-1/(length(statvec)/10);;
end

if ~exist('xvec','var') || isempty(xvec)
  xvec = colvec(linspace(0,2*max(statvec),1000));
end

[pd_dale yvec ylvec yuvec] = daletails(statvec,pl,pu,distrib,xvec);

pd_pareto = paretotails(statvec,0,0.99); % Check what standard parameters are in PALM
%pd_pareto = paretotails(statvec,0,0.999); % Check what standard parameters are in PALM

yvals2 = cat(2,yvec,ylvec,yuvec);
nlog10cdfvec_ecdf = -log10(yvals2);
nlog10cdfvec_dale = -log10(pd_dale.cdf(xvec,'upper'));
nlog10cdfvec_pareto = -log10(fixed_paretotails_cdf(pd_pareto,xvec));

if exist('fh','var') && ~isempty(fh)
  if strcmp(class(fh),'matlab.graphics.axis.Axes')
    axes(fh);
  else
    sfigure(fh); 
  end
  plot(xvec,nlog10cdfvec_ecdf(:,1),xvec,nlog10cdfvec_ecdf(:,2:3),'k:',xvec,nlog10cdfvec_dale,xvec,nlog10cdfvec_pareto,'LineWidth',2);
  h=xlabel('test statistic'); h=ylabel('-log_1_0(1-cdf)'); h=legend({'ecdf' 'flo' 'fup' distrib 'paretotails'},'Location','NW');
  drawnow;
end