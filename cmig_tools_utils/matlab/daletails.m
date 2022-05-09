function [pd yvec ylvec yuvec xvec] = daletails(x,pl,pu,distrib,xvec)

if ~exist('distrib','var') | isempty(distrib)
  distrib = 'weibull';
end

if ~exist('xvec','var') || isempty(xvec)
%  xvec = colvec(linspace(0,5*max(x),1001));
  xvec = colvec(linspace(0,2*max(x),1000));
end

xvec = colvec(xvec);

hc = colvec(histc(x,xvec)); hc = [0; hc(1:end-1)]; chc = cumsum(hc);
[fvec pci] = binofit(chc,sum(hc)); yvec = 1-fvec; ylvec = 1-pci(:,1); yuvec = 1-pci(:,2);

%figure(667); plot(xvec,-log10(yvec),xvec,-log10(ylvec),xvec,-log10(yuvec),'LineWidth',2);


%figure(668); plot(log(xvec(ivec_fit)),log(-log(yvec(ivec_fit))),log(xvec(ivec_fit)),log(-log(ylvec(ivec_fit))),log(xvec(ivec_fit)),log(-log(yuvec(ivec_fit))),'LineWidth',2);

xvec_tmp = colvec(log(xvec)); 
wvec = colvec((log(-log(yuvec))-log(-log(ylvec))).^-2);
xl = prctile(x,100*pl); xu = prctile(x,100*pu);
ivec_fit = find(xvec>=xl&xvec<=xu&isfinite(wvec));
switch lower(distrib) % Should use common code for cases below
  case {'weibull' 'wbl'} % This is broken!!
    yvec_tmp = colvec(-log(yvec)); 
    parvec0 = wblfit(x(x>0));
%    costfun = @(x)sum(wvec(ivec_fit).*(-log(wblcdf(xvec(ivec_fit),x(1),x(2),'upper'))-yvec_tmp(ivec_fit)).^2);
    costfun = @(x)sum(wvec(ivec_fit).*(((yvec_tmp(ivec_fit(1))--log(wblcdf(xvec(ivec_fit(1)),x(1),x(2),'upper')))-log(wblcdf(xvec(ivec_fit),x(1),x(2),'upper')))-yvec_tmp(ivec_fit)).^2); % Must allow for an offset
    parvec = fminsearch(costfun,parvec0,optimset('MaxFunEvals',1000));
    pd_fit = makedist('wbl','A',parvec(1),'B',parvec(2));
  case {'beta' 'logbeta'}
    yvec_tmp = colvec(log(-log(yvec))); 
    parvec0 = betafit(10.^(-x));
    costfun = @(x)sum(wvec(ivec_fit).*(-log(betacdf(10.^(-xvec(ivec_fit)),x(1),x(2)))-yvec_tmp(ivec_fit)).^2);
    parvec = fminsearch(costfun,parvec0,optimset('MaxFunEvals',1000));
%    yvec_fit = betacdf(10.^(-xvec),parvec(1),parvec(2));
    pd_fit = makedist('beta','A',parvec(1),'B',parvec(2));
  case {'sidak' 'minp'}
    yvec_tmp = colvec(log(-log(yvec))); 
    parvec0 = betafit(10.^(-x)); parvec0 = parvec0(2); 
    costfun = @(x)sum(wvec(ivec_fit).*(-log(betacdf(10.^(-xvec(ivec_fit)),1,x(1)))-yvec_tmp(ivec_fit)).^2);
    parvec = fminsearch(costfun,parvec0,optimset('MaxFunEvals',1000));
%    yvec_fit = betacdf(10.^(-xvec),1,parvec(1));
    pd_fit = makedist('beta','A',1,'B',parvec(1));
  case {'gamma'}
    yvec_tmp = colvec(-log(yvec)); 
    parvec0 = gamfit(x(x>0));
%    costfun = @(x)sum(wvec(ivec_fit).*(-log(gamcdf(xvec(ivec_fit),x(1),x(2),'upper'))-yvec_tmp(ivec_fit)).^2);
    costfun = @(x)sum(wvec(ivec_fit).*(((yvec_tmp(ivec_fit(1))--log(gamcdf(xvec(ivec_fit(1)),x(1),x(2),'upper')))-log(gamcdf(xvec(ivec_fit),x(1),x(2),'upper')))-yvec_tmp(ivec_fit)).^2); % Must allow for an offset
    parvec = fminsearch(costfun,parvec0,optimset('MaxFunEvals',1000));
    pd_fit = makedist('gamma','A',parvec(1),'B',parvec(2));
  case {'chi2'}
    yvec_tmp = colvec(-log(yvec)); 
    parvec0 = gamfit(x); parvec0 = parvec0(1);
%    costfun = @(x)sum(wvec(ivec_fit).*(-log(gamcdf(xvec(ivec_fit),x(1),2,'upper'))-yvec_tmp(ivec_fit)).^2); % Must allow for an offset
    costfun = @(x)sum(wvec(ivec_fit).*(((yvec_tmp(ivec_fit(1))--log(gamcdf(xvec(ivec_fit(1)),x(1),2,'upper')))-log(gamcdf(xvec(ivec_fit),x(1),2,'upper')))-yvec_tmp(ivec_fit)).^2); % Must allow for an offset
    parvec = fminsearch(costfun,parvec0,optimset('MaxFunEvals',1000));
    pd_fit = makedist('gamma','A',parvec(1),'B',2);
  case {'extremevalue'} % This is totally broken!
    parvec0 = evfit(x);
    parvec = parvec0;
    pd_fit = makedist('extremevalue','mu',parvec(1),'sigma',parvec(2));
end
yvec_fit = (yvec(ivec_fit(1))/pd_fit.cdf(xvec(ivec_fit(1)),'upper'))*pd_fit.cdf(xvec,'upper');

if 0
  figure(669); plot(xvec(ivec_fit),-log10(yvec(ivec_fit)),xvec(ivec_fit),-log10(yvec_fit(ivec_fit)),'LineWidth',2);
end

% Should update costfun to extrapolate from ecdf @ xl -- check into how paretotails is defined

if 0
  figure(669); plot(log(xvec(ivec_fit)),log(-log(yvec(ivec_fit))),log(xvec(ivec_fit)),log(-log(yvec_fit(ivec_fit))),'LineWidth',2);
end

yvec_pred = yvec_fit;
if pl>0
  blendvec = interp1([xvec(1) xvec(ivec_fit(1)) xvec(ivec_fit(end)) xvec(end)],[0 0 1 1],xvec,'linear','extrap');
  ivec_blend = find(blendvec<1);
  yvec_pred(ivec_blend) = exp(-exp((1-blendvec(ivec_blend)).*log(-log(yvec(ivec_blend)))+blendvec(ivec_blend).*log(-log(yvec_fit(ivec_blend)))));
  yvec_pred(1:(min(find(isfinite(yvec_pred)))-1)) = 1;
end

if 0
  figure(670); plot(xvec,-log10(yvec),xvec,-log10(yvec_pred),'LineWidth',2);
end

pd = struct;
pd.pd_fit = pd_fit;
pd.nlogcdf = @(x)exp(interp1(xvec,log(-log(yvec_pred)),x,'linear','extrap'));
pd.cdf = @(x,varargin)condexp(nargin>1&&strcmp(varargin{1},'upper'),exp(-pd.nlogcdf(x)),1-exp(-pd.nlogcdf(x))); % This should be replaced with analytic form
%pd.cdf = @(x,varargin)pd_fit.cdf(x,varargin{:}); % Use analytic form

if 0
  figure(671); plot(xvec,-log10(yvec),xvec,-log10(pd.cdf(xvec,'upper')),'LineWidth',2);
end

return

% ToDo
%   Make gamma and weibull code consistent, including wvec fix

