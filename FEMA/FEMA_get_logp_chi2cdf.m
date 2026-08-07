function log10p = FEMA_get_logp_chi2cdf(Wald, df, tail, extend)
% Function to compute unsigned -log10p values, given Wald statistics and df
%% Inputs:
% Wald:         numeric     vector or matrix of Wald statistics
% 
% df:           numeric     numerator degrees of freedom
% 
% tail:         character   should be one of the following:
%                               * 'upper' (default)
%                               * 'lower'
% 
% extend:       logical     indicates whether to use extended tail method
%                           set to false to reproduce previous
%                           functionality of using chi2cdf for calculating
%                           p values - this will return Inf for small p
%                           values; specifically, if p is lower than
%                           realmin (2.2251e-308) or 307.6527 on a log10
%                           scale, then these values will be truncated to
%                           Inf by chi2cdf while the extend option is
%                           designed to handle those
%                           
%% Output:
% log10p:       numeric     unsigned -log10(p) values
%
%% Test:
% Points to the right of the dashed line are chi2cdf degrading in the
% denormal range, not the extended method
% df = [1; 100; 1000; 10000];
% for dd = 1:numel(df)
%     Wald = linspace(0, 2*df(dd) + 2000, 20000)';
%     ref  = FEMA_get_logp_chi2cdf(Wald, df(dd), 'upper', false);
%     new  = FEMA_get_logp_chi2cdf(Wald, df(dd), 'upper', true);
% 
%     fin  = isfinite(ref);
%     good = fin & (new < 307.65);
% 
%     subplot(1, numel(df), dd);
%     d = abs(ref(fin) - new(fin));
%     semilogy(new(fin), max(d, 1e-16), '.');
%     xline(307.65, 'k--'); axis tight;
%     xlabel('-log10(p)'); ylabel('abs diff'); title(['df = ' num2str(df(dd))]);
% 
%     disp(['df ', num2str(df(dd)), ...
%         ': max abs diff where chi2cdf is valid: ', num2str(max(abs(ref(good) - new(good)))), ...
%         '; chi2cdf Inf: ', num2str(sum(~fin)), ...
%         '; extended non-finite: ', num2str(sum(~isfinite(new)))]);
% end
% 
% % df = 1 upper tail is the two-tailed normal p, via gammainc instead of erfcx
% z = logspace(0, 2.5, 1000)';
% disp(['df=1 vs normcdf, max relative difference: ', ...
%     num2str(max(abs(FEMA_get_logp_chi2cdf(z.^2, 1, 'upper') - ...
%     FEMA_get_logp_normcdf(z, 'both')) ./ ...
%     FEMA_get_logp_normcdf(z, 'both')))]);

%% Check inputs
if ~exist('Wald', 'var') || isempty(Wald)
    error('Please provide a vector or matrix of z statistics');
else
    if ~isreal(Wald) || any(Wald(:) < 0)
        error('Wald statistics should be real-valued and positive');
    end
end

if ~exist('df', 'var') || isempty(df)
    error('Please provide numerator degrees of freedom');
else
    if ~isscalar(df)
        error('Expected df to be a scalar');
    else
        if ~isreal(df) || df <= 0
            error('df should be real positive scalar');
        end
    end
end

if ~exist('tail', 'var') || isempty(tail)
    tail = 'upper';
else
    if iscell(tail)
        tail = tail{1};
    end
    if ~(ischar(tail) || isstring(tail))
        error('tail should be a character vector');
    end
    tail = lower(char(tail));
    if ~ismember(tail, {'upper', 'lower'})
        error('tail should be either upper or lower');
    end
end

if ~exist('extend', 'var') || isempty(extend)
    extend = true;
else
    if ~isscalar(extend) || ~(islogical(extend) || (isnumeric(extend) && ismember(extend, [0 1])))
        error('extend should be either true or false');
    end
end

%% Compute unsigned -log10(p) values
if extend
    switch tail
        case 'upper'
            y      = Wald/2;
            a      = df/2;
            Q      = gammainc(y, a, 'upper');
            log10p = -log10(Q);
            ext    = ~(Q > 1e-300);
            if any(ext(:))
                log10p(ext) = (y(ext) + gammaln(a+1) - a.*log(y(ext)))/log(10) ...
                                                     - log10(gammainc(y(ext), a, 'scaledupper'));
            end

        case 'lower'
            y      = Wald/2;
            a      = df/2;
            P      = gammainc(y, a, 'lower');
            log10p = -log10(P);
            ext    = ~(P > 1e-300);
            if any(ext(:))
                log10p(ext) = (y(ext) + gammaln(a+1) - a.*log(y(ext)))/log(10) ...
                                                     - log10(gammainc(y(ext), a, 'scaledlower'));
            end
    end
else
    switch tail
        case 'upper'
            log10p = -log10(chi2cdf(Wald, df, 'upper'));

        case 'lower'
            log10p = -log10(chi2cdf(Wald, df));
    end
end