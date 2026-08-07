function log10p = FEMA_get_logp_normcdf(zmat, tail, extend)
% Function to compute unsigned -log10p values, given z statistics
%% Inputs:
% zmat:         numeric     vector or matrix of z statistics (z = beta/SE)
% 
% tail:         character   should be one of the following:
%                               * 'both'
%                               * 'left'
%                               * 'right'
% 
% extend:       logical     indicates whether to use extended tail method
%                           (relevant for |z| > 37; default = true); set to
%                           false to reproduce previous functionality of
%                           using normcdf for calculating p values - this
%                           will return Inf for |z| > 38.475, although it
%                           is recommended to not use normcdf beyond |z| =
%                           37 due to numerical reasons
%                           
%% Output:
% log10p:       numeric     unsigned -log10(p) values
% 
%% Notes:
% In case of one tailed test (left or right), large z values in the
% opposite direction would result in a -log10(p) value of approximately (or
% exactly) zero
% 
% To recover signed log10 p values for two-tailed test, do:
% signed_log10p = -sign(zmat) .* FEMA_get_logp_normcdf(zmat, 'both');
% 
% Note that signed log10 p values are not applicable for one tailed test
%
%% Test: show equality of approaches for |z| <= 37
% zmat  = linspace(-37, 37, 100000)';
% tails = {'both', 'left', 'right'};
% 
% for tt = 1:3
%     ref = FEMA_get_logp_normcdf(zmat, tails{tt}, false);
%     new = FEMA_get_logp_normcdf(zmat, tails{tt}, true);
% 
%     subplot(1,3,tt);
%     scatter(abs(zmat), abs(ref - new), '.');
%     set(gca, 'YScale', 'log');
%     axis tight; xlabel('|z|'); ylabel('abs diff p value'); title(tails{tt});
% 
%     disp([tails{tt} ': max abs. difference: ', num2str(max(abs(ref-new)))]);
% end
%

%% Check inputs
if ~exist('zmat', 'var') || isempty(zmat)
    error('Please provide a vector or matrix of z statistics');
else
    if ~isreal(zmat)
        error('zmat contains complex numbers');
    end
end

if ~exist('tail', 'var') || isempty(tail)
    tail = 'both';
else
    if ~(ischar(tail) || isstring(tail))
        error('tail should be a character vector');
    end
    tail = lower(char(tail));
    if ~ismember(tail, {'both', 'left', 'right'})
        error('tail should be one of: both, left, or right');
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
        case 'both'
            zmat   = abs(zmat);
            log10p = 0.5 * zmat.^2 / log(10) - log10(erfcx(zmat / sqrt(2)));

        case 'left'
            log10p = zeros(size(zmat), 'like', zmat);

            % For negative values
            locNeg         = zmat <= 0;
            log10p(locNeg) = log10(2) + 0.5 * zmat(locNeg).^2 / log(10) - log10(erfcx(-zmat(locNeg) / sqrt(2)));

            % For other values
            log10p(~locNeg) = -log10(normcdf(zmat(~locNeg)));

        case 'right'
            log10p = zeros(size(zmat), 'like', zmat);

            % For positive values
            locPos         = zmat >= 0;
            log10p(locPos) = log10(2) + 0.5 * zmat(locPos).^2 / log(10) - log10(erfcx(zmat(locPos) / sqrt(2)));

            % For other values
            log10p(~locPos) = -log10(normcdf(-zmat(~locPos)));
    end
else
    switch tail
        case 'both'
            log10p = -log10(normcdf(-abs(zmat))*2);

        case 'left'
            log10p = -log10(normcdf(zmat));

        case 'right'
            log10p = -log10(normcdf(-zmat));
    end
end