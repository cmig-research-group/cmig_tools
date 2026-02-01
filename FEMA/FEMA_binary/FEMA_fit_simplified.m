function [sig2tvec,  sig2mat] =                                  ...
          FEMA_fit_simplified(X, iid, eid, fid, ymat_res, sig2tvec,  ...
                   pihatmat, W_1, varargin)

p = inputParser;
% addParamValue(p,'ciflag', false);
addParamValue(p,'MLflag', false);
addParamValue(p,'FamilyStruct', {});
addParamValue(p,'NonnegFlag', true);

parse(p,varargin{:})

MLflag       = p.Results.MLflag;
FamilyStruct = p.Results.FamilyStruct;
NonnegFlag   = p.Results.NonnegFlag;
% ciflag       = p.Results.ciflag;
clusterinfo = FamilyStruct.clusterinfo;
Ss          = FamilyStruct.Ss;
subvec1     = FamilyStruct.subvec1;
subvec2     = FamilyStruct.subvec2;


LHS      = ymat_res(subvec1,:) .* ymat_res(subvec2,:) ./ mean(ymat_res.^2,1); % use normalized residuals
% LHS      = ymat_res(subvec1,:) .* ymat_res(subvec2,:);
sig2mat  = NaN(size(Ss,2), size(ymat_res, 2));
theta_cov = cell(size(ymat_res, 2));

% heterogeneity variance is considered, showed as E ~ N(0, diag(W)^{-1}).
M           = FamilyStruct.M;
subvec_e    = find(M(:,end)); % diagonal position of variance matrix for y 
LHS(subvec_e,:)  = LHS(subvec_e,:) - M(subvec_e,end) .*  W_1 / sig2tvec;
 
% no loop
for coli=1:size(ymat_res, 2)

    % Use new version of lsqnonneg_amd to enfoce non-negative variances
    sig2mat_tmp     = lsqnonneg_amd3(M(:,1:end-1),LHS(:,coli));
    sig2mat(:,coli) = [sig2mat_tmp ; 1/sig2tvec] * sig2tvec;


end

%% Using maximum likelihood solution
% if MLflag % Phenotypes should be pre-normalized! -- now, scale is all over the place
%     options_fmincon = optimoptions('fmincon','Display','off');
%     logLikvec       = nan(1,size(ymat_res,2));
%     [sig2mat_ml, sig2mat_ll, sig2mat_ul] = deal(nan(size(sig2mat)));
%     disp(var(ymat_res));
% 
%     for coli=1:size(ymat_res, 2)
%         f = @(x) (-1 * FEMA_logLik(exp(x), X, ymat_res(:, coli), clusterinfo, Ss));
%         g = @(x) (-1 * FEMA_logLik(x,      X, ymat_res(:, coli), clusterinfo, Ss));
%         sig2vec0 = double(sig2mat(:, coli) * sig2tvec(coli));
%         fprintf(1,'Optimizing using fmincon\n');
%         tic
%         [sig2vec_hat, cost, exitflag, output] = fmincon(g, sig2vec0, [], [], [], [], 0*ones(size(sig2vec0)), [], [], options_fmincon);
%         toc
%         if 1 % exitflag<0
%             fprintf(1,'fmincon exited with exitflag = %d:\n',exitflag);
%             disp(output)
%         end
%         if ciflag % Compute confidence intervals on random effects?
%             fprintf(1,'Computing Confidence Interval\n');
%             tic
%             loglikthresh = chi2inv(1-0.05/2,1)/2;
%             [sig2vec_ll, sig2vec_ul] = deal(nan(size(sig2vec0)));
%             for ri = 1:length(sig2vec_hat)
%                 ivec   = double((1:length(sig2vec_hat))==ri);
%                 tmpfun = @(x) g(fmincon(g, (1-ivec)' .* sig2vec_hat + ivec' * x, [], [], ivec, sig2vec_hat(ri)+x, 0*ones(size(sig2vec0)), [], [], options_fmincon))-cost; % allow other parameters to change
%                 dx0    = 0.01 * sum(sig2vec_hat);
%                 y0     = tmpfun(dx0); % Hack to scale initial step size by variance -- phenotypes should be pre-normalized (unity variance)
%                 if y0 < 0.2 % Increase scale if change in cost is too small
%                     dx0 = dx0 * 4;
%                     y0  = tmpfun(dx0);
%                 elseif y0 > 6 % Decrease scale if change in cost is too large
%                     dx0 = dx0/4;
%                     y0  = tmpfun(dx0);
%                 end
%                 dx1  = dx0 * sqrt(2/y0);
%                 y1   = tmpfun(dx1); % Get scale -- should increase dx0, if y0 is too small (or negative)
%                 x    = [0 max([dx0 dx1]) * [0.5 1]];
%                 y    = [0 tmpfun(x(2)) tmpfun(x(3))];
%                 p    = polyfit(x, y, 2);
%                 xvec = linspace(0, max(x), 101);
%                 yvec = polyval(p, xvec);
% 
%                 figure(coli*10);
%                 subplot(length(sig2vec_hat), 2, (ri-1)*2+2);
%                 plot(xvec, yvec, x, y, '*', 'lineWidth', 2);
%                 drawnow;
%                 ul = sig2vec_hat(ri) + (-p(2)+((p(2)^2-4*p(1)*(p(3)-loglikthresh)))^0.5)/(2*p(1));
%                 ll = 0;
% 
%                 if sig2vec_hat(ri) > 0.01 * sum(sig2vec_hat)
%                     if x(end) > sig2vec_hat(ri)
%                         x = x*sig2vec_hat(ri)/x(end);
%                     end
%                     x     = -x;
%                     y     = [0 tmpfun(x(2)) tmpfun(x(3))];
%                     p     = polyfit(x, y, 2);
%                     xvec  = linspace(min(x), max(x), 101);
%                     yvec  = polyval(p, xvec);
% 
%                     figure(coli*10);
%                     subplot(length(sig2vec_hat), 2, (ri-1)*2+1);
%                     plot(xvec, yvec, x, y, '*', 'lineWidth', 2);
%                     drawnow;
% 
%                     ll = max(0,sig2vec_hat(ri) + (-p(2)-((p(2)^2-4*p(1)*(p(3)-loglikthresh)))^0.5)/(2*p(1)));
%                 end
% 
%                 sig2vec_ll(ri) = ll;
%                 sig2vec_ul(ri) = ul;
%                 fprintf(1,'ri=%d: ll=%f ul=%f (%s)\n',ri, ll, ul, char(datetime));
%                 if ~isreal(ll+ul) || ~isfinite(ll+ul) % Stop if result is imaginary or not finite
%                     fprintf(1,'Invalid confidence interval estimates\n');
%                 end
%             end
%             toc
%         end
%         % disp(num2str(cost,'%0.6e') )
%         sig2mat_ml(:, coli) = sig2vec_hat;
%         if ciflag
%             sig2mat_ll(:, coli) = sig2vec_ll;
%             sig2mat_ul(:, coli) = sig2vec_ul;
%         end
%         % disp(rowvec(sig2mat_ml(:, coli)/sum(sig2mat_ml(:, coli))))
%         logl_ml  = g(sig2mat_ml(:, coli)); % This takes ~0.13s per column
%         logl_mom = g(double(sig2mat(:, coli) * sig2tvec(coli)));
%         logging('pheno %i of %i, perm %i of %i: loglike(MoM)=%.2f, loglike(ML)=%.2f', coli, size(ymat_res, 2), permi, nperms, logl_mom, logl_ml);
%         logLikvec(coli) = -logl_ml;
%     end
%     sig2tvec_ml = sum(sig2mat_ml);
%     sig2mat_ml  = sig2mat_ml ./ sig2tvec_ml;
%     if ciflag
%         sig2mat_ci = cat(3, sig2mat_ll, sig2mat_ul) ./ sig2tvec_ml;
%     end
%     sig2mat  = sig2mat_ml;
%     sig2tvec = sig2tvec_ml;
% end

end