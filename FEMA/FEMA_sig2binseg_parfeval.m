function [beta_hat, beta_se, betacon_hat, betacon_se] = FEMA_sig2binseg_parfeval(X, ymat, contrasts, clusterinfo, binvec, sig2mat, sig2tvec, RandomEffects, GroupByFamType, nfamtypes, famtypevec, OLSflag, SingleOrDouble)

% Initialize
nfam        = length(clusterinfo);
nsig2bins   = max(binvec); 
tvec_bins   = zeros(nsig2bins, 1);
beta_hat    = zeros(size(X, 2), size(ymat,2), class(ymat)); 
beta_se     = zeros(size(beta_hat), class(ymat)); 
[betacon_hat, betacon_se]             = deal(zeros(size(contrasts,1), size(ymat,2), class(ymat))); 
[Vs_famtype, Vis_famtype, Ws_famtype] = deal(cell(1, nfamtypes));
[Vs_fam,     Vis_fam,     Ws_fam]     = deal(cell(1, nfam));

% Check if lsqminnorm exists
if exist('lsqminnorm', 'file')
    doLSQ = true;
else
    doLSQ = false;
end

% Make sure that X is not rank deficient
if rank(X) < size(X, 2)
    lowRank = true;
else
    lowRank = false;
end

for sig2bini = unique(binvec(isfinite(binvec)), 'stable')
  t0         = now; %#ok<*TNOW1>
  ivec_bin   = find(binvec == sig2bini);
  sig2vec    = mean(sig2mat(:, ivec_bin), 2);
  
  if ~isempty(ivec_bin)
      if GroupByFamType % Compute Vs and Vis by family type
          for fi = 1:nfamtypes
              ivec = find(famtypevec==fi);
              Vs_famtype{fi} = 0;
              for ri = 1:length(RandomEffects)
                  Vs_famtype{fi} = Vs_famtype{fi} + sig2vec(ri) * clusterinfo{ivec(1)}.(['V_', RandomEffects{ri}]);
              end
              Vis_famtype{fi} = cast(pinv(double(Vs_famtype{fi})), SingleOrDouble);
              if OLSflag
                  Ws_famtype{fi} = eye(size(Vis_famtype{fi}), class(Vis_famtype{fi}));
              else
                  Ws_famtype{fi} = Vis_famtype{fi};
              end
          end
      else % Compute Vs and Vis for each family
          for fi = 1:nfam
              Vs_fam{fi} = 0;
              for ri = 1:length(RandomEffects)
                  Vs_fam{fi} = Vs_fam{fi} + sig2vec(ri) * clusterinfo{fi}.(['V_', RandomEffects{ri}]);
              end
              Vis_fam{fi} = cast(pinv(double(Vs_fam{fi})), SingleOrDouble);
              if OLSflag
                  Ws_fam{fi} = eye(size(Vis_fam{fi}), class(Vis_fam{fi}));
              else
                  Ws_fam{fi} = Vis_fam{fi};
              end
          end
      end

    XtW     = zeros(fliplr(size(X)), class(X)); 
    % XtWVsWt = zeros(fliplr(size(X)), class(X));

    if GroupByFamType
        for fi = 1:length(clusterinfo)
            XtW(:,clusterinfo{fi}.jvec_fam)     = X(clusterinfo{fi}.jvec_fam,:)'  * Ws_famtype{famtypevec(fi)};
            % XtWVsWt(:,clusterinfo{fi}.jvec_fam) = XtW(:,clusterinfo{fi}.jvec_fam) * Vs_famtype{famtypevec(fi)} * Ws_famtype{famtypevec(fi)}';
        end
    else
        for fi = 1:length(clusterinfo)
            XtW(:,clusterinfo{fi}.jvec_fam)     = X(clusterinfo{fi}.jvec_fam,:)'    * Ws_fam{fi};
            % XtWVsWt(:,clusterinfo{fi}.jvec_fam) = XtW(:,clusterinfo{fi}.jvec_fam)   * Vs_fam{fi} * Ws_fam{fi}';
        end
    end

    % XtWVsWtX = XtWVsWt*X; 
    % B        = XtWVsWtX; 
    B  = XtW * X;
    if lowRank
        if doLSQ
            Bi = lsqminnorm(B, eye(size(B)));
        else
            Bi = pinv(B); % Note that Cov_beta = Bi = pinv(XtWVsWtX)*XtWVsWt*pinv(XtWVsWtX)
        end
    else
        Bi = B \ eye(size(B));
    end

    if OLSflag
        if lowRank
            if doLSQ
                beta_hat_tmp = lsqminnorm(X, ymat(:, ivec_bin));
            else
                beta_hat_tmp = pinv(X) * ymat(:,ivec_bin);
            end
        else
            beta_hat_tmp = X \ ymat(:, ivec_bin);
        end
    else
        beta_hat_tmp = Bi * XtW * ymat(:,ivec_bin); % Should generalize this to work with arbitrary W matrices
    end

    Cov_beta             = nearestSPD(Bi);
    beta_hat(:,ivec_bin) = beta_hat_tmp; % These are some of the most time consuming operations -- use MEX?
    beta_se(:, ivec_bin) = sqrt(diag(Cov_beta) * sig2tvec(ivec_bin));
  end

  for ci = 1:size(contrasts,1)
    betacon_hat(ci,ivec_bin) = contrasts(ci,:) * beta_hat(:,ivec_bin); 
    betacon_se(ci,ivec_bin)  = sqrt(contrasts(ci,:) * Cov_beta * contrasts(ci,:)' * sig2tvec(ivec_bin));
  end
  tvec_bins(sig2bini) = (now-t0) * 3600*24; % Time in seconds
end

% fprintf(1,'Elapsed time %0.2f seconds\n',sum(tvec_bins));

% ToDo
%   Fix bug w. all zeros in beta_hat and beta_se, when GroupByFamType  == false 