function [beta_hat, beta_se, betacon_hat, betacon_se, tvec_bins, nvec_bins] =   ...
          FEMA_sig2binseg_parfeval(X, ymat, contrasts, clusterinfo, binvec,     ...
                                   sig2mat, sig2tvec, RandomEffects,            ...
                                   GroupByFamType, nfamtypes, famtypevec,       ...
                                   OLSflag, SingleOrDouble)

% Initialize
nfam        = length(clusterinfo);
nsig2bins   = max(binvec); 
tvec_bins   = zeros(nsig2bins, 1); 
nvec_bins   = NaN(nsig2bins, 1);
Ws_famtype  = cell(1, nfamtypes);
Ws_fam      = cell(1, nfam);
[betacon_hat, betacon_se] = deal(zeros(size(contrasts,1), size(ymat,2), class(ymat)));
[beta_hat,    beta_se]    = deal(zeros(size(X,2), size(ymat,2), class(ymat)));

% Check if lsqminnorm can be used
if exist('lsqminnorm', 'file')
    useLSQ = true;
else
    useLSQ = false;
end

% Check if X is rank deficient
if rank(X) < size(X, 2)
    lowRank = true;
else
    lowRank = false;
end

% Get ordering of fields in clusterinfo - reasonable to assume that fields
% are always ordered in the same way since clusterinfo is created in the
% same way across all clusters
ff           = fieldnames(clusterinfo{1});
[a, b]       = ismember(ff, strcat('V_', RandomEffects));
[~, RFX_ord] = sort(b(a));
locJVec      = strcmpi(ff, 'jvec_fam');

for sig2bini = unique(binvec(isfinite(binvec)), 'stable')
    t0         = now; %#ok<*TNOW1>
    ivec_bin   = find(binvec==sig2bini);
    nvec_bins(sig2bini) = length(ivec_bin);
    sig2vec    = mean(sig2mat(:, ivec_bin), 2);

    if ~isempty(ivec_bin)
        % Handle the case of OLS
        if OLSflag
            XtX  = X' * X;
            if lowRank
                if useLSQ
                    iXtX = lsqminnorm(XtX, eye(size(XtX)));
                else
                    iXtX = pinv(XtX);
                end
            else
                iXtX = XtX \ eye(size(XtX));
            end
            beta_hat(:, ivec_bin) = iXtX * (X' * ymat(:, ivec_bin));
            beta_se(:,  ivec_bin) = sqrt(diag(iXtX) * sig2tvec(ivec_bin));
            Cov_beta              = iXtX;
        else
            if GroupByFamType
                % Compute Vs and Vis by family type
                for fi = 1:nfamtypes
                    ivec       = find(famtypevec == fi);
                    currClus   = struct2cell(clusterinfo{ivec(1)});
                    tmpSize    = length(currClus{locJVec});
                    Vs_famtype = zeros(tmpSize);

                    % Compute V
                    for ri = 1:length(RandomEffects)
                        Vs_famtype = Vs_famtype + sig2vec(ri) * currClus{RFX_ord(ri)};
                        % Vs_famtype = Vs_famtype + sig2vec(ri) * currClus.(['V_', RandomEffects{ri}]);
                    end

                    % Compute inverse of V
                    Vis_famtype = double(Vs_famtype) \ eye(tmpSize, SingleOrDouble);
                    % Vis_famtype = cast(double(Vs_famtype) \ eye(tmpSize), SingleOrDouble);
                    if any(isnan(Vis_famtype))
                        Vis_famtype = cast(pinv(double(Vs_famtype)), SingleOrDouble);
                    end
                    Ws_famtype{fi} = Vis_famtype;
                    % if OLSflag
                    %     Ws_famtype{fi} = eye(tmpSize, SingleOrDouble);
                    % else
                    %     Ws_famtype{fi} = Vis_famtype;
                    % end
                end
            else
                % Compute Vs and Vis for each family
                for fi = 1:nfam
                    currClus   = struct2cell(clusterinfo{fi});
                    tmpSize    = length(currClus{locJVec});
                    % tmpSize  = length(currClus.jvec_fam);
                    Vs_fam     = zeros(tmpSize);

                    % Compute V
                    for ri = 1:length(RandomEffects)
                        Vs_fam = Vs_fam + sig2vec(ri) * currClus{RFX_ord(ri)};
                        % Vs_fam = Vs_fam + sig2vec(ri) * currClus.(['V_', RandomEffects{ri}]);
                    end

                    % Compute inverse of V
                    Vis_fam = double(Vs_fam) \ eye(tmpSize, SingleOrDouble);
                    % Vis_fam = cast(double(Vs_fam) \ eye(tmpSize), SingleOrDouble);
                    if any(isnan(Vis_fam))
                        Vis_fam = cast(pinv(double(Vs_fam)), SingleOrDouble);
                    end
                    Ws_fam{fi} = Vis_fam;
                    % if OLSflag
                    %     Ws_fam{fi} = eye(tmpSize, SingleOrDouble);
                    % else
                    %     Ws_fam{fi} = Vis_fam;
                    % end
                end
            end

            % Older solution:
            % if ~isempty(ivec_bin)
            %     if GroupByFamType % Compute Vs and Vis by family type
            %         for fi = 1:nfamtypes
            %             ivec = find(famtypevec == fi);
            %             Vs_famtype{fi} = 0;
            %             for ri = 1:length(RandomEffects)
            %                 Vs_famtype{fi} = Vs_famtype{fi} + sig2vec(ri) * clusterinfo{ivec(1)}.(['V_', RandomEffects{ri}]);
            %             end
            %             Vis_famtype{fi} = cast(pinv(double(Vs_famtype{fi})), SingleOrDouble);
            %             if OLSflag
            %                 Ws_famtype{fi} = eye(size(Vis_famtype{fi}), class(Vis_famtype{fi}));
            %             else
            %                 Ws_famtype{fi} = Vis_famtype{fi};
            %             end
            %         end
            %     else % Compute Vs and Vis for each family
            %         for fi = 1:nfam
            %             Vs_fam{fi} = 0;
            %             for ri = 1:length(RandomEffects)
            %                 Vs_fam{fi} = Vs_fam{fi} + sig2vec(ri) * clusterinfo{fi}.(['V_', RandomEffects{ri}]);
            %             end
            %             Vis_fam{fi} = cast(pinv(double(Vs_fam{fi})), SingleOrDouble);
            %             if OLSflag
            %                 Ws_fam{fi} = eye(size(Vis_fam{fi}), class(Vis_fam{fi}));
            %             else
            %                 Ws_fam{fi} = Vis_fam{fi};
            %             end
            %         end
            %     end

            % Compute XtW
            XtW   = zeros(fliplr(size(X)), class(X));
            nClus = length(clusterinfo);

            if GroupByFamType
                for fi = 1:nClus
                    currClus = clusterinfo{fi};
                    XtW(:, currClus.jvec_fam) = X(currClus.jvec_fam,:)' * Ws_famtype{famtypevec(fi)};
                end
            else
                for fi = 1:nClus
                    currClus = clusterinfo{fi};
                    XtW(:, currClus.jvec_fam) = X(currClus.jvec_fam,:)' * Ws_fam{fi};
                end
            end

            % Older solution:
            % Previously, we were calculating XtWVsWt term which was basically:
            % XtW * Vs * Ws
            % Since Ws = inv(Vs), therefore, Vs * Ws = I
            % XtWVsWt = zeros(fliplr(size(X)), class(X));
            % XtW = zeros(fliplr(size(X)), class(X));
            % if GroupByFamType
            %     for fi = 1:length(clusterinfo)
            %         XtW(:,clusterinfo{fi}.jvec_fam) = X(clusterinfo{fi}.jvec_fam,:)' * Ws_famtype{famtypevec(fi)};
            %         XtWVsWt(:,clusterinfo{fi}.jvec_fam) = XtW(:,clusterinfo{fi}.jvec_fam)*Vs_famtype{famtypevec(fi)}*Ws_famtype{famtypevec(fi)}';
            %     end
            % else
            %     for fi = 1:length(clusterinfo)
            %         XtW(:,clusterinfo{fi}.jvec_fam) = X(clusterinfo{fi}.jvec_fam,:)' * Ws_fam{fi};
            %         XtWVsWt(:,clusterinfo{fi}.jvec_fam) = XtW(:,clusterinfo{fi}.jvec_fam)*Vs_fam{fi}*Ws_fam{fi}';
            %     end
            % end

            % XtWVsWtX = XtWVsWt*X;
            % B        = XtWVsWtX;
            % Bi       = pinv(B); % Note that Cov_beta = Bi = pinv(XtWVsWtX)*XtWVsWt*pinv(XtWVsWtX)

            % Compute XtWX
            B  = XtW * X;

            % Calculate inverse of B
            if lowRank
                if useLSQ
                    Bi = lsqminnorm(B, eye(size(B)));
                else
                    Bi = pinv(B);
                end
            else
                Bi = B \ eye(size(B));
            end

            % Calculate beta coefficient
            beta_hat_tmp          = Bi * (XtW * ymat(:, ivec_bin)); % Should generalize this to work with arbitrary W matrices
            Cov_beta              = nearestSPD(Bi);
            beta_hat(:, ivec_bin) = beta_hat_tmp; % These are some of the most time consuming operations -- use MEX?
            beta_se(:,  ivec_bin) = sqrt(diag(Cov_beta) * sig2tvec(ivec_bin));
        end

        % Evaluate contrasts
        for ci = 1:size(contrasts,1)
            betacon_hat(ci, ivec_bin) = contrasts(ci,:) * beta_hat(:,ivec_bin);
            betacon_se(ci,  ivec_bin) = sqrt(contrasts(ci,:) * Cov_beta * contrasts(ci,:)' * sig2tvec(ivec_bin));
        end
        tvec_bins(sig2bini) = (now-t0) * 3600 * 24; % Time in seconds
    end
end
return

% fprintf(1,'Elapsed time %0.2f seconds\n',sum(tvec_bins));

defvec = isfinite(nvec_bins+tvec_bins); %#ok<UNRCH>
P = polyfit(nvec_bins(defvec),tvec_bins(defvec),1);
tvec_bins_pred = polyval(P,nvec_bins);
% Split bins into sets with similar total projected execution time: nansum(tvec_bins_pred)/nsplit each

figure; plot(nvec_bins,tvec_bins,'*')

% ToDo
%   Fix bug w. all zeros in beta_hat and beta_se, when GroupByFamType  == false