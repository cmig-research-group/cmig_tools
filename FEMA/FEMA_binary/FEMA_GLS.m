function [allWsTerms, beta_hat, beta_cov, tvec_bins] = FEMA_GLS(y, X, W_1, sig2mat, RandomEffects, ...
                                                     clusterinfo, nfamtypes, famtypevec, varargin)

p = inputParser;
p.addParameter('ciflag', false);
p.addParameter('GroupByFamType', false);
p.addParameter('allWsTerms', []);
p.addParameter('precision', 'double');
p.addParameter('useLSQ', true);
p.addParameter('GLSflag', true);
p.addParameter('binvec', []);
p.addParameter('max_pages', 100000);
p.addParameter('GLSCache', []);
p.addParameter('nbins', 0);
parse(p, varargin{:});

ciflag         = p.Results.ciflag;
GroupByFamType = p.Results.GroupByFamType;
allWsTerms     = p.Results.allWsTerms;
precision      = p.Results.precision;
useLSQ         = p.Results.useLSQ;
GLSflag        = p.Results.GLSflag;
binvec         = p.Results.binvec;
max_pages      = p.Results.max_pages;
GLSCache       = p.Results.GLSCache;
nbins          = p.Results.nbins;

nobs           = size(X,1);
num_y          = size(y,2);
num_RFX        = numel(RandomEffects);

beta_hat   = zeros(size(X,2), num_y, class(y));
beta_cov   = zeros(size(X,2), size(X,2), num_y, class(y));
if nbins == 0
    tvec_bins = zeros(max(binvec), 1);
else
    tvec_bins = sparse(max(binvec), 1); % only occupied bins need storage
end


%% if allWsTerms is not available, construct the sparsematrix
%% Initialize Precision Matrix Construction
if ~any(strcmp('allWsTerms', varargin))

    allWsTerms = cell(1, num_y);


    if ~isempty(GLSCache)
        nnz_max        = GLSCache.nnz_max;
        nfam           = GLSCache.nfam;
        RFX_ord        = GLSCache.RFX_ord;
        locJVec        = GLSCache.locJVec;
        allR           = GLSCache.allR;
        allC           = GLSCache.allC;
        ivec_cache     = GLSCache.ivec_cache;
        currClus_cache = GLSCache.currClus_cache;
        obs_idx_cache  = GLSCache.obs_idx_cache;
        block_start    = GLSCache.block_start;
    else
        clusterinfo_mat = cell2mat(clusterinfo);
        nnz_max     = sum(cellfun(@length, {clusterinfo_mat.jvec_fam}).^2);
        nfam        = length(clusterinfo);
        
        % Get ordering of fields in clusterinfo - reasonable to assume that fields
        % are always ordered in the same way since clusterinfo is created in the
        % same way across all clusters
        ff           = fieldnames(clusterinfo{1});
        RFX_ord      = zeros(num_RFX,1);
        locJVec      = strcmpi(ff, 'jvec_fam');
        for rfx = 1:num_RFX
            RFX_ord(rfx,1) = find(strcmpi(ff, ['V_', RandomEffects{rfx}]));
        end
    
        % Sparse coordinates and family metadata
        count = 1;
        allR  = zeros(nnz_max, 1);
        allC  = zeros(nnz_max, 1);
    
        if GroupByFamType
            ivec_cache       = cell(nfamtypes, 1);
            currClus_cache   = cell(nfamtypes, 1);
            obs_idx_cache    = cell(nfamtypes, 1);
            block_start      = zeros(nfamtypes, 1);
    
            for fi = 1:nfamtypes
                ivec           = find(famtypevec == fi);
                currClus       = struct2cell(clusterinfo{ivec(1)});
                tmpSize        = length(currClus{locJVec});
                tmpClusterinfo = clusterinfo_mat(ivec(:));
                tmpR           = repmat(vertcat(tmpClusterinfo(:).jvec_fam)', tmpSize, 1);
                tmpC           = repmat(horzcat(tmpClusterinfo(:).jvec_fam),  tmpSize, 1);
                tmp             = numel(tmpR);
                allR(count:count+tmp-1) = tmpR(:);
                allC(count:count+tmp-1) = tmpC(:);
                ivec_cache{fi}       = ivec;
                currClus_cache{fi}   = currClus;
                obs_idx_cache{fi}    = reshape(horzcat(tmpClusterinfo(:).jvec_fam),tmpSize, []);
                block_start(fi)      = count;
                count = count + tmp;
            end
        else
            currClus_cache = cell(nfam, 1);
            obs_idx_cache  = cell(nfam, 1);
            block_start    = zeros(nfam, 1);
            for fi = 1:nfam
                currClus = struct2cell(clusterinfo{fi});
                tmpSize  = length(currClus{locJVec});
                currIDX  = currClus{locJVec};
                tmpR     = repmat(currIDX', tmpSize, 1);
                tmpC     = repmat(currIDX,  tmpSize, 1);
                tmp      = numel(tmpR);
                allR(count:count+tmp-1) = tmpR(:);
                allC(count:count+tmp-1) = tmpC(:);
                currClus_cache{fi} = currClus;
                obs_idx_cache{fi}  = currIDX;
                block_start(fi)    = count;
                count = count + tmp;
            end
        end
    end
end

%% Construct Precision Matrices for the Current Bin
for sig2bini = unique(binvec, 'stable')
    t0       = tic;
    ivec_bin = find(binvec == sig2bini);
    num_bin  = length(ivec_bin);

    if ~any(strcmp('allWsTerms', varargin))
        sig2vec = mean(sig2mat(:,ivec_bin), 2);
        W_1_relative = W_1(:,ivec_bin) .* sig2mat(end,ivec_bin);
        W_1_bin = mean(W_1_relative, 2);
        allV = zeros(nnz_max, 1);

        % Compile terms
        if GroupByFamType

            if nbins == 0

                allV        = zeros(nnz_max, num_bin);
                fam_chunk   = max(1, floor(max_pages / num_bin));
                % Compute Vs and Vis by family type
                for fi = 1:nfamtypes
                    ivec       = ivec_cache{fi};
                    currClus   = currClus_cache{fi};
                    tmpSize    = length(currClus{locJVec});
                    Vs_famtype = zeros(tmpSize);
        
                    % Compute V
                    for ri = 1:(num_RFX-1)
                        Vs_famtype = Vs_famtype + sig2vec(ri) * currClus{RFX_ord(ri)};
                    end
                    
                    W_1_base  = W_1_bin(currClus{locJVec});
                    Vs_famtype = Vs_famtype + currClus{RFX_ord(num_RFX)} .* ...
                      diag(W_1_base); %jvec_fam
        
                    % Compute inverse of V    
                    Vs_double = double(Vs_famtype);
                    if rcond(Vs_double) < eps(class(Vs_double))
                        Vis_famtype = cast(pinv(Vs_double), precision);
                    else
                        Vis_famtype = Vs_double \ eye(tmpSize, precision);
                    end
    
                    % Families are independent blocks. restore W_diag and
                    % compile allWsFam
                    obs_idx_all = obs_idx_cache{fi};
                    for first_fam = 1:fam_chunk:length(ivec)
                        last_fam = min(first_fam + fam_chunk - 1, length(ivec));
                        num_fam  = last_fam - first_fam + 1;
                        obs_idx  = obs_idx_all(:,first_fam:last_fam);
                        W_target = reshape(W_1_relative(obs_idx(:), :), tmpSize, num_fam, num_bin);
                        V_chunk  = restore_W_diag(Vis_famtype, Vs_famtype, W_1_base, W_target, max_pages, num_fam);
                        row_idx  = block_start(fi) + (first_fam-1)*tmpSize^2 : block_start(fi) + last_fam*tmpSize^2 - 1;
                        allV(row_idx,:) = reshape(V_chunk, length(row_idx), num_bin);
                    end
                end
            else
                for fi = 1:nfamtypes
                    ivec       = ivec_cache{fi};
                    currClus   = currClus_cache{fi};
                    obs_idx    = obs_idx_cache{fi};
                    tmpSize    = length(currClus{locJVec});
                    num_famtype = numel(ivec);

                    R_famtype = zeros(tmpSize);
                    for ri = 1:(num_RFX-1)
                        R_famtype = R_famtype + sig2vec(ri) * currClus{RFX_ord(ri)};
                    end
                    error_diag = diag(currClus{RFX_ord(num_RFX)});

                    % The random covariance template is shared by family type
                    fam_chunk = max(1, floor(max_pages));
                    for first_fam = 1:fam_chunk:num_famtype
                        last_fam = min(first_fam + fam_chunk - 1, num_famtype);
                        fam_ids = first_fam:last_fam;
                        num_page = numel(fam_ids);
                        V_page = repmat(double(R_famtype), 1, 1, num_page);
                        W_page = W_1_bin(obs_idx(:, fam_ids));
                        for dd = 1:tmpSize
                            V_page(dd, dd, :) = V_page(dd, dd, :) + ...
                                reshape(double(error_diag(dd) .* W_page(dd, :)), 1, 1, num_page);
                        end
                        Vis_page = pagemldivide(V_page, eye(tmpSize, precision));
                        bad_page = reshape(any(any(~isfinite(Vis_page), 1), 2), 1, []);
                        for pp = find(bad_page)
                            Vis_page(:,:,pp) = pinv(V_page(:,:,pp));
                        end

                        row_idx = block_start(fi) + (first_fam-1)*tmpSize^2 : ...
                                  block_start(fi) + last_fam*tmpSize^2 - 1;
                        allV(row_idx) = reshape(Vis_page, [], 1);
                    end
                end
            end

        else 
            % Compute Vs and Vis for each family
            for fi = 1:nfam
                currClus   = currClus_cache{fi};
                tmpSize    = length(currClus{locJVec});
                Vs_fam     = zeros(tmpSize);
    
                % Compute V
                for ri = 1:(num_RFX-1)
                    Vs_fam = Vs_fam + sig2vec(ri) * currClus{RFX_ord(ri)};
                end
    
                W_1_base = W_1_bin(currClus{locJVec});
                Vs_fam = Vs_fam + currClus{RFX_ord(num_RFX)} .* ...
                 diag(W_1_base); 
    
                % Compute inverse of V    
                Vs_double = double(Vs_fam);
                if rcond(Vs_double) < eps(class(Vs_double))
                    Vis_fam = cast(pinv(Vs_double), precision);
                else
                    Vis_fam = Vs_double \ eye(tmpSize, precision);
                end
    
                % Compile allWsFam
                row_idx = block_start(fi):block_start(fi)+tmpSize^2-1;
                if nbins == 0
                    currIDX = obs_idx_cache{fi};
                    allV(row_idx,:) = restore_W_diag(Vis_fam, Vs_fam, W_1_base, W_1_relative(currIDX,:), max_pages);
                else
                    allV(row_idx) = Vis_fam(:);
                end
            end
    
        end

        % Put together as sparse matrics
        Ws_bin = sparse(allR, allC, allV, nobs, nobs, nnz_max);
        allWsTerms(ivec_bin) = {Ws_bin};
    
    end


    %% perfrom GLS estimation for new beta
    if GLSflag
    
        yy = ivec_bin;

        % Compute XtV
        XtW = double(X)' * allWsTerms{yy(1)};
    
        % Compute XtVX
        B = XtW * X;
    
        % Calculate inverse of XtWX
        if rank(B) < size(B,2)
            if useLSQ
                Bi = lsqminnorm(B, eye(size(B)));
            else
                Bi = pinv(B);
            end
        else
            Bi = B \ eye(size(B));
        end
    
        % Calculate beta: inv(X' * inv(V) * X) * X' * inv(V) * y
        beta_hat(:,yy) = Bi * (XtW * y(:,yy));
        % relative working scale for the coefficient covariance
        cov_bin = cast(nearestSPD(Bi), class(y));
        beta_cov(:,:,yy) = repmat(cov_bin, [1 1 numel(yy)]);
    
    end

    tvec_bins(sig2bini) = toc(t0);
end

end