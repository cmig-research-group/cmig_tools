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
p.addParameter('sig2tvec', []);
parse(p, varargin{:});

ciflag         = p.Results.ciflag;
GroupByFamType = p.Results.GroupByFamType;
allWsTerms     = p.Results.allWsTerms;
precision      = p.Results.precision;
useLSQ         = p.Results.useLSQ;
GLSflag        = p.Results.GLSflag;
binvec         = p.Results.binvec;
sig2tvec       = p.Results.sig2tvec;

nobs           = size(X,1);
num_y          = size(y,2);

beta_hat   = zeros(size(X,2), num_y, class(y));
beta_cov   = zeros(size(X,2), size(X,2), num_y, class(y));
tvec_bins  = zeros(max(binvec), 1); % computation time for every bin


%% if allWsTerms is not available, construct the sparsematrix
%% Initialize Precision Matrix Construction
if ~any(strcmp('allWsTerms', varargin))

    nnz_max     = sum(cellfun(@length, {cell2mat(clusterinfo).jvec_fam}).^2);
    nfam        = length(clusterinfo);
    
    % Get ordering of fields in clusterinfo - reasonable to assume that fields
    % are always ordered in the same way since clusterinfo is created in the
    % same way across all clusters
    ff           = fieldnames(clusterinfo{1});
    RFX_ord      = zeros(length(RandomEffects),1);
    locJVec      = strcmpi(ff, 'jvec_fam');
    for rfx = 1:length(RandomEffects)
        RFX_ord(rfx,1) = find(strcmpi(ff, ['V_', RandomEffects{rfx}]));
    end

    % every outcome with each own precision matrix
    allWsTerms = cell(1, num_y);

end

%% Construct Precision Matrices for the Current Bin
for sig2bini = unique(binvec, 'stable')
    t0       = tic;
    ivec_bin = find(binvec == sig2bini);

    if ~any(strcmp('allWsTerms', varargin))
        sig2vec = mean(sig2mat(:,ivec_bin), 2);
        W_1_relative = W_1(:,ivec_bin) .* sig2mat(end,ivec_bin);
        W_1_bin = mean(W_1_relative, 2);

        % Compile terms
        if GroupByFamType
    
            count       = 1;
            allR        = zeros(nnz_max, 1);
            allC        = zeros(nnz_max, 1);
            allV        = zeros(nnz_max, length(ivec_bin));
            % Compute Vs and Vis by family type
            for fi = 1:nfamtypes
                ivec       = find(famtypevec == fi);
                currClus   = struct2cell(clusterinfo{ivec(1)});
                tmpSize    = length(currClus{locJVec});
                Vs_famtype = zeros(tmpSize);
    
                % Compute V
                for ri = 1:(length(RandomEffects)-1)
                    Vs_famtype = Vs_famtype + sig2vec(ri) * currClus{RFX_ord(ri)};
                end
    
                Vs_famtype = Vs_famtype + currClus{RFX_ord(length(RandomEffects))} .* ...
                  diag(W_1_bin(currClus{locJVec})); %jvec_fam
    
                % Compute inverse of V    
                Vs_double = double(Vs_famtype);
                if rcond(Vs_double) < eps(class(Vs_double))
                    Vis_famtype = cast(pinv(Vs_double), precision);
                else
                    Vis_famtype = Vs_double \ eye(tmpSize, precision);
                end
    
                % Compile allWsFam
                tmpClusterinfo              = cell2mat(clusterinfo(ivec(:)));
                tmpR                        = repmat(vertcat(tmpClusterinfo(:).jvec_fam)', tmpSize, 1);
                tmpC                        = repmat(horzcat(tmpClusterinfo(:).jvec_fam),  tmpSize, 1);
                tmp                         = numel(tmpR);
                allR(count:count+tmp-1)     = tmpR(:);
                allC(count:count+tmp-1)     = tmpC(:);

                for fi2 = 1:length(ivec)
                    currIDX = tmpClusterinfo(fi2).jvec_fam;
                    ivec_c  = count + (fi2-1)*tmpSize^2 : count + fi2*tmpSize^2 - 1;
                    allV(ivec_c,:) = restore_W_diag(Vis_famtype, Vs_famtype, ...
                                                    W_1_bin(currClus{locJVec}), W_1_relative(currIDX,:));
                end
                count                       = count + tmp;
            end
    
        else 
            % Compute Vs and Vis for each family
            count    = 1;
            allR        = zeros(nnz_max, 1);
            allC        = zeros(nnz_max, 1);
            allV        = zeros(nnz_max, length(ivec_bin));
            for fi = 1:nfam
                currClus   = struct2cell(clusterinfo{fi});
                tmpSize    = length(currClus{locJVec});
                Vs_fam     = zeros(tmpSize);
    
                % Compute V
                for ri = 1:(length(RandomEffects)-1)
                    Vs_fam = Vs_fam + sig2vec(ri) * currClus{RFX_ord(ri)};
                end
    
                Vs_fam = Vs_fam + currClus{RFX_ord(length(RandomEffects))} .* ...
                 diag(W_1_bin(currClus{locJVec})); 
    
                % Compute inverse of V    
                Vs_double = double(Vs_fam);
                if rcond(Vs_double) < eps(class(Vs_double))
                    Vis_fam = cast(pinv(Vs_double), precision);
                else
                    Vis_fam = Vs_double \ eye(tmpSize, precision);
                end
    
                % Compile allWsFam
                currIDX                 = currClus{locJVec};
                tmpR                    = repmat(currIDX', tmpSize, 1);
                tmpC                    = repmat(currIDX,  tmpSize, 1);
                tmp                     = numel(tmpR);
                allR(count:count+tmp-1) = tmpR(:);
                allC(count:count+tmp-1) = tmpC(:);
                allV(count:count+tmp-1,:) = restore_W_diag(Vis_fam, Vs_fam, ...
                    W_1_bin(currClus{locJVec}), W_1_relative(currIDX,:));
                count                   = count + tmp;
            end
    
        end

        % Put together as sparse matrics
        for jj = 1:length(ivec_bin)
            allWsTerms{ivec_bin(jj)} = sparse(allR, allC, allV(:,jj), nobs, nobs, nnz_max);
        end
    
    end


    %% perfrom GLS estimation for new beta
    if GLSflag
    
        for jj = 1:length(ivec_bin)
    
            yy = ivec_bin(jj);
    
            % Compute XtV
            XtW = double(X)' * allWsTerms{yy};
        
            % Compute XtVX
            B  = XtW * X;
        
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
            beta_cov(:,:,yy) = cast(nearestSPD(Bi), class(y));
    
        end
    
    end

    tvec_bins(sig2bini) = toc(t0);
end

end