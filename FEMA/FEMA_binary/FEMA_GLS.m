function [allWsTerms, beta_hat, beta_se, beta_cov] = FEMA_GLS(y, X, W_1, sig2mat, RandomEffects, ...
                                                        clusterinfo, nfamtypes, famtypevec, varargin)

p = inputParser;
p.addParameter('ciflag', false);
p.addParameter('GroupByFamType', false);
p.addParameter('allWsTerms', cell(1,size(y,2)));
p.addParameter('SingleOrDouble', 'double');
p.addParameter('useLSQ', true);
p.addParameter('GLSflag', true);
parse(p, varargin{:});

ciflag         = p.Results.ciflag;
GroupByFamType = p.Results.GroupByFamType;
allWsTerms     = p.Results.allWsTerms;
SingleOrDouble = p.Results.SingleOrDouble;
useLSQ         = p.Results.useLSQ;
GLSflag        = p.Results.GLSflag;

nobs           = size(X,1);

[beta_hat,    beta_se]    = deal(zeros(size(X,2), size(y,2), class(y)));
beta_cov                  = cell(1,size(y,2));


%% if allWsTerms is not available, construct the sparsematrix
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

    
    % Compile terms
    if GroupByFamType
        for coli = 1:size(y,2)
            count    = 1;
            allR        = zeros(nnz_max, 1);
            allC        = zeros(nnz_max, 1);
            allV        = zeros(nnz_max, 1);
        % Compute Vs and Vis by family type
            for fi = 1:nfamtypes
                ivec       = find(famtypevec == fi);
                currClus   = struct2cell(clusterinfo{ivec(1)});
                tmpSize    = length(currClus{locJVec});
                Vs_famtype = zeros(tmpSize);
    
                % Compute V
                for ri = 1:(length(RandomEffects)-1)
                    Vs_famtype = Vs_famtype + sig2mat(ri, coli) * currClus{RFX_ord(ri)};
                end
    
                Vs_famtype = Vs_famtype + currClus{RFX_ord(length(RandomEffects))} .* ...
                  diag(W_1(currClus{locJVec}, coli)); %jvec_fam
    
                % Compute inverse of V
                Vis_famtype = double(Vs_famtype) \ eye(tmpSize, SingleOrDouble);
                if any(isnan(Vis_famtype) | isinf(Vis_famtype))
                    Vis_famtype = cast(pinv(double(Vs_famtype)), SingleOrDouble);
                end
    
                % Compile allWsFam
                tmpClusterinfo              = cell2mat(clusterinfo(ivec(:)));
                tmpR                        = repmat(vertcat(tmpClusterinfo(:).jvec_fam)', tmpSize, 1);
                tmpC                        = repmat(horzcat(tmpClusterinfo(:).jvec_fam),  tmpSize, 1);
                temp                        = repmat(Vis_famtype(:), length(ivec), 1);
                tmp                         = numel(tmpR);
                allR(count:count+tmp-1)     = tmpR(:);
                allC(count:count+tmp-1)     = tmpC(:);
                allV(count:count+tmp-1)     = temp;
                count                       = count + tmp;
            end

            % Put together as a sparse matrix
            allWsTerms{coli} = sparse(allR, allC, allV, nobs, nobs, nnz_max);

        end

    else 
        % Compute Vs and Vis for each family
        for coli = 1:size(y, 2)
            count    = 1;
            allR        = zeros(nnz_max, 1);
            allC        = zeros(nnz_max, 1);
            allV        = zeros(nnz_max, 1);
            for fi = 1:nfam
                currClus   = struct2cell(clusterinfo{fi});
                tmpSize    = length(currClus{locJVec});
                Vs_fam     = zeros(tmpSize);
    
                % Compute V
                for ri = 1:(length(RandomEffects)-1)
                    Vs_fam = Vs_fam + sig2mat(ri, coli) * currClus{RFX_ord(ri)};
                end
    
                Vs_fam = Vs_fam + currClus{RFX_ord(length(RandomEffects))} .* ...
                 diag(W_1(currClus{locJVec}, coli)); 
    
                % Compute inverse of V
                Vis_fam = double(Vs_fam) \ eye(tmpSize, SingleOrDouble);
                if any(isnan(Vis_fam) | isinf(Vis_fam))
                    Vis_fam = cast(pinv(double(Vs_fam)), SingleOrDouble);
                end
    
                % Compile allWsFam
                currIDX                 = currClus{locJVec};
                tmpR                    = repmat(currIDX', tmpSize, 1);
                tmpC                    = repmat(currIDX,  tmpSize, 1);
                tmp                     = numel(tmpR);
                allR(count:count+tmp-1) = tmpR(:);
                allC(count:count+tmp-1) = tmpC(:);
                allV(count:count+tmp-1) = Vis_fam(:);
                count                   = count + tmp;
            end

            % Put together as a sparse matrix
            allWsTerms{coli} = sparse(allR, allC, allV, nobs, nobs, nnz_max);

        end

    end

end


%% perfrom GLS estimation for new beta
if GLSflag

    for coli = 1:size(y,2)
        % Compute XtV
        XtW = double(X)' * allWsTerms{coli};
    
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
        beta_hat(:,coli) = Bi * (XtW * y(:,coli));
        beta_cov{coli}   = nearestSPD(Bi);
        beta_se(:,coli)  = sqrt(diag(beta_cov{coli}));
        
    end
end

end