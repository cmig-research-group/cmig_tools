function GLSCache = FEMA_GLS_make_cache(clusterinfo, RandomEffects, nfamtypes, famtypevec, GroupByFamType)
% Precompute the iteration-invariant sparse-assembly coordinates used by FEMA_GLS

    clusterinfo_mat = cell2mat(clusterinfo);
    nnz_max     = sum(cellfun(@length, {clusterinfo_mat.jvec_fam}).^2);
    nfam        = length(clusterinfo);

    % Get ordering of fields in clusterinfo
    ff           = fieldnames(clusterinfo{1});
    RFX_ord      = zeros(length(RandomEffects),1);
    locJVec      = strcmpi(ff, 'jvec_fam');
    for rfx = 1:length(RandomEffects)
        RFX_ord(rfx,1) = find(strcmpi(ff, ['V_', RandomEffects{rfx}]));
    end

    % Sparse coordinates and family metadata
    count = 1;
    allR  = zeros(nnz_max, 1);
    allC  = zeros(nnz_max, 1);

    ivec_cache = [];
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

    GLSCache = struct('nnz_max', nnz_max, 'nfam', nfam,               ...
                      'RFX_ord', RFX_ord, 'locJVec', locJVec,         ...
                      'allR', allR,       'allC', allC,               ...
                      'ivec_cache',     {ivec_cache},                 ...
                      'currClus_cache', {currClus_cache},             ...
                      'obs_idx_cache',  {obs_idx_cache},              ...
                      'block_start',    block_start);
end
