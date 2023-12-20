function allWsTerms = FEMA_compileTerms(clusterinfo,    binvec,  nfamtypes,     ...
                                        famtypevec,     sig2mat, RandomEffects, ...
                                        GroupByFamType, numObs,  SingleOrDouble)
%% Compile Vs and Ws terms for every bin
%% Inputs:
% clusterinfo:      [1 x f]     cell type of f families/clusters, returned
%                               as an output from FEMA_parse_family
% 
% binvec:           [1 x b]     vector of b bin values
% 
% nfamtypes:        [1 x 1]     used if GroupByFamType
% 
% famtypevec:       [1 x f]     vector of family types
% 
% sig2mat:          [r x b]     estimates of r random effects across b bins
% 
% RandomEffects:    [1 x r]     cell having names of the r random effects
% 
% GroupByFamType:   logical     whether to GroupByFamType or not
% 
% SingleOrDouble:   character   single or double precision
%
%% Outputs:
% Ws_famtype and Ws_fam as cell types containing inverses of the "V" term 
% across clusters and bins. 
% 
% Additionally, allWsTerms, a cell type vector is also returned which 
% contains sparse matrices for Ws_famtype or Ws_fam terms in the order of 
% allJVec. However, this is only the case when a shortcut can be used (see
% Notes below)
%
%% Notes:
% Whether or not a shortcut (faster alternate) can be used in
% FEMA_sig2bin_parfeval_GWAS_bf is determined when compiling terms; a
% simple try-catch run is performed, trying to initialize:
%                   zeros(numObs, numObs); 
% If this initialization is successful, then the additional terms
% allWsTerms and allJVec are compiled and useShortcut is returned as true;
% otherwise these output are set to empty and useShortcut is false
%
% The shortcut itself involves a sparse matrix multiplication as opposed to
% performing cluster-wise multiplication - for cases where a full matrix of
% numObs x numObs can be initialized, tests indicate that sparse matrix
% multiplication can be faster

%% Initialize
nnz_max     = sum(cellfun(@length, {cell2mat(clusterinfo).jvec_fam}).^2);
nfam        = length(clusterinfo);
nbins       = length(unique(binvec(isfinite(binvec)), 'stable'));
binLoc      = 1;
allWsTerms  = cell(nbins, 1);
allR        = zeros(nnz_max, 1);
allC        = zeros(nnz_max, 1);
allV        = zeros(nnz_max, 1);

% Get ordering of fields in clusterinfo - reasonable to assume that fields
% are always ordered in the same way since clusterinfo is created in the
% same way across all clusters
ff           = fieldnames(clusterinfo{1});
RFX_ord      = zeros(length(RandomEffects),1);
locJVec      = strcmpi(ff, 'jvec_fam');
for rfx = 1:length(RandomEffects)
    RFX_ord(rfx,1) = find(strcmpi(ff, ['V_', RandomEffects{rfx}]));
end

%% Compile terms
for sig2bini = unique(binvec(isfinite(binvec)), 'stable')
    ivec_bin = find(binvec == sig2bini);
    sig2vec  = mean(sig2mat(:,ivec_bin), 2);
    count    = 1;

    if ~isempty(ivec_bin)
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
                end

                % Compute inverse of V
                Vis_famtype = double(Vs_famtype) \ eye(tmpSize, SingleOrDouble);
                if any(isnan(Vis_famtype))
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
                % allV(count:count+tmp-1)     = nonzeros(blkdiag(temp{:}));
                count                       = count + tmp;

                % if useShortcut
                %     tmpClusterinfo              = cell2mat(clusterinfo);
                %     currIDX                     = horzcat(tmpClusterinfo(ivec(:)).jvec_fam);
                %     temp                        = repmat(Vis_famtype, 1, length(ivec));
                %     allWsFam(currIDX, currIDX)  = blkdiag(temp{:}); %#ok<SPRIX>
                % end
            end

            % Put together as a sparse matrix
            allWsFam = sparse(allR, allC, allV, numObs, numObs, nnz_max);

        else 
            % Compute Vs and Vis for each family
            for fi = 1:nfam
                currClus   = struct2cell(clusterinfo{fi});
                tmpSize    = length(currClus{locJVec});
                Vs_fam     = zeros(tmpSize);

                % Compute V
                for ri = 1:length(RandomEffects)
                    Vs_fam = Vs_fam + sig2vec(ri) * currClus{RFX_ord(ri)};
                end

                % Compute inverse of V
                Vis_fam = double(Vs_fam) \ eye(tmpSize, SingleOrDouble);
                if any(isnan(Vis_fam))
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

                % if useShortcut
                %     currIDX                    = currClus{locJVec};
                %     allWsFam(currIDX, currIDX) = Vis_fam; %#ok<SPRIX>
                % end
            end
            % Put together as a sparse matrix
            allWsFam = sparse(allR, allC, allV, numObs, numObs, nnz_max);
        end
    end

    % Additional terms
    allWsTerms{binLoc, 1} = allWsFam;
    
    % Update bin counter
    binLoc = binLoc + 1;
end