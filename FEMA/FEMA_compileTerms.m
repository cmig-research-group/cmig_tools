function [allWsTerms, timing] = FEMA_compileTerms(clusterinfo, binvec,         ...
                                                  sig2mat,     RandomEffects,  ...
                                                  famtypevec,  GroupByFamType, ...
                                                  CovType,     SingleOrDouble, visitnum)
%% Compile Ws term for every bin
%% Inputs:
% clusterinfo:      [1 x f]     cell type of f families/clusters, returned
%                               as an output from FEMA_parse_family
% 
% binvec:           [1 x b]     vector of b bin values
% 
% sig2mat:          [r x v]     estimates of r random effects across v
%                               voxels in case of analytic covariance type
%                               OR visit * visit * r random effects across
%                               v voxels in case of unstructured covariance
%
% RandomEffects:    [1 x r]     cell having names of the r random effects
%
% famtypevec:       [1 x f]     vector of family types
% 
% GroupByFamType:   logical     whether to GroupByFamType or not
%
% CovType:          character   should be either:
%                                   * 'analytic' (compound symmetry)
%                                   * 'unstructured'
% 
% SingleOrDouble:   character   single or double precision
%
% visitnum:         [e x 1]     event number as output from FEMA_fit
%                               (required if CovType is unstructured)
%
%% Output(s):
% allWsTerms is a cell type containing invers of the "V" term for each bin;
% additionally, timing for compiling terms is output as a structure

%% Start overall timer
tOver = tic;

%% Check if unstructured covariance
if ~exist('CovType', 'var') || isempty(CovType)
    unstructuredCov = false;
else
    if strcmpi(CovType, 'unstructured')
        unstructuredCov = true;
    else
        unstructuredCov = false;
    end
end

% Ensure visitnum is specified for unstructuredcov
if unstructuredCov
    if ~exist('visitnum', 'var') || isempty(visitnum)
        error('visitnum is required if CovType is unstructured');
    end
end

%% Initialize
tInit       = tic;
tempVar     = cell2mat(clusterinfo);
nfamtypes   = length(unique([tempVar.famtype]));
numObs      = length([tempVar.jvec_fam]);
nnz_max     = length(vertcat(tempVar.ivec_fam));
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

% Save timing information for initialization
timing.Initialize = toc(tInit);

%% Compile terms
tCompile = tic;

% Get warning statuses for singular and nearly singular cases;
% temporarily set their display off
statusSingular = warning('off', 'MATLAB:singularMatrix');
statusNearSing = warning('off', 'MATLAB:nearlySingularMatrix');

% Clear last warning
lastwarn('');

% Handle unstructured covariance
if unstructuredCov
    for ivec_bin = 1:length(binvec)
        count = 1;
        for fi = 1:nfam
            currClus   = struct2cell(clusterinfo{fi});
            tmpSize    = length(currClus{locJVec});

            % Compute V
            Vs_fam = zeros(tmpSize);
            b = visitnum(clusterinfo{fi}.jvec_fam);
            for ri = 1:length(RandomEffects)-1
                Vs_fam = Vs_fam + sig2mat(b, b, ri, ivec_bin) .* currClus{RFX_ord(ri)};
            end

            % Compute inverse of V
            Vis_fam = double(Vs_fam) \ eye(tmpSize, SingleOrDouble);
            msg     = lastwarn;
            if ~isempty(msg)
                Vis_fam = cast(pinv(double(Vs_fam)), SingleOrDouble);
                msg = '';
                lastwarn('');
            end
            % if any(isnan(Vis_fam) | isinf(Vis_fam))
            %     Vis_fam = cast(pinv(double(Vs_fam)), SingleOrDouble);
            % end

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
        allWsFam = sparse(allR, allC, allV, numObs, numObs, nnz_max);

        % Additional terms
        allWsTerms{binLoc, 1} = allWsFam;

        % Update bin counter
        binLoc = binLoc + 1;
    end
else
    % Regular compound symmetry
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
                    msg         = lastwarn;
                    if ~isempty(msg)
                        Vis_famtype = cast(pinv(double(Vs_famtype)), SingleOrDouble);
                        msg = '';
                        lastwarn('');
                    end
                    % if any(isnan(Vis_famtype) | isinf(Vis_famtype))
                    %     Vis_famtype = cast(pinv(double(Vs_famtype)), SingleOrDouble);
                    % end
    
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
                    msg     = lastwarn;
                    if ~isempty(msg)
                        Vis_fam = cast(pinv(double(Vs_fam)), SingleOrDouble);
                        msg = '';
                        lastwarn('');
                    end
                    % if any(isnan(Vis_fam) | isinf(Vis_fam))
                    %     Vis_fam = cast(pinv(double(Vs_fam)), SingleOrDouble);
                    % end
    
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
                allWsFam = sparse(allR, allC, allV, numObs, numObs, nnz_max);
            end
        end
    
        % Additional terms
        allWsTerms{binLoc, 1} = allWsFam;
        
        % Update bin counter
        binLoc = binLoc + 1;
    end
end

% Reset the status of warnings
warning(statusSingular);
warning(statusNearSing);

% Save timing information for compiling terms
timing.Compiling = toc(tCompile);

% Save overall timing
timing.Overall = toc(tOver);