function [Ws_famtype, Ws_fam, allWsTerms, allJVec, useShortcut] =       ...
          FEMA_compileTerms(clusterinfo, binvec, nfamtypes, famtypevec, ...
                            sig2mat, RandomEffects,  GroupByFamType,    ...
                            numObs, SingleOrDouble,  OLSflag)
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
% OLSflag:          logical     OLS or GLS (default) solution
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

%% Assign default
% By default, use GLS solution
if ~exist('OLSflag', 'var') || isempty(OLSflag)
    OLSflag = false;
end

%% Determine if shortcut can be used
try
    allWsFam   = zeros(numObs, numObs);
    useShortcut = true;
catch
    try
        allWsFam    = spalloc(numObs, numObs, numObs);
        useShortcut = true;
    catch
        useShortcut = false;
    end
end

%% Initialize
nfam         = length(clusterinfo);
nbins        = length(unique(binvec(isfinite(binvec)), 'stable'));
Vs_famtype   = cell(nbins, nfamtypes); 
Vis_famtype  = cell(nbins, nfamtypes); 
Ws_famtype   = cell(nbins, nfamtypes);
Vs_fam       = cell(nbins, nfam); 
Vis_fam      = cell(nbins, nfam); 
Ws_fam       = cell(nbins, nfam);
binLoc       = 1;

if useShortcut
    allWsTerms   = cell(nbins, 1);
    allJVec      = cell(nbins, 1);
else
    allWsTerms   = [];
    allJVec      = [];
end

for sig2bini = unique(binvec(isfinite(binvec)),'stable')
    ivec_bin = find(binvec==sig2bini);
    sig2vec  = mean(sig2mat(:,ivec_bin), 2);

    if useShortcut
        count    = 1;
        if GroupByFamType
            JVec     = [];
            allWsFam = [];
        else
            JVec     = zeros(numObs, 1);
            allWsFam = zeros(numObs, numObs);
        end
    end

    if ~isempty(ivec_bin)
        if GroupByFamType % Compute Vs and Vis by family type
            for fi = 1:nfamtypes
                ivec = find(famtypevec==fi);
                Vs_famtype{binLoc, fi} = 0;
                for ri = 1:length(RandomEffects)
                    Vs_famtype{binLoc, fi} = Vs_famtype{binLoc, fi} + sig2vec(ri)*clusterinfo{ivec(1)}.(['V_', RandomEffects{ri}]);
                end
                Vis_famtype{binLoc, fi} = cast(double(Vs_famtype{binLoc, fi}) \ eye(size(Vs_famtype{binLoc, fi})), SingleOrDouble);
                % Vis_famtype{binLoc, fi} = cast(pinv(double(Vs_famtype{binLoc, fi})),SingleOrDouble);
                if OLSflag
                    Ws_famtype{binLoc, fi} = eye(size(Vis_famtype{binLoc, fi}), class(Vis_famtype{binLoc, fi}));
                else
                    Ws_famtype{binLoc, fi} = Vis_famtype{binLoc, fi};
                end

                if useShortcut
                    currIDX                     = clusterinfo{ivec(1)}.jvec_fam;
                    tmpLen                      = length(currIDX);
                    JVec(count:count+tmpLen-1)  = currIDX;
                    allWsFam(currIDX, currIDX)  = Ws_famtype{binLoc, fi};
                    count                       = count + tmpLen;
                end
            end
        else % Compute Vs and Vis for each family
            for fi = 1:nfam
                Vs_fam{binLoc, fi} = 0;
                for ri = 1:length(RandomEffects)
                    Vs_fam{binLoc, fi} = Vs_fam{binLoc, fi} + sig2vec(ri)*clusterinfo{fi}.(['V_', RandomEffects{ri}]);
                end
                Vis_fam{binLoc, fi} = cast(double(Vs_fam{binLoc, fi}) \ eye(size(Vs_fam{binLoc, fi})), SingleOrDouble);
                % Vis_fam{binLoc, fi} = cast(pinv(double(Vs_fam{binLoc, fi})),SingleOrDouble);
                if OLSflag
                    Ws_fam{binLoc, fi} = eye(size(Vis_fam{binLoc, fi}), class(Vis_fam{binLoc, fi}));
                else
                    Ws_fam{binLoc, fi} = Vis_fam{binLoc, fi};
                end

                if useShortcut
                    currIDX                     = clusterinfo{fi}.jvec_fam;
                    tmpLen                      = length(currIDX);
                    JVec(count:count+tmpLen-1)  = currIDX;
                    allWsFam(currIDX, currIDX)  = Ws_fam{binLoc, fi};
                    count                       = count + tmpLen;
                end
            end
        end
    end

    % Additional terms
    if useShortcut
        allWsTerms{binLoc, 1} = sparse(allWsFam);
        allJVec{binLoc,    1} = JVec;
    end

    % Update bin counter
    binLoc = binLoc + 1;
end