function residualGenovec = FEMA_residualizeGenotype(genovec,     fixedEffects,   binvec,         ...
                                                    OLSflag,     clusterinfo,    GroupByFamType, ...
                                                    famtypevec,  Ws_famtype,     Ws_fam,         ...
                                                    useShortcut, allJVec,        allWsTerms)
% Residualize a genotype for the fixed effects using OLS or GLS solution
%% Inputs:
% genovec:          [n x m]     matrix of n subjects and m SNPs
% 
% fixedEffects:     [n x p]     matrix of n subjects and p covariates
% 
% binvec:           [1 x b]     vector of bin values on which random
%                               effects are evaluated
% 
% OLSflag:          logical     OLS or GLS solution (default - OLS)
%
% The following inputs are only necessary of GLS solution is desired:
%
% clusterinfo:      [1 x f]     cell type of f families/clusters, returned
%                               as an output from FEMA_parse_family
%
% GroupByFamType:   logical     true/false indicating if GroupByFamType
% 
% famtypevec:       [1 x f]     vector indicating family type
%
% Ws_famtype and Ws_fam:        inverse of the V term returned by 
%                               FEMA_compileTerms (for every bin)
%
% useShortcut:      logical     returned from FEMA_compileTerms; only
%                               relevant for GLS case; the next two
%                               parameters are only used if useShortcut is
%                               true (in this case, the above optional
%                               arguments can be skipped/empty)
%
% allJVec:          cell        cell type containing all locations for
%                               every bin, returned from FEMA_compileTerms
%
% allWsTerms:       cell        inverse of the V term returned by 
%                               FEMA_compileTerms for every bin
%  
%% Output(s):
% residualGenovec   [b x 1]     cell type with each cell entry containing 
%                               [n x m] matrix of n subjects and m SNPs 
%                               with each SNP residualized for the effect 
%                               of the fixedEffects using OLS or GLS
%
%% Notes:
% When using GLS solution, residualization is performed for every bin of 
% phenotype; therefore, every bin of phenotype gets a different 
% residualized genotype matrix. This is because every bin of the phenotype 
% has its own Vs term (and therefore, its own Vis term as well); since the 
% GLS implementation explicitly uses the Vs term, this should be calculated
% for every bin separately. 
%
% In case of an OLS solution, since Vs and Vis are not needed,
% residualization can be performed overall as a single computation
% across all bins
%
% In this residualization step, the "X" variables are the fixed effects and
% the "y" variable is the genotype matrix which needs to be residualized
% for the fixed effects at each bin

%% Determine if lsqminnorm can be used
if exist('lsqminnorm', 'file')
    useLSQ = true;
else
    useLSQ = false;
end

%% Initialize
numBins         = length(unique(binvec(isfinite(binvec)), 'stable'));
residualGenovec = cell(numBins, 1);
for bins        = 1:numBins
    residualGenovec{bins,1} = zeros(size(genovec));
end

%% Handle case of OLS solution
if OLSflag
    XtW = fixedEffects' * fixedEffects;
    if useLSQ
        tmpInv = lsqminnorm(XtW, eye(size(XtW)));
    else
        if rank(XtW) < size(XtW, 2)
            tmpInv = pinv(XtW);
        else
            tmpInv = XtW \ eye(size(XtW));
        end
    end

    beta_hat_tmp = tmpInv * fixedEffects' * genovec;
    computedResd = genovec - (fixedEffects * beta_hat_tmp);
    for bins     = 1:numBins
        residualGenovec{bins,1} = computedResd;
    end
    return;
end

%% GLS solution
binLoc   = 1;
for bins = 1:numBins
    % for sig2bini    = unique(binvec(isfinite(binvec)), 'stable')
    if ~useShortcut
        XtW         = zeros(fliplr(size(fixedEffects)), class(fixedEffects));
        % XtWVsWt     = zeros(fliplr(size(fixedEffects)), class(fixedEffects));

        if GroupByFamType
            for fi = 1:length(clusterinfo)
                XtW(:, clusterinfo{fi}.jvec_fam)    = fixedEffects(clusterinfo{fi}.jvec_fam,:)' * Ws_famtype{binLoc, famtypevec(fi)};
                % XtWVsWt(:,clusterinfo{fi}.jvec_fam) = XtW(:,clusterinfo{fi}.jvec_fam)           * Vs_famtype{binLoc, famtypevec(fi)} * Ws_famtype{binLoc, famtypevec(fi)}';
            end
        else
            for fi = 1:length(clusterinfo)
                XtW(:,clusterinfo{fi}.jvec_fam)     = fixedEffects(clusterinfo{fi}.jvec_fam,:)' * Ws_fam{binLoc, fi};
                % XtWVsWt(:,clusterinfo{fi}.jvec_fam) = XtW(:,clusterinfo{fi}.jvec_fam)           * Vs_fam{binLoc, fi} * Ws_fam{binLoc, fi}';
            end
        end
    else
        XtW = fixedEffects(allJVec{binLoc},:)' * allWsTerms{binLoc};
    end

    % Calculate beta values for all fixedEffects
    XtWVsWtX     = XtW * fixedEffects;
    if useLSQ
        Bi       = lsqminnorm(XtWVsWtX, eye(size(XtWVsWtX)));
    else
        if rank(XtWVsWtX) < size(XtWVsWtX, 2)
            Bi   = pinv(XtWVsWtX);
        else
            Bi   = XtWVsWtX \ eye(size(XtWVsWtX));
        end
    end
    beta_hat_tmp = Bi * (XtW * genovec);
    
    % Save residuals for this bin
    residualGenovec{binLoc,1} = genovec - (fixedEffects * beta_hat_tmp);
    binLoc                    = binLoc + 1;
end