function residualGenovec = FEMA_residualizeGenotype(genovec,    fixedEffects,       clusterinfo,    ...
                                                    binvec,     GroupByFamType,     famtypevec,     ...
                                                    OLSflag,    Ws_famtype,         Ws_fam)
% Residualize a genotype for the fixed effects using OLS or GLS solution
%% Inputs:
% genovec:          [n x m]     matrix of n subjects and m SNPs
% 
% fixedEffects:     [n x p]     matrix of n subjects and p covariates
% 
% clusterinfo:      [1 x f]     cell type of f families/clusters, returned
%                               as an output from FEMA_parse_family
% 
% binvec:           [1 x b]     vector of bin values on which random
%                               effects are evaluated
% 
% GroupByFamType:   logical     true/false indicating if GroupByFamType
% 
% famtypevec:       [1 x f]     vector indicating family type
% 
% OLSflag:          logical     OLS or GLS solution
% 
% Ws_famtype and Ws_fam:        inverse of the V term returned by 
%                               FEMA_compileTerms (for every bin)
% 
%% Output:
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

%% Initialize
numBins         = length(unique(binvec(isfinite(binvec)), 'stable'));
residualGenovec = cell(numBins, 1);
for bins        = 1:numBins
    residualGenovec{bins} = zeros(size(genovec));
end

%% Handle case of OLS solution
if OLSflag
    beta_hat_tmp = pinv(fixedEffects) * genovec;
    computedResd = genovec - (fixedEffects * beta_hat_tmp);
    for bins     = 1:numBins
        residualGenovec{bins} = computedResd;
    end
    return;
end

%% GLS solution
binLoc          = 1;
for sig2bini    = unique(binvec(isfinite(binvec)), 'stable')
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
    
    % Calculate beta values for all fixedEffects
    XtWVsWtX     = XtW * fixedEffects; 
    Bi           = XtWVsWtX \ eye(size(XtWVsWtX)); % pinv(XtWVsWtX);
    beta_hat_tmp = Bi * (XtW * genovec);
    
    % if OLSflag
    %   beta_hat_tmp = pinv(fixedEffects) * genovec;
    % else
    %   beta_hat_tmp = Bi * XtW * genovec;
    % end
    
    % Save residuals for this bin
    residualGenovec{binLoc,1} = genovec - (fixedEffects * beta_hat_tmp);
    binLoc                    = binLoc + 1;
end