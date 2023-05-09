%% Nomenclature / notation
% Variable      Size        Description
% --------      ----        -----------
% n                         number of observations
% p                         number of X variables (covariates)
% v                         number of y variables (phenotypes)
% m                         number of SNPs (genetic data)
% b                         number of bins
% q                         number of basis functions
%
% iid           [1 x n]     list of IDs corresponding to n observations
% fid           [1 x n]     list of family IDs corresponding to n observations
% IA            [k x 1]     list of k unique IDs: [~, IA, c] = unique(iid, 'stable');
%
% X             [n x p]     matrix of covariates (X variables)
% ymat          [n x v]     matrix of phenotypes (y variables)
% pihatmat      [k x k]     genetic relationship matrix [IA x IA]
% clusterinfo   [1 x f]     cell of f clusters

%% Pre-GWAS definitions
% The following variables need to be defined
% Consider using initial parts of FEMA_run_on_synthetic_data
iid             = [];
fid             = [];
pihatmat        = [];
RandomEffects   = {'F', 'A', 'S', 'E'};

% Assume that there is some PLINK file that contains genotyping information
% for the unique list of subjects in iid
bFile = '';

%% Read genetic data
% Read in genetic data for IID subjects and all SNPs
% It is also possible to read data for all subjects or selected SNPs, 
% perform mean centering or standardization of genetic data
SNPList     = [];
onlyCheck   = false;
lowMem      = false;
meanImpute  = true;
roundOff    = true;
transform   = 'none';
stdType     = 'none';
[genomat, Chr, allSNPs, basePair] = FEMA_parse_PLINK(bFile, iid, SNPList, onlyCheck,  ...
                                                     lowMem, meanImpute,  roundOff,   ...
                                                     transform, stdType);

%% Assign 10% of independent SNPs as causal SNPs
listIndepSNPs = [];
nSNPs         = size(genomat, 2);
nIndepSNPs    = length(listIndepSNPs);
causalFrac    = 0.1;
nCausal       = floor(causalFrac * nIndepSNPs);
causalSNPs    = listIndepSNPs(randperm(nIndepSNPs, nCausal));

%% Get clusterinfo
[clusterinfo, Ss, iid, famtypevec, famtypelist, subj_famtypevec] =          ...
 FEMA_parse_family(iid, zeros(length(iid),1), fid, zeros(length(iid), 1),   ...
                   pihatmat, 'RandomEffects', RandomEffects);

%% Simulate phenotype
numObs        = length(iid);
numPhenotypes = 10;
numXVariables = 5;
effSize_SNPs  = [];  % Leave empty for effect sizes to be drawn from standard normal distribution
effSize_FFX   = [];  % Leave empty for effect sizes to be drawn from standard normal distribution
V_FFX         = 1;   % Proportion of variance explained by fixed effects
V_A           = 0.4; % Proportion of variance explained by additive genetic effect
V_F           = 0.2; % Proportion of variance explained by family effects
V_S           = 0.3; % Proportion of variance explained by subject effects

[phenotype,    effSize_SNPs,      locCausal,    XVars,           effSize_FFX,       ...
 effSize_RFX,  propTotVariance,   propTotNames, propResVar_RFX,  propResNames] =    ...
 FEMA_simulateGWAS(numPhenotypes, genomat,      allSNPs,         causalSNPs,        ...
                   effSize_SNPs,  V_A,          numXVariables,   effSize_FFX,       ...
                   V_FFX,         clusterinfo,  {'F', 'S'},      [V_F; V_S], numObs);

%% Fit model
eid             = ones(length(iid),  1);
agevec          = zeros(length(iid), 1);
niter           = 1;
contrasts       = [];
nbins           = 20;
nperms          = 0;
SingleOrDouble  = 'double';
FixedEstType    = 'gls';

[beta_hat,  beta_se, zmat, logpmat, sig2tvec, sig2mat, binvec, ymat_res_bak] =  ...
 FEMA_fit(XVars, iid, eid, fid, agevec, phenotype, niter, contrasts, nbins,     ...
          pihatmat, 'RandomEffects', RandomEffects, 'nperms', nperms,           ...
          'SingleOrDouble', SingleOrDouble, 'FixedEstType', FixedEstType);

%% Settings for GWAS
GroupByFamType  = false;
nfamtypes       = length(famtypelist);
OLSflag         = false;
df              = size(phenotype) - numXVariables; % only required if pValType is 't'
basisFunction   = ones(numObs, 1);                 % equivalent to standard GWAS
pValType        = 'z';                             % use normcdf for calculating p values
outDir          = pwd;

%% First, compile Ws terms
[Ws_famtype, Ws_fam] = FEMA_compileTerms(clusterinfo,    binvec,         nfamtypes,     ...
                                         famtypevec,     sig2mat,        RandomEffects, ...
                                         GroupByFamType, SingleOrDouble, OLSflag);

%% Demo 1 - run GWAS on entire data serially
residualGenovec = FEMA_residualizeGenotype(genomat, XVars,          clusterinfo,    ...
                                           binvec,  GroupByFamType, famtypevec,     ...
                                           OLSflag, Ws_famtype,     Ws_fam);

[snp_beta_hat, snp_beta_se, snp_tStats, snp_logpValues] =                        ...
 FEMA_sig2binseg_parfeval_GWAS_bf(residualGenovec, ymat_res_bak, clusterinfo,    ...
                                  binvec,          sig2tvec,     GroupByFamType, ...
                                  famtypevec,      Ws_famtype,   Ws_fam,         ...
                                  basisFunction,   pValType,     df,             ...
                                  OLSflag,         SingleOrDouble);

%% Demo 2 - run GWAS in chunks (parallel)
numWorkers = 2;
numThreads = 2;
chunkSize  = 10000;
allChunks  = 1:chunkSize:nSNPs;
numChunks  = length(allChunks);
splitInfo  = cell(numChunks, 1);

% Get some information of the genotyping data
% [~, Chr, allSNPs, ~] = FEMA_parse_PLINK(bFile, iid, [], true);

% Divide data into chunks
for chunk = 1:numChunks
    if chunk == numChunks
        tmpLocs = allChunks(chunk):nSNPs;
    else
        tmpLocs = allChunks(chunk):allChunks(chunk)+chunkSize-1;
    end
    splitInfo{chunk}.tmpLocs = tmpLocs;
    splitInfo{chunk}.SNPs    = allSNPs(tmpLocs);
    splitInfo{chunk}.outDir  = outDir;
    splitInfo{chunk}.outName = ['FEMA_GWAS_part', num2str(chunk, '%05d')];
end

% % Instead of dividing into chunks, it is also possible to split by
% % chromosomes
% % Identify chromosomes
% chr      = unique(Chr, 'stable');
% numChr   = length(chr);
% 
% % Split genotype across chromosomes
% splitInfo  = cell(numChr,1);
% for chromo = 1:numChr
%     tmpLocs                      = strcmpi(Chr, chr{chromo});
%     splitInfo{chromo}.tmpLocs    = tmpLocs;
%     splitInfo{chromo}.SNPs       = allSNPs(tmpLocs);
%     splitInfo{chromo}.outDir     = outDir;
%     splitInfo{chromo}.outName    = ['FEMA_GWAS_Chr', num2str(chr{chromo}, '%02d')];
% end

% Create parallel pool, if necessary
currPool         = gcp('nocreate');
local            = parcluster('local');
local.NumThreads = numThreads;
if isempty(currPool)
    local.parpool(numWorkers);
end

% Execute in parallel
% Here we optionally demonstrate reading of genotyping information in parts
% It is possible to read and residualize outside of parallel computing
numParts     = size(splitInfo, 1);
parfor parts = 1:numParts

    % Read data
    tmpGenomat = FEMA_parse_PLINK(bFile,     iid,    splitInfo{parts}.SNPs, ...
                                  onlyCheck, lowMem, meanImpute, roundOff,  transform, stdType);

    % Residualize data
    tmpResdGenomat = FEMA_residualizeGenotype(tmpGenomat, XVars,          clusterinfo,    ...
                                               binvec,    GroupByFamType, famtypevec,     ...
                                               OLSflag,   Ws_famtype,     Ws_fam);

    % Perform analysis - let results for each part be saved
    FEMA_sig2binseg_parfeval_GWAS_bf(tmpResdGenomat, ymat_res_bak,    clusterinfo,    ...
                                     binvec,          sig2tvec,       GroupByFamType, ...
                                     famtypevec,      Ws_famtype,     Ws_fam,         ...
                                     basisFunction,   pValType,       df,             ...
                                     OLSflag,         SingleOrDouble,                 ...
                                     splitInfo{parts}.outDir, splitInfo{parts}.outName);
end

%% Visualize QQ plot
for phen = 1:nYvars
    FEMA_plotQQ('', '', '', '', 10.^-snp_logpValues(:,phen));
end