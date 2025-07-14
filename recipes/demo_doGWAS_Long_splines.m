%% Recipe for longitudinal GWAS with spline interaction
%% Preliminaries
% Make sure that FEMA and FEMA utilities are on the path
addpath('.../cmig_tools/FEMA');
addpath('.../cmig_tools/cmig_tools_utils/matlab');

% Now define paths
dirGenetics = '';
dirOutput   = '';
dirData     = '';

% The name of the PLINK bed/bim/fam format file, extension not required
filePLINK   = '';

% The name of the saved MATLAB file which has covariates and phenotypes
fileData    = '';

% Make output directory, if it does not exist
if ~exist(dirOutput, 'dir')
    mkdir(dirOutput);
end

%% Load design matrix and phenotypic data
% Here, we assume that the data was saved as a mat file
% Hint: to read a table, try readtable; if the extension is uncommon, try
% specifying the 'FileType' as 'text'
load(fullfile(dirData, fileData));

%% Prepare design matrix
% Assuming that there is a table named 'data' which has appropriately named
% variables; edit as needed; if you already have a design matrix, skip
% these steps

% Age: assumes that age is a single column with each visit having the visit
% specific age entry on that row
age = data.age;


% Sex: dummy coding using sex == 0 as reference level
sex = double(data.sex == 1);

% Extract genetic PCs - assuming that the column names start with 'PC'
PCs = data{:, not(cellfun(@isempty, regexpi(data.Properties.VariableNames, 'PC')))};

% Extract genotyping batch - assuming column name is genotyping_batch
% Dummy coding and dropping the last level
batches     = dummyvar(categorical(data.genotyping_batch));
batches     = batches(:,1:end-1);

% Intercept - vector of 1's
intercept   = ones(height(data), 1);

% Phenotype(s) - replace with appropriate column names: here, assuming that
% there are three phenotypes called 'phen1', 'phen2', and 'phen3'
phenotypes  = [data.phen1, data.phen2, data.phen3];

% Family ID
FID = data.FID;

% Individual ID
IID = data.IID;

% Event ID - for example, 1 for first visit, 2 for second visit, so on
eid = data.TimePoint;

% How many visit are there in the data?
numVisits = length(unique(eid));

%% Mean center and standardize the phenotypes (optional)
phenotypes = (phenotypes - mean(phenotypes))./std(phenotypes);

%% Determine where knots should be placed - median of age for each visit
knots = zeros(numVisits, 1);
for v = 1:numVisits
    knots(v,1) = median(age(eid == v));
end

%% Create smooth representation of age
[basisFunction, basisSubset, Xvars, bfRank, settings, timing_basisFunction] = ...
 createBasisFunctions(age, knots, 'nsk', [], 'svd');

%% Put covariates together: intercept, s(age), sex, PCs, batch
X = [intercept, basisFunction, sex, PCs, batches];

% Ensure X is full rank: technically, it is OK to have rank deficient
% design matrix but here we generate an error to force you to think about
% what covariates you are including in your model
if rank(X) ~= size(X,2)
    error('X is not full rank - recheck covariates');
end

%% Now load the GRM that we created in demo_makeGRM
% Get information about ordering of subjects
grmInfo = load(fullfile(dirData, 'GRMData_FEMA.mat'));

% Read binary file
fid = fopen(fullfile(dirData, 'GRM.dat'), 'r');
GRM = fread(fid, [length(grmInfo.uqObservations) length(grmInfo.uqObservations)], 'double');
fclose(fid);

% Make GRM full matrix (we only saved the lower triangle)
GRM = GRM + GRM.';

% Make diagonal equal to 1
tmp = size(GRM,1);
GRM(1:1+tmp:tmp*tmp) = 1;

% Convert to single precision to reduce RAM usage
GRM = single(GRM);

%% Does GRM need re-ordering?
% If the list of subjects is smaller / has changed
[a, b, c] = unique(IID, 'stable');
if length(a) ~= length(grmInfo.uqObservations)
    reorderGRM = true;
else
    if sum(strcmpi(grmInfo.uqObservations, a)) ~= length(a)
        reorderGRM = true;
    else
        reorderGRM = false;
    end
end

% Perform reordering, if required
if reorderGRM
    [~, wch] = ismember(a, grmInfo.uqObservations);
    GRM = GRM(wch, wch);
end

%% FEMA settings
RandomEffects   = {'F', 'A', 'S', 'E'};
returnReusable  = true;
contrasts       = [];
nbins           = 0;
niter           = 1;
CovType         = 'unstructured';

%% GWAS-specific settings
% The chunkSize breaks the genome into chunk sized bits; each chunk is
% estimated separately; a small number will mean many chunks of smaller
% size (good for reduced RAM usage but may take more time) while a large
% number means smaller number of chunks (but more RAM usage); we have found
% 5000 SNPs to be a reasonable balance
chunkSize       = 5000;

% How should the genome be split? By SNPs or by chromosome?
splitBy         = 'snp';

% Ensuring double precision for evaluation; change to single precision if
% the RAM usage is too high
SingleOrDouble  = 'double';

%% Start timer
tinit = tic;

%% FEMA
[beta_hat,      beta_se,        zmat,        logpmat,                   ...
 sig2tvec,      sig2mat,        Hessmat,     logLikvec,                 ...
 beta_hat_perm, beta_se_perm,   zmat_perm,   sig2tvec_perm,             ...
 sig2mat_perm,  logLikvec_perm, binvec_save, nvec_bins,                 ...
 tvec_bins,     FamilyStruct,   coeffCovar,  reusableVars] =            ...
 FEMA_fit(X, IID, eid, FID, age, phenotypes, niter, contrasts, nbins,   ...
          GRM, 'RandomEffects', RandomEffects, 'CovType', CovType,      ...
          'returnReusable', returnReusable, 'SingleOrDouble', SingleOrDouble);

%% Get some info about genetics
[~, Chr, SNPID, BP, check, errMsg, genInfo] = ...
 FEMA_parse_PLINK(fullfile(dirGenetics, filePLINK), IID, [], true);
if ~check
    error(errMsg);
end

%% Compile some variables
[allWsTerms, tCompile] = FEMA_compileTerms(FamilyStruct.clusterinfo, binvec_save,            ...
                                           sig2mat, RandomEffects, FamilyStruct.famtypevec,  ...
                                           reusableVars.GroupByFamType, CovType, SingleOrDouble, reusableVars.visitnum);

%% Divide chromosome into chunks
[splitInfo, timing] = divideSNPs(fullfile(dirGenetics, filePLINK), splitBy, ...
                                 chunkSize, SNPID, Chr, BP, genInfo);

%% Clear up some memory
clear GRM data genInfo

%% Initialize parallel cluster
% If you do not want to do parallel computing, comment this part out
local = parcluster('local');

% How many threads per parallel worker? By default, MATLAB uses 1;
% recommend using at least 2, if not higher
local.NumThreads = 2;

% How many parallel workers? Depends on your cluster configuration
% For example, when running a slurm job, if you had 32 cores with 2
% threads, you could set this to 32 and set number of threads to 2 above
pool = local.parpool(32, 'IdleTimeout', 240);

%% Run FEMA GWAS
% Take residuals from step 1
ymat_res_gls = reusableVars.ymat_res_gls;

% Adding interceptt to the basis functions to capture the longitudinal main
% effect of the SNPs: these functions interact with every SNP
bfSNP        = [intercept, basisFunction];

% Execute each chunk in parallel - switch to regular for loop if you do not
% want parallel computing; output is saved in outDir
parfor parts = 1:length(splitInfo)
    t1 = tic;
    FEMA_fit_GWAS(splitInfo{parts}, ymat_res_gls, binvec_save, X, allWsTerms, ...
                  'outDir', dirOutput, 'outName', splitInfo{parts}.outName,   ...
                  'bfSNP', bfSNP, 'doCoeffCovar', true, 'SingleOrDouble', SingleOrDouble);

    % Also display progress
    disp(['Finished: ', num2str(parts, '%04d'), ' in ', num2str(toc(t1), '%.2f'), 's']);
end

%% End timer
tend = toc(tinit);

%% Delete pool and clear some RAM
delete(pool);
clear local splitInfo Chr SNPID BP ymat_res_gls

%% Save results 
% GWAS results are already saved by each chunk
tmp = whos;
if sum([tmp.bytes]) > 2^31
    save(fullfile(dirOutput, 'allVars.mat'), '-v7.3');
else
    save(fullfile(dirOutput, 'allVars.mat'));
end

%% Generate summary statistics
% Set cleanUp to true if you want chunk-wise summary statistics files to be
% deleted while preparing the full summary statistics
cleanUp = false;
FEMA_gatherGWAS(dirOutput, [], cleanUp, dirOutput);

%% Take a look at FEMA_convert_splines for charting SNP effect over time