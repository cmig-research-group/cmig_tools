%% Recipe for preparing genetic relationship matrix (GRM)
% Generally, creating GRM requires a lot of time and RAM; suggest running
% this on a cluster with large resources

%% Preliminaries
% Make sure that FEMA and FEMA utilities are on the path
addpath('.../cmig_tools/FEMA');
addpath('.../cmig_tools/cmig_tools_utils/matlab');

% Now define paths
dirGenetics = '';
dirOutput   = '';

% The name of the PLINK bed/bim/fam format file, extension not required
filePLINK   = '';

% A pruned list of SNPs to read for preparing GRM
fileSNPs    = '';

% Make output directory, if it does not exist
if ~exist(dirOutput, 'dir')
    mkdir(dirOutput);
end

%% Read in the list of SNPs
fid  = fopen(fullfile(dirGenetics, fileSNPs), 'r');
snps = textscan(fid, '%s');
fclose(fid);
snps = snps{1};

%% Prepare list of unique subjects
% Either read in the design matrix with 'IID' column, or something similar
% Here, we assume that a table type variable 'data' with IID column name
% Hint: to read a table, try readtable; if the extension is uncommon, try
% specifying the 'FileType' as 'text'
%
% FEMA expects the GRM ordering to be the same as unique(IID, 'stable');
[uqObservations, ~, c]  = unique(data.IID,  'stable');

%% Get list of SNPs to work with
bFile = fullfile(dirGenetics, filePLINK);
[~, Chr, SNPID, basePair, check, errMsg] = FEMA_parse_PLINK(bFile, uqObservations, [], true);
[a, b]   = ismember(snps, SNPID);
toSample = SNPID(sort(b));

%% Read standardized full genetic data for calculating GRM
% The GRM is the correlation coefficient between entries in the
% standardized genotyping matrix (0,1,2)
% Here, we read the genotyping matrix, and standardize it (transform ==
% 'std') using the observed allele frequency (stdType == 'emperical')
onlyCheck   = false;
lowMem      = false;
meanImpute  = true;
roundOff    = false;
transform   = 'std';
stdType     = 'emperical';

genomat = FEMA_parse_PLINK(bFile, uqObservations, toSample, onlyCheck, ...
                           lowMem, meanImpute, roundOff, transform, stdType);

%% Compute GRM as the correlation coefficient
% Note that this step requires a significant amount of time and RAM
GRM = corrcoef(genomat');

%% Extract lower triangle, followed by conversion to double precision
GRM = double(tril(GRM));

%% Save as binary file
% You can, alternatively, save this as a mat file (will probably need -v7.3
% switch) or some other file format such as HD5 or parquet (probably
% faster)
fid = fopen(fullfile(dirOutput, 'GRM.dat'), 'w');
fwrite(fid, GRM, 'double');
fclose(fid);

%% Save the order of subjects for later reference
% If previously saving as a mat file, all variables can be saved together
save(fullfile(dirOutput, 'GRMData_FEMA.mat'), 'uqObservations', 'c', '-v7.3');