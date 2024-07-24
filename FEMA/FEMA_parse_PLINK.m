function [genomat, Chr, SNPID, BP, check, errMsg, genInfo, tInfo] =             ...
          FEMA_parse_PLINK(bFile,      iid,      SNPList,   onlyCheck, lowMem,  ...
                           meanImpute, roundOff, transform, stdType,   genInfo)
% Function to parse PLINK files
%% Inputs:
% bFile:            character       full path to a PLINK file (no extension)
%
% iid:              cell type       list of subjects for which to extract 
%                                   genotyping data; leave empty if all
%                                   subjects in the fam file should be read
%
% SNPList:          cell type       list of SNPs to read; leave empty if
%                                   all SNPs in the bim file should be read
%
% onlyCheck:        logical         if true, a check is performed on 
%                                   whether all IIDs exist in the PLINK fam 
%                                   file and whether all SNPs in SNPList 
%                                   exist in the PLINK bim file; output 
%                                   is true if both conditions are met;
%                                   additionally, Chr, SNPID, and basePair 
%                                   information are also returned
%
% lowMem:           logical         if true, genomat is returned as int8 
%                                   data type (output by PlinkRead_binary2);
%                                   no other operations (convert missing to
%                                   NaN, mean imputation, standardization, 
%                                   or rounding) are applied, irrespective
%                                   of the arguments provided by the user
%
% meanImpute:       logical         if true, replace NaN in genotyping 
%                                   matrix with column wise mean value
% 
% roundOff:         logical         if true, mean imputed data is converted
%                                   back to integers (0, 1, or 2); this is
%                                   implemented using the round function
%                                   and any value which exceeds 2, after 
%                                   rounding, is truncated to 2; only 
%                                   meaningful for hardcalls and when 
%                                   impute is true
%
% transform:       	character       type of transformation to perform on
%                                   the genotyping data (see Notes):
%                                       * 'center'
%                                       * 'centre'
%                                       * 'std'
%                                       * 'none'
%
% stdType:          character       determine how the standardization
%                                   should be performed (only relevant when 
%                                   transform is set to std; see Notes):
%                                       * 'emperical'
%                                       * 'gcta'
%
% genInfo:          structure       if genInfo structure containing the
%                                   fields 'locIID', 'locSNPs', and
%                                   'numSubjs' is present, checking is
%                                   skipped (see outputs for details on
%                                   these fields)
%
%% Outputs:
% genomat:          [n x m]         matrix of m SNPs for n subjects (iid)
%
% Chr:              [m x 1]         cell type having chromosome number for
%                                   the m SNPs (character type inside)
%
% SNPID:            [m x 1]         cell type having SNP ID for the m SNPs
%
% BP:               [m x 1]         cell type having base pair coordinates
%                                   for the m SNPs
%
% check:            logical         true if all IDs in IID are present in
%                                   the PLINK fam file AND if all SNPs in 
%                                   SNPList are present in the PLINK bim 
%                                   file; if check fails and onlyCheck is 
%                                   false, an error is generated
%
% errMsg:           character       if check fails, a message indicating
%                                   which of the checks (IID/SNPs) failed;
%                                   empty otherwise
% 
% tInfo:            structure       timing information for different
%                                   aspects of this script; useful for
%                                   identifying bottleneck
%
% genInfo:          structure       contains variables that can be reused
%                                   when reading the full data;
%                                   specifically, contains the following:
%
% locIID:           [n x 1]         vector having the location in the
%                                   genomat for the subjects in IID
%
% locSNPs:          [m x 1]         vector having the location of the SNPs
%                                   within the PLINK files for the SNPs in
%                                   SNPList
%
% numSubjs:         [1 x 1]         number of subjects in the genetics file
%
%% Notes:
% Genotyping matrix can be transformed in the following mways:
% center:           mean centering of the genotyping matrix
% centre:           same as center
% std:              mean center and scaling by standard deviation
% 
% Standardization can be implemented in two ways:
% Emperical:        mean (and standard deviation) of each column of
%                   genotyping matrix is used
% GCTA style:       allele frequency is computed and used
%
% If transformation is set to centering, column wise mean is used
% 
% When using GCTA style standardization, allele frequencies are computed by
% counting integer values of 0, 1, and 2
% 
% If mean imputation is performed, the imputed data will likely contain
% fractional values; these values will not be considered when calculating
% allele frequencies (GCTA style); this behaviour can be modified when
% roundOff is set to true
% 
% For emperical style data transformation, both the mean and the standard
% deviations are calculated by omitting the NaN values; therefore, if mean
% imputation is not done and there are NaN values, the resulting
% transformed genotyping matrix will contain valid numbers, except for
% locations with NaN values
%
% Additionally, note that rounding can make a bit of a difference for 
% emperical transformation as well - when a missing SNP is replaced 
% by the mean of the non-missing values, then during data transformation 
% this SNP will get a value of 0 because the imputed value will be equal 
% to the mean; on the other hand, if rounding to nearest integer is done 
% prior to data transformation then this SNP will get the same value as the
% other SNPs with the same integer value
%
% For calculating allele frequencies, counting them independently seeems to
% be faster than using histc: counts = histc(genomat, [0 1 2]);
% Using histcounts is efficient but it doesn't count across each column
%
%% Defaults:
% iid:              read all
% SNPList:          read all
% onlyCheck:        false
% lowMem:           false
% meanImpute:       false
% roundOff:         false
% transform:        'none'
% stdType:          'emperical'
% locIID:           []
% locSNPs:          []

%% Check inputs and assign defaults
% Initialize timer
tInit         = tic;
ticInputCheck = tic;

% Check for fam, bim, and bed files
[loc, prefix] = fileparts(bFile);
if ~exist([bFile, '.fam'], 'file') || ...
   ~exist([bFile, '.bim'], 'file') || ...
   ~exist([bFile, '.bed'], 'file')
    error(['Unable to find files; check that fam, bim, and bed files with prefix: ', ...
           prefix, ' exist in: ', loc]);
end

% Check IIDs
if ~exist('iid', 'var') || isempty(iid)
    readsubjAll = true;
else
    readsubjAll = false;
end

% Check list of SNPs
if ~exist('SNPList', 'var') || isempty(SNPList)
    readSNPAll = true;
else
    readSNPAll = false;
end

% Check onlyCheck
if ~exist('onlyCheck', 'var') || isempty(onlyCheck)
    onlyCheck = false;
else
    if not(islogical(onlyCheck))
        error('onlyCheck should be either true or false');
    end
end

% Check lowmem
if ~exist('lowMem', 'var') || isempty(lowMem)
    lowMem = false;
else
    if not(islogical(lowMem))
        error('lowmem should be either true or false');
    end
end

% Check if data should be mean imputed
if ~exist('meanImpute', 'var') || isempty(meanImpute)
    meanImpute = false;
else
    if not(islogical(meanImpute))
        error('meanImpute should be either true or false');
    end
end

% Check if rounding up is necessary
if ~exist('roundOff', 'var') || isempty(roundOff)
    roundOff = false;
else
    if not(islogical(roundOff))
        error('roundOff should be either true or false');
    end
end

% Check which data transformation method to apply
if ~exist('transform', 'var') || isempty(transform)
    transform = 'none';
else
    transform = lower(transform);
    if not(ismember(transform, {'center', 'centre', 'std', 'none'}))
        error(['Unrecognized data transform specified: ', transform, ...
               '; should be one of: center, centre, std, or none']);
    end
end

% Data standardization type
if ~exist('stdType', 'var') || isempty(stdType)
    stdType = 'none';
else
    stdType = lower(stdType);
    if not(ismember(stdType, {'gcta', 'emperical', 'none'}))
        error(['Unrecognized transformation type provided: ', stdType, ...
               '; should be one of: gcta, emperical, or none']);
    else
        if strcmpi(transform, 'none') && not(strcmpi(stdType, 'none'))
            warning('Setting standardization type to none');
            stdType = 'none';
        end
    end
end

% Check if reading fam and bim files can be skipped
if ~exist('genInfo', 'var') || isempty(genInfo)
    skipCheck = false;
else
    if ~isstruct(genInfo)
        error('genInfo should be a structure with locIID, locSNPs, and numSubjs fields');
    else
        if ~all(isfield(genInfo, {'locIID', 'locSNPs', 'numSubjs'}))
            warning('genInfo was provided but one of the following was missing: locIID, locSNPs, numSubjs');
            skipCheck = false;
        else
            skipCheck = true;
        end
    end
end

% If skipCheck is true, set readsubjAll to false
if skipCheck
    readsubjAll = false;
end

% End of input checking
tInfo.inputChecking = toc(ticInputCheck);

%% Sanity check
% Initialize timer
ticSanityCheck = tic;

if ~skipCheck
    % Read fam file
    fhandle  = fopen([bFile, '.fam'], 'r');
    fam_file = textscan(fhandle, '%s %s %s %s %s %s');
    fclose(fhandle);
    
    % Read bim file
    fhandle  = fopen([bFile, '.bim'], 'r');
    bim_file = textscan(fhandle, '%s %s %s %s %s %s');
    fclose(fhandle);
    
    % Find the location of the IIDs in the fam file
    if readsubjAll
        iid = fam_file{2};
    end
    [~, locIID] = ismember(iid, fam_file{2});
    if sum(locIID == 0) ~= 0
        check_IID = false;
    else
        check_IID = true;
    end
    
    % Make sure that all SNPs in listSNPs exist in bim file
    if readSNPAll
        SNPList = bim_file{2};
    end
    % Find the location of the SNPs
    [~, locSNPs] = ismember(SNPList, bim_file{2});
    if sum(locSNPs == 0) ~= 0
        check_SNPs = false;
    else
        check_SNPs = true;
    end
    
    % Overall check variable
    check = check_IID & check_SNPs;
    if ~check
        if check_IID && ~check_SNPs
            errMsg = 'One or more SNPs not present in the bim file';
        else
            if ~check_IID && check_SNPs
                errMsg = 'One or more IDs not present in the fam file';
            else
                errMsg = 'One or more IDs and one or more SNPs not present in fam/bim files';
            end
        end
    else
        errMsg = '';
    end

    % Chromosome number is bim{1}, rsID is bim{2}, and basePair is bim{4};
    % number of subjects in the genetics file is the length of fam_file{2}
    Chr      = bim_file{1}(locSNPs);
    SNPID    = bim_file{2}(locSNPs);
    BP       = bim_file{4}(locSNPs);
    numSubjs = length(fam_file{2});

    % Save some information for posteriety
    genInfo.locIID   = locIID;
    genInfo.locSNPs  = locSNPs;
    genInfo.numSubjs = numSubjs;
else
    % Extract values from structure
    locIID   = genInfo.locIID;
    locSNPs  = genInfo.locSNPs;
    numSubjs = genInfo.numSubjs;

    % Assign empty output values
    errMsg  = '';
    Chr     = '';
    SNPID   = '';
    BP      = '';
    check   = true;
end

% End of sanity check
tInfo.sanityCheck = toc(ticSanityCheck);

%% Exit, if onlyCheck
if onlyCheck
    tInfo.overall = toc(tInit);
    genomat       = [];
    return
else
    error(errMsg);
end

%% Read data
ticRead = tic;
if readsubjAll
    genmat  = PlinkRead_binary2(numSubjs, locSNPs, bFile);
else
    genmat  = PlinkRead_binary2_subj(numSubjs, locSNPs, locIID, bFile);
end
tInfo.readPlink = toc(ticRead);

%% If lowMemMode, exit
if lowMem
    if readsubjAll
        genomat = genmat(locIID, :);
    else
        genomat = genmat;
    end
    return;
end

%% Cast as single precision
ticSingle = tic;
if readsubjAll
    genomat = single(genmat(locIID, :));
else
    genomat = single(genmat);
end
tInfo.singlePrecision = toc(ticSingle);

%% Set out of range to NaN
% tmpLocs = genomat ~= 0 & genomat ~= 1 & genomat ~= 2;
ticNaN           = tic;
tmpLocs          = genmat == -1;
genomat(tmpLocs) = NaN;
tInfo.setNaN     = toc(ticNaN);

%% Mean impute NaN, if necessary
if meanImpute
    % Replace with column wise means
    ticMeanImpute        = tic;
    colMeans             = mean(genomat, 'omitnan');
    tmp                  = repmat(colMeans, [size(genomat, 1), 1]);
    toImpute             = isnan(genomat(:));
    genomat(toImpute)    = tmp(toImpute);
    tInfo.meanImpute     = toc(ticMeanImpute);
end

%% Round up mean imputed data, if necessary
if roundOff
    genomat = round(genomat, 0);
    genomat(genomat > 2) = 2;
end

%% Transform data, if necessary
switch transform
    case 'center'
        ticTransform          = tic;
        genomat               = (genomat - mean(genomat, 'omitnan'));
        tInfo.standardize     = toc(ticTransform);
        
    case 'std'
        if strcmpi(stdType, 'emperical')
            ticTransform          = tic;
            genomat               = (genomat - mean(genomat, 'omitnan'))./std(genomat, [], 'omitnan');
            tInfo.standardize     = toc(ticTransform);
        else
            ticTransform          = tic;
            counter0              = sum(genomat == 0);
            counter1              = sum(genomat == 1);
            counter2              = sum(genomat == 2);
            numObs                = counter0 + counter1 + counter2;
            alleleFreq            = (counter1 + counter2*2)./(numObs*2);
            genomat               = (genomat - (2 * alleleFreq)) ./ ...
                                    (sqrt(2 * alleleFreq .* (1 - alleleFreq)));
            tInfo.standardize     = toc(ticTransform);
        end
end

% Overall time taken
tInfo.overall = toc(tInit);