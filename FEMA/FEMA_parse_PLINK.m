function [genomat, Chr, SNPID, basePair, check] = FEMA_parse_PLINK(bFile, iid, onlyCheck, stdType, meanImputeNaN)
% Function to parse PLINK files
%% Inputs:
% bFile:            character       full path to a PLINK file (no extension)
%
% iid:              cell type       list of subjects for which to extract 
%                                   genotyping data
%
% onlyCheck:        logical         if true, only a check is performed on 
%                                   whether all IIDs exist in the PLINK  
%                                   file and check = true is returned
%
% stdType:          character       type of standardization to perform on
%                                   the genotyping data (see Notes):
%                                       * 'gcta'
%                                       * 'emperical'
%                                       * 'none'
%
% meanImputeNan:    logical         whether missing genotyping data should
%                                   be replaced by NaN
% 
%% Outputs:
% genomat:          [n x m]         matrix of m SNPs for n subjects (iid)
%
% Chr:              [m x 1]         cell type having chromosome number for
%                                   the m SNPs (character type inside)
%
% SNPID:            [m x 1]         cell type having SNP ID for the m SNPs
%
% basePair:         [m x 1]         cell type having base pair coordinates
%                                   for the m SNPs
%
% check:            logical         true if all IDs in IID are present in
%                                   the PLINK files; if check fails and
%                                   onlyCheck is false, an error is
%                                   generated
%
%% Notes:
% Standardization of the genotype can be done in two ways:
% GCTA style: allele frequency is computed and used
% emperical:  mean and standard deviation of each column of genomat is used
%
%% ToDo:
% Read data in chunks
% Take a list of SNPs and only read those

% Check for fam, bim, and bed files
[loc, prefix] = fileparts(bFile);
if ~exist([bFile, '.fam'], 'file') || ...
   ~exist([bFile, '.bim'], 'file') || ...
   ~exist([bFile, '.bed'], 'file')
    error(['Unable to find files; check that fam, bim, and bed files with prefix: ', ...
           prefix, ' exist in: ', loc]);
end

if ~exist('stdType', 'var') || isempty(stdType)
    stdType = 'none';
else
    stdType = lower(stdType);
    if not(ismember(stdType, {'gcta', 'emperical', 'none'}))
        error(['Unrecognized standardization type provided: ', stdType, ...
               '; should be one of: gcta, emperical, or none']);
    end
end

if ~exist('onlyCheck', 'var')
    onlyCheck = false;
end

if ~exist('meanImputeNaN', 'var') || isempty(meanImputeNaN)
    meanImputeNaN = false;
end

% Read fam file
fhandle  = fopen([bFile, '.fam'], 'r');
fam_file = textscan(fhandle, '%s %s %s %s %s %s');
fclose(fhandle);

% Make sure that all IDs in iid exist in fam file
if sum(ismember(fam_file{2}, iid)) ~= length(iid)
    check = false;
else
    check = true;
end

if onlyCheck
    genomat  = [];
    Chr      = [];
    SNPID    = [];
    basePair = [];
    return;
end

if ~check
    error('One or more IDs in iid are not present in the fam file');
else
    % Read bim file
    fhandle  = fopen([bFile, '.bim'], 'r');
    bim_file = textscan(fhandle, '%s %s %s %s %s %s');
    fclose(fhandle);
    
    % Chromosome number is bim{1}, rsID is bim{2}, and basePair is bim{4}
    Chr      = bim_file{1};
    SNPID    = bim_file{2};
    basePair = bim_file{4};
    
    % Read bed file
    genomat = PlinkRead_binary2(length(fam_file{2}), 1:length(bim_file{1}), bFile);
    
    % Subset genotype data
    [~, b]  = ismember(iid, fam_file{2});
    genomat = genomat(b, :);
    
    % Cast as single precision
    genomat = single(genomat);
    genomat(genomat ~= 0 & genomat ~= 1 & genomat ~= 2) = NaN;
    % genomat(genomat == -1)  = NaN;
    % gmat     = NaN(size(genomat), 'single');
    % for code = int8([0, 1, 2])
    %     gmat(genomat == code) = single(code);
    % end
    
    % Mean impute NaN, if necessary
    if meanImputeNaN
        % This implementation is not optimal
        tmp = repmat(nanmean(genomat), [size(genomat, 1), 1]);
        genomat(~isfinite(genomat(:))) = tmp(~isfinite(genomat(:)));
        clear tmp
    end
    
    % Standardize
    switch stdType
        case 'gcta'
            % Standardization is done using allele frequency; for
            % calculating allele frequency, counting them independently is
            % still faster than using histc:
            % counts = histc(genomat, [0 1 2]);
            % Using histcounts is more efficient but it doesn't count
            % across each column
            counter0    = sum(genomat == 0);
            counter1    = sum(genomat == 1);
            counter2    = sum(genomat == 2);
            numObs      = counter0 + counter1 + counter2;
            alleleFreq  = (counter1 + counter2*2)./(numObs*2);
            genomat     = (genomat - (2 * alleleFreq)) ./ ...
                          (sqrt(2 * alleleFreq .* (1 - alleleFreq)));

        case 'emperical'
            genomat = (genomat - mean(genomat))./std(genomat);
    end
end