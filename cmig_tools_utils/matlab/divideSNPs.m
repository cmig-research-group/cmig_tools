function [splitInfo, timing] = divideSNPs(bFile, splitBy, chunkSize, SNPID, Chr, BP, genInfo)
% Function that divides SNPs into chunks or by chromosomes
%% Input(s):
% bFile:        character       full path to a bFile without extension
%
%% Optional inputs: 
% splitBy:      character       should be either of the following:
%                                   * 'snp' (default)
%                                   * 'chromosome'
%
% chunkSize:    [1 x 1]         if splitBy is 'snp', then chunkSize
%                               indicates the number of SNPs in each chunk
%                               (defaults to 1000); if splitBy is
%                               'chromosome', then chunkSize is not
%                               relevant
%
%% Other inputs:
% Parsing of bFile can be skipped if these are additionally input:
% Chr:          [m x 1]         cell type having chromosome number for
%                               the m SNPs
%
% SNPID:        [m x 1]         cell type having SNP ID for the m SNPs
%
% BP:           [m x 1]         cell type having base pair coordinates
%                               for the m SNPs
%
% genInfo:      structure       should contain the following fields:
%                                   * 'locSNPs'
%                                   * 'locIID'
%                                   * 'numSubjs'
%
%% Output(s):
% splitInfo:    cell type       each cell contains a structure with the
%                               following fields for chunk/chromosome: 
%                                   * 'fname':   full path to PLINK file
%                                   * 'Locs':    index locations of SNPs
%                                   * 'Chr':     chromosome number of SNPs
%                                   * 'SNPs':    names of the SNPs
%                                   * 'BP':      base pair position
%                                   * 'outName': temporary name for saving
%                                                mat file for every part
%                                   * 'genInfo': contains 'locSNPs',
%                                                'locIID', and 'numSubjs'
%
% timing:       structure       contains timing information

%% Parse inputs and assign defaults
tOver = tic;
tInit = tic;

% Check bFile
if ~exist('bFile', 'var') || isempty(bFile)
    error('Please provide a full path to a PLINK file');
end

% Check splitBy
if ~exist('splitBy', 'var') || isempty(splitBy)
    splitBy = 'snp';
else
    splitBy = lower(splitBy);
    if ~ismember(splitBy, {'snp', 'chromosome'})
        error(['Unknown splitBy value passed: ', splitBy, ' ; should be either snp or chromosome']);
    end
end

% Check chunkSize
if ~exist('chunkSize', 'var') || isempty(chunkSize)
    chunkSize = 1000;
else
    if strcmpi(splitBy, 'snp')
        if ~isnumeric(chunkSize) && ~isscalar(chunkSize)
            error('If splitBy is snp, chunkSize should be a numeric scalar');
        end
    end
end

% Check if bFile reading can be skipped
if (exist('SNPID',   'var') && ~isempty(SNPID))   && ...
   (exist('Chr',     'var') && ~isempty(Chr))     && ...
   (exist('BP',      'var') && ~isempty(BP))      && ...
   (exist('genInfo', 'var') && ~isempty(genInfo)) && ...
    all(isfield(genInfo, {'locIID', 'locSNPs', 'numSubjs'}))
    skipBFile = true;
else
    skipBFile = false;
end

% Save timing for checks
timing.tChecks = toc(tInit);

%% Parse the PLINK file, if required
tPLINK = tic;
if ~skipBFile
    [~, Chr, SNPID, BP, check, errMsg, genInfo] = FEMA_parse_PLINK(bFile, [], [], true);
    if ~check
        error(errMsg);
    end
end

% Save timing for getting info from PLINK files
timing.tPLINK = toc(tPLINK);

%% Divide into chunks
tDivide = tic;
if strcmpi(splitBy, 'snp')
    nSNPs       = length(SNPID);
    allChunks   = 1:chunkSize:nSNPs;
    numChunks   = length(allChunks);
    splitInfo   = cell(numChunks, 1);

    % Divide data into chunks
    for chunk = 1:numChunks
        if chunk == numChunks
            tmpLocs = allChunks(chunk):nSNPs;
        else
            tmpLocs = allChunks(chunk):allChunks(chunk)+chunkSize-1;
        end

        % Assign fields
        splitInfo{chunk}.fname            = bFile;
        splitInfo{chunk}.Locs             = tmpLocs;
        splitInfo{chunk}.Chr              = Chr(tmpLocs);
        splitInfo{chunk}.SNPs             = SNPID(tmpLocs);
        splitInfo{chunk}.BP               = BP(tmpLocs);
        splitInfo{chunk}.outName          = ['FEMA_GWAS_Chunk-', num2str(chunk, '%05d')];
        splitInfo{chunk}.genInfo.locSNPs  = genInfo.locSNPs(tmpLocs);
        splitInfo{chunk}.genInfo.locIID   = genInfo.locIID;
        splitInfo{chunk}.genInfo.numSubjs = genInfo.numSubjs;
    end
else
    % Identify chromosomes
    chr      = unique(Chr, 'stable');
    numChr   = length(chr);

    % Split genotype across chromosomes
    splitInfo  = cell(numChr,1);
    for chromo = 1:numChr
        tmpLocs                            = strcmpi(Chr, chr{chromo});
        splitInfo{chromo}.fname            = bFile;
        splitInfo{chromo}.Locs             = tmpLocs;
        splitInfo{chromo}.Chr              = Chr(tmpLocs);
        splitInfo{chromo}.SNPs             = SNPID(tmpLocs);
        splitInfo{chromo}.BP               = BP(tmpLocs);
        splitInfo{chromo}.outName          = ['FEMA_GWAS_Chr-', num2str(chr{chromo}, '%02d')];
        splitInfo{chromo}.genInfo.locSNPs  = genInfo.locSNPs(tmpLocs);
        splitInfo{chromo}.genInfo.locIID   = genInfo.locIID;
        splitInfo{chromo}.genInfo.numSubjs = genInfo.numSubjs;
    end
end

% Save timing information for dividing data into parts
timing.tDivide = toc(tDivide);

% Overall timing information
timing.tOverall = toc(tOver);