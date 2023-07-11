function splitInfo = divideSNPs(bFile, splitBy, chunkSize, SNPID, Chr)
% Function that divides SNPs into either chunks or by chromosomes
%% Inputs:
% bFile:        full path to a bFile without extension
% splitBy:      should be either of the following:
%                   * 'snp' (default)
%                   * 'chromosome'
% chunkSize:    if splitBy is 'snp', then chunkSize indicates the number of
%               SNPs in each chunk - defaults to 1000
%
% Optionally, it is possible to skip parsing a bFile by providing SNPID and
% (optionally Chr) along with other inputs
%
%% Output(s):
% splitInfo:    cell type; each cell contains a structure with the
%               following fields:
%                   * 'Locs':       locations of the SNPs
%                   * 'SNPs':       names of the SNPs
%                   * 'outName':    temporary name for saving mat file for
%                                   every chromosome/chunk

%% Parse inputs and assign defaults
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
end

% Check bFile and other parameters
if ~exist('bFile', 'var') || isempty(bFile)
    skipBFile = true;
    % Check if SNPID exists
    if ~exist('SNPID', 'var') || isempty(SNPID)
        error('Both bFile and SNPID not provided; unable to proceed');
    else
        % Check if splitBy is chromosome; if so, check if Chr exists
        if strcmpi(splitBy, 'chromosome')
            if ~exist('Chr', 'var') || isempty(Chr)
                error('Both bFile and Chr not provided; unable to proceed');
            end
        end
    end
else
    skipBFile = false;
end

%% Parse the PLINK file, if required
if ~skipBFile
    [~, Chr, SNPID, ~, check, errMsg] = FEMA_parse_PLINK(bFile, [], [], true);
    if ~check
        error(errMsg);
    end
end

%% Divide into chunks
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
        splitInfo{chunk}.Locs    = tmpLocs;
        splitInfo{chunk}.SNPs    = SNPID(tmpLocs);
        splitInfo{chunk}.outName = ['FEMA_GWAS_Chunk-', num2str(chunk, '%05d')];
    end
else
    % Identify chromosomes
    chr      = unique(Chr, 'stable');
    numChr   = length(chr);

    % Split genotype across chromosomes
    splitInfo  = cell(numChr,1);
    for chromo = 1:numChr
        tmpLocs                   = strcmpi(Chr, chr{chromo});
        splitInfo{chromo}.Locs    = tmpLocs;
        splitInfo{chromo}.SNPs    = SNPID(tmpLocs);
        splitInfo{chromo}.outName = ['FEMA_GWAS_Chr-', num2str(chr{chromo}, '%02d')];
    end
end