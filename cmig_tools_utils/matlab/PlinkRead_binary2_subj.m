function genomat = PlinkRead_binary2_subj(nsubj, snps, subjLoc, fileprefix)
% Function to read PLINK bed file and return int8 type genotyping data,
% given a list of SNPs and a list of subjects to parse
%% Inputs:
% nsubj:        total number of subjects in the file
% snps:         a sorted list of SNP locations to fetch
% subjLoc:      a list of subject locations to fetch
% fileprefix:   full path to a PLINK bed file (without extension)
%
%% Output(s):
% genomat:      a subjLoc x snps matrix containing genotyping information,
%               such that -1 indicates missing data, 0 indicates homozygous 
%               recessive, 1 indicates heterozygous, and 2 indicates 
%               homozygous dominant
%
%% Notes:
% Based on PlinkRead_binary2 written by Chun 2015
% The code expects a sorted list of SNPs - this restriction could be
% circumvented as all fseek calls are with reference to the beginning of
% the file

persistent geno_values
 
if ~issorted(snps), error('PlinkRead_binary2 expect sorted list of snps'); end

nsnp = length(snps);

% bit shift to generate the genovalue matrix
bedprefix = sprintf('%s.bed', fileprefix);

if isempty(geno_values)
    geno_values = [
     2    2    2    2
    -1    2    2    2
     1    2    2    2
     0    2    2    2
     2   -1    2    2
    -1   -1    2    2
     1   -1    2    2
     0   -1    2    2
     2    1    2    2
    -1    1    2    2
     1    1    2    2
     0    1    2    2
     2    0    2    2
    -1    0    2    2
     1    0    2    2
     0    0    2    2
     2    2   -1    2
    -1    2   -1    2
     1    2   -1    2
     0    2   -1    2
     2   -1   -1    2
    -1   -1   -1    2
     1   -1   -1    2
     0   -1   -1    2
     2    1   -1    2
    -1    1   -1    2
     1    1   -1    2
     0    1   -1    2
     2    0   -1    2
    -1    0   -1    2
     1    0   -1    2
     0    0   -1    2
     2    2    1    2
    -1    2    1    2
     1    2    1    2
     0    2    1    2
     2   -1    1    2
    -1   -1    1    2
     1   -1    1    2
     0   -1    1    2
     2    1    1    2
    -1    1    1    2
     1    1    1    2
     0    1    1    2
     2    0    1    2
    -1    0    1    2
     1    0    1    2
     0    0    1    2
     2    2    0    2
    -1    2    0    2
     1    2    0    2
     0    2    0    2
     2   -1    0    2
    -1   -1    0    2
     1   -1    0    2
     0   -1    0    2
     2    1    0    2
    -1    1    0    2
     1    1    0    2
     0    1    0    2
     2    0    0    2
    -1    0    0    2
     1    0    0    2
     0    0    0    2
     2    2    2   -1
    -1    2    2   -1
     1    2    2   -1
     0    2    2   -1
     2   -1    2   -1
    -1   -1    2   -1
     1   -1    2   -1
     0   -1    2   -1
     2    1    2   -1
    -1    1    2   -1
     1    1    2   -1
     0    1    2   -1
     2    0    2   -1
    -1    0    2   -1
     1    0    2   -1
     0    0    2   -1
     2    2   -1   -1
    -1    2   -1   -1
     1    2   -1   -1
     0    2   -1   -1
     2   -1   -1   -1
    -1   -1   -1   -1
     1   -1   -1   -1
     0   -1   -1   -1
     2    1   -1   -1
    -1    1   -1   -1
     1    1   -1   -1
     0    1   -1   -1
     2    0   -1   -1
    -1    0   -1   -1
     1    0   -1   -1
     0    0   -1   -1
     2    2    1   -1
    -1    2    1   -1
     1    2    1   -1
     0    2    1   -1
     2   -1    1   -1
    -1   -1    1   -1
     1   -1    1   -1
     0   -1    1   -1
     2    1    1   -1
    -1    1    1   -1
     1    1    1   -1
     0    1    1   -1
     2    0    1   -1
    -1    0    1   -1
     1    0    1   -1
     0    0    1   -1
     2    2    0   -1
    -1    2    0   -1
     1    2    0   -1
     0    2    0   -1
     2   -1    0   -1
    -1   -1    0   -1
     1   -1    0   -1
     0   -1    0   -1
     2    1    0   -1
    -1    1    0   -1
     1    1    0   -1
     0    1    0   -1
     2    0    0   -1
    -1    0    0   -1
     1    0    0   -1
     0    0    0   -1
     2    2    2    1
    -1    2    2    1
     1    2    2    1
     0    2    2    1
     2   -1    2    1
    -1   -1    2    1
     1   -1    2    1
     0   -1    2    1
     2    1    2    1
    -1    1    2    1
     1    1    2    1
     0    1    2    1
     2    0    2    1
    -1    0    2    1
     1    0    2    1
     0    0    2    1
     2    2   -1    1
    -1    2   -1    1
     1    2   -1    1
     0    2   -1    1
     2   -1   -1    1
    -1   -1   -1    1
     1   -1   -1    1
     0   -1   -1    1
     2    1   -1    1
    -1    1   -1    1
     1    1   -1    1
     0    1   -1    1
     2    0   -1    1
    -1    0   -1    1
     1    0   -1    1
     0    0   -1    1
     2    2    1    1
    -1    2    1    1
     1    2    1    1
     0    2    1    1
     2   -1    1    1
    -1   -1    1    1
     1   -1    1    1
     0   -1    1    1
     2    1    1    1
    -1    1    1    1
     1    1    1    1
     0    1    1    1
     2    0    1    1
    -1    0    1    1
     1    0    1    1
     0    0    1    1
     2    2    0    1
    -1    2    0    1
     1    2    0    1
     0    2    0    1
     2   -1    0    1
    -1   -1    0    1
     1   -1    0    1
     0   -1    0    1
     2    1    0    1
    -1    1    0    1
     1    1    0    1
     0    1    0    1
     2    0    0    1
    -1    0    0    1
     1    0    0    1
     0    0    0    1
     2    2    2    0
    -1    2    2    0
     1    2    2    0
     0    2    2    0
     2   -1    2    0
    -1   -1    2    0
     1   -1    2    0
     0   -1    2    0
     2    1    2    0
    -1    1    2    0
     1    1    2    0
     0    1    2    0
     2    0    2    0
    -1    0    2    0
     1    0    2    0
     0    0    2    0
     2    2   -1    0
    -1    2   -1    0
     1    2   -1    0
     0    2   -1    0
     2   -1   -1    0
    -1   -1   -1    0
     1   -1   -1    0
     0   -1   -1    0
     2    1   -1    0
    -1    1   -1    0
     1    1   -1    0
     0    1   -1    0
     2    0   -1    0
    -1    0   -1    0
     1    0   -1    0
     0    0   -1    0
     2    2    1    0
    -1    2    1    0
     1    2    1    0
     0    2    1    0
     2   -1    1    0
    -1   -1    1    0
     1   -1    1    0
     0   -1    1    0
     2    1    1    0
    -1    1    1    0
     1    1    1    0
     0    1    1    0
     2    0    1    0
    -1    0    1    0
     1    0    1    0
     0    0    1    0
     2    2    0    0
    -1    2    0    0
     1    2    0    0
     0    2    0    0
     2   -1    0    0
    -1   -1    0    0
     1   -1    0    0
     0   -1    0    0
     2    1    0    0
    -1    1    0    0
     1    1    0    0
     0    1    0    0
     2    0    0    0
    -1    0    0    0
     1    0    0    0
     0    0    0    0
    ];
    geno_values = cast(geno_values, 'int8');
end

if isempty(geno_values)
    geno_values = zeros(256,4,'int8');
    geno_code   = [-1,1,0,2];
    shiftind    = [0,2,4,6];
    indvec      = zeros(1,4);

    for i = 1:256
        ishift  = int16(i-1);
        for j   = 1:4
            indvec(j) = bitand(bitsra(ishift,shiftind(j)),3) ;
        end
        indvec(indvec == 0) = 4;
        geno_values(i,:) = geno_code(indvec);
    end
end

% Read in the binary file
bedid   = fopen(bedprefix);
genobin = uint16(fread(bedid, 3));

% Check magic number
if genobin(1) ~= 108
	error('- Not a valid Plink BED file \r\n');
elseif genobin(2) ~= 27
	error('- Not a valid Plink BED file \r\n');
elseif genobin(3) ~= 1
	error('- Not in SNP-major format \r\n');
end

% Determine block length - number of subjects/4 rounded up
n_bytes = ceil(nsubj/4);

% Initialize to number of subjects and SNPs to subset - setting to int8
% Seems to be more efficient to subset in one go at the end, although this
% will increase memory requirement as all subjects need to be read
genomat = zeros(length(subjLoc), nsnp, 'int8');
% genmat = zeros(nsubj, nsnp, 'int8');

% Go over every SNP location
for i = 1:nsnp

    % Skip the first three bytes; move from beginning of file to the start
    % of the SNP specific line which contains data for all subjects stored
    % as two bytes
    fseek(bedid, 3 + (snps(i) - 1) * n_bytes, 'bof');
    genobin = uint16(fread(bedid, n_bytes));
    if length(genobin) ~= n_bytes
        error('-- Invalid number of entries from the bed \r\n'); 
    end
    tmp_values   = geno_values(genobin + 1, :)';
    tmp_values   = tmp_values(:);

    genomat(:,i) = tmp_values(subjLoc);
    % genmat(:,i) = tmp_values(1:nsubj);
end
fclose(bedid);

% Only retain subject locations that are relevant
% genomat = genmat(subjLoc,:);

end