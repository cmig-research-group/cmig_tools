function FEMA_writeGWAS(Chr, SNPID, basePair, sampleSize, beta, SE, stats, log10pValues, df, outDir, fname)
% Function to write out a formatted GWAS table for every phenotype
%% Inputs:
% Chr:              [k x 1]     vector of chromosome numbers for k SNPs
% 
% SNPID:            [k x 1]     vector of SNP ID for k SNPs
% 
% basePair:         [k x 1]     vector of base pair positions
% 
% sampleSize:       [1 x 1]     sample size used for analysis
% 
% beta:             [k x v]     matrix of beta values for k SNPs and v
%                               outcome variables
% 
% SE:               [k x v]     matrix of standard error values for k SNPs
%                               and v outcome variables
% 
% stats:            [k x v]     matrix of T/Z (or other) statistics values
%                               for k SNPs and v outcome variables
% 
% log10pValues:     [k x v]     matrix of -log10 p values for k SNPs and v
%                               outcome variables
% 
% df:               [1 x 1]     degrees of freedom for the analysis
% 
% outDir:           character   full path to where the results are written
% 
% fname:            [1 x v]     cell type having file names for v outcome
%                               variables (optional)
%
%% Outputs:
% GWAS summary files are written out in outDir; if fname is not specified,
% then the files are named 'GWAS_Summary_Phenotype_xx', where xx is the
% phenotype number (zero prefixed to four decimal places). The summary
% statistics csv files are tab separated

%% Specify foramt for writing files
fmt          = '%s \t %s \t %s \t %d \t %g \t %g \t %g \t %g \t %d \n';
GWASvarNames = {'Chromosome', 'SNPID', 'BasePair', 'SampleSize', 'Beta', 'SE', 'TStatistics', 'log10pValue', 'DF'};

%% Some initialization
numYvars = size(beta, 2);
numSNPs  = size(beta, 1);

%% Prepare file names, if not specified
if not(exist('fname', 'var')) || isempty(fname)
    fname = cellstr(strcat('GWAS_Summary_Phenotype_', num2str((1:numYvars)', '%04d'), '.csv'));
end

%% Loop over phenotypes and write out
for phen = 1:numYvars
    % Prepare table
    GWASTable           = cell(numSNPs, length(GWASvarNames));
    GWASTable(:, 1)     = Chr;
    GWASTable(:, 2)     = SNPID;
    GWASTable(:, 3)     = basePair;
    GWASTable(:, 4:9)   = num2cell([repmat(sampleSize, numSNPs, 1), beta(:,phen), SE(:,phen), stats(:,phen), log10pValues(:,phen), repmat(df, numSNPs, 1)]);
    GWASTable           = GWASTable';
    
    % Write out the formatted GWAS table
    fid = fopen(fullfile(outDir, fname{phen}), 'w');
    fprintf(fid, [repmat('%s \t ', 1, length(GWASvarNames)-1), ' %s \n'], GWASvarNames{:});
    fprintf(fid, fmt, GWASTable{:});
    fclose(fid);
end