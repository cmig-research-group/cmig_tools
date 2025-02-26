function FEMA_gatherGWAS(inDir, lookupName, cleanUp, outDir, outName)
% Function to compile split GWAS measures into one
%% Inputs:
% inDir:        full path to a directory having part-wise GWAS results
%
% lookupName:   common part of the file names that need to be looked up
%               (is looked up using regexpi)
%
% cleanUp:      true or false indicating if part-wise files should be
%               deleted afterwards (default: false)
%
% outDir:       full path to a directory where results should be saved
%               (else results are written to inDir)
%
% outName:      if specified, results are saved in outDir using this name;
%               else results are saved in outDir with the name
%               FEMA_GWAS_Aggregate.mat
%
%% Initialize timer
tInit = tic;

%% Check inputs
if ~exist('inDir', 'var') || isempty(inDir)
    error('Please specify input directory');
else
    if ~exist(inDir, 'dir')
        error(['Unable to find directory: ', inDir]);
    end
end

if ~exist('lookupName', 'var') || isempty(lookupName)
    lookupName = 'FEMA_GWAS_Ch';
end

if ~exist('cleanUp', 'var') || isempty(cleanUp)
    cleanUp = false;
else
    if ~islogical(cleanUp)
        error('cleanUp should be either true or false');
    end
end

if ~exist('outDir', 'var') || isempty(outDir)
    outDir = inDir;
else
    if ~exist(outDir, 'dir')
        mkdir(outDir);
    end
end

if ~exist('outName', 'var') || isempty(outName)
    outName = 'FEMA_GWAS_Aggregate';
end

%% Search for files
ls = dir(fullfile(inDir, '*.mat'));
ls = {ls(:).name}';
ls = ls(not(cellfun(@isempty, regexpi(ls, lookupName))));
disp(['Found ', num2str(length(ls)), ' files']);

% Convert these to full files
ls = fullfile(inDir, ls);

%% Get some info for initialization
info1 = load(ls{1},   'genStruct');
infoL = load(ls{end}, 'genStruct', 'beta_hat', 'Wald_F', 'coeffCovar');

% Number of SNPs
numSNPs = length(info1.genStruct.Chr) * (length(ls)-1) + length(infoL.genStruct.Chr);

% Number of coefficients and number of phenotypes
if ndims(infoL.beta_hat) == 3
    is3D                   = true;
    [~, numCoeff, numPhen] = size(infoL.beta_hat);
else
    is3D         = false;
    numCoeff     = 1;
    [~, numPhen] = size(infoL.beta_hat);
end

% Number of contrasts
if isempty(infoL.Wald_F)
    doWald = false;
else
    doWald = true;
    numCon = size(infoL.Wald_F, 3);
end

%% Initialize
if is3D
    [beta_hat, beta_se, tStats, logpValues] = deal(zeros(numSNPs, numCoeff, numPhen));
    if doWald
        [Wald_F, Wald_p] = deal(zeros(numSNPs, numPhen, numCon));
    end
else
    [beta_hat, beta_se, tStats, logpValues] = deal(zeros(numSNPs, numPhen));
    Wald_F = [];
    Wald_p = [];
end

% Also initialize chromosome, base pair location, and rsID
[Chr, BP, rsID] = deal(cell(numSNPs, 1));

% Do we need coefficient covariance?
if isempty(infoL.coeffCovar)
    coeffCovar = [];
else
    coeffCovar = zeros(numSNPs, numCoeff, numCoeff, numPhen);
end

% Error flag
errFlag = false(numSNPs, numPhen);

%% Loop over files, load, and aggregate
init = 1;
if is3D
    if doWald
        for files = 1:length(ls)
            % Load data
            tmp = load(ls{files}, 'beta_hat', 'beta_se', 'tStats', 'logpValues', 'Wald_F', 'Wald_p', 'genStruct', 'coeffCovar', 'errFlag');

            % How many rows?
            chk = size(tmp.beta_hat, 1);

            % Assign
            beta_hat(init:init+chk-1,   :, :)    = tmp.beta_hat;
            beta_se(init:init+chk-1,    :, :)    = tmp.beta_se;
            tStats(init:init+chk-1,     :, :)    = tmp.tStats;
            logpValues(init:init+chk-1, :, :)    = tmp.logpValues;
            Wald_F(init:init+chk-1,     :, :)    = tmp.Wald_F;
            Wald_p(init:init+chk-1,     :, :)    = tmp.Wald_p;
            Chr(init:init+chk-1)                 = tmp.genStruct.Chr;
            BP(init:init+chk-1)                  = tmp.genStruct.BP;
            rsID(init:init+chk-1)                = tmp.genStruct.SNPs;
            coeffCovar(init:init+chk-1, :, :, :) = tmp.coeffCovar;
            errFlag(init:init+chk-1, :)          = tmp.errFlag;

            % Update init
            init = init + chk;
        end
    else
        for files = 1:length(ls)
            % Load data
            tmp = load(ls{files}, 'beta_hat', 'beta_se', 'tStats', 'logpValues', 'genStruct', 'errFlag');

            % How many rows?
            chk = size(tmp.beta_hat, 1);

            % Assign
            beta_hat(init:init+chk-1,   :, :) = tmp.beta_hat;
            beta_se(init:init+chk-1,    :, :) = tmp.beta_se;
            tStats(init:init+chk-1,     :, :) = tmp.tStats;
            logpValues(init:init+chk-1, :, :) = tmp.logpValues;
            Chr(init:init+chk-1)              = tmp.genStruct.Chr;
            BP(init:init+chk-1)               = tmp.genStruct.BP;
            rsID(init:init+chk-1)             = tmp.genStruct.SNPs;
            errFlag(init:init+chk-1, :)       = tmp.errFlag;

            % Update init
            init = init + chk;
        end
    end
else
    for files = 1:length(ls)
        % Load data
        tmp = load(ls{files}, 'beta_hat', 'beta_se', 'tStats', 'logpValues', 'genStruct', 'errFlag');

        % How many rows?
        chk = size(tmp.beta_hat, 1);

        % Assign
        beta_hat(init:init+chk-1,   :) = tmp.beta_hat;
        beta_se(init:init+chk-1,    :) = tmp.beta_se;
        tStats(init:init+chk-1,     :) = tmp.tStats;
        logpValues(init:init+chk-1, :) = tmp.logpValues;
        Chr(init:init+chk-1)           = tmp.genStruct.Chr;
        BP(init:init+chk-1)            = tmp.genStruct.BP;
        rsID(init:init+chk-1)          = tmp.genStruct.SNPs;
        errFlag(init:init+chk-1, :)    = tmp.errFlag;

        % Update init
        init = init + chk;
    end
end

%% Convert Chr and BP to numbers
Chr = cellfun(@str2double, Chr);
BP  = cellfun(@str2double, BP);

%% Save
% Save name - make sure .mat is not already part of outName
saveName = fullfile(outDir, [strrep(outName, '.mat', ''), '.mat']);

% Variables to save
toSave = {'beta_hat', 'beta_se', 'tStats', 'logpValues', 'Wald_F', 'Wald_p', 'Chr', 'BP', 'rsID', 'errFlag', 'coeffCovar'};

% Get an estimate of size of the variables
tmpInfo = whos;
if sum([tmpInfo(ismember({tmpInfo(:).name}', toSave)).bytes]) > 2^31
    save(saveName, 'beta_hat', 'beta_se', 'tStats', 'logpValues', 'Wald_F', 'Wald_p', 'Chr', 'BP', 'rsID', 'errFlag', 'coeffCovar', '-v7.3');
else
    save(saveName, 'beta_hat', 'beta_se', 'tStats', 'logpValues', 'Wald_F', 'Wald_p', 'Chr', 'BP', 'rsID', 'errFlag', 'coeffCovar');
end

%% Clean up, if required
if cleanUp
    for files = 1:length(ls)
        delete(ls{files});
    end
end

%% End timer
disp(['Finished aggregating statistics in:, ', num2str(toc(tInit)), 's']);