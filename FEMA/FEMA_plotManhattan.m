function FEMA_plotManhattan(fileName, outDir, strTitle, fileFlag, GWASThresh, Chr, bpLoc, logpVals)
% Function that reads a FEMA output GWAS summary file and prepares a 
% Manhattan plot
%% Inputs:
% fileName:     full path to a GWAS summary file
%
% outDir:       full path to where the figure should be saved
%
% strTitle:     character indicating the title of the figure (if left
%               empty, no title will be put on top of the figure); this is
%               also used for deciding output name; if empty, the output
%               plot name is 'ManhattanPlot.png'
%
% fileFlag:     indicate which software was used to generate the GWAS 
%               summary file (not case sensitive):
%                   * 'FEMA'    (default)
%                   * 'GCTA'    (assumes fastGWA output)
%                   * 'Regenie'
%                   * 'PLINK'   (assumes PLINK 2)
% 
% GWSAThresh:   number indicating the GWAS significance threshold level;
%               assumes -log10 p value
%               (default: -log10(5 x 10^-8))
%
% The following three inputs can be used instead of a fileName to create a
% plot from some custom/pre-read GWAS summary file (ignored if fileName is 
% provided):
% Chr:          vector of chromosome numbers for m SNPs
% 
% bpLoc:        vector of base pair locations for m SNPs
%
% logpVals:     vector of -log10(p) values for m SNPs
%
%% TO-DO:
% Detect and handle sex chromosomes
% Remove fixed number of chromosomes
% Handle situation where one or more chromosome(s) may be missing
% Add space between chromosomes

%% Check inputs
if ~exist('fileName', 'var') || isempty(fileName)
    % Check if Chr, bpLoc, and pValues exist
    if (~exist('Chr',      'var') || isempty(Chr))    || ...
       (~exist('bpLoc',    'var') || isempty(bpLoc))  || ...
       (~exist('logpVals', 'var') || isempty(logpVals))
            warning('No relevant inputs provided;, nothing to do'); 
            return
    else
        % Ensure Chr, bpLoc, and pValues have same sizes
        numSNPs = size(Chr,1);
        if length(bpLoc) ~= numSNPs || length(logpVals) ~= numSNPs
            error('Mismatch in size of Chr, bpLoc, and pValue');
        else
            readFile = false;
        end
    end
else
    % Make sure file exists
    if ~exist(fileName, 'file')
        error(['Unable to find file: ', fileName]);
    else
        readFile = true;
        % Check fileFlag
        if ~exist('fileFlag', 'var') || isempty(fileFlag)
            fileFlag = 'fema';
        else
            fileFlag = lower(fileFlag);
            if ~ismember(fileFlag, {'fema', 'gcta', 'regenie', 'plink'})
                error('Unknown fileFlag provided');
            end
        end
    end
end

% Check outDir
if ~exist('outDir', 'var') || isempty(outDir)
    error('Please specify where the plot should be saved');
else
    if ~exist(outDir, 'dir')
        mkdir(outDir);
    end
end

% Check strTitle
if ~exist('strTitle', 'var') || isempty(strTitle)
    doTitle  = false;
    strTitle = 'ManhattanPlot';
else
    doTitle = true;
end

% % Check keepNaN
% if ~exist('keepNaN', 'var') || isempty(keepNaN)
%     keepNaN = true;
% else
%     if ~islogical(keepNaN)
%         error('keepNaN should be either true or false');
%     end
% end

% Check GWASThresh
if ~exist('GWASThresh', 'var') || isempty(GWASThresh)
    GWASThresh = -log10(5e-08);
end

%% Read file, if required
if readFile
    fileData = readtable(fileName, 'FileType', 'text');

    % Extract relevant information
    switch fileFlag
        case 'fema'
            Chr      = fileData.Chromosome;
            bpLoc    = fileData.BasePair;
            logpVals = fileData.log10pValue;
            
        case 'gcta'
            Chr      = fileData.CHR;
            bpLoc    = fileData.POS;
            logpVals = -log10(fileData.P);
            
        case 'regenie'
            Chr      = fileData.CHROM;
            bpLoc    = fileData.GENPOS;
            logpVals = fileData.LOG10P;
            
        case 'plink'
            % Restrict to rows which are additive tests
            fileData(~strcmpi(fileData.TEST, 'ADD'), :) = [];
            Chr      = fileData.x_CHROM;
            bpLoc    = fileData.POS;
            logpVals = -log10(fileData.P);
    end
    clear fileData
end

% %% Remove NaN values, if required
% if ~keepNaN
%     locNaN           = isnan(logpVals);
%     [a, b]           = histcounts(categorical(Chr(locNaN)));
%     Chr(locNaN)      = [];
%     bpLoc(locNaN)    = [];
%     logpVals(locNaN) = [];
%     disp(['Removing ', num2str(sum(locNaN)), ' NaN values']);
%     info = strcat('Removed', {' '}, cellstr(num2str(a')), ' SNPs from chromosome', {' '}, b');
%     fprintf('%s \n', info{:});
% end

%% Create figure
fig = figure('Units', 'centimeters', 'Position', [10 10 16 14]);
if doTitle
    tight_subplot(1, 1, 0, [0.045 0.04], [0.1 0.02]);
else
   tight_subplot(1, 1, 0, [0.045 0.01], [0.1 0.02]); 
end
ax       = gca;
offset   = 0;
numChr   = 22;
tickLoc  = zeros(numChr,1);
% spacer  = max(bpLoc(Chr == 1))/5;
hold(ax, 'on');

% Define colours - 3-class BuPu [colorbrewer2.org]
colOdd  = [136, 86,  167]./255;
colEven = [158, 188, 218]./255;

% Do plot
for chr = 1:numChr
    % Find locations to plot
    locs = Chr == chr;
    
    % Plot
    if logical(mod(chr,2))
        scatter(bpLoc(locs) + offset, logpVals(locs), 'Marker', '.', 'MarkerEdgeColor', colEven, 'SizeData', 45);
    else
        scatter(bpLoc(locs) + offset, logpVals(locs), 'Marker', '.', 'MarkerEdgeColor', colOdd,  'SizeData', 45);
    end
    
    % Figure out where to place the tick
    tickLoc(chr) = offset + ((offset + max(bpLoc(locs))) - (offset + min(bpLoc(locs))))/2;

    % Update offset to the last base pair
    offset = offset + max(bpLoc(locs));
end

% Add GWAS line
plot([0, offset - max(bpLoc(Chr == numChr))], [GWASThresh, GWASThresh], 'Color', 'k', 'LineStyle', ':', 'LineWidth', 1.5);

% Adjust X axis settings
xlim([0, offset - max(bpLoc(Chr == numChr))]);
xticks(tickLoc);
xticklabels(1:numChr);
xtickangle(90);
ax.XAxis.FontSize = 8;

% Adjust Y axis settings
ylim([min(logpVals) max(logpVals)+1]);
ax.YAxis.TickValues = sort([ax.YAxis.TickValues, round(GWASThresh, 2)]);
ax.YAxis.TickLabelsMode = 'auto';
ylabel('-log_{10} \itp', 'Interpreter', 'tex', 'FontSize', 12);

% Put a title, if user wants
if doTitle
    title(strTitle, 'FontSize', 12);
end

%% Print figure
print(fullfile(outDir, [strTitle, '.png']), '-dpng', '-r600');
close(fig);