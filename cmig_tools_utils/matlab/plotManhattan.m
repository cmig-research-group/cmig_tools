function fig = plotManhattan(Chr, bpLoc, logpVals, strTitle, markers)
% Function that reads a GWAS summary file and prepares a Manhattan plot
%% Inputs:
% fileName:     full path to a GWAS summary file
%
% outDir:       full path to where the figure should be saved (if empty,
%               figure handle is returned and image is not saved)
%
% strTitle:     character indicating the title of the figure (if left
%               empty, no title will be put on top of the figure); this is
%               also used for deciding output name; if empty, the output
%               plot name is 'ManhattanPlot.png'; else 'Manhattan_strTitle'
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
% Check strTitle
if ~exist('strTitle', 'var') || isempty(strTitle)
    doTitle  = false;
else
    doTitle = true;
end

% Check GWASThresh
GWASThresh = -log10(5e-08);

%% Create figure
fig = figure('Units', 'centimeters', 'Position', [10 10 19 16]);
if doTitle
    tight_subplot(1, 1, 0, [0.045 0.04], [0.1 0.02]);
else
   tight_subplot(1, 1, 0, [0.045 0.01], [0.1 0.02]); 
end
ax       = gca;
offset   = 0;
numChr   = 22;
tickLoc  = zeros(numChr,1);
hold(ax, 'on');

% Define colours - 3-class BuPu [colorbrewer2.org]
colOdd  = [216, 179, 101]./255;
colEven = [90,  180, 172]./255;

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
plot([0, offset + max(bpLoc(Chr == numChr))], [GWASThresh, GWASThresh], 'Color', 'k', 'LineStyle', ':', 'LineWidth', 1.5);

% Adjust X axis settings
xlim([0, offset + max(bpLoc(Chr == numChr))]);
xticks(tickLoc);
xticklabels(1:numChr);
xtickangle(90);
ax.XAxis.FontSize = 8;
ax.XAxis.FontName = 'Courier';

% Adjust Y axis settings
if not(exist('markers', 'var'))
    ylim([min(logpVals) max(logpVals)+1]);
    ax.YAxis.TickValues = sort([ax.YAxis.TickValues, round(GWASThresh, 2)]);
    ax.YAxis.TickLabelsMode = 'auto';
else
    ylim([min(markers) max(markers)]);
    ax.YAxis.TickValues = markers;
    ax.YAxis.TickLabelsMode = 'auto';
end
ylabel('-log_{10} \itp', 'Interpreter', 'tex', 'FontSize', 12, 'FontName', 'Courier');
ax.YAxis.FontName = 'Courier';

% Put a title, if user wants
if doTitle
    title(strTitle, 'FontSize', 12, 'FontName', 'Courier');
end