function fig = plotManhattan(Chr, bpLoc, logpVals, strTitle, markers, GWASThresh, style)
% Function to create Manhattan plot
%% Inputs:
% Chr:          [g x 1]     numeric vector of values indicating chromosome
%                           number for every genetic variant g
%
% bpLoc:        [g x 1]     numeric vector of values indicating the base
%                           pair location for every genetic variant g
%
% logpVals:     [g x 1]     numeric vector of values indicating the -log10
%                           p values for every genetic variant g
%
% strTitle:     <character> text to use as a title to the plot (optional)
%
% markers:      <numeric>   vector of values to use as tick values for y
%                           axis (optional)
%
% GWSAThresh:   <numeric>   number indicating the GWAS significance
%                           threshold level; assumes -log10 p value; 
%                           (default: -log10(5 x 10^-8))
%
% style:        <character> should be one of the following (see Notes):
%                               * 'mono'
%                               * 'dark'
%                               * 'dark-diff' (default)
%                               * 'diverge'
%                               * 'diverge-diff'
%
%% TO-DO:
% Reintroduce the reading of summary statistics file
% Detect and handle sex chromosomes
% Remove fixed number of chromosomes
% Handle situation where one or more chromosome(s) may be missing
% Add space between chromosomes
% Detect p values instead of -log10 p values
% Allow an axis handle to be passed in

%% Check inputs
% Check strTitle
if ~exist('strTitle', 'var') || isempty(strTitle)
    doTitle  = false;
else
    doTitle = true;
end

% Check GWASThresh
if ~exist('GWASThresh', 'var') || isempty(GWASThresh)
    GWASThresh = -log10(5e-08);
end

% Check style
if ~exist('style', 'var') || isempty(style)
    style = 'dark-diff';
else
    style = lower(style);
    if ~ismember(style, {'mono', 'dark', 'diverge', 'dark-diff', 'diverge-diff'})
        error('Unknown style type provided; should be one of: mono, dark, diverge, dark-diff, diverge-diff');
    end
end

%% Some sanity checks on p values
% Delete any NaN, Inf, or complex values
locNaN  = isnan(logpVals);
locInf  = isinf(logpVals);
locImag = imag(logpVals) ~= 0;
locDel  = locNaN | locInf | locImag;

if sum(locDel) > 0
    warning(['Removing ', num2str(sum(locNaN)),  ' NaN values, ',     ...
                          num2str(sum(locInf)),  ' Inf values, and ', ...
                          num2str(sum(locImag)), ' complex valued p values']);
    logpVals(locDel) = [];
    Chr(locDel)      = [];
    bpLoc(locDel)    = [];
end

% Ensure we are working with unsigned p values
logpVals = abs(logpVals);

%% Style configuration
switch(style)
    case 'mono'
        % Single colour: 3-class Dark2 [colorbrewer2.org]
        colOdd      = [117,112,179]./255;
        colEven     = [117,112,179]./255;
        colOdd_ns   = [117,112,179]./255;
        colEven_ns  = [117,112,179]./255;        
        col_thresh  = [217,95,2]./255;

    case 'dark'
        % Two colours: 3-class Dark2 [colorbrewer2.org]
        colOdd      = [27,158,119]./255;
        colEven     = [117,112,179]./255;
        colOdd_ns   = [27,158,119]./255;
        colEven_ns  = [117,112,179]./255;
        col_thresh  = [217,95,2]./255;

    case 'diverge'
        % Define colours - 3-class BrBG [colorbrewer2.org]
        colOdd      = [216,179,101]./255;
        colEven     = [90,180,172]./255;        
        colOdd_ns   = [216,179,101]./255;
        colEven_ns  = [90,180,172]./255;
        col_thresh  = [0 0 0];

    case 'dark-diff'
        % Define colours - 3-class Dark2 [colorbrewer2.org]
        colOdd  = [27,158,119]./255;
        colEven = [117,112,179]./255;
        
        % Define colours - 3-class Set2 [colorbrewer2.org]
        colOdd_ns  = [102,194,165]./255;
        colEven_ns = [141,160,203]./255;

        % Color threshold
        col_thresh = [217,95,2]./255;

    case 'diverge-diff'
        % Define colours: 4-class BrBG
        colOdd      = [166,97,26]./255;
        colEven     = [1,133,113]./255;        
        colOdd_ns   = [223,194,125]./255;
        colEven_ns  = [128,205,193]./255;
        col_thresh  = [0 0 0];
end

%% Create figure
fig  = figure('Units', 'centimeters', 'Position', [10 10 18 12]);
if doTitle
    allH = tight_subplot(1, 1, 0, [0.05 0.04], [0.08 0.01]);
else
   allH = tight_subplot(1, 1, 0, [0.05 0.0], [0.08 0.01]); 
end
hold(allH(1), 'on');

%% Some settings
offset   = 0;
numChr   = 22;
tickLoc  = zeros(numChr,1);
sizeData = 60;

%% Do plot
for chr = 1:numChr
    % Find locations to plot - non significant ones
    locs = Chr == chr & logpVals <= GWASThresh;

    % Plot
    if logical(mod(chr,2))
        scatter(allH(1), bpLoc(locs) + offset, logpVals(locs), 'Marker', '.', 'MarkerEdgeColor', colEven_ns, 'SizeData', sizeData);
    else
        scatter(allH(1), bpLoc(locs) + offset, logpVals(locs), 'Marker', '.', 'MarkerEdgeColor', colOdd_ns,  'SizeData', sizeData);
    end

    % Find locations to plot - non significant ones
    locs = Chr == chr & logpVals > GWASThresh;

    % Plot
    if logical(mod(chr,2))
        scatter(allH(1), bpLoc(locs) + offset, logpVals(locs), 'Marker', '.', 'MarkerEdgeColor', colEven, 'SizeData', sizeData);
    else
        scatter(allH(1), bpLoc(locs) + offset, logpVals(locs), 'Marker', '.', 'MarkerEdgeColor', colOdd,  'SizeData', sizeData);
    end

    % Figure out where to place the tick
    locs = Chr == chr;
    tickLoc(chr) = offset + mean(bpLoc(locs));
    % tickLoc(chr) = offset + ((offset + max(BP(locs))) - (offset + min(BP(locs))))/2;

    % Update offset to the last base pair
    offset = offset + max(bpLoc(locs));
end

%% Add GWAS line
plot(allH(1), [0, offset + max(bpLoc(Chr == numChr))], [GWASThresh, GWASThresh], 'Color', col_thresh, 'LineStyle', ':', 'LineWidth', 2);

%% Adjust X axis settings
% xlim([0, offset + max(BP(Chr == numChr))]);
xlim(allH(1), [0, offset]);
xticks(tickLoc);
xticklabels(num2str((1:numChr)'));
xtickangle(90);
allH(1).XAxis.FontSize   = 10;
allH(1).XAxis.TickLength = [0 0];

%% Adjust Y axis
if ~exist('markers', 'var') || isempty(markers)
    ylim(allH(1), [min(logpVals) max(logpVals)+1]);
    allH(1).YAxis.TickValues     = sort([allH(1).YAxis.TickValues, round(GWASThresh, 2)]);
    allH(1).YAxis.TickLabelsMode = 'auto';
    allH(1).YAxis.FontSize       = 10;
else
    ylim(allH(1), [min(markers) max(markers)]);
    allH(1).YAxis.TickValues     = markers;
    allH(1).YAxis.TickLabelsMode = 'auto';
    allH(1).YAxis.FontSize       = 10;
end

ylabel(allH(1), '-log_{10} \itp', 'Interpreter', 'tex', 'FontSize', 12);

%% Put a title, if user wants
if doTitle
    title(allH(1), strTitle, 'FontSize', 12);
end