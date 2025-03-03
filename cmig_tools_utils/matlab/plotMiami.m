function fig = plotMiami(Chr, bpLoc, logpVals1, logpVals2, strTitle1, strTitle2, GWASThresh, style)
% Function to create a Miami plot (aka mirrored Manhattan plot)
%% Inputs:
% Chr:          [g x 1]     numeric vector of values indicating chromosome
%                           number for every genetic variant g
%
% bpLoc:        [g x 1]     numeric vector of values indicating the base
%                           pair location for every genetic variant g
%
% logpVals1:    [g x 1]     numeric vector of values indicating the first
%                           set of -log10 p values for every genetic variant
%
% logpVals2:    [g x 1]     numeric vector of values indicating the second
%                           set of -log10 p values for every genetic variant
%
% strTitle1:    <character> text to use as a title for the top plot (optional)
%
% strTitle2:    <character> text to use as a title for the bottom plot (optional)
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
% Handle missing / NaN / Inf values
% Detect p values instead of -log10 p values
% Allow an axis handle to be passed in
% Allow for custom y axis ticks
% Add chromosome labels

%% Check inputs
% Check strTitle1
if ~exist('strTitle1', 'var') || isempty(strTitle1)
    doTitle1  = false;
else
    doTitle1 = true;
end

% Check strTitle2
if ~exist('strTitle2', 'var') || isempty(strTitle2)
    doTitle2  = false;
else
    doTitle2 = true;
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
locNaN  = isnan(logpVals1)      | isnan(logpVals2);
locInf  = isinf(logpVals1)      | isinf(logpVals2);
locImag = imag(logpVals1) ~= 0  | imag(logpVals2) ~= 0;
locDel  = locNaN | locInf | locImag;

if sum(locDel) > 0
    warning(['Removing ', num2str(sum(locNaN)),  ' NaN values, ',     ...
                          num2str(sum(locInf)),  ' Inf values, and ', ...
                          num2str(sum(locImag)), ' complex valued p values']);
    logpVals1(locDel) = [];
    logpVals2(locDel) = [];
    Chr(locDel)       = [];
    bpLoc(locDel)     = [];
end

% Ensure we are working with unsigned p values
logpVals1 = abs(logpVals1);
logpVals2 = abs(logpVals2);

%% Overall max
ovMax = max([logpVals1; logpVals2], [], 'all');

%% Convert second set of p values to negative values
logpVals2 = -logpVals2;

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
fig = figure('Units', 'centimeters', 'Position', [10 10 18 12]);
if doTitle1 || doTitle2
    allH = tight_subplot(2, 1, [0.06 0.00], [0.01 0.04], [0.1 0.01]);
else
    allH = tight_subplot(2, 1, [0.02 0.00], [0.01 0.01], [0.1 0.01]);
end
hold(allH(:), 'on');

%% Some settings
offset   = 0;
numChr   = 22;
tickLoc  = zeros(numChr,1);
sizeData = 60;

%% Plot the first set of p values
for chr = 1:numChr
    % Find locations to plot - non significant ones
    locs = Chr == chr & logpVals1 <= GWASThresh;

    % Plot
    if logical(mod(chr,2))
        scatter(allH(1), bpLoc(locs) + offset, logpVals1(locs), 'Marker', '.', 'MarkerEdgeColor', colEven_ns, 'SizeData', sizeData);
    else
        scatter(allH(1), bpLoc(locs) + offset, logpVals1(locs), 'Marker', '.', 'MarkerEdgeColor', colOdd_ns,  'SizeData', sizeData);
    end

    % Find locations to plot - non significant ones
    locs = Chr == chr & logpVals1 > GWASThresh;

    % Plot
    if logical(mod(chr,2))
        scatter(allH(1), bpLoc(locs) + offset, logpVals1(locs), 'Marker', '.', 'MarkerEdgeColor', colEven, 'SizeData', sizeData);
    else
        scatter(allH(1), bpLoc(locs) + offset, logpVals1(locs), 'Marker', '.', 'MarkerEdgeColor', colOdd,  'SizeData', sizeData);
    end

    % Figure out where to place the tick
    locs = Chr == chr;
    tickLoc(chr) = offset + mean(bpLoc(locs));

    % Update offset to the last base pair
    offset = offset + max(bpLoc(locs));
end

%% Plot second set of p values
offset  = 0;
for chr = 1:numChr
    % Find locations to plot - non significant ones
    locs = Chr == chr & abs(logpVals2) <= GWASThresh;

    % Plot
    if logical(mod(chr,2))
        scatter(allH(2), bpLoc(locs) + offset, logpVals2(locs), 'Marker', '.', 'MarkerEdgeColor', colEven_ns, 'SizeData', sizeData);
    else
        scatter(allH(2), bpLoc(locs) + offset, logpVals2(locs), 'Marker', '.', 'MarkerEdgeColor', colOdd_ns,  'SizeData', sizeData);
    end

    % Find locations to plot - non significant ones
    locs = Chr == chr & abs(logpVals2) > GWASThresh;

    % Plot
    if logical(mod(chr,2))
        scatter(allH(2), bpLoc(locs) + offset, logpVals2(locs), 'Marker', '.', 'MarkerEdgeColor', colEven, 'SizeData', sizeData);
    else
        scatter(allH(2), bpLoc(locs) + offset, logpVals2(locs), 'Marker', '.', 'MarkerEdgeColor', colOdd,  'SizeData', sizeData);
    end

    % Update offset to the last base pair
    locs = Chr == chr;
    offset = offset + max(bpLoc(locs));
end

%% Add GWAS line
plot(allH(1), [0, offset + max(bpLoc(Chr == numChr))], [GWASThresh,   GWASThresh], 'Color', col_thresh, 'LineStyle', ':', 'LineWidth', 2);
plot(allH(2), [0, offset + max(bpLoc(Chr == numChr))], [-GWASThresh, -GWASThresh], 'Color', col_thresh, 'LineStyle', ':', 'LineWidth', 2);

%% Adjust X axis settings
xlim(allH(1), [0, offset]);
xlim(allH(2), [0, offset]);

allH(1).XAxis.TickLength = [0 0];

allH(2).XAxis.Visible = 'off';

% For labeling chromosomes - gets in the way of title2; either move title2
% to the bottom of the figure or do not show chromosome labels
% xticks(allH(1), tickLoc);
% xticklabels(allH(1), num2str((1:numChr)'));
% xtickangle(allH(1), 90);
% allH(1).XAxis.FontSize   = 10;

%% Adjust Y axis
ylim(allH(1), [0 ovMax+1]);
ylim(allH(2), [-ovMax-1 0]);

allH(1).YAxis.TickValues     = setdiff(sort([allH(1).YAxis.TickValues, round(GWASThresh, 2)]), 0);
allH(2).YAxis.TickValues     = setdiff(sort([allH(2).YAxis.TickValues, round(-GWASThresh, 2)]), 0);

allH(1).YAxis.TickLabelsMode = 'auto';
allH(2).YAxis.TickLabelsMode = 'auto';
allH(2).YAxis.TickLabels     = abs(allH(2).YAxis.TickValues);

allH(1).YAxis.FontSize       = 10;
allH(2).YAxis.FontSize       = 10;

ylabel(allH(1), '-log_{10} \itp', 'Interpreter', 'tex', 'FontSize', 12);
ylabel(allH(2), '-log_{10} \itp', 'Interpreter', 'tex', 'FontSize', 12);

%% Put a title, if user wants
if doTitle1
    title(allH(1), strTitle1, 'FontSize', 12);
end

if doTitle2
    title(allH(2), strTitle2, 'FontSize', 12);
end