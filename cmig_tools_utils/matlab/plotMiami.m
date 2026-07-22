function [fig, pieceWise] = plotMiami(Chr, bpLoc, logpVals1, logpVals2, strTitle1, strTitle2, GWASThresh, style, filterThresh, tickMarks, toMarkPos)
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
%                               * 'mono2'
%                               * 'dark'
%                               * 'dark-diff' (default)
%                               * 'diverge'
%                               * 'diverge-diff'
%                               * 'diverge-diff2'
%                               * 'stark'
%
% filterThresh: <numeric>   -log10 threshold that should be applied to
%                           filter out variants
% 
% tickMarks:    <numeric>   vector of -log10 values to mark on y axes
%
% toMarkPos:    <struct>    structure type having the following fields;
%                           these SNPs are highlighted by a large marker
%                           around it; fields should include (only for the
%                           top panel): p positions, c chromosomes
%                               * 'Chr'         1 x c
%                               * 'bpLoc'       p x c
%                               * 'logpVals'    p x c
%                               * 'rsID'        p x c
%
%% TO-DO:
% Reintroduce the reading of summary statistics file
% Detect and handle sex chromosomes
% Remove fixed number of chromosomes
% Handle situation where one or more chromosome(s) may be missing
% Add space between chromosomes
% Allow an axis handle to be passed in
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
    if ~ismember(style, {'mono', 'mono2', 'dark', 'diverge', 'dark-diff', 'diverge-diff', 'diverge-diff2', 'stark'})
        error('Unknown style type provided; should be one of: mono, mono2, dark, diverge, dark-diff, diverge-diff, diverge-diff2, stark');
    end
end

% Check filterThresh and get rid of some p values if user wants
if exist('filterThresh', 'var') && ~isempty(filterThresh)
    toDelete            = logpVals1 < filterThresh & logpVals2 < filterThresh;
    logpVals1(toDelete) = [];
    logpVals2(toDelete) = [];
    Chr(toDelete)       = [];
    bpLoc(toDelete)     = [];
    ovMin               = filter_thresh;
else
    ovMin = 0;
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

% Make sure these are logp values
if min(logpVals1) >= 0 && max(logpVals1) <= 1
    warning('Detected p values instead of log p values for logpVals1; converting to -log10(p) values');
    logpVals1 = -log10(logpVals1);
end

if min(logpVals2) >= 0 && max(logpVals2) <= 1
    warning('Detected p values instead of log p values for logpVals2; converting to -log10(p) values');
    logpVals2 = -log10(logpVals2);
end


%% Overall max
ovMax = max([logpVals1; logpVals2], [], 'all');

%% Convert second set of p values to negative values
logpVals2 = -logpVals2;

%% Style configuration
switch(style)
    case 'mono'
        % Single colour: 3-class Dark2 [colorbrewer2.org]
        colOdd_ax1      = [117,112,179]./255;
        colEven_ax1     = [117,112,179]./255;
        colOdd_ns_ax1   = [117,112,179]./255;
        colEven_ns_ax1  = [117,112,179]./255;
        col_thresh      = [217,95,2]./255;

        % Same colours for axis 2
        colOdd_ax2      = [117,112,179]./255;
        colEven_ax2     = [117,112,179]./255;
        colOdd_ns_ax2   = [117,112,179]./255;
        colEven_ns_ax2  = [117,112,179]./255;

    case 'mono2'
        % Single colour: 3-class Dark2 [colorbrewer2.org]
        colOdd_ax1      = [117,112,179]./255;
        colEven_ax1     = [117,112,179]./255;
        colOdd_ns_ax1   = [117,112,179]./255;
        colEven_ns_ax1  = [117,112,179]./255;
        col_thresh      = [217,95,2]./255;

        % Single colour: 3-class Dark2 [colorbrewer2.org]
        colOdd_ax2      = [27,158,119]./255;
        colEven_ax2     = [27,158,119]./255;
        colOdd_ns_ax2   = [27,158,119]./255;
        colEven_ns_ax2  = [27,158,119]./255;

    case 'dark'
        % Two colours: 3-class Dark2 [colorbrewer2.org]
        colOdd_ax1      = [27,158,119]./255;
        colEven_ax1     = [117,112,179]./255;
        colOdd_ns_ax1   = [27,158,119]./255;
        colEven_ns_ax1  = [117,112,179]./255;
        col_thresh      = [217,95,2]./255;

        % Same colours for axis 2
        colOdd_ax2      = [27,158,119]./255;
        colEven_ax2     = [117,112,179]./255;
        colOdd_ns_ax2   = [27,158,119]./255;
        colEven_ns_ax2  = [117,112,179]./255;

    case 'diverge'
        % Define colours - 3-class BrBG [colorbrewer2.org]
        colOdd_ax1      = [216,179,101]./255;
        colEven_ax1     = [90,180,172]./255;
        colOdd_ns_ax1   = [216,179,101]./255;
        colEven_ns_ax1  = [90,180,172]./255;
        col_thresh      = [0 0 0];

        % Same colours for axis 2
        colOdd_ax2      = [216,179,101]./255;
        colEven_ax2     = [90,180,172]./255;
        colOdd_ns_ax2   = [216,179,101]./255;
        colEven_ns_ax2  = [90,180,172]./255;

    case 'dark-diff'
        % Define colours - 3-class Dark2 [colorbrewer2.org]
        colOdd_ax1  = [27,158,119]./255;
        colEven_ax1 = [117,112,179]./255;

        % Define colours - 3-class Set2 [colorbrewer2.org]
        colOdd_ns_ax1  = [102,194,165]./255;
        colEven_ns_ax1 = [141,160,203]./255;

        % Color threshold
        col_thresh = [217,95,2]./255;

        % Same colours for axis 2
        colOdd_ax2      = [27,158,119]./255;
        colEven_ax2     = [117,112,179]./255;
        colOdd_ns_ax2   = [102,194,165]./255;
        colEven_ns_ax2  = [141,160,203]./255;

    case 'diverge-diff'
        % Define colours: 4-class BrBG
        colOdd_ax1      = [166,97,26]./255;
        colEven_ax1     = [1,133,113]./255;
        colOdd_ns_ax1   = [223,194,125]./255;
        colEven_ns_ax1  = [128,205,193]./255;
        col_thresh      = [0 0 0];

        % Same colours for axis 2
        colOdd_ax2      = [166,97,26]./255;
        colEven_ax2     = [1,133,113]./255;
        colOdd_ns_ax2   = [223,194,125]./255;
        colEven_ns_ax2  = [128,205,193]./255;

    case 'diverge-diff2'
        % Define colours: 4-class BrBG
        colOdd_ax1      = [166,97,26]./255;
        colEven_ax1     = [1,133,113]./255;
        colOdd_ns_ax1   = [223,194,125]./255;
        colEven_ns_ax1  = [128,205,193]./255;
        col_thresh      = [0 0 0];

        % Define colours: 4-class PRGn
        colOdd_ax2      = [123,50,148]./255;
        colEven_ax2     = [0,136,55]./255;
        colOdd_ns_ax2   = [194,165,207]./255;
        colEven_ns_ax2  = [166,219,160]./255;

    case 'stark'
        % Define colours: 6-class BrBG
        colOdd_ax1      = [140,81,10]./255;
        colEven_ax1     = [1,102,94]./255;
        colOdd_ns_ax1   = [246,232,195]./255;
        colEven_ns_ax1  = [199,234,229]./255;
        col_thresh      = [0 0 0];

        % Same colours for axis 2
        colOdd_ax2      = [140,81,10]./255;
        colEven_ax2     = [1,102,94]./255;
        colOdd_ns_ax2   = [246,232,195]./255;
        colEven_ns_ax2  = [199,234,229]./255;
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
sizeData = 30;

%% Plot the first set of p values
% Some variables to hold axis elements
offset_chr     = zeros(numChr,1);
ax1_ns_points  = cell(numChr, 1);
ax1_sig_points = cell(numChr, 1);

for chr = 1:numChr
    try
        % Find locations to plot - non significant ones
        locs = Chr == chr & logpVals1 <= GWASThresh;
    
        % Plot
        if logical(mod(chr,2))
            ax1_ns_points{chr} = scatter(allH(1), bpLoc(locs) + offset, logpVals1(locs), 'Marker', '.', 'MarkerEdgeColor', colEven_ns_ax1, 'SizeData', sizeData);
        else
            ax1_ns_points{chr} = scatter(allH(1), bpLoc(locs) + offset, logpVals1(locs), 'Marker', '.', 'MarkerEdgeColor', colOdd_ns_ax1,  'SizeData', sizeData);
        end
    catch
    end

    % Find locations to plot - significant ones
    try
        locs = Chr == chr & logpVals1 > GWASThresh;
    
        % Plot
        if logical(mod(chr,2))
            ax1_sig_points{chr} = scatter(allH(1), bpLoc(locs) + offset, logpVals1(locs), 'Marker', '.', 'MarkerEdgeColor', colEven_ax1, 'SizeData', sizeData);
        else
            ax1_sig_points{chr} = scatter(allH(1), bpLoc(locs) + offset, logpVals1(locs), 'Marker', '.', 'MarkerEdgeColor', colOdd_ax1,  'SizeData', sizeData);
        end
    catch
    end

    % Figure out where to place the tick
    locs = Chr == chr;
    tickLoc(chr) = offset + mean(bpLoc(locs));

    % Update offset to the last base pair
    offset_chr(chr) = offset;
    offset          = offset + max(bpLoc(locs));
end

%% Plot second set of p values
offset         = 0;
ax2_ns_points  = cell(numChr, 1);
ax2_sig_points = cell(numChr, 1);

for chr = 1:numChr
    try
        % Find locations to plot - non significant ones
        locs = Chr == chr & abs(logpVals2) <= GWASThresh;
    
        % Plot
        if logical(mod(chr,2))
            ax2_ns_points{chr} = scatter(allH(2), bpLoc(locs) + offset, logpVals2(locs), 'Marker', '.', 'MarkerEdgeColor', colEven_ns_ax2, 'SizeData', sizeData);
        else
            ax2_ns_points{chr} = scatter(allH(2), bpLoc(locs) + offset, logpVals2(locs), 'Marker', '.', 'MarkerEdgeColor', colOdd_ns_ax2,  'SizeData', sizeData);
        end
    catch
    end

    try
        % Find locations to plot - significant ones
        locs = Chr == chr & abs(logpVals2) > GWASThresh;
    
        % Plot
        if logical(mod(chr,2))
            ax2_sig_points{chr} = scatter(allH(2), bpLoc(locs) + offset, logpVals2(locs), 'Marker', '.', 'MarkerEdgeColor', colEven_ax2, 'SizeData', sizeData);
        else
            ax2_sig_points{chr} = scatter(allH(2), bpLoc(locs) + offset, logpVals2(locs), 'Marker', '.', 'MarkerEdgeColor', colOdd_ax2,  'SizeData', sizeData);
        end
    catch
    end

    % Update offset to the last base pair
    locs = Chr == chr;
    offset = offset + max(bpLoc(locs));
end

%% Add GWAS line
GWASline{1} = plot(allH(1), [0, offset + max(bpLoc(Chr == numChr))], [GWASThresh,   GWASThresh], 'Color', col_thresh, 'LineStyle', ':', 'LineWidth', 2);
GWASline{2} = plot(allH(2), [0, offset + max(bpLoc(Chr == numChr))], [-GWASThresh, -GWASThresh], 'Color', col_thresh, 'LineStyle', ':', 'LineWidth', 2);

%% Adjust X axis settings
xlim(allH(1), [0, offset]);
xlim(allH(2), [0, offset]);

allH(1).XAxis.TickLength = [0 0];
allH(2).XAxis.Visible    = 'off';

% For labeling chromosomes - gets in the way of title2; either move title2
% to the bottom of the figure or do not show chromosome labels
% xticks(allH(1), tickLoc);
% xticklabels(allH(1), num2str((1:numChr)'));
% xtickangle(allH(1), 90);
% allH(1).XAxis.FontSize   = 10;

%% Adjust Y axis
ylim(allH(1), [ovMin ovMax+1]);
ylim(allH(2), [-ovMax-1 -ovMin]);

allH(1).YAxis.TickValues     = setdiff(sort([allH(1).YAxis.TickValues, round(GWASThresh, 2)]),  0);
allH(2).YAxis.TickValues     = setdiff(sort([allH(2).YAxis.TickValues, round(-GWASThresh, 2)]), 0);

% Mark y axes based on user choice
if exist('tickMarks', 'var') && ~isempty(tickMarks)
    allH(1).YAxis.TickValues = tickMarks;
    allH(2).YAxis.TickValues = sort(-1 .* tickMarks);
end

allH(1).YAxis.TickLabelsMode = 'auto';
allH(2).YAxis.TickLabelsMode = 'auto';
allH(2).YAxis.TickLabels     = abs(allH(2).YAxis.TickValues);
allH(1).YAxis.FontSize       = 10;
allH(2).YAxis.FontSize       = 10;

try
    allH(1).YAxis.TickLabels = cellstr(allH(1).YAxis.TickLabels);
    allH(2).YAxis.TickLabels = cellstr(allH(2).YAxis.TickLabels);
catch
end

ylabel(allH(1), '-log_{10} \itp', 'Interpreter', 'tex', 'FontSize', 12);
ylabel(allH(2), '-log_{10} \itp', 'Interpreter', 'tex', 'FontSize', 12);

%% Put a title, if user wants
if doTitle1
    plotTitle1 = title(allH(1), strTitle1, 'FontSize', 12);
end

if doTitle2
    plotTitle2 = title(allH(2), strTitle2, 'FontSize', 12);
end

%% Make some markers prominent
if exist('toMarkPos', 'var') && ~isempty(toMarkPos)
    allPoints_txt = cell(1,1);
    allPoints     = cell(1,1);
    for chr = 1:length(toMarkPos.Chr)
        if logical(mod(toMarkPos.Chr(chr), 2))
            allPoints{chr} = scatter(allH(1), toMarkPos.bpLoc(:,chr) + offset_chr(toMarkPos.Chr(chr)), toMarkPos.logpVals(:,chr), 'Marker', 'o', 'MarkerEdgeColor', colEven_ax1, 'SizeData', 50, 'LineWidth', 1.25);
        else
            allPoints{chr} = scatter(allH(1), toMarkPos.bpLoc(:,chr) + offset_chr(toMarkPos.Chr(chr)), toMarkPos.logpVals(:,chr), 'Marker', 'o', 'MarkerEdgeColor', colOdd_ax1,  'SizeData', 50, 'LineWidth', 1.25);
        end

        % Put a text next to it
        for tt = 1:length(toMarkPos.rsID(:,chr))
            allPoints_txt{chr,tt} = text(allH(1), toMarkPos.bpLoc(tt,chr) + offset_chr(toMarkPos.Chr(chr)), toMarkPos.logpVals(tt,chr), toMarkPos.rsID(tt,chr), 'FontSize', 6, 'VerticalAlignment', 'middle', 'HorizontalAlignment', 'left', 'FontWeight', 'bold', 'Color', 'k');
        end
    end
end

% Additionally, return parts of the plot for potential selective
% rasterization, if user wants
if nargout > 1
    pieceWise.allH           = allH;
    pieceWise.ax1_ns_points  = ax1_ns_points;
    pieceWise.ax2_ns_points  = ax2_ns_points;
    pieceWise.ax1_sig_points = ax1_sig_points;
    pieceWise.ax2_sig_points = ax2_sig_points;
    pieceWise.GWASline       = GWASline;
    pieceWise.offset_chr     = offset_chr;
    try
        pieceWise.allPoints     = allPoints;
        pieceWise.allPoints_txt = allPoints_txt;
    catch
    end
    try
        pieceWise.plotTitle1 = plotTitle1;
    catch
    end
    try
        pieceWise.plotTitle2 = plotTitle2;
    catch
    end
end