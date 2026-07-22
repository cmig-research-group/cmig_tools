function [fig, pieceWise] = plotManhattan(Chr, bpLoc, logpVals, strTitle, markers, GWASThresh, style, filterThresh, toMarkPos)
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
%                               * 'stark'
%
% filterThresh: <numeric>   -log10 threshold that should be applied to
%                           filter out variants
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
    if ~ismember(style, {'mono', 'dark', 'diverge', 'dark-diff', 'diverge-diff', 'stark'})
        error('Unknown style type provided; should be one of: mono, dark, diverge, dark-diff, diverge-diff, stark');
    end
end

% Check filterThresh and get rid of some p values if user wants
if exist('filterThresh', 'var') && ~isempty(filterThresh)
    toDelete           = logpVals < filterThresh;
    logpVals(toDelete) = [];
    Chr(toDelete)      = [];
    bpLoc(toDelete)    = [];
    ovMin              = filter_thresh;
else
    ovMin = 0;
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

% Make sure these are logp values
if min(logpVals) >= 0 && max(logpVals) <= 1
    warning('Detected p values instead of log p values; converting to -log10(p) values');
    logpVals = -log10(logpVals);
end

%% Overall max
ovMax = max(logpVals, [], 'all');

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

    case 'stark'
        % Define colours: 6-class BrBG
        colOdd      = [140,81,10]./255;
        colEven     = [1,102,94]./255;
        colOdd_ns   = [246,232,195]./255;
        colEven_ns  = [199,234,229]./255;
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
sizeData = 30;

%% Do plot
% Some variables to hold axis elements
offset_chr = zeros(numChr,1);
ns_points  = cell(numChr, 1);
sig_points = cell(numChr, 1);

for chr = 1:numChr
    try
        % Find locations to plot - non significant ones
        locs = Chr == chr & logpVals <= GWASThresh;
    
        % Plot
        if logical(mod(chr,2))
            ns_points{chr} = scatter(allH(1), bpLoc(locs) + offset, logpVals(locs), 'Marker', '.', 'MarkerEdgeColor', colEven_ns, 'SizeData', sizeData);
        else
            ns_points{chr} = scatter(allH(1), bpLoc(locs) + offset, logpVals(locs), 'Marker', '.', 'MarkerEdgeColor', colOdd_ns,  'SizeData', sizeData);
        end
    catch
    end

    try
        % Find locations to plot - significant ones
        locs = Chr == chr & logpVals > GWASThresh;
    
        % Plot
        if logical(mod(chr,2))
            sig_points{chr} = scatter(allH(1), bpLoc(locs) + offset, logpVals(locs), 'Marker', '.', 'MarkerEdgeColor', colEven, 'SizeData', sizeData);
        else
            sig_points{chr} = scatter(allH(1), bpLoc(locs) + offset, logpVals(locs), 'Marker', '.', 'MarkerEdgeColor', colOdd,  'SizeData', sizeData);
        end
    catch
    end

    % Figure out where to place the tick
    locs = Chr == chr;
    tickLoc(chr) = offset + mean(bpLoc(locs));
    % tickLoc(chr) = offset + ((offset + max(BP(locs))) - (offset + min(BP(locs))))/2;

    % Update offset to the last base pair
    offset_chr(chr) = offset;
    offset          = offset + max(bpLoc(locs));
end

%% Add GWAS line
GWASline{1} = plot(allH(1), [0, offset + max(bpLoc(Chr == numChr))], [GWASThresh, GWASThresh], 'Color', col_thresh, 'LineStyle', ':', 'LineWidth', 2);

%% Adjust X axis settings
% xlim([0, offset + max(BP(Chr == numChr))]);
xlim(allH(1), [0, offset]);
xticks(tickLoc);
xticklabels(num2str((1:numChr)'));
xtickangle(90);
allH(1).XAxis.FontSize   = 10;
allH(1).XAxis.TickLength = [0 0];

%% Adjust Y axis
ylim(allH(1), [ovMin ovMax+1]);
allH(1).YAxis.TickValues = setdiff(sort([allH(1).YAxis.TickValues, round(GWASThresh, 2)]),  0);

% Mark y axes based on user choice
if exist('markers', 'var') && ~isempty(markers)
    allH(1).YAxis.TickValues = markers;
end
allH(1).YAxis.TickLabelsMode = 'auto';
allH(1).YAxis.FontSize       = 10;

% if ~exist('markers', 'var') || isempty(markers)
%     ylim(allH(1), [ovMin max(logpVals)+1]);
%     allH(1).YAxis.TickValues     = sort([allH(1).YAxis.TickValues, round(GWASThresh, 2)]);
%     allH(1).YAxis.TickLabelsMode = 'auto';
%     allH(1).YAxis.FontSize       = 10;
% else
%     ylim(allH(1), [min(markers) max(markers)]);
%     allH(1).YAxis.TickValues     = markers;
%     allH(1).YAxis.TickLabelsMode = 'auto';
%     allH(1).YAxis.FontSize       = 10;
% end

try
    allH(1).YAxis.TickLabels = cellstr(allH(1).YAxis.TickLabels);
catch
end

ylabel(allH(1), '-log_{10} \itp', 'Interpreter', 'tex', 'FontSize', 12);

%% Put a title, if user wants
if doTitle
    plotTitle = title(allH(1), strTitle, 'FontSize', 12);
end

%% Make some markers prominent
if exist('toMarkPos', 'var') && ~isempty(toMarkPos)
    allPoints_txt = cell(1,1);
    allPoints     = cell(1,1);
    for chr = 1:length(toMarkPos.Chr)
        if logical(mod(toMarkPos.Chr(chr), 2))
            allPoints{chr} = scatter(allH(1), toMarkPos.bpLoc(:,chr) + offset_chr(toMarkPos.Chr(chr)), toMarkPos.logpVals(:,chr), 'Marker', 'o', 'MarkerEdgeColor', colEven, 'SizeData', 50, 'LineWidth', 1.25);
        else
            allPoints{chr} = scatter(allH(1), toMarkPos.bpLoc(:,chr) + offset_chr(toMarkPos.Chr(chr)), toMarkPos.logpVals(:,chr), 'Marker', 'o', 'MarkerEdgeColor', colOdd,  'SizeData', 50, 'LineWidth', 1.25);
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
    pieceWise.allH          = allH;
    pieceWise.ns_points     = ns_points;
    pieceWise.sig_points    = sig_points;
    pieceWise.GWASline      = GWASline;
    pieceWise.offset_chr    = offset_chr;
    try
        pieceWise.allPoints     = allPoints;
        pieceWise.allPoints_txt = allPoints_txt;
    catch
    end
    try
        pieceWise.plotTitle = plotTitle;
    catch
    end
end