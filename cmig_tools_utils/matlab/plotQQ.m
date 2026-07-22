function [fig, h_ax, ll] = plotQQ(h_ax, pVals, strTitle, colOver, axLabel)
% Function that reads a GWAS summary file and creates a QQ plot
%% Inputs:
% h_ax:         handle to an existing figure axis that should be used for
%               plotting (leave empty for creating a new figure); use this 
%               option for creating plots that should be added to the same
%               figure 
%
% pVals:        vector of p values (not -log10 values)
%
% strTitle:     character indicating the title of the figure (if left
%               empty, no title will be put on top of the figure)
%
% colOver:      color to be used for plotting
%
% axLabel:      true or false indicating if x and y axes should be labeled
%               (default: true)
%
%% Outputs:
% fig:          handle to the figure where the plot exists
%
% h_ax:         handle to the axis within the figure
%
% ll:           handle to the QQ plot actual line
%
%% TO-DO:
% Reintroduce the reading of summary statistics file
% Add binomial confidence interval

%% Check inputs
% Check h_ax
if ~exist('h_ax', 'var') || isempty(h_ax)
    createFig = true;
else
    if ishandle(h_ax)
        createFig = false;
    else
        createFig = true;
        warning('Invalid axis handle provided; ignoring h_ax as input');
    end
end

% Check p values
if ~exist('pVals', 'var') || isempty(pVals)
    error('Please provide a vector of p values');
else
    if size(pVals,2) > 1
        if size(pVals,1) == 1
            pVals = pVals';
        else
            error('Multiple sets of plotting not yet supported');
        end
    else
        % Delete any NaN, Inf, or complex values
        locNaN  = isnan(pVals);
        locInf  = isinf(pVals);
        locImag = imag(pVals) ~= 0;
        locDel  = locNaN | locInf | locImag;
        if sum(locDel) > 0
            warning(['Removing ', num2str(sum(locNaN)),  ' NaN values, ',     ...
                                  num2str(sum(locInf)),  ' Inf values, and ', ...
                                  num2str(sum(locImag)), ' complex valued p values']);
            pVals(locDel) = [];
        end
        if max(abs(pVals)) > 1
            warning('Detected -log10 p values; converting to p values');
            pVals = 10.^-pVals;
        else
            % Ensure we are working with unsigned p values
            pVals = abs(pVals);
        end
    end
end

% Check strTitle
if ~exist('strTitle', 'var') || isempty(strTitle)
    doTitle  = false;
else
    doTitle = true;
end

% Check colOver
% Define colour - 3-class Dark2 [colorbrewer2.org]
if ~exist('colOver', 'var') || isempty(colOver)
    colTrue  = [217, 95, 2]./255;
else
    colTrue = colOver;
end

% Check axLabel
if ~exist('axLabel', 'var') || isempty(axLabel)
    axLabel = true;
else
    if isnumeric(axLabel)
        axLabel = logical(axLabel);
    else
        if ~islogical(axLabel)
            error('axLabel should be logical');
        end
    end
end

%% Create figure
if createFig
    fig = figure('Units', 'centimeters', 'Position', [10 10 10 10]);
    if doTitle
        tight_subplot(1, 1, 0, [0.11 0.06], [0.1 0.02]);
    else
        tight_subplot(1, 1, 0, [0.11 0.01], [0.1 0.02]);
    end
    h_ax = gca;
else
    try
        fig = h_ax.Parent;
        hold(h_ax, 'on');
    catch
        % Maybe this is a figure handle?
        fig  = h_ax;
        h_ax = gca;
    end
end

hold(h_ax, 'on');

% Create sorted data
% Sort p values from lowest to highest
p_observed = sort(pVals);

% Largest -log10 p value
maxVal = max(-log10(p_observed));

% Create expected p values
p_expected = (1:length(p_observed)) / length(p_observed);

% Plot p values
ll = plot(h_ax, -log10(p_expected), -log10(p_observed), '.', 'Color', colTrue);

% Plot null
plot(h_ax, [0, maxVal], [0, maxVal], 'Color', [0 0 0], 'LineStyle', '--', 'LineWidth', 1);

% Round up maxVal
maxVal = ceil(maxVal);

% Customize axes limits
xlim(h_ax, [0  maxVal+0.1]);
ylim(h_ax, [0  maxVal+0.1]);

% Place ticks and tick labels
% This is a bit arbitrary
if maxVal > 18
    xticks(h_ax,      0:4:maxVal);
    yticks(h_ax,      0:4:maxVal);
    xticklabels(h_ax, 0:4:maxVal);
    yticklabels(h_ax, 0:4:maxVal);
else    
    if maxVal > 6
        xticks(h_ax,      0:2:maxVal);
        yticks(h_ax,      0:2:maxVal);
        xticklabels(h_ax, 0:2:maxVal);
        yticklabels(h_ax, 0:2:maxVal);
    else
        xticks(h_ax,      0:1:maxVal);
        yticks(h_ax,      0:1:maxVal);
        xticklabels(h_ax, 0:1:maxVal);
        yticklabels(h_ax, 0:1:maxVal);
    end
end

% Change font size for the tick labels
h_ax.FontSize = 8;

% Label the axes
if axLabel
    xlabel(h_ax, 'Expected -log_{10} \itp', 'Interpreter', 'tex', 'FontSize', 10);
    ylabel(h_ax, 'Observed -log_{10} \itp', 'Interpreter', 'tex', 'FontSize', 10);
end

% Put a title, if user wants
if doTitle
    title(h_ax, strTitle, 'FontSize', 11);
end