function fig = FEMA_plotQQ(fileName, outDir, strTitle, fileFlag, h_ax, pVals)
% Function that reads a GWAS summary file and creates a QQ plot
%% Inputs:
% fileName:     full path to a GWAS summary file
%
% outDir:       full path to where the figure should be saved (if empty, 
%               figure handle is returned and image is not saved)
%
% strTitle:     character indicating the title of the figure (if left
%               empty, no title will be put on top of the figure); this is
%               also used for deciding output name; if empty, the output
%               plot name is 'QQPlot.png'; else 'QQ_strTitle'
%
% fileFlag:     indicate which software was used to generate the GWAS 
%               summary file (not case sensitive):
%                   * 'FEMA'    (default)
%                   * 'GCTA'    (assumes fastGWA output)
%                   * 'Regenie'
%                   * 'PLINK'   (assumes PLINK 2)
%
% h_ax:         handle to a figure axis that should be used for plotting
%               (optional input; use this for creating plots that should 
%               be added to the same figure)
% 
% The following input can be used instead of a fileName to create a
% plot from some custom/pre-read GWAS summary file (ignored if fileName is 
% provided):
% pVals:        vector of p values (not -log10 values)
%
%% Check inputs
if ~exist('fileName', 'var') || isempty(fileName)
    % Check if logpVals exist
    if (~exist('pVals', 'var') || isempty(pVals))
        warning('No relevant inputs provided;, nothing to do');
        return
    else
        % Sanity check p values
        minP = min(pVals);
        maxP = max(pVals);
        if minP < 0 || maxP > 1
            error('Possibly incorrect p values; please check data');
        end
        readFile = false;
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
    writeFile = false;
else
    writeFile = true;
    if ~exist(outDir, 'dir')
        mkdir(outDir);
    end
end

% Check strTitle
if ~exist('strTitle', 'var') || isempty(strTitle)
    doTitle  = false;
else
    doTitle = true;
end

% Output name, if required
if writeFile
    if doTitle
        outName = ['QQ_', strrep(strTitle, ' ', ''), '.png'];
    else
        outName = 'QQPlot.png';
    end
end

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

%% Read file, if required
if readFile
    fileData = readtable(fileName, 'FileType', 'text');

    % Extract relevant information
    switch fileFlag
        case 'fema'
            log10pVals  = fileData.log10pValue;
            pVals       = 10.^-log10pVals;
            
        case 'gcta'
            log10pVals  = -log10(fileData.P);
            pVals       = 10.^-log10pVals;
            
        case 'regenie'
            log10pVals  = fileData.LOG10P;
            pVals       = 10.^-log10pVals;
            
        case 'plink'
            % Restrict to rows which are additive tests
            fileData(~strcmpi(fileData.TEST, 'ADD'), :) = [];
            pVals       = fileData.P;
    end
    clear fileData
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
    fig = h_ax.Parent;
end

hold(h_ax, 'on');

% Define colour - 3-class Dark2 [colorbrewer2.org]
colTrue  = [217, 95, 2]./255;

% Create sorted data
% Sort p values from lowest to highest
p_observed = sort(pVals);

% Largest -log10 p value
maxVal = max(-log10(p_observed));

% Create expected p values
p_expected = (1:length(p_observed)) / length(p_observed);

% Plot p values
plot(h_ax, -log10(p_expected), -log10(p_observed), '.', 'Color', colTrue);

% Plot null
plot(h_ax, [0, maxVal], [0, maxVal], 'Color', [0 0 0], 'LineStyle', '--', 'LineWidth', 1);

% Round up maxVal
maxVal = ceil(maxVal);

% Customize axes limits
xlim(h_ax, [0  maxVal+0.1]);
ylim(h_ax, [0  maxVal+0.1]);

% Place ticks and tick labels
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

% Change font size for the tick labels
h_ax.FontSize = 8;

% Label the axes
xlabel(h_ax, 'Expected -log_{10} \itp', 'Interpreter', 'tex', 'FontSize', 10);
ylabel(h_ax, 'Observed -log_{10} \itp', 'Interpreter', 'tex', 'FontSize', 10);

% Put a title, if user wants
if doTitle
    title(h_ax, strTitle, 'FontSize', 11);
end

%% Print figure, if required
if writeFile
    print(fullfile(outDir, outName), '-dpng', '-r600');
    if createFig
        close(fig);
    end
end