function fh = showSurfPlot(vertvals, str, cmap, polarity, climMin, climMax, ico, varargin)
%	 Function to plot vertexwise results of a FEMA analysis.
%	 Inputs:
%	 vertvals: vertexwise statistics for plotting
%	
%	 str: string for colorbar label
%	
%	 cmap: colormap for plotting
%	
%	 polarity: polarity of the statistic 
%	 			set to 1 for unipolar values e.g. values that range 0 to 1
%	 			set to 2 for bipolar values e.g. values that range -1 to 1	
%	
%	 climMin: minimum value for colorbar
%	
%	 climMax: maximum value for colorbar
%
%	 ico: icosahedral order for plotting
%
%	 Optional Inputs:
%	 fhtitle: title for the plot
%
%	 legendPosition: position of the colorbar legend, can be 'South' (default) or 'East'
%
%	 viewname: cell array of view names for each subplot, e.g. {'left', 'right', 'right', 'left'}
%	 
% 	 includemat: matrix indicating which hemisphere to plot in each subplot, e.g. [1 0; 0 1; 1 0; 0 1]
%	 
% 	 subplotDims: dimensions of the subplot grid, e.g. [2 2]
%	 
%    figPosition: position of the figure in centimeters, e.g. [0 0 16 10]
%
%	 Output:
%	 fh: figure handle
%	
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

	p = inputParser;
	addParamValue(p, 'fhtitle', '');
	addParamValue(p, 'legendPosition', 'South');
    addParamValue(p, 'viewname', {'left' 'right' 'right' 'left'}, @iscellstr);
    addParamValue(p, 'includemat', [1 0; 0 1; 1 0; 0 1], @ismatrix);
	addParamValue(p, 'subplotDims', [2 2], @(x) isnumeric(x) && numel(x) == 2);  
	addParamValue(p, 'figPosition', [0 0 16 10], @(x) isnumeric(x) && numel(x) == 4);

	parse(p, varargin{:});
	fhtitle = p.Results.fhtitle;
	legendPosition = p.Results.legendPosition;
    viewname = p.Results.viewname;
    includemat = p.Results.includemat;
	subplotDims = p.Results.subplotDims;
	figPosition = p.Results.figPosition;

	numSubplots = prod(subplotDims);
	if length(viewname) > numSubplots
		warning('Number of views exceeds available subplot slots. Truncating view list.');
		viewname = viewname(1:numSubplots);
		includemat = includemat(1:numSubplots, :);
	end

	load SurfView_surfs.mat

	icnum = ico + 1; % index for ico number (icnum = ico + 1)
	% check icnum matches the number of vertices
	if ~isequal(size(icsurfs{icnum}.faces, 1), length(vertvals))
		% they dont match, find which icsurfs{icnum} does match the number of vertices
		try
			for i = 1:length(icsurfs)
				if isequal(size(icsurfs{i}.faces, 1), length(vertvals))
					icnum = i;
					break;
				end
			end
		catch
			error('Number of vertices does not match the icosahedral order');
		end 
	end

	icnvert = size(icsurfs{icnum}.vertices,1); % indices of vertices for specified icosahedral order per hemisphere 
	vertvals_lh = vertvals(1:icnvert); % divide statistics by hemisphere for plotting
	vertvals_rh = vertvals(icnvert+[1:icnvert]);

	% specify limits for plot based on vertvals
	fmax = climMax; % max limit for plotting purposes
	fmin = 0.0; % min limit for plotting purposes
	fmid = fmax/2; % middle value for plotting purposes
	fvals = [fmin fmid fmax]; % this will be passed to the SurfView_show_new function

	curvcontrast = [0.2 0.2]; % contrast of gyri/sulci
	bgcol        = [0 0 0]; % change to [1 1 1] for white background

	% set colorbar limits - usually [fmin fmax] or [-fmax fmax]
	clim = [climMin fmax];

	% Define space
	hvgap = [0.02 0.02];
	%lrgap = [0.02 0.02];
	%btgap = [0.12 0.01];
	if strcmpi(legendPosition, 'south')
		lrgap = [0.02 0.02]; %controls the space from the left of the figure and the right
		if ~isempty(fhtitle)
			btgap = [0.12 0.08]; %controls the space from the bottom of the figure and the top 
		else
			btgap = [0.12 0.01];
		end
	else
		if strcmpi(legendPosition, 'east')
			lrgap = [0.02 0.138];
			if ~isempty(fhtitle)
				btgap = [0.018 0.08];
			else
				btgap = [0.018 0.018];[colnames_model{zi} group]
			end
		end
	end

	% Background - black
	bgcol = [0 0 0];

	% Create figure
	fh = figure('Units', 'centimeters', 'Position', p.Results.figPosition, 'Color', bgcol, 'InvertHardcopy', 'off');

	% Create axes
	allH = tight_subplot(subplotDims(1), subplotDims(2), hvgap, btgap, lrgap);
	hold(allH(:), 'on');

	if ischar(cmap) || isstring(cmap)
		cm = eval(cmap);
	else
		cm = cmap;
	end

    for vi=1:length(viewname)
        view = viewname{vi};
        includevec = includemat(vi,:); % location to plot to
	    axes(allH(vi)); % selecting subplot axes
        SurfView_show_new(surf_lh_pial,surf_rh_pial,vertvals_lh,vertvals_rh,fvals,cm,view, includevec ,curvvec_lh,curvvec_rh,icsurfs{icnum},polarity,curvcontrast,bgcol); 
        set(gca,'visible','off'); 
        axis tight;
    end 

	if ~isempty(fhtitle)
		titleAx = axes;
		set(titleAx,'position',[0 0 1 1],'units','normalized');axis off;
		text(titleAx, 0.5,1,fhtitle,'color','w','fontweight','bold','interpreter','none','verticalalignment','top','horizontalalignment','center','fontsize',14)
	end

	% Set colorbar	
	colormap(cmap);
	cb = colorbar('color', 'w');
	cb.Label.String = str; 
	cb.Label.FontSize = 12;
	cb.FontSize = 10;
	cb.Box = 'off';
	cb.Label.FontWeight = 'bold';   
	cb.Location = 'south';
	caxis(clim);

	if strcmpi(legendPosition, 'south')
		cb.Location = 'south';
		if ~isempty(fhtitle)
			%cb.Position(1) = allH(1).Position(1);
			cb.Position(1) = allH(1).Position(1) + hvgap(1);
			cb.Position(2) = cb.Position(2) - hvgap(1);
			%cb.Position(3) = allH(1).Position(3)*2 + hvgap(1);
			cb.Position(3) = allH(1).Position(3)*2 - hvgap(1);
		else
			%cb.Position(1) = allH(1).Position(1);
			cb.Position(1) = allH(1).Position(1) + hvgap(1);
			cb.Position(2) = cb.Position(2) - btgap(1);
			%cb.Position(3) = allH(1).Position(3)*2 + hvgap(1);
			cb.Position(3) = allH(1).Position(3)*2 - hvgap(1);
		end
	else
		if strcmpi(legendPosition, 'east')
			cb.Location = 'eastoutside';
			if ~isempty(fhtitle)
				cb.Position(1) = allH(4).Position(1) + allH(4).Position(3) + 0.01;
				cb.Position(2) = allH(3).Position(2);
				cb.Position(4) = allH(1).Position(4)*2 + hvgap(1);
			else
				cb.Position(1) = allH(4).Position(1) + allH(4).Position(3) + 0.16;
				cb.Position(2) = allH(3).Position(2);
				cb.Position(4) = allH(1).Position(4)*2 + hvgap(1);
			end
		end
	end
end