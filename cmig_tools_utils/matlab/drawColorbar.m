function drawColorbar(ax, colormap, limits, labelColor)
% drawColorbar  draw an imaging colorbar
%
%   drawColorbar(ax, colormap, limits, [labelColor])
%
%   ax:         axis handle to draw in
%   colormap:   Nx3 or Nx4 (with alpha) colormap
%   limits:     [min max] or [min max threshold] labels for colorbar
%   labelColor: optional. 'white' (default), 'black' or other latex color name (see latexcolor.com)
%

%TODO create axis if none is specified, e.g. ax='lrv' vertical axis in lower right of current figure, etc

if ~exist('labelColor','var')
  labelColor = 'white';
end

cla(ax)
cbar = colormap(:,1:3);
if size(colormap,2)==4
  cbar = cbar .* colormap(:,4); %apply the alpha masking
end

image(ax, permute(cbar,[1 3 2]))
axis(ax,'xy')
box(ax,'off')

yl = ylim(ax);
if numel(limits) == 2
  ylabelopts = {'YtickLabel',{['\color{' labelColor '}\bf' num2str(limits(1),3)] ...
    ['\color{' labelColor '}\bf' num2str(mean(limits),3)] ...
    ['\color{' labelColor '}\bf' num2str(limits(2),3)]},...
    'YTick',[yl(1) range(yl)/2 yl(2)]};
elseif numel(limits) == 3
  %display threshold points
  tr = limits(3)/limits(2); % threshold wrt max
  mid = range(yl)/2;
  ut = mid + mid*tr;
  lt = mid - mid*tr;
  ylabelopts = {'YtickLabel',{['\color{' labelColor '}\bf' num2str(limits(1),3)]...
    ['\color{' labelColor '}\bf' num2str(-limits(3),3)]...
    ['\color{' labelColor '}\bf' num2str(limits(3),3)]...
    ['\color{' labelColor '}\bf' num2str(limits(2),3)]},...
    'YTick',[yl(1) lt ut yl(2)]};
end

set(ax,'XTickLabel','','XTick',-.5,'YColor','k',ylabelopts{:})

