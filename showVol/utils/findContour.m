function varargout = findContour(surfStruct, vertVals, isoVal, varargin)
% Solve the isocontour problem on a triangular mesh surface.
%
% findContour(surfStruct, vertVals, isoVal, [contourColor], [lineWidth])
%
% Outputs:
% -------
% 1) lh = findContour(...) returns the line handles for the contour
%
% 2) [startPts, endPts] = findContour(...) returns the points on the
%                         ends of the contour line segments
%
% Note:
% ----
% 1) contourColor is an RGB triple, defaults to [1 1 0] (yellow)
% 2) lineWidth defaults to 3

contourColor = [1 1 0];
if nargin >= 4
  contourColor = varargin{1};
end

lineWidth = 3;
if nargin >= 5
  lineWidth = varargin{2};
end

vertVals = vertVals - isoVal;

[Xs, Ys, Zs, Xe, Ye, Ze] = findContourMEX(surfStruct.faces, ...
					  surfStruct.vertices, ...
					  vertVals);

if (nargout == 0) | (nargout == 1)
  
  numSegments = length(Xs);
  lh = zeros(numSegments, 1);
  for ss = 1:numSegments
    lh(ss) = line([Xs(ss); Xe(ss)], [Ys(ss); Ye(ss)], [Zs(ss); Ze(ss)]);
    set(lh(ss), 'LineWidth', lineWidth);
    set(lh(ss), 'Color', contourColor);
  end

  varargout{1} = lh;
  
end
  
if nargout == 2
  varargout{1} = [Xs Ys Zs];
  varargout{2} = [Xe Ye Ze];
end



