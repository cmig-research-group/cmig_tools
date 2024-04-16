% r,c,s coords of UI remain in scale of first volume passed in; use this to convert to
% current volume's coords
function x = scaleCoord(x,sf)
%x = floor(sf*(x-1)+1); %this is for one-based indexing
x = floor(sf*x); %but new 1mm space is centered as if 0-based indexing
