%% ------------------------------------------------------------
%% --- Fly to row/col/slice number stopIndex
%% ------------------------------------------------------------

function handles = fly(handles, stopIndex, varargin)

if nargin < 3
  pauseTime = 0;
else
  pauseTime = varargin{1};
end

switch handles.ORIENTATION
  case 1
    handles = fly_ss(handles, stopIndex, pauseTime);
  case 2
    handles = fly_cc(handles, stopIndex, pauseTime);
  case 3
    handles = fly_rr(handles, stopIndex, pauseTime);
end

%% end of function handles = fly(handles, stopIndex, varargin)
%% ------------------------------------------------------------