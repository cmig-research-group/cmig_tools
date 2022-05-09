function [] = logging(text, varargin)
%logging - Standard logging utility
%
% Syntax: [] = logging(text)
%
% Logs in the following format: 04-Apr-2022 18:28:01 - mfilename - text
stack = dbstack;
if length(stack) == 1
    current_file = 'running script';
else
    current_file = stack(2).name;
end

if length(varargin) == 0
    fprintf(1, '%s in %s - %s\n', datestr(now()), current_file, text);
else
    fprintf(1, '%s in %s - %s\n', datestr(now()), current_file, sprintf(text, varargin{:}));
end

end