function [rows, cols] = calc_rows_cols_subplot(num_plots, balance_value)
% Function to calculate the number of rows and columns in a subplot given
% the number of plots to be made
%% Inputs:
% num_plots:     the actual number of plots to be made
% balance_value: difference between the number of rows and columns
% 
%% Outputs:
% rows:          the number of rows in the subplot
% cols:          the number of columns in the subplot
% 
%% Notes:
% Depending on the value of balance_value, the resulting row and column
% combination can lead to many more plots than are actually needed
% 
% If balance_value = 0, that is, only same number of rows and columns
% are expected. The largest value of num_plots in this case is 10000
% 
% If num_plots <= 2, rows = 1, cols = num_plots
% 
% Also see numSubplots function on FileExchange: 
% https://www.mathworks.com/matlabcentral/fileexchange/26310-numsubplots-neatly-arrange-subplots
% 
%% Default:
% balance_value: 3
% 
%% Author(s)
% Parekh, Pravesh
% August 23, 2018
% MBIAL

%% Check num_plots
if ~exist('num_plots', 'var') || isempty(num_plots)
    error('Number of plots should be provided');
else
    % Return as columns if num_plots <= 2
    if num_plots <= 2
        rows = 1;
        cols = num_plots;
        return
    end
end

% Check balance_value
if ~exist('balance_value', 'var') || isempty(balance_value)
    balance_value = 3;
else
    if balance_value < 0
        warning('Balance value cannot be less than 0; using default of 3');
        balance_value = 3;
    else
        if balance_value == 0 && num_plots > 10000
            error('Cannot compute symmetrical subplots for values greater than 10000');
        end
    end
end

%% Attempt to figure out rows and cols
% Check special case of balance_value = 0
if balance_value == 0
    % Find squares of numbers till a number is found who's square is larger
    % than the input value of num_plots; only check the first 100 numbers;
    rows = find(((1:100).^2)>=num_plots,1);
    cols = rows;
    return
end

% Check if odd number
if mod(num_plots,2) == 1
    num_plots = num_plots+1;
end

% Check for squares
if mod(sqrt(num_plots),1) == 0
    rows = sqrt(num_plots);
    cols = rows;
    return;
end
    
% Get factors
facs = factor(num_plots);
% If two factors, 
% rows = factor(1), 
% cols = factor(2)
if length(facs) == 2
    rows = facs(1);
    cols = facs(2);
else
    % If three factors, 
    % rows = factor(1)*factor(2), 
    % cols = factor(3)
    if length(facs) == 3
        rows = facs(1)*facs(2);
        cols = facs(3);
    else
        % If factors>3, 
        % rows = multiplication of factors(1:end-2), 
        % cols = multiplication of last two factors
        tmp = cumprod(facs(1:end-2));
        rows = tmp(end);
        cols = facs(end-1)*facs(end);
    end
end

% Recompute if difference between the number of rows and columns exceeds
% the balance_value
if abs(rows - cols) > balance_value
    num_plots = num_plots + 1;
    [rows, cols] = calc_rows_cols_subplot(num_plots, balance_value);
end