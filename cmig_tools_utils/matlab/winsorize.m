function output =  winsorize(inputVars, varargin);

    % Winsorize 

    p = inputParser;
    addRequired(p, 'inputVars', @(x) validateattributes(x, {'numeric'}, {'2d'}));
    addOptional(p, 'percentile', [1, 99], @(x) validateattributes(x, {'numeric'}, {'vector', 'numel', 2}));

    parse(p, inputVars, varargin{:});
    inputVars = p.Results.inputVars;
    percentile = p.Results.percentile;

    % Compute the lower and upper bounds for winsorization
    lowerBound = prctile(inputVars, percentile(1), 1);
    upperBound = prctile(inputVars, percentile(2), 1);

    % Apply winsorization
    output = inputVars;
    output(output < lowerBound) = lowerBound;
    output(output > upperBound) = upperBound;

end 
