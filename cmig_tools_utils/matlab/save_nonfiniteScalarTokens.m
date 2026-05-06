function v = save_nonfiniteScalarTokens(v)
    if isnumeric(v) && isscalar(v) && ~isinteger(v)
        if isnan(v)
            v = 'NaN';
        elseif isinf(v) && v > 0
            v = 'Inf';
        elseif isinf(v)
            v = '-Inf';
        end
    end
end