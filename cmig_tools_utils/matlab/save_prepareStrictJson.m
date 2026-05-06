function out = save_prepareStrictJson(data)
    if isnumeric(data) && ~isinteger(data)
        if isempty(data)
            out = data;
            return
        end
        if isscalar(data)
            if isnan(data)
                out = 'NaN';
            elseif isinf(data) && data > 0
                out = 'Inf';
            elseif isinf(data)
                out = '-Inf';
            else
                out = data;
            end
            return
        end
        if all(isfinite(data(:)))
            out = data;
            return
        end
        out = save_numericToJsonCells(data);
        return
    end
    if isstruct(data)
        out = data;
        fn = fieldnames(out);
        for i = 1:numel(out)
            for j = 1:numel(fn)
                out(i).(fn{j}) = save_prepareStrictJson(out(i).(fn{j}));
            end
        end
        return
    end
    if iscell(data)
        out = cell(size(data));
        for k = 1:numel(data)
            out{k} = save_prepareStrictJson(data{k});
        end
        return
    end
    if istable(data)
        out = save_prepareStrictJson(table2struct(data));
        return
    end
    out = data;
end
