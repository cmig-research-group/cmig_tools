function out = save_numericToJsonCells(A)
    % Encode numeric array as nested cells so jsonencode preserves rows; non-finite -> char tokens.
    A = squeeze(A);
    if isempty(A)
        out = A;
        return
    end
    if ndims(A) > 2
        szA = size(A);
        C = cell(szA(1), 1);
        idxCell = repmat({':'}, 1, ndims(A));
        for i = 1:szA(1)
            idxCell{1} = i;
            C{i} = save_numericToJsonCells(A(idxCell{:}));
        end
        out = C;
        return
    end
    if isvector(A)
        row = num2cell(A(:).');
        for j = 1:numel(row)
            row{j} = save_nonfiniteScalarTokens(row{j});
        end
        out = row;
        return
    end
    [m, n] = size(A);
    C = cell(m, 1);
    for i = 1:m
        row = num2cell(A(i, :));
        for j = 1:n
            row{j} = save_nonfiniteScalarTokens(row{j});
        end
        C{i} = row;
    end
    out = C;
end