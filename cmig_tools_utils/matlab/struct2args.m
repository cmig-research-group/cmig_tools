function nvPairs = struct2args(structure, prefix)
%   struct2args  Convert nested structure to name-value cell array
%
%   Recursively flattens nested structures and cell arrays of structs.
%   Produces dot-separated field names suitable for use in name-value lists.

    if nargin < 2
        prefix = '';
    end

    nvPairs = {};
    f = fieldnames(structure);

    for i = 1:numel(f)
        fieldName = f{i};
        value = structure.(fieldName);

        % build dotted name
        if isempty(prefix)
            fullName = fieldName;
        else
            fullName = [prefix '.' fieldName];
        end

        % handle value types
        if isstruct(value)
            % recursively expand structs
            nvPairs = [nvPairs, struct2args(value, fullName)];

        elseif iscell(value)
            % handle cell arrays (e.g. {struct, struct, ...})
            for j = 1:numel(value)
                elem = value{j};
                if isstruct(elem)
                    nvPairs = [nvPairs, struct2args(elem, sprintf('%s{%d}', fullName, j))];
                else
                    nvPairs = [nvPairs, {sprintf('%s{%d}', fullName, j), elem}];
                end
            end

        elseif isempty(value)
            % optionally skip empty fields
            continue;

        else
            % base case: regular name–value pair
            nvPairs = [nvPairs, {fullName, value}];
        end
    end
end