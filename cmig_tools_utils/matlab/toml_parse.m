function s = toml_parse(fname)
% TOML_PARSE  Parse a TOML config file into a MATLAB struct.
%
% Supports the subset of TOML needed for FEMA non-DEAP config files:
%   - Bare sections:       [section] and [section.subsection]
%   - Arrays of tables:    [[section.array]]
%   - Key = value pairs:   strings, numbers, booleans, inline arrays
%   - Arrays may span multiple lines (e.g. fstem = [ "a", "b" ])
%
% Usage:
%   s = toml_parse('/path/to/config.toml')
%
% Returns a nested MATLAB struct.

    lines = strsplit(fileread(fname), '\n');

    s           = struct();
    currentPath = {};     % cell of keys for current section (plain [x])
    arrayPath   = {};     % cell of keys for current array-of-tables ([[x]])
    inArray     = false;  % true when we are inside an [[...]] block

    iLine = 1;
    while iLine <= numel(lines)
        line = strtrim(lines{iLine});

        % Strip inline comments (# not inside a string)
        line = stripComment(line);
        if isempty(line)
            iLine = iLine + 1;
            continue
        end

        % ---- Array of tables: [[a.b.c]] ----
        if numel(line) >= 4 && strcmp(line(1:2), '[[') && strcmp(line(end-1:end), ']]')
            keyStr   = strtrim(line(3:end-2));
            arrayPath = strsplit(keyStr, '.');
            currentPath = {};
            inArray   = true;

            % Grow the array at that path
            s = appendToArray(s, arrayPath);
            iLine = iLine + 1;
            continue
        end

        % ---- Plain section: [a.b.c] ----
        if numel(line) >= 2 && line(1) == '[' && line(end) == ']' && ...
                ~(numel(line) >= 3 && line(2) == '[')
            keyStr      = strtrim(line(2:end-1));
            currentPath = strsplit(keyStr, '.');
            arrayPath   = {};
            inArray     = false;
            iLine = iLine + 1;
            continue
        end

        % ---- Key = value ----
        eqIdx = find(line == '=', 1);
        if isempty(eqIdx)
            iLine = iLine + 1;
            continue
        end
        key   = strtrim(line(1:eqIdx-1));
        value = strtrim(line(eqIdx+1:end));

        if isempty(key)
            iLine = iLine + 1;
            continue
        end

        % Multi-line array: first line is e.g. key = [ or key = [ "a",
        if ~isempty(value) && value(1) == '[' && ~isBracketArrayComplete(value)
            [value, iLine] = accumulateMultilineArray(lines, iLine, value);
        end

        parsedValue = parseValue(value);

        if inArray
            % Set field on the LAST element of the array
            s = setInLastArrayElement(s, arrayPath, key, parsedValue);
        else
            s = setNested(s, [currentPath, {key}], parsedValue);
        end

        iLine = iLine + 1;
    end
end

% -------------------------------------------------------------------------
function [value, iEnd] = accumulateMultilineArray(lines, iFirst, value)
% Append subsequent physical lines until [...] has balanced brackets.
    i = iFirst;
    v = value;
    while ~isBracketArrayComplete(v)
        i = i + 1;
        if i > numel(lines)
            break
        end
        nextLine = strtrim(lines{i});
        nextLine = stripComment(nextLine);
        v = [v, ' ', nextLine]; %#ok<AGROW>
    end
    value = v;
    iEnd = i;
end

% -------------------------------------------------------------------------
function tf = isBracketArrayComplete(str)
% True when square brackets are balanced (quotes ignored), for array literals.
    str = strtrim(str);
    if isempty(str) || str(1) ~= '['
        tf = true;
        return
    end
    depth = 0;
    inDouble = false;
    inSingle = false;
    for k = 1:numel(str)
        c = str(k);
        if c == '"' && ~inSingle
            inDouble = ~inDouble;
        elseif c == '''' && ~inDouble
            inSingle = ~inSingle;
        elseif ~inDouble && ~inSingle
            if c == '['
                depth = depth + 1;
            elseif c == ']'
                depth = depth - 1;
            end
        end
    end
    tf = (depth == 0);
end

% -------------------------------------------------------------------------
function line = stripComment(line)
% Remove everything after an unquoted '#'
    inSingle = false;
    inDouble = false;
    for k = 1:numel(line)
        c = line(k);
        if c == '''' && ~inDouble
            inSingle = ~inSingle;
        elseif c == '"' && ~inSingle
            inDouble = ~inDouble;
        elseif c == '#' && ~inSingle && ~inDouble
            line = strtrim(line(1:k-1));
            return
        end
    end
end

% -------------------------------------------------------------------------
function v = parseValue(str)
% Parse a TOML value string into a MATLAB value.
    str = strtrim(str);

    % Inline array: [...]
    if ~isempty(str) && str(1) == '[' && str(end) == ']'
        v = parseArray(str(2:end-1));
        return
    end

    % Double-quoted string
    if numel(str) >= 2 && str(1) == '"' && str(end) == '"'
        v = str(2:end-1);
        return
    end

    % Single-quoted string
    if numel(str) >= 2 && str(1) == '''' && str(end) == ''''
        v = str(2:end-1);
        return
    end

    % Boolean
    if strcmpi(str, 'true')
        v = true;
        return
    end
    if strcmpi(str, 'false')
        v = false;
        return
    end

    % Numeric
    num = str2double(str);
    if ~isnan(num)
        v = num;
        return
    end

    % Bare string (unquoted — treat as string)
    v = str;
end

% -------------------------------------------------------------------------
function arr = parseArray(inner)
% Parse the contents of [...] into a MATLAB cell array or numeric array.
    inner   = strtrim(inner);
    if isempty(inner)
        arr = {};
        return
    end

    elements = splitArrayElements(inner);
    elements = cellfun(@strtrim, elements, 'UniformOutput', false);
    elements = elements(~cellfun('isempty', elements));

    if isempty(elements)
        arr = {};
        return
    end

    % Parse each element
    parsed = cellfun(@parseValue, elements, 'UniformOutput', false);

    % If all elements are numeric scalars, return numeric vector
    if all(cellfun(@(x) isnumeric(x) && isscalar(x), parsed))
        arr = cell2mat(parsed);
    else
        arr = parsed;
    end
end

% -------------------------------------------------------------------------
function parts = splitArrayElements(str)
% Split comma-separated array elements, respecting quoted strings.
    parts    = {};
    depth    = 0;
    inDouble = false;
    inSingle = false;
    start    = 1;
    for k = 1:numel(str)
        c = str(k);
        if c == '"' && ~inSingle
            inDouble = ~inDouble;
        elseif c == '''' && ~inDouble
            inSingle = ~inSingle;
        elseif c == '[' && ~inDouble && ~inSingle
            depth = depth + 1;
        elseif c == ']' && ~inDouble && ~inSingle
            depth = depth - 1;
        elseif c == ',' && depth == 0 && ~inDouble && ~inSingle
            parts{end+1} = str(start:k-1); %#ok<AGROW>
            start = k + 1;
        end
    end
    parts{end+1} = str(start:end);
end

% -------------------------------------------------------------------------
function s = setNested(s, pathParts, value)
% Set s.a.b.c = value given pathParts = {'a','b','c'}.
    if isempty(pathParts)
        return
    end
    key = matlab.lang.makeValidName(pathParts{1});
    if numel(pathParts) == 1
        s.(key) = value;
    else
        if ~isfield(s, key) || ~isstruct(s.(key))
            s.(key) = struct();
        end
        s.(key) = setNested(s.(key), pathParts(2:end), value);
    end
end

% -------------------------------------------------------------------------
function s = appendToArray(s, pathParts)
% Append an empty struct to the cell array of structs at s.a.b.
% Arrays-of-tables ([[...]]) are stored as cell arrays to avoid MATLAB's
% requirement that struct array elements all have identical fields.
    key = matlab.lang.makeValidName(pathParts{1});
    if numel(pathParts) == 1
        if ~isfield(s, key)
            s.(key) = {struct()};
        else
            existing = s.(key);
            if iscell(existing)
                s.(key){end+1} = struct();
            else
                % First repeat: convert single struct to cell array, then grow
                s.(key) = {existing, struct()};
            end
        end
    else
        if ~isfield(s, key) || ~isstruct(s.(key))
            s.(key) = struct();
        end
        s.(key) = appendToArray(s.(key), pathParts(2:end));
    end
end

% -------------------------------------------------------------------------
function s = setInLastArrayElement(s, arrayPath, key, value)
% Set a field on the last element of the cell array of structs at arrayPath.
    key = matlab.lang.makeValidName(key);
    if numel(arrayPath) == 1
        akey = matlab.lang.makeValidName(arrayPath{1});
        arr  = s.(akey);
        if iscell(arr)
            arr{end}.(key) = value;
        else
            arr(end).(key) = value;
        end
        s.(akey) = arr;
    else
        akey = matlab.lang.makeValidName(arrayPath{1});
        s.(akey) = setInLastArrayElement(s.(akey), arrayPath(2:end), key, value);
    end
end
