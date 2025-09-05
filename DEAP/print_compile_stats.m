function print_compile_stats()
%PRINT_COMPILE_STATS Display build information from compile_stats.json
%
% This function reads the compile statistics JSON file and displays
% the build information including host, user, git details, and timestamp.

    % Try multiple possible locations for the compile stats file
    % Priority: current dir (where executable runs), /app (container), relative paths
    possible_files = {'./compile_stats.json', 'compile_stats.json', '/app/compile_stats.json', '../compile_stats.json'};
    compile_stats_file = '';
    
    for i = 1:length(possible_files)
        if exist(possible_files{i}, 'file')
            compile_stats_file = possible_files{i};
            break;
        end
    end
    
    % Check if the compile stats file exists
    if isempty(compile_stats_file)
        fprintf('=== BUILD INFO ===\n');
        fprintf('Compile stats file not found in any expected location.\n');
        fprintf('Build information unavailable.\n');
        fprintf('==================\n\n');
        return;
    end
    
    try
        % Read the JSON file
        fid = fopen(compile_stats_file, 'r');
        if fid == -1
            error('Could not open compile stats file');
        end
        
        raw_json = fread(fid, inf, 'char=>char')';
        fclose(fid);
        
        % Parse JSON (using jsondecode if available, otherwise manual parsing)
        if exist('jsondecode', 'builtin') || exist('jsondecode', 'file')
            stats = jsondecode(raw_json);
        else
            % Fallback for older MATLAB versions - simple manual parsing
            stats = parse_simple_json(raw_json);
        end
        
        % Display the build information
        fprintf('=== BUILD INFO ===\n');
        fprintf('Host: %s\n', get_field_safe(stats, 'build_host', 'unknown'));
        fprintf('User: %s\n', get_field_safe(stats, 'build_user', 'unknown'));
        fprintf('Git Branch: %s\n', get_field_safe(stats, 'git_branch', 'unknown'));
        fprintf('Git Hash: %s\n', get_field_safe(stats, 'git_hash', 'unknown'));
        fprintf('Git Status: %s\n', get_field_safe(stats, 'git_status', 'unknown'));
        fprintf('Build Type: %s\n', get_field_safe(stats, 'build_type', 'unknown'));
        fprintf('Build Time: %s\n', get_field_safe(stats, 'build_timestamp', 'unknown'));
        fprintf('==================\n\n');
        
    catch ME
        fprintf('=== BUILD INFO ===\n');
        fprintf('Error reading compile stats: %s\n', ME.message);
        fprintf('Build information unavailable.\n');
        fprintf('==================\n\n');
    end
end

function value = get_field_safe(struct_data, field_name, default_value)
%GET_FIELD_SAFE Safely get field value from struct with default fallback
    if isfield(struct_data, field_name)
        value = struct_data.(field_name);
        if isempty(value)
            value = default_value;
        end
    else
        value = default_value;
    end
end

function stats = parse_simple_json(json_str)
%PARSE_SIMPLE_JSON Simple JSON parser for basic key-value pairs
% This is a fallback for older MATLAB versions without jsondecode
    
    stats = struct();
    
    % Remove whitespace and braces
    json_str = strrep(json_str, sprintf('\n'), '');
    json_str = strrep(json_str, sprintf('\r'), '');
    json_str = strrep(json_str, sprintf('\t'), '');
    json_str = strrep(json_str, ' ', '');
    json_str = strrep(json_str, '{', '');
    json_str = strrep(json_str, '}', '');
    
    % Split by commas
    pairs = strsplit(json_str, ',');
    
    for i = 1:length(pairs)
        pair = pairs{i};
        if contains(pair, ':')
            parts = strsplit(pair, ':');
            if length(parts) == 2
                key = strrep(parts{1}, '"', '');
                value = strrep(parts{2}, '"', '');
                stats.(key) = value;
            end
        end
    end
end
