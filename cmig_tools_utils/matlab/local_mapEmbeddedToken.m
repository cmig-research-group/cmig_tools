function [name_out, tf] = local_mapEmbeddedToken(name_in, map_use, direction)
    name_out = '';
    tf = false;
    keys_map = keys(map_use);
    % Longest-first so e.g. lobfrt__lh wins over lh inside ...lobfrt__lh_mean
    [~, ord] = sort(cellfun(@length, keys_map), 'descend');
    keys_map = keys_map(ord);
    for k = 1:numel(keys_map)
        key_name = keys_map{k};
        patt = ['(^|__)', regexptranslate('escape', key_name), '(?=_|$)'];
        [st, ~] = regexp(name_in, patt, 'start', 'end', 'once');
        if isempty(st)
            continue
        end
        if st > 1 && st+1 <= length(name_in) && strcmp(name_in(st:st+1), '__')
            key_start = st + 2;
        else
            key_start = st;
        end
        key_end = key_start + length(key_name) - 1;
        map_val = map_use(key_name);
        if strcmp(direction, 'tab2atlas')
            name_out = map_val;
        else
            name_out = [name_in(1:key_start-1), map_val, name_in(key_end+1:end)];
        end
        tf = true;
        return
    end
end