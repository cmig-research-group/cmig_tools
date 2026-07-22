function vararginOut = FEMA_mergeArgs(fname_json, vararginIn)
% FEMA_MERGEARGS  Apply params.makeDesign from JSON before inputParser.
%
% Non-DEAP TOML writes an optional [makeDesign] table via FEMA_createInputJSON into
% cfg.params.makeDesign. Name-value pairs listed here match FEMA_makeDesign addParameter
% names. Order is [defaults from JSON, vararginIn] so explicit MATLAB arguments win.
%
% Valid field names: iid, eid, dataFile, dirTabulated, dropMissing, outDir,
% outName, outType, study, fname_fam

    vararginOut = vararginIn;
    if isempty(fname_json) || ~ischar(fname_json)
        return
    end
    if ~exist(fname_json, 'file')
        return
    end

    cfg = jsondecode(fileread(fname_json));
    if ~isfield(cfg, 'params') || ~isfield(cfg.params, 'makeDesign')
        return
    end

    md = cfg.params.makeDesign;
    if isempty(md)
        return
    end
    if ~isstruct(md)
        return
    end
    if numel(md) ~= 1
        warning('FEMA_mergeArgs: params.makeDesign must be a single struct; using first element.');
        md = md(1);
    end

    known = {'iid', 'eid', 'dataFile', 'dirTabulated', ...
             'dropMissing', 'outDir', 'outName', 'outType', 'study', 'fname_fam'};

    ff = fieldnames(md);
    userKeys = vararginNameKeysLower(vararginIn);

    jsonNV = {};
    for k = 1:numel(ff)
        rawName = ff{k};
        ix = find(strcmp(lower(rawName), lower(known)), 1);
        if isempty(ix)
            warning('FEMA_mergeArgs: ignoring unknown [makeDesign] field "%s".', rawName);
            continue
        end
        pname = known{ix};
        if ismember(lower(pname), userKeys)
            continue
        end
        val = md.(rawName);
        val = coerceMakeDesignValue(pname, val);
        jsonNV = [jsonNV, {pname, val}]; %#ok<AGROW>
    end

    % JSON defaults first; explicit varargin wins on duplicate (inputParser last-wins)
    vararginOut = [jsonNV, vararginIn];
end


function keys = vararginNameKeysLower(args)
    keys = {};
    for k = 1:2:numel(args) - 1
        if ischar(args{k}) || isstring(args{k})
            keys{end+1} = lower(char(args{k})); %#ok<AGROW>
        end
    end
end


function v = coerceMakeDesignValue(pname, v)
    switch lower(pname)
        case {'dropmissing'}
            if isnumeric(v) || islogical(v)
                v = logical(v);
            end
        otherwise
            % char, string, string array -> leave; json may decode string array as cell
    end
end
