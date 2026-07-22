function txt = save_jsonencode(data, varargin)
    %JSONENCODE with strict JSON: non-finite floating-point values become string tokens.
    txt = jsonencode(save_prepareStrictJson(data), varargin{:});
end