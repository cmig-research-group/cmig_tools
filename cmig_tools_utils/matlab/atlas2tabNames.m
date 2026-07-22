function atlasName = atlas2tabNames(roiAtlas, direction)

    % atlas2tabNames - Convert atlas names to tabulated names or vice versa

    % Inputs:
    %   roiAtlas - cell array of atlas names
    %   direction - mapping direction (optional):
    %               'atlas2tab' or 'tab2atlas' (default) 

    % Mapping 
    % Tabulated atlas names | Atlas names
    % --------------------- | --------------------
    % dsk                   | aparc
    % dst                   | aparc_a2009s
    % at                    | fiber
    % aseg                  | aseg
    % fzy                   | fuzzy
    % gp                    | gordon
    % --------------------- | --------------------


    if nargin < 2 || isempty(direction)
        direction = 'tab2atlas';
    end
    validatestring(direction, {'atlas2tab', 'tab2atlas'}, mfilename, 'direction', 2);

    key = lower(roiAtlas);

    switch direction
        case 'atlas2tab'
            switch key
                case {'aparc', 'aparcaseg'}
                    atlasName = 'dsk';
                case {'aparc_a2009s'}
                    atlasName = 'dst';
                case {'fiber'}
                    atlasName = 'at';
                case {'aseg'}
                    atlasName = 'aseg';
                case {'fuzzy'}
                    atlasName = 'fzy';
                otherwise
                    atlasName = key;
            end
        case 'tab2atlas'
            switch key
                case {'dsk', 'aparc'}
                    atlasName = 'aparc';
                case {'dst', 'aparc_a2009s'}
                    atlasName = 'aparc_a2009s';
                case {'at'}
                    atlasName = 'fiber';
                case {'aseg'}
                    atlasName = 'aseg';
                case {'fzy'}
                    atlasName = 'fuzzy';
                otherwise
                    atlasName = key;
            end
    end
end