function atlasName = atlas2tabNames(roiAtlas, direction)
    if nargin < 2 || isempty(direction)
        direction = 'tab2atlas';
    end
    validatestring(direction, {'atlas2tab', 'tab2atlas'}, mfilename, 'direction', 2);

    key = lower(roiAtlas);

    switch direction
        case 'atlas2tab'
            switch key
                case {'dsk', 'aparc'}
                    atlasName = 'aparc';
                case {'dst', 'aparc_a2009s'}
                    atlasName = 'aparc_a2009s';
                case {'at'}
                    atlasName = 'fiber';
                case {'aseg'}
                    atlasName = 'aseg';
                otherwise
                    atlasName = key;
            end
        case 'tab2atlas'
            switch key
                case {'aparc'}
                    atlasName = 'dsk';
                case {'aparc_a2009s'}
                    atlasName = 'dst';
                case {'fiber'}
                    atlasName = 'at';
                case {'aseg'}
                    atlasName = 'aseg';
                otherwise
                    atlasName = key;
            end
    end
end