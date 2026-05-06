function roi2atlas = roi2FreesurferAtlas(icoNum, fname_out, varargin)
    % roi2FreesurferAtlas - Convert ROIs to surface atlas format
    % works with surface atlases aparc, aparc, aparc.a2009s, yeo17 and yeo7
    % create a logical matrix to map roi data to surface atlases
    %
    % Required input:
    %   icoNum        - ico-number for vertexwise analyses (0-based)
    %   fname_out - output filename to save the roi2atlas structure
    % 
    % Optional inputs:
    %   parc_in - parcellation to process, can be a string or cell array of strings
    %             default: {'aparc', 'aparc.a2009s', 'yeo7', 'yeo17'}
    %   dirname_fs - Freesurfer home directory
    %               default: getenv('FREESURFER_HOME')
    %   splitLR - whether to split left and right hemispheres
    %             default: false
    % 
    % Outputs 
    %   roi2atlas - a structure with fields for each parcellation

    p = inputParser;
    addParameter(p, 'parc_in', {'aparc', 'aparc.a2009s'}, ...
                               @(x) ischar(x) || isstring(x) || iscell(x) && ...
                                ismember(x, {'aparc', 'aparc.a2009s'}));
    addParameter(p, 'dirname_fs', getenv('FREESURFER_HOME'), ...
                               @(x) ischar(x) || isstring(x));
    addParameter(p, 'splitLR', false, @(x) islogical(x) || ismember(x, [0 1]));
    addParameter(p, 'save_flag', true, @(x) islogical(x));


    parse(p, varargin{:});
    parc_in = p.Results.parc_in;  
    dirname_fs = char(p.Results.dirname_fs);
    splitLR = p.Results.splitLR;
    save_flag = p.Results.save_flag;

    % check if fname_out is empty if so set save_flag to false
    if isempty(fname_out)
        save_flag = false;
    end

    % make parc_in a cell array for uniform processing
    if ischar(parc_in) || isstring(parc_in)
        parc_in = {char(parc_in)};
    elseif iscell(parc_in)
        parc_in = cellfun(@char, parc_in, 'UniformOutput', false);
    end

    for pi = 1:length(parc_in)
        parc_name = parc_in{pi};
        switch parc_name
            case {'aparc', 'dsk'}
                annot = 'aparc';
            case {'aparc.a2009s', 'dst'}
                annot = 'aparc.a2009s';
            otherwise
                error('Unknown parcellation name: %s', parc_name);
        end
        
        data  = read_fs_annot(dirname_fs, icoNum, annot, splitLR);
                
        % roicodes, roinames, roirgb
        roicodes = data.labels.key;
        roinames = data.labels.name;
        roirgb = data.labels.rgba;

        % create logical matrix of roi vertices
        nroi = length(roicodes); % assuming both hemispheres have same number of rois
        numVertices = numIcoVertices(icoNum)*2; % this is for a single hemisphere
        roimat = false(nroi, numVertices);
        % left comes first 
        for r = 1:nroi
            % check key and name is same for both hemispheres
            roi_code = roicodes(r);
            roi_name = roinames{r};
            %find vertices
            idx = find(data.cdata == roi_code);
            roimat(r, idx) = true;
        end
        %add 'lh', 'rh' and 'all' as rois to match tabulated data
        roicodes(end+1) = -1; % all lh
        roicodes(end+1) = -2; % all rh
        roicodes(end+1) = -3; % all vertices
        roinames{end+1} = 'lh';
        roinames{end+1} = 'rh';
        roinames{end+1} = 'all';
        roimat(end+1, 1:numVertices/2) = true; % all lh
        roimat(end+1, numVertices/2+1:end) = true; % all rh
        roimat(end+1, :) = true; % all vertices
        

        % Store in output structure with parcellation name as field
        if splitLR
            % find left and right hemisphere roicodes, roinames, roirgb
            idx_rh = find(contains(roinames, 'rh'));
            idx_lh = find(contains(roinames, 'lh'));
            idx_other = setdiff(1:length(roinames), [idx_lh; idx_rh]);
            roi2atlas.(matlab.lang.makeValidName(parc_name)).roimat_lh = roimat(idx_lh, 1:numVertices/2);
            roi2atlas.(matlab.lang.makeValidName(parc_name)).roimat_rh = roimat(idx_rh, numVertices/2+1:end);
            roi2atlas.(matlab.lang.makeValidName(parc_name)).roinames_lh = roinames(idx_lh); 
            roi2atlas.(matlab.lang.makeValidName(parc_name)).roinames_rh = roinames(idx_rh);
            roi2atlas.(matlab.lang.makeValidName(parc_name)).roicodes_lh = roicodes(idx_lh);
            roi2atlas.(matlab.lang.makeValidName(parc_name)).roicodes_rh = roicodes(idx_rh);
            roi2atlas.(matlab.lang.makeValidName(parc_name)).numVertices_lh = numVertices/2;
            roi2atlas.(matlab.lang.makeValidName(parc_name)).numVertices_rh = numVertices/2;
        else 
            roi2atlas.(matlab.lang.makeValidName(parc_name)).roimat = roimat;
            roi2atlas.(matlab.lang.makeValidName(parc_name)).roinames = roinames;
            roi2atlas.(matlab.lang.makeValidName(parc_name)).roicodes = roicodes;
            roi2atlas.(matlab.lang.makeValidName(parc_name)).numVertices = numVertices;
        end 
        roi2atlas.(matlab.lang.makeValidName(parc_name)).icoNum = icoNum;
        roi2atlas.(matlab.lang.makeValidName(parc_name)).splitLR = splitLR;

        % make sure that all roinames, codes have the same length
        roi2atlas.(matlab.lang.makeValidName(parc_name)).roicodes = reshape(roi2atlas.(matlab.lang.makeValidName(parc_name)).roicodes, length(roi2atlas.(matlab.lang.makeValidName(parc_name)).roicodes), 1);
        roi2atlas.(matlab.lang.makeValidName(parc_name)).roinames = reshape(roi2atlas.(matlab.lang.makeValidName(parc_name)).roinames, length(roi2atlas.(matlab.lang.makeValidName(parc_name)).roinames), 1);

    end
    
    % save output structure
    if save_flag
        save(fname_out, 'roi2atlas');
    end

end

 
