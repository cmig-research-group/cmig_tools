function atlasFile = showVolAtlasFile(atlasStr)
% showVolAtlasFile  return default showVol atlas file for an atlasVersion
%
% used internally by showVol

atlasStr = validateAtlasVersion(atlasStr); %canonicalize

cfg = abcdConfig('showVol');

atlasFile = fullfile(cfg.data.showVolData,'Atlas',atlasStr, ['showVolAtlases_' atlasStr '.mat']);
        
end
