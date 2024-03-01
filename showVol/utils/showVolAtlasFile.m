function atlasFile = showVolAtlasFile(atlasVersion)
% showVolAtlasFile  return default showVol atlas file for an atlasVersion
%
% used internally by showVol. Needs to be updated when new atlas versions are created

atlasVersion = validateAtlasVersion(atlasVersion); %canonicalize

cfg = abcdConfig('showVol');

switch atlasVersion
  
  case '3.0_ABCD1_cor10'
    atlasFile = fullfile(cfg.data.showVolData,'Atlas','showVolAtlases_ABCD1_cor10.mat');
    
  case '4.0_ABCD2_cor10'
    atlasFile = fullfile(cfg.data.showVolData,'Atlas','showVolAtlases_ABCD2_cor10.mat');
    
  case '5.0_ABCD3_cor10'
    atlasFile = fullfile(cfg.data.showVolData,'Atlas','5.0_ABCD3_cor10', 'showVolAtlases_5.0_ABCD3_cor10.mat');
    
  otherwise
    error('Unrecognized atlas version (%s)',atlasVersion)
        
end
