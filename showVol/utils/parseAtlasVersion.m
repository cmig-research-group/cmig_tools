function [abcdRelease, atlasVersion, atlasSuffix] = parseAtlasVersion(atlasVersion)
% parseAtlasVersion companion to validateAtlasVersion splits atlas into abcdVersion, atlasVersion and atlasSuffix
%
%  [abcdRelease, atlasVersion, atlasSuffix] = parseAtlasVersion(atlasVersion)
%
%  e.g. 5.0_ABCD3_cor10 -> '5.0' 'ABCD3' 'cor10'

parts = regexp(atlasVersion,'(?<ver>\d\.?\d)?_?(?<atlas>ABCD\d)?_?(?<suffix>.+)?','tokens');
abcdRelease = parts{1}{1};
atlasVersion = parts{1}{2};
atlasSuffix = parts{1}{3};