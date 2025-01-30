function atlasVersionFull = validateAtlasVersion(atlasVersion)
% validateAtlasVersion Create valid atlas version string (x.y_ATLASz_cor10) from various abbreviations of inputs
%
%   atlasVersion = validateAtlasVersion(atlasVersionPartial)
%
% This function will take in various abbreviations and return a full atlas Version in the form x.y_ATLASz_cor10
% where x.y is the ABCD release and ATLASz is the atlas version
%
% * NOTE *: when a new release/atlas appears you will have to do a few things
%       1) Create a new prepareAtlases_... function (e.g. basing from prepareAtlases_50_ABCD3_cor10.m) to generate showVol atlas files
%       2) If the atlas as a different center point, update Mvxl2lph.m and showVol (this was the case for 5.0_ABCD3)
%
%
% examples of partial inputs and full atlasVersions:
% NOTE: ABCD1 and ABCD2 atlases are one-to-one with releases 3.0 and 4.0 and thus did not use the release version to identify
% from release 5.0 and onwards (Feb 2024), the release version is included to disambiguate atlas versions

% 'ABCD1'               -> 'ABCD1_cor10'
% 'ABCD2'               -> 'ABCD2_cor10'
% '5.0' or '5.0_ABCD3'  -> '5.0_ABCD3_cor10'
% '5.1' or '5.1_ABCD3'  -> '5.0_ABCD3_cor10' % NOTE--sub-releases do not change atlas. If that ever becomes untrue, this will need to change
% '6.0' or '6.0_ABCD3'  -> '6.0_ABCD3_cor10'
% note unlike before, 'ABCD3' alone will give an error, since there are now release-dependent versions
%
% JRI (jiversen@mcmaster.ca) Feb 2024

%see what parts of the version string we have
[abcdRelease, atlasVersion, atlasSuffix] = parseAtlasVersion(atlasVersion);

if isempty(abcdRelease) && isempty(atlasVersion)
  error('The  input must contain at least an atlas name (ABCD#) or an ABCD release version (x.y)')
end

%fill in what's missing
if isempty(abcdRelease)
  abcdRelease = atlas2rel(atlasVersion);
end

if isempty(atlasVersion)
  atlasVersion = rel2atlas(abcdRelease);
end

if isempty(atlasSuffix)
  atlasSuffix = 'cor10'; %this is default and seems never to change
end

%reconstitute the full version string. Since used in filenames (including mfiles)
% NOTE: ABCD1 and ABCD2 atlases are one-to-one with releases 3.0 and 4.0 and thus did not use the release version to identify
% from release 5.0 and onwards, the release version is included to disambiguate atlas versions
switch atlasVersion
  case {'ABCD1','ABCD2'}
    atlasVersionFull = [atlasVersion '_' atlasSuffix];
  otherwise
    atlasVersionFull = [abcdRelease(1) '.0_' atlasVersion '_' atlasSuffix];
end

% Helper Funcs
% --------------------------------------------------------------
function atlasVersion = rel2atlas(abcdRelease)
% rel2atlas  map ABCD release to atlas version
%
%   atlasVersion = rel2atlas(abcdRelease)
%
% from Feb 2024 will have version-specific atlases. This function will create a complete atlasVersion string
%  in the form  5.0_ABCD3_cor10 -- i.e. the version prepended to the atlas
%
% abcdRelease is a string in the form x.y
%
% questions--will atlas ever be different between x.0 and x.y releases?

switch(abcdRelease)
  
  case {'6.0'}
    atlasVersion = 'ABCD3';
  
  case {'5.0', '5.1'}
    atlasVersion = 'ABCD3';
    
  case {'4.0'}
    atlasVersion = 'ABCD2';
    
  case {'3.0'}
    atlasVersion = 'ABCD1';
    
  otherwise
    error('%s not a recognized ABCD release',abcdRelease)
    
end

% --------------------------------------------------------------
function abcdRelease = atlas2rel(atlasVersion)
% rel2atlas  map ABCD atlas version to release (when possible)
%
%   abcdRelease = atlas2rel(atlasVersion)
%
% from Feb 2024 will have version-specific atlases. This function will create a complete atlasVersion string
%  in the form  5.0_ABCD3_cor10 -- i.e. the version prepended to the atlas
%
% this is for backwards compatability, to catch any uses of the non-release-qualified atlasVersions
%

switch(atlasVersion(1:5))
  
  case {'ABCD3'}
    %error('ABCD3 atlas has different release-specific versions. Update your code to specify the release in the atlas name, e.g. ''50_ABCD3''')
    abcdRelease = '6.0';

  case {'ABCD2'}
    abcdRelease = '4.0';
    
  case {'ABCD1'}
    abcdRelease = '3.0'; 
    
  otherwise
    error('%s not a recognized atlas version',atlasVersion)
    
end