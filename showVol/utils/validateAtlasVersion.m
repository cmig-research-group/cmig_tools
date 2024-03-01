function atlasVersionFull = validateAtlasVersion(atlasVersion)
% validateAtlasVersion create valid atlas version string (x.y_ATLASz_cor10) from various abbreviations of inputs
%
%   atlasVersion = validateAtlasVersion(atlasVersionPartial)
%
% From Feb 2024 we are adding release to atlas version identifier as there are now multiple release-specific atlases. 
% This function will take in various abbreviations and return a full atlas Version in the form xy_ATLASz_cor10
% we'll have to one-time rename a bunch of files in showVolData. Note no decimal in version 5.0 -> 50 for file naming purposes.
%
% NOTE: This file, but only this, will need updating when new releases/atlases appear
%         Everything else using atlases will call this function
%
%
% examples of partial inputs and full atlasVersions:
% 'ABCD1'               -> '30_ABCD1_cor10'
% 'ABCD2'               -> '40_ABCD2_cor10'
% '5.0' or '5.0_ABCD3'  -> '50_ABCD3_cor10'
% '5.1' or '5.1_ABCD3'  -> '50_ABCD3_cor10' % NOTE--sub-releases do not change atlas. If that ever becomes untrue, this will need to change
% '6.0' or '6.0_ABCD3'  -> '60_ABCD3_cor10'
% note unlike before, 'ABCD3' alone will give an error, since there are now release-dependent versions
%
% JRI (jiversen@mcmaster.ca) Feb 2024

%see what parts of the version string we have
parts = regexp(atlasVersion,'(?<ver>\d\.?\d)?_?(?<atlas>ABCD\d)?_?(?<suffix>.+)?','tokens');
abcdRelease = parts{1}{1};
atlasVersion = parts{1}{2};
atlasSuffix = parts{1}{3};

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

%reconstitute the full version string. Since used in filenames (including mfiles), we render 5.0 as 50 as matlab gets funny with additional '.' in filenames
atlasVersionFull = [abcdRelease(1) '.0_' atlasVersion '_' atlasSuffix];


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
% rel2atlas  map ABCD release to atlas version
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
    error('ABCD3 atlas has different release-specific versions. Update your code to specify the release in the atlas name, e.g. ''50_ABCD3''')
    
  case {'ABCD2'}
    abcdRelease = '4.0';
    
  case {'ABCD1'}
    abcdRelease = '3.0'; 
    
  otherwise
    error('%s not a recognized atlas version',atlasVersion)
    
end