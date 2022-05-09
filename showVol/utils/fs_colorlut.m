function [roicodes, roinames, rgbv] = fs_colorlut(fname)
% [roicodes roinames rgb] = fs_colorlut(fname)
%
% Purpose: reads a freesurfer color lookup table
%
% Optional input:
%   fname: name of color LUT file
%     {default = $FREESURFER_HOME/FreeSurferColorLUT.txt}
%
% Output:
%  roicodes: vector of numeric ROI codes
%  roinames: cell array of ROI names
%  rgbv: matrix of rgb values for each ROI
%
% Created:  06/02/05 by Doug Greve
% Last Mod: 10/15/12 by Don Hagler
%

% $Id: read_fscolorlut.m,v 1.1 2005/06/02 22:42:34 greve Exp $

roicodes = [];
roinames = [];
rgbv = [];

if ~exist('fname','var') || isempty(fname)
  fname = sprintf('%s/FreeSurferColorLUT.txt',getenv('FREESURFER_HOME'));
end;

if ~exist(fname,'file'), error('file %s not found',fname); end;
fp = fopen(fname,'r');
if(fp == -1) ,error('failed to open file %s',fname); end;

tmp_roinames = '';
nthitem = 1;
while(1)
  tline = fgetl(fp);
  if ~isempty(tline) && tline(1)==-1, break; end; % end of file
  if isempty(deblank(tline)), continue; end; % skip blank lines
  if tline(1)=='#', continue; end; % skip comments

  c = sscanf(tline,'%d',1);
  n = sscanf(tline,'%*d %s',1);
  r = sscanf(tline,'%*d %*s %d',1);
  g = sscanf(tline,'%*d %*s %*d %d',1);
  b = sscanf(tline,'%*d %*s %*d %*d %d',1);
  v = sscanf(tline,'%*d %*s %*d %*d %*d %d',1);
  roicodes(nthitem,1) = c;

  n = char(n);
  if size(n,1)>size(n,2), n=n'; end;
  tmp_roinames = strvcat(tmp_roinames,n);
  if isempty(v), v = 0; end;
  rgbv(nthitem,:) = [r g b v];

  nthitem = nthitem + 1;
end

fclose(fp);

nrois = length(roicodes);
roinames = cell(nrois,1);
for i=1:nrois
  roinames{i} = deblank(tmp_roinames(i,:));
end;
