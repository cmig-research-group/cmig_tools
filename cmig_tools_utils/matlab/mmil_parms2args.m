function args = mmil_parms2args(parms,tags)
%function args = mmil_parms2args(parms,[tags])
% 
% Purpose:
%   Convert from mmil parms object to 
%   a named_argument-value pair list.
%
% Required Parameters:
%   parms: mmil parms object to convert
%    This is just a matlab structure with arbitrary field names
%
% Optional Parameters:
%   tags: cell array of field names to extract from parms
%    If empty, take all fields from parms
%     {default = []}
% 
% See also: mmil_args2parms
%
% Created:  08/01/07 by Ben Cipollini
% Last Mod: 12/31/09 by Don Hagler
%

mmil_check_nargs(nargin, 1);
if ~exist('tags','var'), tags=[]; end;

% initialize args
args   = {};

% get the field names from parms
fields = fieldnames(parms);

% reduce to tags only
if ~isempty(tags)
  fields = intersect(fields,tags);
end;

% append field names and values to args list
for i=1:length(fields)
  args{end+1} = fields{i};
  args{end+1} = getfield(parms,fields{i});
end;

