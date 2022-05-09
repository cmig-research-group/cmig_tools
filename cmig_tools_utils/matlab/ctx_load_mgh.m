function [ctx_vol,mr_parms,volsz] = ctx_load_mgh(fname)
%function [ctx_vol,mr_parms,volsz] = ctx_load_mgh(fname)
%
% early mod: 02/29/08 by Don Hagler
% last mod:  12/11/10 by Don Hagler
%

if (~mmil_check_nargs(nargin,1)) return; end;

[vol,M,mr_parms,volsz] = fs_load_mgh(fname,[],1); % first frame only
ctx_vol = ctx_mgh2ctx(vol,M);

