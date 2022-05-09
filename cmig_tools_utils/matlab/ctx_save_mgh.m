function status = ctx_save_mgh(ctx_vol,fname,mr_parms)
%function status = ctx_save_mgh(ctx_vol,fname,mr_parms)
%
% early mod: 12/11/10 by Don Hagler
% last mod:  12/11/10 by Don Hagler 
%

status = 1;
if (~mmil_check_nargs(nargin,2)) return; end;
if ~exist('mr_parms','var'), mr_parms=[]; end;

[vol,M] = ctx_ctx2mgh(ctx_vol);
fs_save_mgh(vol,fname,M,mr_parms);

status = 0;
