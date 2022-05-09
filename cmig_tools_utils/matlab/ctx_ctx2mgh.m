function [vol,M] = ctx_ctx2mgh(vol_ctx)

vol = single(vol_ctx.imgs);
M = M_LPH_TO_RAS * vol_ctx.Mvxl2lph;
