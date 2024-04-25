%% ------------------------------------------------------------
%% -- Go to the given LPH point
%% ------------------------------------------------------------

function handles = gotoPointLPH(handles, targetLPH)
Mlph2vxl = inv(handles.Mvxl2lph_atlas);
targetRCS = Mlph2vxl * targetLPH;
handles = gotoPointRCS(handles, targetRCS);

%% ------------------------------------------------------------