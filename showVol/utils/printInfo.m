%% ------------------------------------------------------------
%% --- Print info into text boxes
%% ------------------------------------------------------------
function printInfo(handles)

rr = handles.rr; cc = handles.cc; ss = handles.ss;

if isfield(handles.vols{handles.currentVol},'imgs')
  sf = calculateScale(handles);
  voxVal = handles.vols{handles.currentVol}.imgs(scaleCoord(rr, sf.r), scaleCoord(cc,sf.c), scaleCoord(ss,sf.s), 1);
  CLim1 = handles.CLim(handles.currentVol, 1);
  CLim2 = handles.CLim(handles.currentVol, 2);
  scrVal = (voxVal - CLim1)/(CLim2 - CLim1);
  scrVal = min(scrVal, 1);
  scrVal = max(scrVal, 0);
else %pre-rendered, e.g. FOD, vox val unintersting
  voxVal = 0;
  scrVal = 0;
end

% Convert RCS to LPH
tmp = [rr; cc; ss; 1];
tmp = handles.Mvxl2lph_atlas * tmp;
ll = tmp(1); pp = tmp(2); hh = tmp(3);

set(handles.edit_rr, 'String', sprintf('%d', rr));
set(handles.edit_cc, 'String', sprintf('%d', cc));
set(handles.edit_ss, 'String', sprintf('%d', ss));

set(handles.text_ll, 'String', sprintf('%.3f', ll));
set(handles.text_pp, 'String', sprintf('%.3f', pp));
set(handles.text_hh, 'String', sprintf('%.3f', hh));

set(handles.text_voxval, 'String', sprintf('%.2f', voxVal));
set(handles.text_scrval, 'String', sprintf('%.2f', 100*scrVal));


%% end of function printInfo(handles)
%% ------------------------------------------------------------