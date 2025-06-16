%% ------------------------------------------------------------
%% -- Go to the given RCS point
%% ------------------------------------------------------------

function handles = gotoPointRCS(handles, targetRCS)
  rr = round(targetRCS(1));
  cc = round(targetRCS(2));
  ss = round(targetRCS(3));
  if ((rr >= 1) & (rr <= handles.rrMax) & ...
      (cc >= 1) & (cc <= handles.ccMax) & ...
      (ss >= 1) & (ss <= handles.ssMax))
    handles.rr = rr;
    handles.cc = cc;
    handles.ss = ss;
    handles = updateDisplay_newvol(handles);
    broadcastLPH(handles)
    printInfo(handles);
  else
    set(handles.edit_command, 'String', ...
      'Outside the volume!');
  end
  
  %% ------------------------------------------------------------