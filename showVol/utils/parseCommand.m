%% ------------------------------------------------------------
%% --- Parse string given in command edit box
%% ------------------------------------------------------------

function handles = parseCommand(handles, str)

%% Empty -- do nothing
str = deblank(str);
if isempty(str)
  return
end

%% Jump to slice/col/row?
num = str2double(str);
if isfinite(num)
  switch handles.ORIENTATION
    case 1,
      ss = round(num);
      if (ss >= 1) & (ss <= handles.ssMax)
        handles.ss = ss;
        handles = updateDisplay_ss(handles);
        printInfo(handles);
      end
    case 2,
      cc = round(num);
      if (cc >= 1) & (cc <= handles.ccMax)
        handles.cc = cc;
        handles = updateDisplay_cc(handles);
        printInfo(handles);
      end
    case 3,
      rr = round(num);
      if (rr >= 1) & (rr <= handles.rrMax)
        handles.rr = rr;
        handles = updateDisplay_rr(handles);
        printInfo(handles);
      end
  end
  return;
end

%% Possible command in str
[strFirst, strRest] = strtok(str);
switch strFirst
  
  %% ------------------------------------------------------------
  case 'o'
    %% ------------------------------------------------------------
    
    handles = cycleOrientationsForward(handles);
    
    %% ------------------------------------------------------------
  case {'f', 'fly'}
    %% ------------------------------------------------------------
    
    WRONG_FORMAT = 1;
    if ~isempty(strRest)
      [strStopIndex, strRest] = strtok(strRest);
      stopIndex = round(str2double(strStopIndex));
      if isfinite(stopIndex)
        if isempty(strRest)
          handles = fly(handles, stopIndex);
          WRONG_FORMAT = 0;
        else
          [strPauseTime, strRest] = strtok(strRest);
          pauseTime = str2double(strPauseTime);
          if isfinite(pauseTime)
            if pauseTime >= 0
              handles = fly(handles, stopIndex, pauseTime);
              WRONG_FORMAT = 0;
            end
          end
        end
      end
    end
    if WRONG_FORMAT
      set(handles.edit_command, 'String', ...
        'fly stopIndex [pauseTime]');
    end
    
    %% ------------------------------------------------------------
  case {'c', 'cycle'}
    %% ------------------------------------------------------------
    
    WRONG_FORMAT = 1;
    numLoops = 10;
    pauseTime = 0.1;
    if isempty(strRest)
      WRONG_FORMAT = 0;
      handles = cycleVolumesMovie(handles, numLoops, pauseTime);
    else
      [strNumLoops, strRest] = strtok(strRest);
      numLoops = round(str2double(strNumLoops));
      if numLoops > 0
        if isempty(strRest)
          WRONG_FORMAT = 0;
          handles = cycleVolumesMovie(handles, numLoops, pauseTime);
        else
          [strPauseTime, strRest] = strtok(strRest);
          pauseTime = str2double(strPauseTime);
          if pauseTime >= 0
            WRONG_FORMAT = 0;
            handles = cycleVolumesMovie(handles, numLoops, pauseTime);
          end
        end
      end
    end
    if WRONG_FORMAT
      set(handles.edit_command, 'String', ...
        'cycle [numLoops] [pauseTime]');
    end
    
    %% ------------------------------------------------------------
  case 'lph'
    %% ------------------------------------------------------------
    
    WRONG_FORMAT = 1;
    if ~isempty(strRest)
      [strL, strRest] = strtok(strRest);
      targetL = str2double(strL);
      if isfinite(targetL)
        if ~isempty(strRest)
          [strP, strRest] = strtok(strRest);
          targetP = str2double(strP);
          if isfinite(targetP)
            [strH, strRest] = strtok(strRest);
            targetH = str2double(strH);
            if isfinite(targetH)
              targetLPH = [targetL; targetP; targetH; 1];
              handles = gotoPointLPH(handles, targetLPH);
              WRONG_FORMAT = 0;
            end
          end
        end
      end
    end
    if WRONG_FORMAT
      set(handles.edit_command, 'String', ...
        'lph gotoL gotoP gotoH');
    end
    
    %% ------------------------------------------------------------
  case 'rcs'
    %% ------------------------------------------------------------
    
    WRONG_FORMAT = 1;
    if ~isempty(strRest)
      [strR, strRest] = strtok(strRest);
      targetR = str2double(strR);
      if isfinite(targetR)
        if ~isempty(strRest)
          [strC, strRest] = strtok(strRest);
          targetC = str2double(strC);
          if isfinite(targetC)
            [strS, strRest] = strtok(strRest);
            targetS = str2double(strS);
            if isfinite(targetS)
              targetRCS = [targetR; targetC; targetS; 1];
              handles = gotoPointRCS(handles, targetRCS);
              WRONG_FORMAT = 0;
            end
          end
        end
      end
    end
    if WRONG_FORMAT
      set(handles.edit_command, 'String', ...
        'rcs gotoR gotoC gotoS');
    end
    
end % switch strFirst

%% end of function handles = parseCommand(handles, str)
%% ------------------------------------------------------------