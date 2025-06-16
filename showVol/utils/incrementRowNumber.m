function handles = incrementRowNumber(handles)
    if handles.rr < handles.rrMax
      handles.rr = handles.rr + 1;
    end
    handles = updateDisplay_rr(handles);
    broadcastLPH(handles)
    printInfo(handles);