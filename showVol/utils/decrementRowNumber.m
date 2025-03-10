function handles = decrementRowNumber(handles)
    if handles.rr > 1
      handles.rr = handles.rr - 1;
    end
    handles = updateDisplay_rr(handles);
    broadcastLPH(handles)
    printInfo(handles);