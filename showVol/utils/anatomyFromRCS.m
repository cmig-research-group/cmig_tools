
%% ------------------------------------------------------------
%% --- Anatomical ROI Lookup from RCS coordinates
%% ------------------------------------------------------------

function roi = anatomyFromRCS(handles, type)

if ~handles.hasABCDBrain
  return
end

sf = calculateScale(handles,handles.anat.(type));

prob = squeeze( handles.anat.(type).prob( scaleCoord(handles.rr, sf.r), scaleCoord(handles.cc, sf.c), scaleCoord(handles.ss, sf.s), :) );
[sp,si] = sort(prob,'descend');

switch type
  case 'aseg'
    probThreshold = 0.5;
  case 'fiber'
    probThreshold = 0.1; %fiber probs in general much lower
  case 'aparcaseg'
    probThreshold = 0.5;
end
show = find(sp >= probThreshold);
if isempty(show), show = find(sp >= probThreshold/5); end %FIXME - show at least one ROI as long as it has non-negligible probability
sp = sp(show);
si = si(show);

str = '';
if ~isempty(sp)
  for i = 1:length(sp)
    %FIXME: this has changed a lot
    try
      roiIdx   = find(handles.anat.(type).indlist(si(i)) == handles.anat.(type).roicodes); %as of May 2021, ABCD2
    catch
      try
      roiIdx   = find(handles.anat.(type).filist(si(i)) == handles.anat.(type).roicodes);
      catch
        roiIdx = si(i); %as of Sept 2021
      end
    end
    if isempty(roiIdx), continue, end
    roiProb(i)  = sp(i);
    roinames{i} = handles.anat.(type).roinames{roiIdx};
    roicode(i)  = handles.anat.(type).roicodes(roiIdx);
    %rgb(i,:)    = handles.anat.(type).roicolors(roiIdx,1:3);
    str{i}      = sprintf('  %.2f: %s',roiProb(i), roinames{i});
  end
end
if isempty(str)
  roiProb = 0;
  roinames{1} = 'none';
  roicode = nan;
  %rgb = [0 0 0];
  str = {'  none'};
end
roi.roinames = roinames;
roi.roicode = roicode;
roi.prob = roiProb;
%roi.rgb = rgb;
roi.str = str;