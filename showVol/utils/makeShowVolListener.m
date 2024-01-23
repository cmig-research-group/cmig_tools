function makeShowVolListener(figh, funch, context)
% turn a figure intoo a showVol listener
%
%    makeShowVolListener(figh, funch, context)
%
% figure figh will receive [l p h] coordinates whenever the user changes the cursor
%   in any showVol window
%
% funch is pointer @func to a function that will be alerted of coordinate changes
%   it can do whatever it wants with the coordinate, such as plot a voxel timeseries, 
%   group stats for that voxel, seed to voxel correlation, etc. Default calls a test function
%   embedded in the code below
%
% context is a struct with whatever information the callback needs to do its thing. For example, if showVol
%      was displaying a seed-wholebrain correlation maps, the context might include the whole brain and seed
%      timeseries and perhaps some labels. The callback would get the RCS, pick the voxel timeseries and
%      then plot the seed and voxel timeseries or a scatterplot to inspect the data underlying the correlation map
%
%
% basic outline of a callback: 
% function func(~,event)
%
% destFig = event.AffectedObject
% handles = guidata(destFig); %calling makeShowVolListener adds context to guidata of this figure
% context = handles.context;
%
% ud = destFig.UserData;        %showVol passes coordinates by setting this figure's UserData
% rcs = floor(inv(ud.M)*[ud.lph(:); 1]);
% rr = rcs(1); cc = rcs(2); ss = rcs(3);
% fprintf('Figure %d Received RCS coords from showVol: [%d %d %d]\n',destFig.Number,rr,cc,ss);


if ~exist('funch','var') | isempty(funch)
  funch = @testfunc;
end

if ~exist('context','var')
  context = [];
end

handles = guidata(figh);
handles.figure = figh;
handles.CoordinateListener = addlistener(handles.figure,'UserData', 'PostSet', funch);
handles.context = context;
guidata(handles.figure, handles)
figh.Tag = 'showVol_listener';      %this tag is how showVol knows to send this figure the coordinate data


%=======================================================
% demo callback function: simply reports received RCS

function testfunc(~,event)
destFig=event.AffectedObject;
ud=destFig.UserData;
rcs = floor(inv(ud.M)*[ud.lph(:); 1]);
rr = rcs(1); cc = rcs(2); ss = rcs(3);
fprintf('Figure %d Received RCS coords from showVol: [%d %d %d]\n',destFig.Number,rr,cc,ss);
