%% ------------------------------------------------------------
%% --- Draw all Annotations
%% ------------------------------------------------------------
function handles = annotationDraw(handles,scale)

if ~isfield(handles,'ORIENTATION'), return, end %this gets called once before init is complete

% draw labels on big axis, dits on smaller ones
delete(findobj(handles.figure,'tag','annotation'))

if isempty(handles.annotation) || ~handles.annotation.visible || isempty(handles.annotation.A), return, end

%find unhidden annotations
A = handles.annotation.A;
handles.annotation.uuid.inView = {}; %keep track of on-screen visible handles

for iView = 1:3 %loop through three views
  switch iView
    case 1
      drawax = handles.axesCR;
      rr = [1 handles.rrMax];
      cc = [1 handles.ccMax];
      ss = handles.ss;
      scale = size(handles.imCR.CData,1)/handles.rrMax;
    case 2
      drawax = handles.axesSR;
      rr = [1 handles.rrMax];
      cc = handles.cc;
      ss = [1 handles.ssMax];
      scale = size(handles.imSR.CData,1)/handles.rrMax;
    case 3
      drawax = handles.axesCS;
      rr = handles.rr;
      cc = [1 handles.ccMax];
      ss = [1 handles.ssMax];
      scale = size(handles.imCS.CData,2)/handles.ccMax;
  end
  doLabel = iView == handles.ORIENTATION;
  
  inView = annotationsAtPoint(A,rr,cc,ss);
  if ~any(inView), continue, end
  As = A(inView,:);
  
  handles.annotation.uuid.inView = union(handles.annotation.uuid.inView, As.uuid);
  
  [~,iVisible] = setdiff(As.uuid, handles.annotation.uuid.hidden);
  As = As(iVisible,:);
  
  authentic = annotationsAtPoint(As,rr,cc,ss,0); %points actually in this slice (not extended)
  
  %draw the annotations
  len = 3 * scale / 2^handles.zoom;
  inset = .25;
  nA = height(As);
  displayed = false(nA,1);
  
  % loop over annotations in view
  for iA = 1:nA
    
    if displayed(iA), continue, end %skip if label already was displayed
    
    r0 = As.r(iA); c0 = As.c(iA); s0 = As.s(iA);
    switch iView
      case 1 %CR
        x = scale*c0 + [inset  len];
        y = scale*r0 + [inset  len];
      case 2 %SR
        x = scale*s0 + [inset  len];
        y = scale*r0 + [inset  len];
      case 3 %CS
        x = scale*c0 + [inset  len];
        y = scale*s0 + [inset  -len];
    end
    
    isFOD = isfield(handles.vols{handles.currentVol},'imgs1');
    
    if doLabel && ~isFOD% DRAW LABEL in main axis
      %prepare label
      label = As.abbrev{iA};
      if isempty(label)
        label = As.label{iA}(1:min(length(label),4));
      end
      if handles.annotation.showAuthor, label = sprintf('{%s}_{%s}', label, As.author{iA}(1:2)); end
      %add any additional annotations at this same point to the label
      atCoord = find(annotationsAtPoint(As, r0, c0, s0) );
      atCoord(atCoord == iA) = [];
      if ~isempty(atCoord)
        for iC = atCoord'
          tmp = As.abbrev{iC};
          if isempty(tmp)
            tmp = As.label{iC}(1:min(length(label),4));
          end
          if handles.annotation.showAuthor, tmp = sprintf('{%s}_{%s}', tmp, As.author{iC}(1:2)); end
          label = [label ';' tmp];
          displayed(iC)=true;
        end
      end
      
      if authentic(iA)
        color = 'w';
      else
        color = [.65 .65 .65];
      end
      
      h = [ outlinetext(handles.axes1, x(2), y(2), label, color, 'fontsize',11,'fontweight','bold','horizontalalignment','left','tag','annotation') ...
        line(handles.axes1, x, y, 'color','k','linewidth',3,'tag','annotation')...
        line(handles.axes1, x, y, 'color',color,'linewidth',1,'tag','annotation') ];
      
      % TRY: add invisible ui element with tooltip containing description of label
      % cool: replace button with screengrab of image below it, making it invisible
      % But it is very slow, kills interactivity
      %     tipstr = sprintf('%s (%s)  [%s]\n\n%s',As.label, As.abbrev, As.author(1:2), As.note);
      %     pixpos = get(h(5),'position'); %text label
      %     ax0 = getpixelposition(handles.axes1,true);
      %     axpos = [pixpos(1) pixpos(2)-10 25 20];
      %     figpos = axpos + [ax0(1) ax0(2) 0 0];
      %     shot = getframe(handles.axes1, axpos);
      %     h2 = uicontrol(handles.axes1.Parent,'style','pushbutton','String','','foregroundColor',[1 0 0],...
      %         'Tooltip',tipstr,'visible','on','unit','pixels','position',figpos,'tag','annotation',...
      %         'CData',shot.cdata);
      
    else % DRAW A DOT in small axes
      atCoord = find( annotationsAtPoint(As, r0, c0, s0) );
      atCoord(atCoord == iA) = [];
      displayed(atCoord) = true; %so only draw the first marker at a point with multiple
      
      if authentic(iA)
        color = 'w';
      else
        color = [.75 .75 .75];
      end
      
      h =  [ line(drawax, [x(1) x(1)]-inset, [y(1) y(1)]-inset,'marker','o', ...
        'markerfacecolor','k','markeredgecolor','k', 'markersize',6,'tag','annotation'), ...
        line(drawax, [x(1) x(1)]-inset, [y(1) y(1)]-inset,'marker','o', ...
        'markerfacecolor',color,'markeredgecolor','none', 'markersize',3,'tag','annotation') ];
    end
  end %annotations in slice
end %view loop

% notify annotationUI of changes
updateAnnotationUI(handles)

%% ------------------------------------------------------------