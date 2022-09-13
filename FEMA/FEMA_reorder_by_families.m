function [X,iid,eid,fid,agevec,ymat,pihatmat]=FEMA_reorder_by_families(X,iid,eid,fid,agevec,ymat,pihatmat)
  % change the order of observations so that members of the same fid form a
  % consecutive block; this is done to turn all covariance matrices into
  % block-diagonal shape.
  % 
  % First dimension of X, iid, eid, fid, agevec, ymat corresponds to
  % observations, i.e. shape(X, 1) == nobs (number of observations)
  %
  % pihatmat is different - it's a square matrix of size nsubj x nsubj,
  % potentially sparse. The order of subjects in pihatmat corresponds to
  % the order in which subjects appear in iid vec, that is to IA from
  %   [~, IA] = unique(iid, 'stable');
  % As such after shuffling observations to group families we also need to 
  % re-order pihatmat accordingly.

  % step1 - find how to reorder observations
  jvec=zeros(length(fid), 1);
  [fid_list, ~, IC_fam] = unique(fid,'stable'); nfam = length(fid_list);
  i=1;
  for fi = 1:nfam
    jvec_fam = rowvec(find(IC_fam==fi)); % Identify all observations for a given family
    jvec(i:(i+length(jvec_fam) - 1)) = jvec_fam;
    i=i+length(jvec_fam);
  end

  % step2 - reordering pihatmat is a tricky stuff, see comment above
  % note this step needs to be done before step 3 because that changes ordering in iid.
  [~, IA] = ismember( unique(iid(jvec), 'stable'), unique(iid, 'stable'));
  pihatmat = pihatmat(IA, IA);

  % step 3 - reordering observations (X, iid, eid, fid, agevec and ymat)
  if ~isempty(X), X = X(jvec, :); end
  if ~isempty(iid), iid = iid(jvec); end
  if ~isempty(eid), eid = eid(jvec); end
  if ~isempty(fid), fid = fid(jvec); end
  if ~isempty(agevec), agevec = agevec(jvec); end
  if ~isempty(ymat), ymat = ymat(jvec, :); end
end
