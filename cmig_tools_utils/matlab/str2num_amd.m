function X = str2num_amd(S)

if ~iscell(S)
  if ischar(S)
    [X OK] = str2num(S);
    if ~OK, X = NaN; end
  else
    X = S;
  end
else
  X = NaN(size(S));
  for i = 1:prod(size(X));
    if ischar(S{i})
      X(i) = str2num_amd(S{i});
    else
      X(i) = S{i};
    end
  end
end

