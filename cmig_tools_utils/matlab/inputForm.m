function s = inputForm(A)

if length(A)==0
  s = '[]';
  return;
end

if isnumeric(A)
  s = num2str(A);
  if length(A)>1
      s = ['[' s ']'];
  end
  return;
end

if ischar(A)
  s = sprintf('''%s''',A);
  return;
end

if iscell(A)
  s = '{';
elseif isnumeric(A)
  s = '[';
end

if length(A)==1 & iscell(A)
      s = ['{' sprintf('''%s''',char(A)) '}'];
      return
elseif length(A) > 1 & iscell(A)
  for fi = 1:length(A)
    s = [s ' ' inputForm(A{fi})];
  end
end

if iscell(A)
  s = [s '}'];
elseif isnumeric(A)
  s = [s ']'];
end
