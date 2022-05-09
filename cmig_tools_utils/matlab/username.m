function name = username
% username  return current username

if ~ispc
  name = getenv('USER');
else
  name = getenv('USERNAME'); % FIXME: untested on PC
end