function home = userhome
% userhome  user's home directory

if ~ispc
  home = getenv('HOME'); %char(java.lang.System.getProperty('user.home')); %supposed to be x-platform but buggy?
else
  home = getenv('USERPROFILE'); % FIXME: untested on PC
end