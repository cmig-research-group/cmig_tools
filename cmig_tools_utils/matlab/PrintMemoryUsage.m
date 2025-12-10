pid = feature('getpid');
cmd = sprintf('grep VmPeak /proc/%d/status',pid); % Also get CPU time, I/O, etc.?
[s r]  = jsystem(cmd);
ind = findstr('VmPeak:',r);
VmPeak = NaN;
if ~isempty(ind)
  r = r(ind:end); 
  [dummy remain] = strtok(r); [VmPeak remain] = strtok(remain); [units remain] = strtok(remain);
  VmPeak = str2num(VmPeak);
end
if strcmpi('kb',units), VmPeak = VmPeak/1e6; end
if strcmpi('mb',units), VmPeak = VmPeak/1e3; end
cmd = sprintf('grep VmSize /proc/%d/status',pid); % Also get CPU time, I/O, etc.?
[s r]  = jsystem(cmd);
ind = findstr('VmSize:',r);
VmSize = NaN;
if ~isempty(ind)
  r = r(ind:end); 
  [dummy remain] = strtok(r); [VmSize remain] = strtok(remain); [units remain] = strtok(remain);
  VmSize = str2num(VmSize); ind = findstr('VmSize:',r);
end
if strcmpi('kb',units), VmSize = VmSize/1e6; end
if strcmpi('mb',units), VmSize = VmSize/1e3; end
fprintf(1,'VmSize: %0.1fGB VmPeak: %0.1fGB\n', VmSize,VmPeak);

% ToDos
%   Add support for macOS: using top -l 1 ?

