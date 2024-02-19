function img = maskCSF(img, prob, atlas)
% maskCSF  mask out CSF portions of image volume
%
%   img = maskCSF(img, prob, atlas)
%
%  regions within CSF aseg ROIs > prob are replaced with nan
%  atlas is ABCD1_cor10,  ABCD2_cor10 or ABCD3_cor10 (or just ABCD1, ABCD2 or ABCD3) - As of Dec 2023, ABCD3 by default if not specified
%  prob is optional, defines the threshold level for being 'in' CSF. default is 0.7
%
%   if input img is [] just returns the mask

persistent CSFmask_prob %cache the CSF probability

if ~exist('atlas','var')
  atlas = 'ABCD3_cor10';
  warning('maskCSF; Using ABCD3 atlas by default. Specify atlas if you want something else.')
end
if ~exist('prob','var') || isempty(prob)
  prob = 0.7;
end

% cache mask first time, since involves a file read
if isempty(CSFmask_prob)
  cfg = abcdConfig('showVol');
  disp('Loading aseg to cache CSFmask')
  load(fullfile(cfg.data.showVolData,'Atlas',['showVolAtlases_' atlas '.mat']),'aseg')
  
  CSF_list = {'Left-Lateral-Ventricle','Left-Inf-Lat-Vent','3rd-Ventricle','4th-Ventricle','CSF','Right-Lateral-Ventricle','Right-Inf-Lat-Vent','5th-Ventricle'};
  CSF_idx = contains(aseg.uiNames,CSF_list);
  c = aseg.prob;
  if iscell(c) %decompress
    aseg.prob = zeros(c{1}, 'single');
    aseg.prob(c{2}) = c{3};
  end
  CSFmask_prob = aseg.prob(:,:,:,aseg.uiRoiIdx(CSF_idx));
  CSFmask_prob = sum(CSFmask_prob, 4);
end

CSFmask = (CSFmask_prob>prob);

if ~exist('img','var') || isempty(img)
  img = CSFmask;
else
  img(CSFmask) = nan;
end
