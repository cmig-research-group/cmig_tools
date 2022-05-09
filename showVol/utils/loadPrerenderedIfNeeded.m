function loadPrerenderedIfNeeded(atlas, fodtype)
% loadPrerenderedIfNeeded  cache prerendered images used by showVol
%
%  loadPrerenderedIfNeeded(atlas)
%
%  loads prerendered volumes into a global variable PRI
%
%     by default loads ABCD2 images
%
% PRI.ABCD2
%               FOD:  % voxelwise normalized FOD (emphasizing gray matter FOD)
%           FOD_p80:  % globally normalized (80%)
%                T1:  % mean T1 1mm
%                T2:  % mean T2 1mm
%                FA:  % Fractional anisotropy
%                CO:  % mean fiber orientation rgb
%          aseg_rgb:  % aseg atlas rgb
%         fiber_rgb:  % fiber atlas rgb
%     aparcaseg_rgb:  % aparcaseg atlas rgb
%
% FOD types
%   ABCD2 has different FOD normalizations: voxelwise, and global
%   this can be specified by a second argument
%
%       function loadPrerenderedIfNeeded(atlas, fodtype)
%         fodtype = 'v';   %voxelwise, emphasizes gray matter FOD
%         fodtype = 'p80'; %global 80% normlization, emphasizes white matter FOD
%         fodtype = 'p50'; %global 50% normlization 
%   multiple types can be specified in a cell array
%
%
% Backwards compatibility with older atlases
%
%     can also load ABCD1 atlas, or both
%   loadPrerenderedIfNeeded('ABCD1')
%   loadPrerenderedIfNeeded({'ABCD1','ABCD2'})
%
% PRI.ABCD1
%       FOD_1mm    % normalized FOD (emphasizing gray matter FOD)
%      aseg_rgb    % aseg atlas 1mm
%     fiber_rgb    % fiber atlas 1mm
%            CO    % mean fiber orientation rgb
%            T1    % mean T1 1mm
%            T2    % mean T2 1mm
%


cfg = abcdConfig('showVol');

if ~exist('atlas','var') || isempty(atlas)
  atlas = 'ABCD2';
  fprintf('%s: Using ABCD2_cor10 atlas by default. Specify atlas explicitly if you want something else.\n',mfilename)
end

if ~exist('fodtype','var')
  fodtype = {''}; %no suffix means voxel-wise normalized
  fodsuffix = {''};
else
  if ~iscell(fodtype), fodtype = {fodtype}; end
  fodsuffix = cellfun(@(x) ['_' x], fodtype,'UniformOutput',false);
  fodsuffix = strrep(fodsuffix,'_v','');
end

global PRI

% Load if needed

if contains(atlas,'ABCD1')
  if ~exist('PRI','var') || isempty(PRI) || ~isfield(PRI,'vol_fiber_rgb')
    fprintf('%s -- %s.m: Cacheing pre-rendered images for %s (slow, but only have to do it once per session)\r\n', datestr(now), mfilename, atlas);
    
    
    pre_rendered_path = fullfile(cfg.data.showVolData,'atlas_dspace','rgbsnap_ABCD1_cor10.mat');
    FOD_path = fullfile(cfg.data.showVolData,'FOD','FOD_1mm_c.mat');
    if ~exist(FOD_path,'file')
      FOD_path = fullfile(cfg.data.showVolData,'FOD','FOD_ABCD1_cor10_c.mat');
    end
    
    fprintf('%s -- %s.m: Loading FOD images (May take a few minutes)\r\n', datestr(now), mfilename);
    
    load(FOD_path,'vol_FOD_1mm')
    
    vol_FOD_1mm.imgs3 = vol_FOD_1mm.imgsRC; %quick hack--rename
    vol_FOD_1mm.imgs2 = vol_FOD_1mm.imgsRS;
    vol_FOD_1mm.imgs1 = vol_FOD_1mm.imgsSC;
    vol_FOD_1mm = rmfield(vol_FOD_1mm,{'imgsRC','imgsRS','imgsSC'});
    vol_FOD_1mm.name = 'FOD 1mm (normalized)';
    PRI.vol_FOD_1mm = vol_FOD_1mm;
    
    
    fprintf('%s -- %s.m: Loading prerendered images \r\n', datestr(now), mfilename);
    tmp = load(pre_rendered_path);
    
    %convert to 0-based centering
    Mvxl2lph = M_RAS_TO_LPH*tmp.M_atl; %M_atl is in RAS, apparently
    Mvxl2lph(1,4) = -100;
    Mvxl2lph(2,4) = 130;
    Mvxl2lph(3,4) = 100;
    
    
    PRI.vol_aseg_rgb = tmp.vol_aseg_rgb;
    PRI.vol_aseg_rgb.Mvxl2lph = Mvxl2lph;
    PRI.vol_aseg_rgb.name = 'aseg';
    
    PRI.vol_fiber_rgb = tmp.vol_fiber_rgb;
    PRI.vol_fiber_rgb.Mvxl2lph = Mvxl2lph;
    PRI.vol_fiber_rgb.name = 'fiber atlas';
    
    PRI.vol_CO = ctx_mgh2ctx(tmp.vol_CO_mu, M_LPH_TO_RAS*Mvxl2lph); %this assumes M is in RAS
    PRI.vol_CO.name = 'diffusion direction';
    
    PRI.vol_T1 = ctx_mgh2ctx(tmp.vol_T1_mu, M_LPH_TO_RAS*Mvxl2lph);
    PRI.vol_T1.name = 'T1';
    
    PRI.vol_T2 = ctx_mgh2ctx(tmp.vol_T2_mu, M_LPH_TO_RAS*Mvxl2lph);
    PRI.vol_T2.name = 'T2';
    
    %JRI preference--add ghost of T1 to fiber atlas image for
    % orientation help
    maxFV = max(PRI.vol_fiber_rgb.imgs(:));
    maxT1 = max(PRI.vol_T1.imgs(:));
    PRI.vol_fiber_rgb.imgs = PRI.vol_fiber_rgb.imgs + PRI.vol_T1.imgs * (0.15 * maxFV / maxT1);
    
    % for forward compatibility, (also) put them in ABCD1 field, ditch the 'vol'
    if 1
      PRI.ABCD1.FOD = PRI.vol_FOD_1mm;
      PRI.ABCD1.aseg_rgb = PRI.vol_aseg_rgb;
      PRI.ABCD1.fiber_rgb = PRI.vol_fiber_rgb;
      PRI.ABCD1.CO = PRI.vol_CO;
      PRI.ABCD1.T1 = PRI.vol_T1;
      PRI.ABCD1.T2 = PRI.vol_T2;
      
      PRI = rmfield(PRI,{'vol_FOD_1mm','vol_aseg_rgb','vol_fiber_rgb','vol_CO','vol_T1','vol_T2'});
    end
    
  else
    fprintf('%s -- %s.m: Using cached prerendered images for atlas %s (''clear global PRI'' to force reloading)\r\n', datestr(now), mfilename, atlas);
  end
end

if contains(atlas,'ABCD2')
  if ~exist('PRI','var') || isempty(PRI) || ~isfield(PRI,'ABCD2')
    fprintf('%s -- %s.m: Cacheing pre-rendered images for %s (slow, but only have to do it once per session)\r\n', datestr(now), mfilename, atlas);
    
    for iF = 1:length(fodsuffix)
      FOD_path = fullfile(cfg.data.showVolData,'FOD',['FOD_ABCD2_cor10' fodsuffix{iF} '_c.mat']);
      if ~exist(FOD_path,'file'), disp(['Unable to load FOD. Missing ' FOD_path '!']), continue, end
      
      fprintf('%s -- %s.m: Loading FOD%s images (May take a few minutes)\r\n', datestr(now), mfilename,fodsuffix{iF});
      load(FOD_path) %
      vol_FOD_1mm.imgs3 = vol_FOD_1mm.imgsRC; %quick hack--rename
      vol_FOD_1mm.imgs2 = vol_FOD_1mm.imgsRS;
      vol_FOD_1mm.imgs1 = vol_FOD_1mm.imgsSC;
      vol_FOD_1mm = rmfield(vol_FOD_1mm,{'imgsRC','imgsRS','imgsSC'});
      FODname = ['FOD' fodsuffix{iF}];
      disp(FODname)
      if isempty(fodsuffix{iF})
        vol_FOD_1mm.name = 'FOD [ABCD2] voxelwise normalized'; %temporary FIXME
      end
      PRI.ABCD2.(FODname) = vol_FOD_1mm;
    end
    
    fprintf('%s -- %s.m: Loading prerendered images \r\n', datestr(now), mfilename);

    vols = {'T1','T2','FA','CO','aseg_rgb','fiber_rgb','aparcaseg_rgb'}; %FIXME: removed wmparc
    
    for iV = 1:length(vols)
      varname = [vols{iV} '_ABCD2_cor10'];
      fname = fullfile(cfg.data.showVolData,'Atlas',[varname '.mat']);
      try
        tmp = load(fname);
        fieldname = [vols{iV} '_ABCD2_cor10'];
        PRI.ABCD2.(vols{iV}) = tmp.(fieldname);
      catch
        disp(['problem with ' vols{iV}])
      end
    end
    
    %JRI preference--add ghost of T1 to fiber atlas image for
    % orientation help
    maxFV = max(PRI.ABCD2.fiber_rgb.imgs(:));
    maxT1 = max(PRI.ABCD2.T1.imgs(:));
    PRI.ABCD2.fiber_rgb.imgs = PRI.ABCD2.fiber_rgb.imgs + PRI.ABCD2.T1.imgs * (0.15 * maxFV / maxT1);
    
  else
      fprintf('%s -- %s.m: Using cached prerendered images for atlas %s (''clear global PRI'' to force reloading)\r\n', datestr(now), mfilename, atlas);
  end
end