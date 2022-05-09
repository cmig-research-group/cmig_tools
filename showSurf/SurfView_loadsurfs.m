%SurfView_loadsurfs  Helper script to load SurfView_surfs.mat
%
% This software is Copyright (c) 2022 The Regents of the University of California. All Rights Reserved.
% See LICENSE.

%by default, included in showSurf repository
mfilepath = fileparts(which(mfilename));
fname_cache = fullfile(mfilepath, 'SurfView_surfs.mat');

if exist(fname_cache,'file')
    load(fname_cache);
else
  error('%s: surf file (%s) not found. Make sure you have the latest showSurf from github.',mfilename,  fname_cache)
  
  % the following is included to document how SurfView_surfs.mat is created, but a default 
  surf_lh_inflated = fs_load_subj_amd('VETSA_average','lh','inflated');
  surf_rh_inflated = fs_load_subj_amd('VETSA_average','rh','inflated');
  surf_lh_inflated.vertices(:,1) = -0.5+surf_lh_inflated.vertices(:,1)-max(surf_lh_inflated.vertices(:,1));
  surf_rh_inflated.vertices(:,1) = 0.5+surf_rh_inflated.vertices(:,1)-min(surf_rh_inflated.vertices(:,1));
  surf_lh_pial = fs_load_subj_amd('VETSA_average','lh','pial');
  surf_rh_pial = fs_load_subj_amd('VETSA_average','rh','pial');
  curvvec_lh = read_curv('~/subjects/VETSA_average/surf/lh.sulc');
  curvvec_rh = read_curv('~/subjects/VETSA_average/surf/rh.sulc');
  label_lh = fs_read_label('~/subjects/VETSA_average/label/lh.cort.label');
  label_rh = fs_read_label('~/subjects/VETSA_average/label/rh.cort.label');
  surfmask_lh = zeros(size(curvvec_lh)); surfmask_lh(label_lh) = 1;
  surfmask_rh = zeros(size(curvvec_rh)); surfmask_rh(label_rh) = 1;
  for i = 0:7 % ic0==ic1==ic2==ic3==ic4, ic5==ic6==ic7 -- 
    tmp = fs_read_trisurf(sprintf('/Applications/freesurfer/lib/bem/ic%d.tri',i));
%    icsurf = struct(); icsurf.vertices = tmp.coords; icsurf.faces = tmp.faces;
   icsurf = tmp;
    icsurfs_orig{i+1} = icsurf;
  end
  icsurfs = cell(size(icsurfs_orig)); icsurfs{end} = icsurfs_orig{end};
  vertindvecs = cell(size(icsurfs_orig)); vertindvecs{end} = [1:size(icsurfs{end}.vertices,1)];
  for i = length(icsurfs):-1:2 % Create lower-order icsurfs by recursive subsampling of ic7
    [icsurfs{i-1} vindvec] = subsample_icsurf(icsurfs{i});
    vertindvecs{i-1} = vertindvecs{i}(vindvec);
  end
  vertind_orig2new = cell(1,length(icsurfs));
  vertind_new2orig = cell(1,length(icsurfs));
  for i = 1:length(icsurfs)
    lookupvec_orig2new{i} = NaN(1,size(icsurfs{i}.vertices,1)); 
    for vi = 1:size(icsurfs{i}.vertices,1)
      dvec = sum((icsurfs_orig{i}.vertices-icsurfs{i}.vertices(vi,:)).^2,2);
      [mv mi] = min(dvec);
      vertind_orig2new{i}(vi) = mi; 
      vertind_new2orig{i}(mi) = vi; 
    end
    fprintf(1,'i=%d: sse=%f\n',i,sum(sum(abs(icsurfs{i}.vertices-icsurfs_orig{i}.vertices(vertind_orig2new{i},:)).^2,2),1));
    fprintf(1,'i=%d: sse=%f\n',i,sum(sum(abs(icsurfs{i}.vertices(vertind_new2orig{i},:)-icsurfs_orig{i}.vertices).^2,2),1));
  end
  save(fname_cache,'surf_lh_inflated','surf_rh_inflated','surf_lh_pial','surf_rh_pial','curvvec_lh','curvvec_rh','surfmask_lh','surfmask_rh','icsurfs','icsurfs_orig','vertindvecs','vertind_orig2new','vertind_new2orig');
end

% ToDo
%   Update, and make sure it runs in current environment

