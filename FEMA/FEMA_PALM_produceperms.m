function [perms, VG, EB, iid_EB, nperms] = FEMA_PALM_produceperms(iid, fam_id, agevec, nperms, varargin)

      %addpath('~/matlab/ABCD_PALM/PALM-master')

      p = inputParser;
      addParamValue(p,'pihatmat',[]);
      addParamValue(p,'save_perms',1);
      parse(p,varargin{:})
      pihatmat = p.Results.pihatmat;
      save_perms = p.Results.save_perms;


      [clusterinfo, EB, iid_EB]=FEMA_parse_family(iid,fam_id,agevec,pihatmat);

      EB = palm_reindex(EB,'fixleaves');
      
      %To view EB tree --> need to use subsample because full ABCD_sample is too large
      %M = (1:size(EB,1))';
      %Ptree = palm_tree(EB,M');
      %palm_ptree2dot(Ptree,'ptree.dot'); %for visualisation of EBs


      %Use this to generate permutation shufflings

      nperms=nperms+1;

      [perms, VG]=palm_quickperms([],EB,nperms);

      if save_perms==1
            if ~exist('~/matlab/FEMA_perms')
                  mkdir('~/matlab/FEMA_perms')
            end
            save(sprintf('~/matlab/FEMA_perms/perms_nsubj%d_nperms%d.mat',size(perms,1),nperms),'perms','VG','EB','iid_EB');
      end


end