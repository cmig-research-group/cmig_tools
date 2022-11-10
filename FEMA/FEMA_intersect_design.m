function [X,iid, eid, fid,agevec,ymat,contrasts, colnames_model, pihatmat, PregID, HomeID] = FEMA_intersect_design(fname_design, ymat, iid_concat, eid_concat, varargin)
%
% FEMA_intersect_design intersect imaging data and design matrix - used internally by FEMA_wrapper
%
% INPUTS
%   fname_design <cell>        :  cell array with path to file with design matrix saved (readable by readtable) --> if want to batch can add multiple filepaths as separate rows within fname_design as a cell array
%   ymat                       :  matrix of imaging data (n x v)
%   iid_concat                 :  src_subject_id
%   eid_concat                 :  eventname
%
% Optional inputs:
%   contrasts <path>           :  contrast matrix, or path to file containing contrast matrix (readable by readtable)
%   pihat                      :  pihat structure (output from FEMA_process_data)
%   preg                       :  pregnancy IDs
%   address                    :  address IDs
%
% OUTPUTS (inputs for FEMA_fit)
%   X <num>                    :  design matrix (n x p)
%   iid <char>                 :  subject IDs to match imaging data                           
%   eid <char>                 :  eventname
%   fid <num>                  :  family ID (members of the same family unit have same value)
%   agevec <num>               :  participants age
%   ymat <num>                 :  matrix of imaging data (n x v)
%   contrasts <num>            :  contrast vector
%   pihatmat <num>             :  matrix of genetic relatedness --> already intersected to match X and Y sample (empty if no matrix given)
%

%
% This software is Copyright (c) 2021 The Regents of the University of California. All Rights Reserved.
% See LICENSE.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

starttime = now();
logging('***Start***');

p = inputParser;
addParamValue(p,'contrasts',[]);
addParamValue(p,'pihat',[]);
addParamValue(p,'preg',[]);
addParamValue(p,'address',[]);

parse(p,varargin{:})
contrasts = str2num_amd(p.Results.contrasts);
pihat = p.Results.pihat;
preg = p.Results.preg;
address = p.Results.address;

  % Read in design matrix
  logging('Reading design matrix from %s',fname_design);
  tic
  design_mat = readtable(fname_design);
  
  %design matrix must have first four columns: src_subject_id, eventname, rel_family_id, age
  % these are NOT used as IVs, which are in columns 5 onwards
  iid_design = design_mat{:,1};
  eid_design = design_mat{:,2};
  fid_design = design_mat{:,3};
  agevec_design = design_mat{:,4};
  cids = strcat(iid_design,'_',eid_design);
  
  colnames_design = design_mat.Properties.VariableNames;  
  % if one of the modeled IVs is named 'age', that clashes with column 4 and thus will be renamed 'age_1' when loaded. Fix this
  ageIdx = strmatch('age_1',colnames_design,'exact');
  if ~isempty(ageIdx)
    colnames_design{ageIdx} = 'age';
  end
  
  colnums_design = [5:length(colnames_design)]; %remove first four, unmodeled ID columns
  Cmat_design = design_mat{:,colnums_design};
  colnames_Cmat = colnames_design(colnums_design);
  toc
  
  % Merge design matrix with imaging data
  [~, idi, idc] = intersect(strcat(iid_concat(:),'_',eid_concat(:)),cids,'stable'); % IDS_match is loaded from the concatenated .mat
  ymat = ymat(idi,:);
  Cmat = Cmat_design(idc,:);
  iid = iid_design(idc);
  fid = fid_design(idc);
  agevec = agevec_design(idc);
  eid = eid_design(idc,:);
  
  % Make design matrix, add contrasts
  X = double(Cmat);
  colnames_model = colnames_Cmat;
  if ~isempty(contrasts) && size(contrasts,2)<size(X,2)
    contrasts = cat(2,contrasts,zeros(size(contrasts,1),(size(X,2)-size(contrasts,2))));
  end
  %add contrast names to colnames
  if ~isempty(contrasts)
    for c=1:size(contrasts, 1)
      contrast_names{c} = sprintf('%s_%s', colnames_model{find(contrasts(c,:)>0)}, colnames_model{find(contrasts(c,:)<0)});
    end
    colnames_model = cat(2, contrast_names, colnames_model);
  end
  
  
  % Get rid of rows with missing data
  defvec = isfinite(sum(ymat,2)) & isfinite(sum(X,2)); %potentially remove from code and add warning instead - should be NO NaNs?
  
  if any(defvec==0)
    warning(sprintf('NaNs detected in ymat or X (%d rows removed).',length(find(defvec==0))))
  end

if size(ymat,1) > sum(defvec), ymat = ymat(defvec,:); end
if size(X,1) > sum(defvec),   X = X(defvec,:); end
if length(iid) > sum(defvec), iid = iid(defvec); end
if length(fid) > sum(defvec), fid = fid(defvec); end
if length(agevec) > sum(defvec), agevec = agevec(defvec); end
if length(eid) > sum(defvec), eid = eid(defvec); end

      if ~isempty(pihat)
            [iid_list, IA, IC_subj] = unique(iid,'stable'); nsubj = length(iid_list);
            %[fid_list IA IC_fam] = unique(fids,'stable'); nfam = length(fid_list);
            [~, IA, IB_acs] = intersect(iid_list,pihat.iid_list,'stable'); % Why is setdiff(iid_list,tmp_pihat.iid_list) not empty?
            pihatmat = NaN(nsubj,nsubj); pihatmat(IA,IA) = pihat.GRM(IB_acs,IB_acs); % Make genetic relatedness matrix consistent with iid_list
      elseif isempty(pihat)
            pihatmat=[];
      end

% get final list of pregnancy and address IDs
if ~isempty(preg)
  [iid_list, IA, IC_subj] = unique(iid,'stable'); nsubj = length(iid_list);
  [~, IA, IB_acs] = intersect(iid_list,preg.pguid,'stable');
  PregID = NaN([nsubj 1]); PregID(IA) = preg.pregnancyID(IB_acs); % creates PregID in same order as iid_list
elseif isempty(preg)
  PregID = [];
end

if ~isempty(address)
  [iid_list, IA, IC_subj] = unique(iid,'stable'); nsubj = length(iid_list);
  [~, IA, IB_acs] = intersect(iid_list,address.pguid,'stable');
  HomeID(IA) = address.address_id(IB_acs); % creates HomeID in same order as iid_list
  HomeID = transpose(HomeID); 
elseif isempty(address)
  HomeID = [];
end

logging('Final sample for analysis: %d observations',sum(defvec));

end
      
