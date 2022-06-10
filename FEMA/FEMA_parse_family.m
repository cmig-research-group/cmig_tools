function [clusterinfo, Ss, iid, famtypevec, famtypelist, subj_famtypevec]=FEMA_parse_family(iid,eid,fid,agevec,pihatmat,varargin)

% Function called from within FEMA_fit
%
% Parse family relatedness to create working covariance for estimation of random effects
% If outputEB==1 will also create exchangeability blocks for permutation testing using PALM (PALM code needed in path) --> code specific to ABCD data family structure

% Fan et al., (2021) - FEMA: Fast and efficient mixed-effects algorithm for population-scale whole brain imaging data, BioRxiv (https://doi.org/10.1101/2021.10.27.466202)

% INPUTS
%   iid <char>                 :  subject IDs to match imaging data                           
%   eid <char>                 :  eventname
%   fid <num>                  :  family ID (members of the same family unit have same value)
%   agevec <num>               :  participants age
%   ymat <num>                 :  matrix of imaging data (n x v)
%   niter <num>                :  number of iterations (default 1)
%   contrasts <num> OR <path>  :  contrast matrix, or path to file containing contrast matrix (readable by readtable)
%   nbins <num>                :  number of bins across Y for estimating random effects (default 20)
%   pihatmat <num>             :  matrix of genetic relatedness --> already intersected to match X and Y sample
%   outputEB <boolean>         :  default 0 --> will output exchangeability blocks for permutation testing


% supported random effects
% V_F - family - extended family effect
% V_S - subject (longitudinal stability)
% V_A - additive genetic similarity (pihat)
% V_D - dominant (quadradic) genetic effects, V_A .^ 2
% V_M - mama effect (same father ID)
% V_P - papa effect (same mother ID)
% V_H - home effect (fraction of time spend at the same address, currently defined as  V_M+V_P/2)
% V_T - twins effect (same pregnancy ID)
% V_E - environment

      if ~exist('pihatmat','var')
            pihatmat = [];
      end

      p = inputParser;
addParamValue(p,'SingleOrDouble','single');
addParamValue(p,'RandomEffects',{'F' 'S' 'E'}); % Default to Family, Subject, and eps
addParamValue(p,'FatherID',{}); % Father ID, ordered same as pihatmat
addParamValue(p,'MotherID',{}); % Mother ID, ordered same as pihatmat
addParamValue(p,'PregID',{}); % Pregnancy effect (same ID means twins), ordered same as pihatmat
parse(p,varargin{:})
SingleOrDouble = p.Results.SingleOrDouble;
RandomEffects = p.Results.RandomEffects;

[iid_list IA IC_subj] = unique(iid,'stable'); nsubj = length(iid_list); nobs = length(iid);
[fid_list IA IC_fam] = unique(fid,'stable'); nfam = length(fid_list);

iid_repeated_count = NaN(nobs,1);
for i=1:nobs, iid_repeated_count(i)=sum(IC_subj==IC_subj(i)); end

subj_famtypevec=NaN(length(iid),1);

% Find all timepoints for each subject
jvecs_subj = cell(1,nsubj);
for subji = 1:length(iid_list)
  jvec_subj = rowvec(find(IC_subj==subji));
  [mv mi] = min(agevec(jvec_subj));
  jvecs_subj{subji} = jvec_subj;
end

% Should classify each family as one of multiple types: # of subjects, # of observations each: sort subjects & observations in consistent manner
logging('Parsing family structure');
%tic
clusterinfo = cell(1,nfam);
famtypelist = {}; famtypevec = NaN(1,nfam);
nfmemvec = NaN(1,nfam);
for fi = 1:nfam
  jvec_fam = rowvec(find(IC_fam==fi)); % Identify all observations (rows) for a given family
  subj_fam = IC_subj(jvec_fam); % Subject number for each observation
  subj_unique = unique(IC_subj(jvec_fam)); % List of unque subjects
  freq_jvec = NaN(size(jvec_fam)); freq_unique = NaN(size(subj_unique));
  for j = 1:length(subj_unique)
    freq_unique(j) = sum(IC_subj(jvec_fam)==subj_unique(j));
  end
  [sv si] = sort(freq_unique);
  freq_unique = freq_unique(si);
  subj_unique = subj_unique(si); % Re-order subjects to canonical form
  subji_jvec = NaN(size(subj_fam));
  for subjii = 1:length(subj_unique)
    ivec_tmp = find(IC_subj(jvec_fam)==subj_unique(subjii));
    subji_jvec(ivec_tmp) = subjii;
  end
  [sv si] = sortrows([subji_jvec agevec(jvec_fam)]);
  jvec_fam = jvec_fam(si); % Re-order to cannonical form
  ivec = find(cellfun(@(x)isequal(x,freq_unique),famtypelist));
  if length(ivec)==0
    famtypelist = cat(2,famtypelist,{freq_unique});
    famtypevec(fi) = length(famtypelist);
  elseif length(ivec)==1
    famtypevec(fi) = ivec(1);
  else
    keyboard % This shouldn't happen
  end

  subj_famtypevec(jvec_fam)=famtypevec(fi);
  nfmemvec(fi) = length(jvec_fam);

  V_E = eye(length(jvec_fam),length(jvec_fam),SingleOrDouble);
  V_F = true(length(jvec_fam),length(jvec_fam));
  V_S = false(length(jvec_fam),length(jvec_fam));
  if ~isempty(pihatmat)
    V_A = pihatmat(subj_fam(si),subj_fam(si)); % Reorder to canonical form
    if any(~isfinite(V_A(:))) % Impute missing values within families (assume all proper siblings => phat = 0.5)
      V_A(find(~isfinite(V_A))) = 0.5;
    end
%    V_A = V_A - 0.5; % Re-center around 0.5 -- probably need to adjust the non-negativity constraint to make this work
    V_D = V_A.^2; % Should perhaps residualize linear component?
  else
    V_A = []; V_D = [];
  end

  V_P = []; V_M = []; V_H = []; V_T = [];
  if ~isempty(p.Results.FatherID), tmp = repmat(rowvec(p.Results.FatherID(subj_fam(si))), [length(jvec_fam), 1]); V_P = (tmp == tmp'); end
  if ~isempty(p.Results.MotherID), tmp = repmat(rowvec(p.Results.MotherID(subj_fam(si))), [length(jvec_fam), 1]); V_M = (tmp == tmp'); end
  if ~isempty(p.Results.FatherID) && ~isempty(p.Results.MotherID), V_H = (V_P  + V_M) / 2; end
  if ~isempty(p.Results.PregID), tmp = repmat(rowvec(p.Results.PregID(subj_fam(si))), [length(jvec_fam), 1]); V_T = (tmp == tmp'); end

  for ji = 1:length(jvec_fam)
    jvec_tmp = jvecs_subj{IC_subj(jvec_fam(ji))};
    ivec_tmp = ismember(jvec_fam,jvec_tmp);
    V_S(ji,ivec_tmp) = true;
  end
  clusterinfo{fi} = struct('V_E',V_E,'V_S',V_S,'V_F',V_F,'V_A',V_A,'V_D',V_D, ...
                           'V_P',V_P,'V_M',V_M,'V_H',V_H,'V_T',V_T, ...
                           'jvec_fam',jvec_fam);
  clusterinfo{fi}.famtype = famtypevec(fi);
end
nfamtypes = length(famtypelist);
%toc

% Should make list of random effects an optiuonal argument  -- change variable names to be consistent with equation (coordinate with Chun)

% Should modify code above and below to use arbitray list of clusterinfo{:}.Vs
nnz_max = sum(nfmemvec.^2);
Ss = repmat({spalloc(nobs,nobs,nnz_max)},[1 length(RandomEffects)]);
%tic
for fi = 1:nfam
  for ri = 1:length(RandomEffects)
    Ss{ri}(clusterinfo{fi}.jvec_fam,clusterinfo{fi}.jvec_fam) = getfield(clusterinfo{fi},sprintf('V_%s',RandomEffects{ri}));
  end
end
%toc


end
