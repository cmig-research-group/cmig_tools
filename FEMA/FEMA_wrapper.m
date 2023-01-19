function [fpaths_out beta_hat beta_se zmat logpmat sig2tvec sig2mat beta_hat_perm beta_se_perm zmat_perm sig2tvec_perm sig2mat_perm inputs mask tfce_perm colnames_interest save_params logLikvec Hessmat] = FEMA_wrapper(fstem_imaging,fname_design,dirname_out,dirname_tabulated,dirname_imaging,datatype,varargin)
%
% Wrapper function to run whole FEMA pipeline:
%     1) To load and process imaging data (FEMA_process_data)
%     2) To intersect with design matrices (FEMA_intersect_design)
%     3) To run linear mixed effects models (FEMA_fit)
%     4) To save outputs in specified format
%
% USAGE: FEMA_wrapper(fstem_imaging,fname_design,dirname_out,dirname_tabulated,dirname_imaging,datatype,varargin)
%
% INPUTS
%   fstem_imaging <char>       :  name of vertex/voxel-mapped phenotype (e.g., 'thickness-sm16', 'FA')
%   fname_design <cell>        :  cell array with path to file with design matrix saved (readable by readtable) --> if want to batch can add multiple filepaths as separate rows within fname_design as a cell array
%   dirname_out <cell>         :  cell array with path to output directory --> if batching design matrices must include separate output directory as rows within dirname_out as a cell array
%   dirname_tabulated <char>   :  path to tabulated data directory with NDA downloaded txt files
%   dirname_imaging <char>     :  path to imaging data directory
%   datatype <char>            :  'voxel','vertex','external', 'corrmat'
%                                   NB: Other than 'external' all code is written to expect ABCD data
%
% Optional input arguments:
%   contrasts <num> OR <path>  :  contrast matrix, or path to file containing contrast matrix (readable by readtable)
%   ico <num>                  :  ico-number for vertexwise analyses (0-based, default 5)
%   ranknorm <boolean>         :  rank normalise imaging data (default 0)
%   output <string>            :  'mat' (default) or 'nifti' or 'deap' or concatenations to write multiple formats.
%   ivnames <string>           :  comma-separated list of IVs to write [this is used only for DEAP]
%   RandomEffects <cell>       :  list of random effects to estimate (default {'F','S','E'}):
%                                   family relatedness (F)
%                                   subject - required for longitudinal analyses (S)
%                                   error (E) - always required
%                                   additive genetic relatedness (A) - must include file path to genetic relatedness data (pihat) for this option
%   pihat_file <char>          :  path to genetic relatedness data (pihat) - default [] - only required if A random effect specified
%   preg_file <char>           :  path to pregnancy data - default [] - only required if T random effect specified
%   address_file <char>        :  path to address data - default [] - only required if H random effect specified
%   nperms <num>               :  default 0 --> if >0 will run and output permuted effects
%   mediation <num>            :  default 0 --> if 1 will ensure same seed used for resampling of the two models used for a mediation analysis
%   niter <num>                :  input for FEMA_fit - default 1
%   nbins <num>                :  input for FEMA_fit - default 20 - number of bins across Y for estimating random effects
%   CovType <char>             :  input for FEMA_fit - default 'analytic' --> no other options currently available
%   FixedEstType <char>        :  input for FEMA_fit - default 'GLS' --> other option: 'OLS'
%   GroupByFamType <boolean>   :  input for FEMA_fit - default true
%   Parallelize <boolean>      :  input for FEMA_fit - default false
%   NonnegFlag <blooean>       :  input for FEMA_fit - default true - non-negativity constraint on random effects estimation
%   SingleOrDouble <char>      :  input for FEMA_fit - default 'double' --> other option: 'single' - for precision
%   logLikflag <boolean>       :  input for FEMA_fit - default 0
%   permtype <char>            :  input for FEMA_fit - options: 
%                                   'wildbootstrap' - residual boostrap --> creates null distribution by randomly flipping the sign of each observation 
%                                   'wildbootstrap-nn' - non-null boostrap --> estimates distribution around effect of interest using sign flipping (used for sobel test)
%   tfce <num>                 :  default 0 --> if 1 will run TFCE
%   colsinterest <num>         :  used to specify IVs of interest in design matrix (cols in X) for resampling output and tfce (default 1, i.e. 1st column of X) - only used if nperms>0
%
% OUTPUTS
%   fpaths_out                 :  results will be saved here in the specified format
%   All other outputs are optional if user wants results to be output in the MATLAB workspace
%   All permutation outputs will be empty is nperms==0
%


% This software is Copyright (c) 2021 The Regents of the University of California. All Rights Reserved.
% See LICENSE.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

logging('***Start FEMA v2.0, 4/5/2022)');
starttime = now();
rng shuffle %Set random number generator so different every time

%PARSING INPUTS

if nargin < 6
      logging('Usage: FEMA_wrapper(fstem_imaging,fname_design,dirname_out,dirname_tabulated,dirname_imaging,varargin)');
      error('Incorrect number of input arguments')
end

if isdeployed
      logging('***FEMA_wrapper_app Compiled 11/8/2021, FEMA v2.0***\n\n'); %TODO: remember to change this before compiling
end

inputs = inputParser;
addParamValue(inputs,'ranknorm',0);
addParamValue(inputs,'ico',5);
addParamValue(inputs,'contrasts',[]);
addParamValue(inputs,'output','mat');
addParamValue(inputs,'synth',0); % AMD - put back synth option
addParamValue(inputs,'ivnames','');
addParamValue(inputs,'RandomEffects',{'F' 'S' 'E'}); % Default to Family, Subject, and eps
addParamValue(inputs,'pihat_file',[]);
addParamValue(inputs,'preg_file',[]);
addParamValue(inputs,'address_file',[]);
addParamValue(inputs,'nperms',0);
addParamValue(inputs,'mediation',0);
addParamValue(inputs,'tfce',0);
addParamValue(inputs,'colsinterest',1);

%FEMA_fit variable inputs
addParamValue(inputs,'niter',1);
addParamValue(inputs,'nbins',20);
addParamValue(inputs,'CovType','analytic');
addParamValue(inputs,'FixedEstType','GLS');
addParamValue(inputs,'RandomEstType','MoM');
addParamValue(inputs,'GroupByFamType',true);
addParamValue(inputs,'Parallelize',false);
addParamValue(inputs,'NonnegFlag',true); % Perform lsqnonneg on random effects estimation
addParamValue(inputs,'SingleOrDouble','double');
addParamValue(inputs,'logLikflag',0);
addParamValue(inputs,'Hessflag',false);
addParamValue(inputs,'ciflag',false);
addParamValue(inputs,'permtype','wildbootstrap');

addParamValue(inputs,'reverse_cols',1); % AMD in development
addParamValue(inputs,'reverseinferenceflag',0); % AMD in development

parse(inputs,varargin{:})
% Display input arguments for log 
disp(inputs.Results)
niter = str2num_amd(inputs.Results.niter);
ico = str2num_amd(inputs.Results.ico);
contrasts = str2num_amd(inputs.Results.contrasts);
if ~isfinite(contrasts)
      fname_contrasts = inputs.Results.contrasts;
      logging('Reading contrast matrix from %s',fname_contrasts);
      contrasts = readtable(fname_contrasts);
end
ranknorm = str2num_amd(inputs.Results.ranknorm);
outputFormat = inputs.Results.output;
nbins = str2num_amd(inputs.Results.nbins);
if ~(contains(outputFormat,'mat') || contains(outputFormat,'nifti'))
      error('Incorrect output format (%s)',outputFormat)
end

if ~isempty(inputs.Results.ivnames)
      ivnames = split(inputs.Results.ivnames,',');
else
      ivnames = {};
end

if ~isempty(ivnames) || isdeployed
      logging('%d IVs specified (%s)',length(ivnames), inputs.Results.ivnames);
end

RandomEffects = inputs.Results.RandomEffects;
fname_pihat = inputs.Results.pihat_file;
fname_address = inputs.Results.address_file;
fname_pregnancy = inputs.Results.preg_file;
CovType = inputs.Results.CovType;
FixedEstType = inputs.Results.FixedEstType;
RandomEstType = inputs.Results.RandomEstType;
GroupByFamType = inputs.Results.GroupByFamType;
Parallelize = inputs.Results.Parallelize;
NonnegFlag = inputs.Results.NonnegFlag;
SingleOrDouble = inputs.Results.SingleOrDouble;
OLSflag = ismember(lower(FixedEstType),{'ols'});
GLSflag = ismember(lower(FixedEstType),{'gls'});
logLikflag = inputs.Results.logLikflag;
Hessflag = inputs.Results.Hessflag;
ciflag = inputs.Results.ciflag;
nperms = inputs.Results.nperms;
permtype = inputs.Results.permtype;
mediation=inputs.Results.mediation;
synth=inputs.Results.synth;
tfce=inputs.Results.tfce;
colsinterest=inputs.Results.colsinterest;

reverse_cols=inputs.Results.reverse_cols; % AMD -- should replace this with colsinterest
reverseinferenceflag=inputs.Results.reverseinferenceflag; % AMD -- should rename this swapdirectionflag


if ~iscell(fname_design)
      fname_design = {fname_design};
end

if ~iscell(dirname_out)
      dirname_out = {dirname_out};
end

if length(fname_design)~=length(dirname_out)
      error('if using cell array specifying multiple designs, fname_design and dirname_out must both have an equal number of items.')
end

if ~ismember(lower(datatype),{'voxel' 'vertex' 'external' 'corrmat'})
      error('Input error: invalid datatype')
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% LOAD AND PROCESS IMAGING DATA FOR ANALYSIS - ABCD specific function unless datatype='external'
[ymat, iid_concat, eid_concat, ivec_mask, mask, colnames_imaging, pihat, preg, address] = FEMA_process_data(fstem_imaging,dirname_tabulated,dirname_imaging,datatype,'ranknorm',ranknorm,'ico',ico,'pihat_file',fname_pihat,'preg_file',fname_pregnancy,'address_file',fname_address);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% INTERSECT WITH DESIGN MATRIX

%Loops over multiple design matrices (rows in fname_design cell array) to run several models with the same imaging data

pihat_bak=pihat;
ymat_bak=ymat;
cont_bak=contrasts;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% IF RUNNING MEDIATION ANALYSIS

if mediation==1

      seed=rng; % Save rng in order to ensure resampling is the same for both models

      if ~ismember(lower(permtype),{'wildbootstrap-nn'})
            error('Change permtype input to wildboostrap-nn to perform non-null WB needed for mediation analysis OR set mediation=0.')
      end

      if length(fname_design)~=2
            error('fname_design must be a cell array containing 2 nested design matrices to test for mediation.')
      end

      if nperms==0
            error('Change nperms>0 to run mediation analysis OR set mediation=0.')
      end

      if tfce==1
            error('Cannot run mediation and TFCE as different resampling scheme required. Set tfce=0.')
      end

end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

fpaths_out = {};
for des=1:length(fname_design)
            
      [X,iid,eid,fid,agevec,ymat,contrasts,colnames_model,pihatmat,PregID,HomeID] = FEMA_intersect_design(fname_design{des}, ymat_bak, iid_concat, eid_concat, 'contrasts',cont_bak,'pihat',pihat_bak,'preg',preg,'address',address);
      if synth==1 % Make synthesized data
            [ymat sig2tvec_true sig2mat_true] = FEMA_synthesize(X,iid,eid,fid,agevec,ymat,pihatmat,'nbins',nbins,'RandomEffects',RandomEffects); % Make pihatmat and zygmat optional arguments? % Need to update SSE_synthesize_dev to accept list of random effects to include, and range of values

%            sig2mat_true(length(RandomEffects),:) = 1-sum(sig2mat_true(1:length(RandomEffects)-1,:),1); % This shouldn't be needed, if RandomEffects include 'E'
      else
        sig2tvec_true = []; sig2mat_true = [];
      end

      %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
      
      % IF RUNNING MODELS FOR MEDIATION ANALYSIS

      if mediation==1 

            rng(seed); %set same seed for both runs of FEMA_fit

            if des==1
                  X_bak=X;
            end

            if des==2
                  if size(X_bak,1)~=size(X,1)
                        error('Design matrices must have the same number of observations for mediation analysis.') % Check after intersection
                  end
            end

      end

      %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

      synthstruct = struct('sig2tvec_true',sig2tvec_true,'sig2mat_true',sig2mat_true);

      % FIT MODEL
      [beta_hat beta_se zmat logpmat sig2tvec sig2mat Hessmat logLikvec beta_hat_perm beta_se_perm zmat_perm sig2tvec_perm sig2mat_perm] = FEMA_fit(X,iid,eid,fid,agevec,ymat,niter,contrasts,nbins, pihatmat,'RandomEffects',RandomEffects,...
            'nperms',nperms,'CovType',CovType,'FixedEstType',FixedEstType,'RandomEstType',RandomEstType,'GroupByFamType',GroupByFamType,'Parallelize',Parallelize,'NonnegFlag',NonnegFlag,'SingleOrDouble',SingleOrDouble,'logLikflag',logLikflag,'Hessflag',Hessflag,'ciflag',ciflag,...
            'permtype',permtype,'PregID',PregID,'HomeID',HomeID,'synthstruct',synthstruct);

            if sum(~mask)>0

                  z_tmp=zeros(size(zmat,1),size(mask,2));
                  p_tmp=zeros(size(logpmat,1),size(mask,2));
                  beta_tmp=zeros(size(beta_hat,1),size(mask,2));
                  betase_tmp=zeros(size(beta_se,1),size(mask,2));
                  sig2mat_tmp=zeros(size(sig2mat,1),size(mask,2));
                  sig2tvec_tmp=zeros(size(sig2tvec,1),size(mask,2));

                  z_tmp(:,ivec_mask)=zmat;
                  p_tmp(:,ivec_mask)=logpmat;
                  beta_tmp(:,ivec_mask)=beta_hat;
                  betase_tmp(:,ivec_mask)=beta_se;
                  sig2mat_tmp(:,ivec_mask)=sig2mat;
                  sig2tvec_tmp(:,ivec_mask)=sig2tvec;

                  zmat=z_tmp;
                  logpmat=p_tmp;
                  beta_hat=beta_tmp;
                  beta_se=betase_tmp;
                  sig2mat=sig2mat_tmp;
                  sig2tvec=sig2tvec_tmp;

                  if nperms>0

                        zperm_tmp=zeros(size(zmat_perm,1),size(mask,2),size(zmat_perm,3));
                        betaperm_tmp=zeros(size(beta_hat_perm,1),size(mask,2),size(beta_hat_perm,3));
                        betaseperm_tmp=zeros(size(beta_se_perm,1),size(mask,2),size(beta_se_perm,3));
                        sig2matperm_tmp=zeros(size(sig2mat_perm,1),size(mask,2),size(sig2mat_perm,3));
                        sig2tvecperm_tmp=zeros(size(sig2tvec_perm,1),size(mask,2),size(sig2tvec_perm,3));

                        zperm_tmp(:,ivec_mask,:)=zmat_perm;
                        betaperm_tmp(:,ivec_mask,:)=beta_hat_perm;
                        betaseperm_tmp(:,ivec_mask,:)=beta_se_perm;
                        sig2matperm_tmp(:,ivec_mask,:)=sig2mat_perm;
                        sig2tvecperm_tmp(:,ivec_mask,:)=sig2tvec_perm;

                        zmat_perm=zperm_tmp;
                        beta_hat_perm=betaperm_tmp;
                        beta_se_perm=betaseperm_tmp;
                        sig2mat_perm=sig2matperm_tmp;
                        sig2tvec_perm=sig2tvecperm_tmp;

                  end

            end

      %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

      % IF RUNNING TFCE

      if tfce==1 && nperms>0
            tfce_perm = FEMA_tfce(zmat_perm,colsinterest,datatype,mask);
      elseif tfce==1 && nperms==0
            error('Change nperms>0 to run TFCE OR set tfce=0.')
      else
            tfce_perm=[];
      end

      %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

      % IF RUNNING RESAMPLING

      % Only save permutations for IVs of interest to reduce size of output
      % If wanting to save permutations for all columns of X make colsinterest=[1:size(X,2)]

      if nperms>0
            beta_hat_perm=beta_hat_perm(colsinterest,:,:);
            beta_se_perm=beta_se_perm(colsinterest,:,:);
            zmat_perm=zmat_perm(colsinterest,:,:);
            colnames_interest=colnames_model(colsinterest);
      elseif nperms==0
            colnames_interest=colnames_model;
      end

      %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            
      % SAVE OUTPUT

      if nperms>0 && strcmp(permtype,'wildbootstrap')
            dirname_out{des}=sprintf('%s/nullWB_%dperms',dirname_out{des},nperms);
      elseif nperms>0 && strcmp(permtype,'wildbootstrap-nn')
            dirname_out{des}=sprintf('%s/nonnullWB_%dperms',dirname_out{des},nperms);
      end

      save_params = struct('fstem_imaging',fstem_imaging,'datatype',datatype,'outdir',dirname_out{des},'synth',synth);
      base_variables_to_save = {'X','iid','eid','colnames_model','contrasts','datatype','inputs','zmat','logpmat','beta_hat','beta_se','sig2mat','sig2tvec','save_params','mask'};

      if ~exist(dirname_out{des},'dir'), mkdir(dirname_out{des}); end

      if synth==0
            fpath_out = sprintf('%s/FEMA_wrapper_output_%s_%s.mat',dirname_out{des},datatype,fstem_imaging);
      elseif synth==1
            fpath_out = sprintf('%s/FEMA_wrapper_output_%s_%s_synth.mat',dirname_out{des},datatype,fstem_imaging);
      end

      %write column names to json for DEAP
      fname_col = sprintf('%s/FEMA_results_colnames.json',dirname_out{des});
      out = struct('colnames_model',{colnames_model},'RandomEffects',{RandomEffects});
      jsonStr = jsonencode(out);
      fid = fopen(fname_col,'w');
      fprintf(fid,'%s\n',jsonStr);
      fclose(fid);

      % =========================================================================
      % Write VOXEL results (mat, nifti, or deap)
      % =========================================================================
      if strcmpi(datatype,'voxel')
            
            vol_z = zeros([size(mask) size(zmat,1)]);
            vol_logp = zeros([size(mask) size(zmat,1)]);
            vol_beta_hat = zeros([size(mask) size(zmat,1)]);
            vol_beta_se = zeros([size(mask) size(zmat,1)]);
            for j = 1:size(zmat,1)
                  vol_z(:,:,:,j) = single(fullvol(zmat(j,:),mask));
                  vol_logp(:,:,:,j) = single(fullvol(logpmat(j,:),mask));
                  vol_beta_hat(:,:,:,j) = single(fullvol(beta_hat(j,:),mask));
                  vol_beta_se(:,:,:,j) = single(fullvol(beta_se(j,:),mask));
            end
            
            vol_sig2t = zeros([size(mask) 1]);
            vol_sig2t(ivec_mask) = single(sig2tvec);
            vol_sig2 = zeros([size(mask) size(sig2mat,1)]);
            for j = 1:size(sig2mat,1)
                  vol_sig2(:,:,:,j) = single(fullvol(sig2mat(j,:),mask));
            end

            
            % ============================================================================================================================
            % == MAT Output ==
            if contains(outputFormat, 'mat')

                  if nperms>0 & tfce==0
                        save(fpath_out,base_variables_to_save{:},'vol_z','vol_beta_hat','zmat_perm','beta_hat_perm','colnames_interest','colsinterest','-v7.3');
                  elseif nperms>0 & tfce==1
                        save(fpath_out,base_variables_to_save{:},'vol_z','vol_beta_hat','zmat_perm','beta_hat_perm','tfce_perm','colnames_interest','colsinterest','-v7.3');
                  elseif nperms==0
                        save(fpath_out,base_variables_to_save{:},'vol_z','vol_beta_hat','logpmat','vol_sig2','vol_sig2t','-v7.3');
                  end
                  logging('Results written to %s',fpath_out);

            end
            
            % ============================================================================================================================
            % == NIFTI Output == FIXME: no longer used for DEAP
            if contains(outputFormat, 'nifti')
                  results = struct('beta_hat',vol_beta_hat,'beta_se',vol_beta_se,'zmat',vol_z,'logpmat',vol_logp,'sig2tvec',vol_sig2t,'sig2mat',vol_sig2);
                  writeNIFTI(results, dirname_out{des}, fstem_imaging, ivnames, colnames_model); 
            end
            
      % =========================================================================
      % Write VERTEX results
      % =========================================================================
      elseif strcmpi(datatype, 'vertex')
            
            if contains(outputFormat,'mat')

                  if nperms>0 & tfce==0
                        save(fpath_out,base_variables_to_save{:},'zmat_perm','beta_hat_perm','colnames_interest','colsinterest','-v7.3');
                  elseif nperms>0 & tfce==1
                        save(fpath_out,base_variables_to_save{:},'zmat_perm','beta_hat_perm','tfce_perm','colnames_interest','colsinterest','-v7.3');
                  elseif nperms==0
                        save(fpath_out,base_variables_to_save{:},'-v7.3');
                  end

                  logging('Results written to %s',fpath_out);
            end
            
            if contains(outputFormat, 'nifti') %FIXME: these are much smaller, so haven't added the same optimization as for voxelwise

                  randomFields = {'sig2tvec', 'sig2mat'};

                  results = struct('beta_hat',beta_hat,'beta_se',beta_se,'zmat',zmat,'logpmat',logpmat,'sig2tvec',sig2tvec,'sig2mat',sig2mat);
                  fieldnamelist = fieldnames(results);
                  for fi = 1:length(fieldnamelist)
                        fieldname = fieldnamelist{fi};
                        fname_tmp = sprintf('%s/FEMA_results_vertexwise_%s_%s.nii',dirname_out{des},fstem_imaging,fieldname);
                        vol_nifti = reshape2nifti(getfield(results,fieldname)');
                        niftiwrite(vol_nifti,fname_tmp,'Compressed',true);
                        fprintf(1,'file %s written (dims = [%s])\n',fname_tmp,num2str(size(vol_nifti),'%d '));
                  end

            end
            
      % =========================================================================
      % Write EXTERNAL results (mat)
      % =========================================================================
      elseif strcmpi(datatype, 'external')
            
            if contains(outputFormat, 'mat')
                  warning('Only outputFormat "mat" currently supported for external csv data')
            end
            
            if nperms>0
                  save(fpath_out,base_variables_to_save{:},'colnames_imaging','zmat_perm','beta_hat_perm','colnames_interest','colsinterest','-v7.3');
            elseif nperms==0
                  save(fpath_out,base_variables_to_save{:},'colnames_imaging','-v7.3');
            end

            logging('Results written to %s',fpath_out);
            
      % =========================================================================
      % Write CORRMAT results FIXME: saving is not implemented
      % =========================================================================
      elseif strcmpi(datatype, 'corrmat')
            
            if 0
                  figure; im = reshape(beta_hat(1,:),dims(2:end)); imagesc(im,max(abs(im(:)))*[-1 1]); colormap(blueblackred); axis equal tight;
                  figure; im = reshape(beta_hat(2,:),dims(2:end)); imagesc(im,max(abs(im(:)))*[-1 1]); colormap(blueblackred); axis equal tight;
                  figure; im = reshape(zmat(1,:),dims(2:end)); imagesc(im,max(abs(im(:)))*[-1 1]); colormap(blueblackred); axis equal tight;
                  figure; im = reshape(zmat(2,:),dims(2:end)); imagesc(im,max(abs(im(:)))*[-1 1]); colormap(blueblackred); axis equal tight;
            end
            %beta_hat = reshape(beta_hat,[size(beta_hat,1) dims(2:end)]);
            %beta_se = reshape(beta_se,[size(beta_hat,1) dims(2:end)]);
            %zmat = reshape(zmat,[size(beta_hat,1) dims(2:end)]);
            %logpmat = reshape(logpmat,[size(beta_hat,1) dims(2:end)]);
            %sig2tvec = reshape(sig2tvec,[size(sig2tvec,1) dims(2:end)]);
            %sig2mat = reshape(sig2mat,[size(sig2mat,1) dims(2:end)]);
            
            if nperms>0
                  save(fpath_out,base_variables_to_save{:},'colnames_imaging','zmat_perm','beta_hat_perm','colnames_interest','colsinterest','-v7.3');
            else
                  save(fpath_out,base_variables_to_save{:},'colnames_imaging','-v7.3');
            end
            logging('Results written to %s',fpath_out);
            
      end
      
      if ~isempty(fpath_out)
            fpaths_out = cat(2,fpaths_out,fpath_out);
      end
      
end  % LOOP over design matrices

endtime = now();
logging('***Done*** (%0.2f seconds)',(endtime-starttime)*3600*24);
