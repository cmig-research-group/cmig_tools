function [fpath_out, mediation_tstat, mediation_pvals, mediation_effect, mediation_se, mediation_effect_CInorm, mediation_effect_CIprctile, save_params, mask] = FEMA_sobel_test(FEMA_outfile_reduced,FEMA_outfile_full,IVname,alpha, varargin)

% Runs a mediation analysis (Sobel's test) on nested models using output from FEMA_wrapper.m (Sobel, 1990)
% The mediation effect is calculated as (τ – τ'). This is the change in the magnitude of the effect that
% the independent variable (IV) has on the dependent variable (DV) after controlling for the mediator.
% `IVname` specifies the IV of interest.
% The DV is always the imaging phenotype.
% The mediator is the variable that is present in the full model but not the reduced model.

% NOTE: Resampling must be done using an identical resampling scheme for both the reduced and full models
% in order for the mediation analysis to be accurate. Make sure 'mediation' is set to 1 when running
% FEMA_wrapper.m with both design matrices.

% To save the output from FEMA_sobel_test.m include optional argument 'save_params'
% with the save_params structure output from FEMA_wrapper. This will save the p-values in
% a file within the FEMA directory structure and naming conventions.

% USAGE: FEMA_sobel_test(FEMA_outfile_reduced,FEMA_outfile_full,outdir,IVname,fstem_imaging,datatype)

% INPUTS
%     FEMA_outfile_reduced <char>      :  path with filename for FEMA output for reduced model
%     FEMA_outfile_full <char>         :  path with filename for FEMA output for full model
%     IVname <char>                    :  variable name of independent variable of interest for mediation model
%     alpha <num>                      :  alpha level for calculation of confidence intervals (e.g. alpha=0.05 ~ 95% CI)

% Optional input arguments
%   save_params <struct>               :  structure output from FEMA_wrapper with parameters to save file with FEMA naming and directory conventions
%                                           N.B. - outputs will only be saved if save_params provided as an input
%   mask <num>                         :  imaging mask

% OUTPUTS
%     fpath_out <char>                 :  output filepath (empty if save_params not provided)
%     mediation_tstat <num>            :  sobel t-statistic for the mediation effect
%     mediation_pvals <num>            :  p-value for sobel test statistic --> estimated from fitting a chi-squared distribution to the empirical resampled distribution
%     mediation_effect <num>           :  difference between beta coefficients for IVname in reduced - full model (i.e., tau - tau')
%     mediation_se <num>               :  empirical standard error estimated from the resampled distribution of mediation effects
%     mediation_effect_CInorm <num>    :  upper and lower CI at input alpha level estimated using the normal distribution approach
%     mediation_effect_CIprctile <num> :  upper and lower CI at input alpha level estimated using the percentile approach


% This software is Copyright (c) 2021 The Regents of the University of California. All Rights Reserved.
% See LICENSE.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

sobel_inputs = inputParser;
addParamValue(sobel_inputs,'save_params',[]);
addParamValue(sobel_inputs,'mask',[]);

parse(sobel_inputs,varargin{:})
save_params = sobel_inputs.Results.save_params;
mask = sobel_inputs.Results.mask;

      logging('Loading data...');
      tic
      model_red=load(FEMA_outfile_reduced);
      model_full=load(FEMA_outfile_full);
      toc

      % Check if FEMA_wrapper was run with correct input
      if model_red.inputs.Results.mediation~=1 | model_full.inputs.Results.mediation~=1
            error('FEMA_wrapper was not run with mediation==1.  Models will not necessarily have same resampling scheme.  Mediation analysis will be inaccurate.  Please re-run models using FEMA_wrapper with mediation==1.')
      end
      
      % Check sizes of design matrices to specify full and reduced models
      if size(model_red.X,2)==size(model_full.X,2)
            error('Design matrices for mediation analysis have same number of columns.');
      elseif size(model_full.X,2)<size(model_red.X,2)
            warning('Model full has fewer columns than model reduced. Switching models to calculate mediation effect correctly (reduced - full).')
            model_tmp=model_red;
            model_red=model_full;
            model_full=model_tmp;
            model_tmp=[];
      end   

      logging('Calculating sobel test statistic...');

      if isempty(IVname) && length(model_full.colnames_interest)==1
        IVname = model_full.colnames_interest{1};
      end

      if isfield(model_red, 'colnames_interest')
            DVind_red=find(~cellfun('isempty',regexp(IVname,model_red.colnames_interest,'match')));
      else
            DVind_red=find(~cellfun('isempty',regexp(IVname,model_red.colnames_model,'match')));
      end

      if isfield(model_full, 'colnames_interest')
            DVind_full=find(~cellfun('isempty',regexp(IVname,model_full.colnames_interest,'match')));
      else
            DVind_full=find(~cellfun('isempty',regexp(IVname,model_full.colnames_model,'match')));
      end

      if isempty(DVind_red) | isempty(DVind_full)
            error('IVname is not in one or both of your models.')
      end

      tau = permute(model_red.beta_hat_perm(DVind_red,:,:),[3,2,1]);
      tau_prime = permute(model_full.beta_hat_perm(DVind_full,:,:),[3,2,1]);

      %Calculate mediation effect and std from bootstrapped distribution
      mediation_effect_perms = tau - tau_prime;
      mediation_effect = mediation_effect_perms(1,:);
      mediation_se = std(mediation_effect_perms(2:end,:)); % AMD: include permutations only

      %Calculate confidence intervals using:
      %1) normal approach
      mediation_effect_CInorm(1,:) = mediation_effect-(abs(norminv(alpha/2)).*mediation_se);
      mediation_effect_CInorm(2,:) = mediation_effect+(abs(norminv(alpha/2)).*mediation_se);
      %2) percentile approach
      mediation_effect_CIprctile = prctile(mediation_effect_perms, [alpha/2,1-(alpha/2)],1); 

      %Calculate sobel T statistic for every voxel and boostrap
      mediation_tstat = mediation_effect_perms./mediation_se; 

      mediation_tstat = mediation_tstat - mean(mediation_tstat(2:end,:)); % Remove mean

      statvec = colvec(mediation_tstat(2:end,:)).^2;
      xvec =  linspace(0,100,10001);
      hc = hist(statvec,xvec); chc = cumsum(hc)/sum(hc);  

      %pd_fit = makedist('gamma','a',0.5,'b',2); % Hard-code Chi-Square distribution
      distname = 'gamma';
%      distname = 'weibull';
      statvec=double(statvec);
      pd_fit = fitdist(statvec,distname); % Fit to specific distribution 
      chc_fit = pd_fit.cdf(xvec);
      %figure; plot(xvec,-log10(1-chc),xvec,-log10(1-chc_fit),'LineWidth',2); 

      mediation_pvals = pd_fit.cdf(mediation_tstat.^2,'upper');

      fields={'beta_hat_perm','zmat_perm'};
      model_red = rmfield(model_red, fields);
      model_full = rmfield(model_full, fields);

      mediation_tstat = mediation_tstat(1,:);
      mediation_pvals = mediation_pvals(1,:);

      if ~isempty(save_params)

            outdir=save_params.outdir;
            fstem_imaging=save_params.fstem_imaging;
            datatype=save_params.datatype;

            if save_params.synth==0
                  fpath_out = sprintf('%s/FEMA_sobel_results_%s_%s_%s.mat',outdir,IVname,datatype,fstem_imaging);
            elseif save_params.synth==1
                  fpath_out = sprintf('%s/FEMA_sobel_results_%s_%s_%s_synth.mat',outdir,IVname,datatype,fstem_imaging);
            end

            if strcmpi(datatype,'voxel')

                  mask=model_red.mask;
                  ivec_mask = find(mask>0.5);

                  vol_med_tstat = zeros([size(mask) 1]);
                  vol_med_tstat(ivec_mask) = mediation_tstat;
                  vol_med_effect = zeros([size(mask) 1]);
                  vol_med_effect(ivec_mask) = mediation_effect;
                  vol_med_pval = zeros([size(mask) 1]);
                  vol_med_pval(ivec_mask) = mediation_pvals;

                  save(fpath_out,'model_red','model_full','mediation_effect', 'mediation_se', 'mediation_effect_CInorm', 'mediation_effect_CIprctile', 'mediation_tstat', 'mediation_pvals','vol_med_tstat','vol_med_effect','vol_med_pval','mask','IVname','save_params','-v7.3');
      
            elseif strcmpi(datatype, 'vertex') | strcmpi(datatype, 'external') | strcmpi(datatype, 'corrmat')

                  save(fpath_out,'model_red','model_full','mediation_effect', 'mediation_se', 'mediation_effect_CInorm', 'mediation_effect_CIprctile', 'mediation_tstat', 'mediation_pvals','IVname','save_params','mask','-v7.3');

            end

      elseif isempty(save_params)
            fpath_out=[];
      end

keyboard
save('~/FEMA_sobel_test_snap.mat');

end

