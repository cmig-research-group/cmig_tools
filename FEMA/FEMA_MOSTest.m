function [fpath_out,rank_log10pval_mostest, extrap_log10pval_mostest, save_params] = FEMA_MOSTest(zmat_perm, colnames_interest, varargin)

% FEMA_MOSTest - Computes the MOSTest statistics on the permuted z-scores.
%
% Syntax: [rank_log10pval_mostest, extrap_log10pval_mostest, mostvec] = FEMA_MOSTest(zmat_perm, colnames_interest, outdir, fstem_imaging, datatype, varargin)
%
% INPUTS:
%   zmat_perm          : (n_variables x n_pixels x (n_permutations+1))
%   colnames_interest  : cf. FEMA_wrapper

% Optional input arguments
%   save_params <struct>       :  structure output from FEMA_wrapper with parameters to save file with FEMA naming and directory conventions
%                                   N.B. - outputs will only be saved if save_params provided as an input
%   k: (double) optimised regularisation parameter

% OUTPUTS:
%   fpath_out <char>                 :  output filepath (empty if save_params not provided)
%   rank_log10pval_mostest    : log10 p-value based on rank from empirical distribution
%   extrap_log10pval_mostest  : log10 p-value based on `FEMA_tails` extrapolation

    most_inputs = inputParser;
    addParamValue(most_inputs,'k','auto');
    addParamValue(most_inputs,'save_params',[]);
    parse(most_inputs,varargin{:});
    save_params = most_inputs.Results.save_params;

    % Setting the MOSTest regularization parameter
    if (isempty(most_inputs.Results.k) | strcmp(most_inputs.Results.k, 'auto')) & ~isempty(save_params)
        if strcmp(save_params.datatype, 'vertex')
            if contains(save_params.fstem_imaging, 'area')
                    k = 79;
            elseif contains(save_params.fstem_imaging, 'thickness')
                    k = 47;
            else
                    warning("The MOSTest regularization parameter hasn't been optimized for this set of parameters yet. Defaulting to 100.")
                    k = 100;
            end
        else
            warning("MOSTest is not yet implemented for analyses other than vertex.")
            return
        end
    elseif isempty(most_inputs.Results.k) & isempty(save_params)
        warning('Missing parameters. Need a value for `k` or for `save_params`.')
    else
        warning("Setting your own regularization parameter. NB: MOSTest has not yet been validated for analyses other than vertex.")
        k = most_inputs.Results.k
    end

    starttime = now();

    nperms = size(zmat_perm,3) - 1; %TODO should we do -1?

    rank_log10pval_mostest=NaN(size(zmat_perm,1),1);
    extrap_log10pval_mostest=NaN(size(zmat_perm,1),1);
    mostest_stat=NaN(size(zmat_perm,1),1);
    mostest_mat=NaN(size(zmat_perm,1),size(zmat_perm,3));

    for coli=1:size(zmat_perm,1) 

        % Extract z-scores of interest
        zmat = double(squeeze(zmat_perm(coli,:,:)));
        zmat(~isfinite(zmat)) = 0;

        % Compute covariance matrix (excluding the non-permuted column)
        logging('Estimating covariance matrix');
        C0 = cov(zmat(:,2:end)');
        % Compute inverse of C0, with regularization
        logging('Inverting covariance matrix');
        [U S V] = svd(C0); s = diag(S);
        s_reg = max(s(min(k,length(s))),s);
        C0_inv = U*diag(s_reg.^-1)*U';

        % Partial saving (since computing C0_inv is compute intensive)
        if ~isempty(save_params)
            outdir=save_params.outdir;
            fstem_imaging=save_params.fstem_imaging;
            datatype=save_params.datatype;
            if save_params.synth==0
                  fpath_out = sprintf('%s/FEMA_mostest_partial_results_%s_%s.mat',outdir,datatype,fstem_imaging);
            elseif save_params.synth==1
                  fpath_out = sprintf('%s/FEMA_mostest_partial_results_%s_%s_synth.mat',outdir,datatype,fstem_imaging);
            end
            save(fpath_out,'C0_inv','zmat','coli','nperms','save_params','-v7.3')
        end

        % Compute MOSTest (using Mahalanobis distance)
        logging('Computing MOSTest statistic');
        mostvec = sum(zmat'*C0_inv.*zmat',2);
        mostest_stat(coli)=mostvec(1);
        mostest_mat(coli,:)=mostvec;
        mostvec_sorted = sortrows(mostvec(2:end));
        rank_log10pval_mostest(coli) = -log10((length(find(mostvec_sorted>=mostvec(1)))+1)/nperms);
        % Extrapolate p-values (using different distributions, but using dale tails later on)
        pl = 0.9; pu = 0.999;
        [nlog10cdfvec_ecdf nlog10cdfvec_dale nlog10cdfvec_pareto pd_dale pd_pareto xvals] = FEMA_tails(mostvec(2:end),pl,pu,[],'gamma');

        % Estimate p-value from fitted distribution
        extrap_log10pval_mostest(coli) = interp1(xvals,nlog10cdfvec_dale,mostest_stat(coli),'linear','extrap');
    end

    if ~isempty(save_params)

        if save_params.synth==0
            fpath_out = sprintf('%s/FEMA_mostest_results_%s_%s.mat',outdir,datatype,fstem_imaging);
        elseif save_params.synth==1
            fpath_out = sprintf('%s/FEMA_mostest_results_%s_%s_synth.mat',outdir,datatype,fstem_imaging);
        end
        save(fpath_out,'zmat_perm','rank_log10pval_mostest', 'extrap_log10pval_mostest', 'mostest_stat','mostest_mat','colnames_interest','save_params','-v7.3');
        logging('MOSTest results saved to %s', fpath_out)
    elseif isempty(save_params)
      fpath_out=[];
    end

    for kk=1:length(colnames_interest)
        logging('MOSTest results (%s): log10 p-value=%.3f, log10 extrapolated p-value=%.3f.',colnames_interest{kk}, rank_log10pval_mostest(kk), extrap_log10pval_mostest(kk));
    end
    return
end