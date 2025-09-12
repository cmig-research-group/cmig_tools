function FEMA_classify(fstem_imaging,fname_design,dirname_out,dirname_imaging,datatype,target,varargin)
%
% Wrapper function to run whole FEMA pipeline:
%     1) To load and process imaging data (FEMA_process_data)
%     2) To intersect with design matrices (FEMA_intersect_design)
%     3) To run linear mixed effects models (FEMA_fit)
%     4) To save outputs in specified format
%
% USAGE: FEMA_wrapper(fstem_imaging,fname_design,dirname_out,dirname_imaging,datatype,varargin)
%
% INPUTS
%   fstem_imaging <char>       :  name of vertex/voxel-mapped phenotype (e.g., 'thickness-sm16', 'FA')
%   fname_design <cell>        :  cell array with path to file with design matrix saved (readable by readtable) --> if want to batch can add multiple filepaths as separate rows within fname_design as a cell array
%   dirname_out <cell>         :  cell array with path to output directory --> if batching design matrices must include separate output directory as rows within dirname_out as a cell array
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

% logging('***Start FEMA v2.0, 4/5/2022)');
logging('***Start***');
disp(FEMA_info);

starttime = now();
rng shuffle %Set random number generator so different every time

%PARSING INPUTS

if nargin < 6
    logging('Usage: FEMA_wrapper(fstem_imaging,fname_design,dirname_out,dirname_imaging,varargin)');
    error('Incorrect number of input arguments')
end

if isdeployed
    logging('***FEMA_wrapper_app Compiled 11/8/2021, FEMA v2.0***\n\n'); %TODO: remember to change this before compiling
end

inputs = inputParser;
addParamValue(inputs,'ranknorm',0);
addParamValue(inputs,'varnorm',0);
addParamValue(inputs,'ico',5);
addParamValue(inputs,'contrasts',[]);
addParamValue(inputs,'synth',0); % AMD - put back synth option
addParamValue(inputs,'pihat_file',[]);
addParamValue(inputs,'preg_file',[]);
addParamValue(inputs,'address_file',[]);
addParamValue(inputs,'binsz',2.5);
addParamValue(inputs,'do_agebins', 1); % split by age bins or not 

parse(inputs,varargin{:})
% Display input arguments for log
disp(inputs.Results)
ranknorm = str2num_amd(inputs.Results.ranknorm);
varnorm = str2num_amd(inputs.Results.varnorm);
ico = inputs.Results.ico;
contrasts = inputs.Results.contrasts; 
synth = str2num_amd(inputs.Results.synth);
fname_pihat = inputs.Results.pihat_file;
fname_address = inputs.Results.address_file;
fname_pregnancy = inputs.Results.preg_file;
binsz = inputs.Results.binsz;
do_agebins = inputs.Results.do_agebins;

% make it possible to submit multiple targets same as design matrices
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
if ~iscell(target)
    target = {target};
end 
if length(target)~=length(fname_design)
	target = repmat(target,1,length(fname_design));
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% LOAD AND PROCESS IMAGING DATA FOR ANALYSIS - ABCD specific function unless datatype='external'

[ymat, iid_concat, eid_concat, ivec_mask, mask, colnames_imaging, pihat, preg, address] = FEMA_process_data(fstem_imaging,dirname_imaging,datatype,'ranknorm',ranknorm,'varnorm',varnorm,'ico',ico,'pihat_file',fname_pihat,'preg_file',fname_pregnancy,'address_file',fname_address); % removed infostruct to match standard FEMA code

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% INTERSECT WITH DESIGN MATRIX

%Loops over multiple design matrices (rows in fname_design cell array) to run several models with the same imaging data

pihat_bak=pihat;
ymat_bak=ymat;
cont_bak=contrasts;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

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
    % Try crossvalidated multivariate classification

    colnum = find(strcmp(colnames_model, target{des}));
    if isempty(colnum)
        warning(sprintf('target %s not found in design matrix', target{des}));
        continue;
    end

    k = 10;
    fidlist = unique(fid,'stable'); 
    iidlist = unique(iid,'stable');

    Y = X(:,colnum)==1;

    if 1
        M = X(:,setdiff(1:size(X,2),colnum)); 
        ymat_resid = ymat - M*pinv(M)*ymat; % Residualize imaging by covariates, assuming first column of X is the variable of interest
    else
        ymat_resid = ymat;
    end

    fprintf(1,'Performing SVD\n');
    tic
    [U S V] = svd(ymat_resid,'econ');
    toc

    % Perform multiple repetitions of fit across ages, and by age bins

    %nbins = 5; nr = 2; nc = 3; agebinlist = linspace(9,19,nbins+1);
    %agebinlist = [8.5 10.8 13.2 15.5 18]; nbins = length(agebinlist)-1; nr = 1; nc = 4;
    %agebinlist = [8.5 11.0 13.5 16.0 18.5]; nbins = length(agebinlist)-1; nr = 1; nc = 4; % Probably makes more sense to align bins to visits
    
	% check if agevec is in months or years
    if max(agevec)>25
        agevec = X(:,2)/12; 
    end 

    % the min and max age at each time point
	%[agerange_00A]  = [min(agevec(find(strcmp(eid, 'ses-0A')))), max(agevec(find(strcmp(eid, 'ses-00A'))))];
	%[agerange_02A]  = [min(agevec(find(strcmp(eid, 'ses-02A')))), max(agevec(find(strcmp(eid, 'ses-02A'))))];
	%[agerange_04A]  = [min(agevec(find(strcmp(eid, 'ses-04A')))), max(agevec(find(strcmp(eid, 'ses-04A'))))];
	%[agerange_06A]  = [min(agevec(find(strcmp(eid, 'ses-06A')))), max(agevec(find(strcmp(eid, 'ses-06A'))))];

    %agebinlist = [8.0 11.0 13.0 15.0];
	agevec_min =  floor(min(agevec));
	agevec_max =  ceil(max(agevec));
	agebinlist = agevec_min:binsz:agevec_max;
    nbins = length(agebinlist)-1; 
    nr = 1; 
    nc = nbins; % Probably makes more sense to align bins to visits

	for bini = 1:nbins % Should perhaps crossvalidate based on only families/subjects in bin?
		if bini == nbins
			agecat{bini} = strcat(num2str(agebinlist(bini)), '-', num2str(agebinlist(bini+1)));
		else 
			agecat{bini} = strcat(num2str(agebinlist(bini)), '-', num2str(agebinlist(bini+1)-0.1));
		end 
	end 
    
    nrep = 10;
    scorevecs_cat = repmat({[]},[1 nbins]); 
    scorevec_cat = []; 
    AUCvec = NaN(1,nrep); 
    AUCmat = NaN(nbins,nrep); 
    Y_cat = [];

    for repi = 1:nrep
        fprintf(1,'repi=%d/%d (%s)\n',repi,nrep,datestr(now));
        fcvind = crossvalind('Kfold',length(fidlist),k); 
        cvind = NaN(size(fid));
        for ki = 1:k
            cvind(ismember(fid,fidlist(fcvind==ki))) = ki;
        end
        scorevec = NaN(size(cvind)); 
        ncomp = min(size(U,2),250);
        for ki = 1:k
            fprintf(1,'  ki=%d/%d (%s)\n',ki,k,datestr(now));
            ivec_disc = find(cvind~=ki); 
            ivec_repl = find(cvind==ki);
            U_resid = U - M*pinv(M(ivec_disc,:))*U(ivec_disc,:);
            if 0 % Re-compute svd based on ivec_disc only
                [U_tmp S_tmp V_tmp] = svd(U_resid(ivec_disc,:)*S,'econ'); 
                D = U_resid*S*V_tmp;  % Note that in previous version we used U, not scaled by S -- should perhaps instead scale by diag(sqrt(diag(S).*diag(S_tmp))?)
            else
                D = U_resid;
            end
            mdl = fitcsvm(D(ivec_disc,1:ncomp),Y(ivec_disc),'Standardize',true,'KernelFunction','RBF','KernelScale','auto');
            [label,score] = predict(mdl,D(ivec_repl,1:ncomp));
            scorevec(ivec_repl) = score(:,2);
        end
        if 0 % Whould check that D is very similar to ymat_resid -- i.e., that the U, S, and V computed based on the
            figure; 
            imagesc(U_resid(ivec_disc,:)); 
            colormap(blueblackred); 
            colorbar; 
            xlim([1 50]); 
            ylim([1 500])
            figure; 
            imagesc(D(ivec_disc,:)); 
            colormap(blueblackred); 
            colorbar; 
            xlim([1 50]); 
            ylim([1 500])
        end
        Y_cat = cat(1,Y_cat,Y); 
        scorevec_cat = cat(1,scorevec_cat,scorevec);
        %  sfigure(2); [AUC SPEC SENS ACC F80] = plot_ROC(colvec(scorevec(Y==1)),colvec(scorevec(Y==0))); title(sprintf('%s predicting %s AUC=%0.2f',fstem_imaging,colnames_model{colnum},AUC),'Interpreter','none'); axis equal tight;
        sfigure(2); 
        [AUC SPEC SENS ACC F80] = plot_ROC(colvec(scorevec_cat(Y_cat==1)),colvec(scorevec_cat(Y_cat==0))); 
        title(sprintf('%s predicting %s AUC=%0.3f',fstem_imaging,colnames_model{colnum},AUC),'Interpreter','none'); axis equal tight; % Why is the performance so poor?
        drawnow;
        AUCvec(repi) = AUC;
        scorevec_cat = cat(1,scorevec_cat,scorevec);

        if do_agebins % ignore age bins for now
            for bini = 1:nbins % Should perhaps crossvalidate based on only families/subjects in bin?
                ivec_bin = agevec>=agebinlist(bini) & agevec<agebinlist(bini+1);
                scorevec(:) = NaN;
                for ki = 1:k % This loop is suspiciously fast, AUC curves suspiciously jagged -- check sum(isfinite(scorevec)) scorevec_cat and Y_cat
                    fprintf(1,'  bini=%d/%d ki=%d/%d (%s)\n',bini,nbins,ki,k,datestr(now));
                    ivec_disc = find(cvind~=ki&ivec_bin); 
                    ivec_repl = find(cvind==ki&ivec_bin);
                    U_resid = U - M*pinv(M(ivec_disc,:))*U(ivec_disc,:);
                    if 1
                        [U_tmp S_tmp V_tmp] = svd(U_resid(ivec_disc,:)*S,'econ'); 
                        D = U_resid*S*V_tmp;
                    else
                        D = U_resid;
                    end
                    mdl = fitcsvm(D(ivec_disc,1:ncomp),Y(ivec_disc),'Standardize',true,'KernelFunction','RBF','KernelScale','auto');
                    [label,score] = predict(mdl,D(ivec_repl,1:ncomp));
                    scorevec(ivec_repl) = score(:,2);
                end
                scorevecs_cat{bini} = cat(1,scorevecs_cat{bini},scorevec);
                %    sfigure(20); subplot(nr,nc,bini); [AUC SPEC SENS ACC F80] = plot_ROC(colvec(scorevec(Y==1)),colvec(scorevec(Y==0))); title(sprintf('%s predicting %s AUC=%0.2f',fstem_imaging,colnames_model{colnum},AUC),'Interpreter','none'); axis equal tight; drawnow;
                sfigure(20); 
                subplot(nr,nc,bini); 
                [AUC SPEC SENS ACC F80] = plot_ROC(colvec(scorevecs_cat{bini}(Y_cat==1)),colvec(scorevecs_cat{bini}(Y_cat==0))); 
                title(sprintf('%s predicting %s AUC=%0.3f',fstem_imaging,colnames_model{colnum},AUC),sprintf('Age [%0.1f %0.1f>',agebinlist(bini),agebinlist(bini+1)),'Interpreter','none'); 
                axis equal tight;
                drawnow;
                corrval = nancorr(scorevec,Y);
                cohensdval = cmig_tools_cohensd(nanmean(scorevec(Y==1)),nanmean(scorevec(Y==0)),nanstd(scorevec(Y==1)),nanstd(scorevec(Y==0)));
                AUCmat(bini,repi) = AUC;
                corrmat(bini,repi) = corrval;
                cohensdmat(bini,repi) = cohensdval;
                %    disp(nancorr(scorevec,Y));
                %    disp(cohensd(nanmean(scorevec(Y==1)),nanmean(scorevec(Y==0)),nanstd(scorevec(Y==1)),nanstd(scorevec(Y==0)))) 
            end 
        end
    end

    AUC_mean_allage=round(mean(AUCvec)',2);
    AUC_std_allage=round(std(AUCvec)/sqrt(1-1/k),3);

    AUC_mean_perage=round(mean(AUCmat'),2);
    AUC_std_perage=round(std(AUCmat')/sqrt(1-1/k),3);

    fontScaleFactor=1.3;
    figure(2); 
    [AUC SPEC SENS ACC F80] = plot_ROC(colvec(scorevec_cat(Y_cat==1)),colvec(scorevec_cat(Y_cat==0)));
    title({'All Ages'; 
    sprintf('AUC=%0.2f (SD=%0.3f)',AUC_mean_allage,AUC_std_allage)}); 
    ax=gca; 
    set(ax, 'FontSize', ax.FontSize * fontScaleFactor); 
    axis equal tight;

	% check output directory exist 
	if ~exist(dirname_out{des},'dir')
		mkdir(dirname_out{des});
	end

    f2=figure(2);
    saveas(f2,sprintf('%s/%s_%s_predicting_%s_across_age.png',dirname_out{des},datatype,fstem_imaging,colnames_model{colnum}));

    if do_agebins % ignore age bins for now 
        fontScaleFactor=0.8;
        for bini = 1:nbins
            figure(20); 
            subplot(nr,nc,bini); 
            [AUC SPEC SENS ACC F80] = plot_ROC(colvec(scorevecs_cat{bini}(Y_cat==1)),colvec(scorevecs_cat{bini}(Y_cat==0))); 
            title({sprintf('%s',agecat{bini}); 
            sprintf('AUC=%0.2f (SD=%0.3f)',AUC,AUC_std_perage(bini))}); 
            xlabel('');
            ylabel('');
            ax=gca; 
            set(ax, 'FontSize', ax.FontSize * fontScaleFactor); 
            axis equal tight;
            han=axes(figure(20),'visible','off');
            han.XLabel.Visible='on';
            han.YLabel.Visible='on';
            ylabel(han,'True Positive Rate (Sensitivity)');
            xlabel(han,'False Positive Rate (1-Specificity)');
        end

        f20=figure(20);
        saveas(f20,sprintf('%s/%s_%s_predicting_%s_100rep_SVMout.png',dirname_out{des},datatype,fstem_imaging,target{des}));
    end 

    fpath_out=sprintf('%s/%s_%s_predicting_%s_100rep_SVMout.mat',dirname_out{des},datatype,fstem_imaging,target{des});
    if do_agebins
        base_variables_to_save={'AUCvec','AUCmat','Y','scorevec','scorevec_cat','scorevecs_cat','X','cvind','U','fid','fidlist','corrmat','cohensdmat'};
    else
        base_variables_to_save={'AUCvec','AUCmat','Y','scorevec','scorevec_cat','scorevecs_cat','X','cvind','U','fid','fidlist'};
    end 
    save(fpath_out,base_variables_to_save{:},'-v7.3');
    
    PrintMemoryUsage
end
%keyboard


% ToDos
%   Look at weighting of svd by association with effect of interest in discovery sample
%   Use likelihood with random effects from FEMA to perform classification from trajectories and age bins (see code below)
%   Explore changing ncomp and fitcsvm options
%   Explore other "built-in" classifier methods
%   Explore mixture distributions


