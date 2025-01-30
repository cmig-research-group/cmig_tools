function prepareAtlases__5060_ABCD3_cor10
	% prepareAtlases__5060_ABCD3_cor10 prepare 5.0 or 6.0_ABCD3_cor10 atlas files for use with showVol
	%
	% This saves ShowVoslAtlases_*.mat containing a data structure 'anat' with atlas ROIs:
	%   anat.{aseg, aparc, fiber, thalamus, pauli} 
	% and saves rgb images of the atlas ROIs and other anatomical images {T1, T2, C0...}.
	%
	% Nov 2024, update for 6.0 (Note, with new file naming consistencies, this code should work with new release versions moving forward)
	%
	% Feb 2024, now identifying atlas by atlas version AND ABCD release
	%   
	% Dec 2023, relevant ROI files are now in 
	%   /space/syn65/1/data/abcd-sync/5.0/imaging_concat/voxelwise/roi (atlases) and /smri (T1, etc)
	%
	% ShowVolData still in
	%   /space/amdale/1/tmp/ABCD_cache/ABCD-Atlas/showVolData
	%
	% atlas_dspace in
	%   /home/dale/atlas_dspace/atlas_dspace_ABCD3_cor10.mat 
	%   (also seems to be present at /space/amdale/1/data/atlas_dspace/ -- identical?)
	% 
	% no atlas_dspace_ABCD3_cor10_FOD/
	%
	% OUTPUTS
	%   outputs to cfg.data.showVolData (/space/amdale/1/tmp/ABCD_cache/ABCD-Atlas/Atlas/showVolData)
	%   ROIs: showVolAtlases_6.0_ABCD3_cor10.mat
	%   rgb_<atlasName>_atlas,for atlasName in {aseg, aparc, fiber, thalamus, pauli} as of May 2021
	%   and various T1_,T2_, C0 as available"
	%

	% optimization:remove all-zero volumes and compress volumes
	% this is for ABCD2, building on atlas_dspace_ABCD2_cor10 and using
	%   code from TractQuant/TractQuant_ABCD_makejobs.m to generate rgb and labelprobs data structures
	%   relies on fs_colorlut, now in /utils

	% 9/28/2021: return to practice of normalized ROIs for Pauli/Thalamus atlases (which are at base binary ROIs)
	% 12/2/2023: update for ABCD3
	% 9/29/2024: update for new file names
	% 11/08/2024: update for 6.0 (pauli/thalamus file format changed)

	%% CHANGE HERE for new version
	atlasVersion = '6.0_ABCD3_cor10'; % '6.0_ABCD3_cor10' or '5.0_ABCD3_cor10'
	processDate = datestr(now,'yyyy/mm/dd');

	adjustCenter = true; %false;
	if adjustCenter
  		disp('Adjusting volume center')
  		rcscenter_1mm = [100 99 130]; % ** Emprically, the  ABCD3 atlas_dspace T1 has its center at 100 99 130
  	% this doesn't match current 101 101 131 in M_atl
	end

	%% CHANGE HERE if locations differ
	source = '/space/syn65/1/data/abcd-sync/6.0/imaging_concat/voxelwise/roi';
	sourceAnat = '/space/syn65/1/data/abcd-sync/6.0/imaging_concat/voxelwise/smri';
	atlas_dspace = '/home/dale/tmp/atlas_dspace';
	%%

	% ===============================================================
	persistent AD

	atlasVersion = validateAtlasVersion(atlasVersion);
	atlasVersionVar = atlasVersion(atlasVersion ~= '.'); 
	[abcdRelease, atlasStr, atlasSuffix] = parseAtlasVersion(atlasVersion);
	atlasID = [abcdRelease ' ' atlasStr];

	useAtlasDspace = true; % previously for T1, T2 etc we used 'means' calculated by Diliana--for now take directly from atlas_dspace
	meandir = '';

	cfg = abcdConfig('showVol');
	outdir = fullfile(cfg.data.showVolData, 'Atlas', atlasVersion,'');
	if ~exist(outdir,'dir')
  	mkdir(outdir)
	end
	%outdir = fullfile('~/Images', 'Atlas','') %for TESTING

	atlases = {'aseg','aparcaseg','fiber','thalamus','pauli'};

	files = {fullfile(source,meandir,'aseg.mat'), ... %FIXME: 09/2024: filenames simplified to e.g. 'aseg.mat'
    		fullfile(source,meandir,'aparcaseg.mat'), ...
    		fullfile(source,meandir,'fiber.mat'), ...
    		fullfile(source, 'thalamus_abcd3_cor10.mat'), ... %NOTE these do not depend on the ABCD release
    		fullfile(source, 'pauli_subcortnuclei_abcd3_cor10.mat')};
	
	notes = {[atlasVersion ' ASEG ROIs, ' processDate], ...
		  	[atlasVersion ' APARCASEG ROIs, ' processDate], ...
		  	[atlasVersion ' FIBER ROIs, ' processDate], ...
		  	[atlasVersion ' Thalamus atlas, ' processDate], ...
		  	[atlasVersion ' Pauli atlas, ' processDate]};


	%whether to normalize probabilities to max within each ROI. The reason for this is that ROI outlines
	% in showVol are drawn at a fixed probability threshold. Some ROIs (typically smlaller or more variable ones)
	% may have peak probabilities substantially les than one so they may not be drawn. 
	% In general, most ROIs have peaks >0.85, with many close to 1, so not an issue
	
	%doNormalize = [0 0 0 1 1]; %Prior to Nov 2024 we normalized thalamus & pauli
	doNormalize = [0 0 0 0 0]; %Nov 24 on, no normalization--add as switch in showVol if needed


	overlayColors = {[1 0 1], [1 0 .5], [0 1 .5],[],[]};
	defaultThreshold = [0.7 0.7 0.25 0.5 0.5];

	atlasFname = fullfile(outdir, ['showVolAtlases_' atlasVersion '.mat']);

	if isempty(AD)
  	disp('Loading atlas_dspace...')
  	load(fullfile(atlas_dspace,'atlas_dspace_ABCD3_cor10.mat')); %NOTE, this does not depend on release version?
  	AD = atlas_dspace; %for convenience
	else
  	disp('Using Cached atlas_dspace')
	end

	% volume anatomical parameters
	Mvxl2lph = M_RAS_TO_LPH * AD.M_atl; %M_atl is in RAS

	if adjustCenter
  	%shift to actual anatomical center of new 1mm atlas
  	Mvxl2lph(1,4) = -rcscenter_1mm(2);
  	Mvxl2lph(2,4) = rcscenter_1mm(3);
  	Mvxl2lph(3,4) = rcscenter_1mm(1);
	else %use center from M_atl
  	rcscenter_1mm = [Mvxl2lph(3,4) -Mvxl2lph(1,4) Mvxl2lph(2,4)];
	end

	voxSize = 1;

	%% ROI Probabilities %%

	%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
	% loop over ABCD ROIs
	for iA = 1:3
  	atlas = atlases{iA};
  	disp(atlas)
	
  	A = load(files{iA});
	
  	%remove empty volumes
  	try
    	probMax = max(A.volmat_prob,[],1:3);
  	catch
    	probMax = max(max(max(A.volmat_prob,[],1),[],2),[],3);
  	end
  	goodvols = squeeze(probMax)>0;
  	%A.volmat_prob = A.volmat_prob ./ probMax; %no normalization any longer
  	A.roiabbreviations = [];

  	%remove derivative & redundant ROIs
  	if strcmp(atlas,'fiber')
    	%remove derivative ROIs from fiber atlas:
    	% 38-42 are all fibers, R/L all fibers (no cc), and R/L all fibers
    	% ROIs 36-37 are R/L frontal fornix--keep now
    	goodvols(38:end)=false;
    	disp('Removing allFibers ROIs from fiber atlas')
  	end

  	%remove all of the aseg rois from aparcaseg. In principle, does fitting jointly change them from the aseg-only data?
  	%  Checking the probMax for aseg ROIs, they are identical in the aseg and aparcaseg atlases, from
  	%  which we assume the probs are the same (except for CC_All and Ventricles_All...)
  	if strcmp(atlas,'aparcaseg')
    	ctxRois = contains(A.roinames,'ctx-');
    	goodvols(~ctxRois) = false;
  	end

  	%remove empty/unneded volumes
  	if sum(~goodvols)
    	fprintf(' Removing %d empty/unneeded vols:\n',sum(~goodvols))
    	disp(A.roinames(~goodvols(:)))
    	A.volmat_prob = A.volmat_prob(:,:,:,goodvols);
    	A.roicodes    = A.roicodes(goodvols);
    	A.roinames    = A.roinames(goodvols);
    	A.roirgb      = A.roirgb(goodvols,:);
  	end  

  	%re-sort ROIs, add abbreviations
  	if strcmp(atlas,'aparcaseg')
    	disp('re-sorting aparc cortical ROIs')
    	aparcAbbrev = {'bSTS','cACC','cMFC','Cun','Ent','Fus','IPC','ITC','isthCing','lOcc','lOFC',...
      	'Ling','mOFC','mTem','paraHipp','paraCen','POperc','POrb','PTri','periCalc','postCen','PCC',...
      	'preCen','preCun','rACC','rMFC','SFC','SPC','STC','SMar','FP','TP','tTem','Ins'};
    	roiabbreviations = [aparcAbbrev aparcAbbrev];
    	order = [1:34;35:68];
    	order = order(:);
    	A.volmat_prob = A.volmat_prob(:,:,:,order);
    	A.roicodes    = A.roicodes(order);
    	A.roinames    = A.roinames(order);
    	A.roirgb      = A.roirgb(order,:);
    	A.roiabbreviations = roiabbreviations(order);
  	end

  	prob = single(A.volmat_prob);
  	if doNormalize(iA)
    	disp('Normalizing Probabilities')
    	prob = prob ./ max(prob,[],1:3);
  	end

  	%set up anat struct representing an atlas
  
  	atl=[];
  	atl.overlayColor = overlayColors{iA};
  	atl.prob = prob;
  	atl.roicodes = A.roicodes;
  	atl.roinames = A.roinames;
  	atl.roiabbreviations = A.roiabbreviations;
  	atl.roicolors = A.roirgb;
  	atl.subjlist = A.subjlist;
  	atl.size = size(atl.prob(:,:,:,1));
  	atl.roi_text = [];
  	atl.file = files{iA};
  	atl.probabilityThreshold = defaultThreshold(iA);
  	atl.showNames = false;
  	atl.outlineColor = 'w';

  	%make a list of ROI names for use in UI for highlighting ROIs
  	atl.uiNames = {'None'};
  	atl.uiRoiIdx = 0;
  	for iV = 1:size(atl.prob,4)
    	atl.uiNames{end+1} = atl.roinames{iV};
    	atl.uiRoiIdx(end+1) = iV;
  	end
  	%for UI state
  	atl.uiAtlasName = [atlas ' [' atlasVersion ']'];
  	atl.uiRoiOverlayIdx = []; %indices into roi volume of selected ROIs to overlay
  	atl.uiRoiOverlayImg = [];
  	atl.uiRoiOverlaySelected = []; %indices into uiNames

  	%add additional shape info
  	atl.Mvxl2lph = Mvxl2lph;
  	atl.voxelSize = voxSize;
  	atl.dimr=size(atl.prob,1); atl.dimc = size(atl.prob,2); atl.dimd = size(atl.prob,3);
  	atl.vx=atl.voxelSize; atl.vy=atl.voxelSize; atl.vz=atl.voxelSize;
  	atl.rcscent = rcscenter_1mm(:);
  	lphcent = atl.Mvxl2lph * [atl.rcscent(:); 1];
  	atl.lphcent = lphcent(1:3);

  	atl.datePrepared = datestr(now);
  	atl.notes = notes{iA};
  	atl.atlasVersion = atlasVersion;

  	anat.(atlas) = atl;

  	% write out RGB image
  	rgb_vol = ctx_mgh2ctx(single(zeros(atl.size)), M_LPH_TO_RAS*Mvxl2lph); %starter volume
  	rgb_vol.imgs = A.volmat_prob_rgbsum;

  	%add additional shape info
  	rgb_vol.Mvxl2lph = Mvxl2lph;
  	rgb_vol.voxelSize = voxSize;
  	rgb_vol.dimr=size(rgb_vol.imgs,1); rgb_vol.dimc = size(rgb_vol.imgs,2); rgb_vol.dimd = size(rgb_vol.imgs,3);
  	rgb_vol.vx=rgb_vol.voxelSize; rgb_vol.vy=rgb_vol.voxelSize; rgb_vol.vz=rgb_vol.voxelSize;
  	rgb_vol.rcscent = rcscenter_1mm(:);
  	lphcent = rgb_vol.Mvxl2lph * [rgb_vol.rcscent(:); 1];
  	rgb_vol.lphcent = lphcent(1:3);

  	rgb_vol.name = [atlas ' rgb [' atlasVersion ']'];
  	rgb_vol.datePrepared = datestr(now);
  	rgb_vol.notes = [notes{iA} ' RGB'];
  	rgb_vol.atlasVersion = atlasVersion;
  	varName = [atlas '_rgb_' atlasVersionVar];
  	eval([varName ' = rgb_vol;']); %assign to atlas-specific name
  	save(fullfile(outdir,[atlas '_rgb_' atlasVersion '.mat']), varName)

	end %loop over ABCD ROI atlases

	%% save ABCD PRERENDERED IMAGES - T1, T2, FA, C0 %%

	if useAtlasDspace
	
		disp(' save T1 (atlas_dspace)')
		clear vol
		vol.imgs = AD.muvols_T1;
		vol.Mvxl2lph = Mvxl2lph;
		vol.dimc = 200; vol.dimr=200; vol.dimd=260;
		vol.vx=1;vol.vy=1;vol.vz=1;
		vol.lphcent=[0 0 0]';
		vol.atlasVersion = atlasVersion;
		vol.name = ['T1 [' atlasID ']'];
		varName = ['T1_' atlasVersionVar];
		eval([varName '= vol;'])
		save(fullfile(outdir,['T1_' atlasVersion '.mat']), varName)
		
		disp(' save T2 (atlas_dspace)')
		clear vol
		vol.imgs = AD.muvols_T2;
		vol.Mvxl2lph = Mvxl2lph;
		vol.dimc = 200; vol.dimr=200; vol.dimd=260;
		vol.vx=1;vol.vy=1;vol.vz=1;
		vol.lphcent=[0 0 0]';
		vol.atlasVersion = atlasVersion;
		vol.name = ['T2 [' atlasID ']'];
		varName = ['T2_' atlasVersionVar];
		eval([varName '= vol;'])
		save(fullfile(outdir,['T2_' atlasVersion '.mat']), varName)
		
		disp(' save CO (atlas_dspace)')
		vol = AD.vol_CO_rgb;
		vol.Mvxl2lph = Mvxl2lph;
		vol.atlasVersion = atlasVersion;
		vol.name = ['diffusion direction [' atlasID ']'];
		varName = ['CO_' atlasVersionVar];
		eval([varName '= vol;'])
		save(fullfile(outdir,['CO_' atlasVersion '.mat']), varName)
	
	else
  		error('This non atlas_dspace code hasn''t been updated')
  		%T1
  		disp(' save T1')
  		fname = fullfile(sourceAnat,'volmat_nu.mat');
  		tmp = load(fname);
  		vol = ctx_mgh2ctx(single(tmp.vol_mean), M_LPH_TO_RAS*Mvxl2lph);
  		vol.lphcent = [0 0 0]'; %have to do this, because we are using 0-based M, while the function assumes 1-based
  		vol.atlasVersion = atlasVersion;
  		vol.name = 'T1 [ABCD3, volmat]';
  		T1_50_ABCD3_cor10 = vol;
  		save(fullfile(outdir,['T1_' atlasVersion '.mat']), 'T1_50_ABCD3_cor10')
		
  		%T2
  		disp(' save T2')
  		fname = fullfile(sourceAnat,'volmat_T2.mat');
  		tmp = load(fname);
  		vol = ctx_mgh2ctx(single(tmp.vol_mean), M_LPH_TO_RAS*Mvxl2lph);
  		vol.lphcent = [0 0 0]'; %have to do this, because we are using 0-based M, while the function assumes 1-based
  		vol.atlasVersion = atlasVersion;
  		vol.name = 'T2 [ABCD3, volmat]';
  		T2_50_ABCD3_cor10  = vol;
  		save(fullfile(outdir,['T2_' atlasVersion '.mat']), 'T2_50_ABCD3_cor10')
		
  		%FA
  		disp(' save FA')
  		fname = fullfile(sourceAnat,'volmat_FA.mat');
  		tmp = load(fname);
  		vol = ctx_mgh2ctx(single(tmp.vol_mean), M_LPH_TO_RAS*Mvxl2lph);
  		vol.lphcent = [0 0 0]'; %have to do this, because we are using 0-based M, while the function assumes 1-based
  		vol.atlasVersion = atlasVersion;
  		vol.name = 'FA [ABCD3, volmat]';
  		FA_50_ABCD3_cor10  = vol;
  		save(fullfile(outdir,['FA_' atlasVersion '.mat']), 'FA_50_ABCD3_cor10')
	
	end


	%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
	%% Thalamus & Pauli (MNI Atlases from other studies)

	%hand-code abbreviations for Thalamus. 
	thalamusAbbrev = {'PUL','A','MD','VLD','CllpmPUL','VA','VLV','PUL','A','MD','VLD','CllpmPUL','VA','VLV'}';

	for iA = 4:5
  		if ~exist(files{iA},'file')
			continue
		end

    	atlas = atlases{iA};
    	disp(atlas)
    	A = load(files{iA});
	
    	nroi = size(A.vol_labels_reg,4);
	
    	prob = single(A.vol_labels_reg);
    	if doNormalize(iA)
    		disp('Normalizing Probabilities')
    		prob = prob ./ max(prob,[],1:3);
    	end
	
    	atl=[];
    	atl.prob = prob;
	
    	LUT = A.([atlas '_lut']);
		
    		%there are some differences between pauli and thalamus naming 
    	if strcmp(atlas,'thalamus')
    	  	name = LUT.LabelName;
    	  	atl.roiabbreviations = thalamusAbbrev;
    	else
    	  	name = LUT.Name;
    	  	atl.roiabbreviations = LUT.Abbreviation;
    	end
    	atl.roinames = name;
	
    	atl.roicolors = [LUT.R LUT.G LUT.B];
    	m = max(atl.roicolors(:));
    	if m > 1
    		atl.roicolors = atl.roicolors / 255; %rescale to matlab [0 1]
    	end
    	atl.size = size(atl.prob(:,:,:,1));
    	atl.roi_text = []; %for text display of current ROI names in figure
    	atl.file = files{iA};
    	atl.overlayColor = []; %unused
    	atl.probabilityThreshold = 0.5;
    	atl.showNames = true;

    	atl.uiNames = {'None'}; atl.uiRoiIdx = 0;
    	for iV = 1:nroi
    		atl.uiNames{end+1} = sprintf('%s (%s)',atl.roinames{iV}, atl.roiabbreviations{iV}) ;
    		atl.uiRoiIdx(end+1) = iV;
    	end
    	atl.uiRoiOverlayIdx = []; %indices into roi volume of selected ROIs to overlay
    	atl.uiRoiOverlayImg = [];
    	atl.uiRoiOverlaySelected = []; %indices into uiNames
	
    	%add additional shape info
    	atl.Mvxl2lph = Mvxl2lph;
    	atl.voxelSize = voxSize;
    	atl.dimr=size(atl.prob,1); atl.dimc = size(atl.prob,2); atl.dimd = size(atl.prob,3);
    	atl.vx=atl.voxelSize; atl.vy=atl.voxelSize; atl.vz=atl.voxelSize;
    	atl.rcscent = rcscenter_1mm(:);
    	lphcent = atl.Mvxl2lph * [atl.rcscent(:); 1];
    	atl.lphcent = lphcent(1:3);
		
    	atl.notes = notes{iA};
    	atl.dataPrepared = datestr(now);
    	atl.atlasVersion = atlasVersion;
	
    	anat.(atlas) = atl;
	
    	%% also save some high-res T1 and T2 images and rgb atlas if present
    	if isfield(A,'vol_T1_reg')
    	    disp(' save T1')
    	    vol = ctx_mgh2ctx(single(A.vol_T1_reg), M_LPH_TO_RAS*Mvxl2lph);
    	    vol.lphcent = [0 0 0]'; %have to do this, because we are using 0-based M, while the function assumes 1-based
    	    vol.atlasVersion = atlasVersion;
    	    vol.name = [atlas ' atlas T1 [' atlasID ']'];
			vol.imgs(isnan(vol.imgs)) = 0; %remove NaNs
    	    varName = [atlas '_T1'];
    	    eval([varName ' = vol;']); %assign to atlas-specific name
    	    save(fullfile(outdir,[atlas '_T1_' atlasVersion '.mat']), varName);
    	end
    	if isfield(A,'vol_T2_reg')
    		disp(' save T2')
    		vol = ctx_mgh2ctx(single(A.vol_T2_reg), M_LPH_TO_RAS*Mvxl2lph);
    		vol.lphcent = [0 0 0]';
    		vol.atlasVersion = atlasVersion;
    		vol.name = [atlas ' atlas T2 [' atlasID ']'];
			vol.imgs(isnan(vol.imgs)) = 0; %remove NaNs
    		varName = [atlas '_T2_' atlasVersionVar];
    		eval([varName ' = vol;']); %assign to atlas-specific name
    		save(fullfile(outdir,[atlas '_T2_' atlasVersion '.mat']), varName);
    	end
    	if isfield(A,'vol_labels_reg_rgb')
    		disp(' save rgb')
    		vol = A.vol_labels_reg_rgb;
    		vol.Mvxl2lph = Mvxl2lph;
    		vol.lphcent = [0 0 0]';
    		vol.atlasVersion = atlasVersion;
    		vol.name = [atlas ' atlas rgb [' atlasID ']'];
			vol.imgs(isnan(vol.imgs)) = 0; %remove NaNs
    		varName = [atlas '_rgb_' atlasVersionVar];
    		eval([varName ' = vol;']); %assign to atlas-specific name
    		save(fullfile(outdir,[atlas '_rgb_' atlasVersion '.mat']),varName);
    	end
    	if isfield(A,'vol_labels_T1_reg_rgb')
    		disp(' save T1 + rgb')
    		vol = A.vol_labels_T1_reg_rgb;
    		vol.Mvxl2lph = Mvxl2lph;
    		vol.lphcent = [0 0 0]';
			vol.imgs(isnan(vol.imgs)) = 0; %remove NaNs
    		vol_labels_rgb.atlasVersion = atlasVersion;
    		vol_labels_rgb.name = [atlas ' atlas rgb+T1 [' atlasID ']'];
    		varName = [atlas '_rgb_T1_' atlasVersionVar];
    		eval([varName ' = vol;']); %assign to atlas-specific name
    		save(fullfile(outdir,[atlas '_rgb_T1_' atlasVersion '.mat']),varName);
    	end
    	if isfield(A,'vol_labels_T2_reg_rgb')
    		disp(' save T2 + rgb')
    		vol = A.vol_labels_T2_reg_rgb;
    		vol.Mvxl2lph = Mvxl2lph;
    		vol.lphcent = [0 0 0]';
			vol.imgs(isnan(vol.imgs)) = 0; %remove NaNs
    		vol_labels_rgb.atlasVersion = atlasVersion;
    		vol_labels_rgb.name = [atlas ' atlas rgb+T2 [' atlasID ']'];
    		varName = [atlas '_rgb_T2_' atlasVersionVar];
    		eval([varName ' = vol;']); %assign to atlas-specific name
    		save(fullfile(outdir,[atlas '_rgb_T2_' atlasVersion '.mat']),varName);
    	end
		
	end %pauli and thalamus

	%% compress probabilities and unpack atlases into individual variables (so we can load as needed)
	% FIXME: could also save with '-struct' argument in newer matlab
	for iA = 1:length(atlases)
  		fprintf('%s: ', atlases{iA})
  		anat.(atlases{iA}).prob = compressVol(anat.(atlases{iA}).prob);
  		anat.(atlases{iA}).atlasVersion = atlasVersion;
  		eval([atlases{iA} ' = anat.' atlases{iA} ';'])
	end

	%% set up ROI overlay defaults
	params.showOverlay = true; %toggle for overlay
	params.binarizeOverlay = false; %entire ROI > threshold rendered same color, false: show probability
	params.showOutline = true; %outline ROIs
	params.outlineColor = [1 1 1];
	params.outlineWidth = 1;
	params.overlayAlpha = 0.4;
	params.imageFade = 1; %start with no fade
	params.showRoiLabel = true;
	params.currentAtlas = '';
	params.generatedDate = datestr(now);
	params.notes = [atlasVersion ', '  processDate];

	%% SAVE
	disp('save atlases')
	save(atlasFname,'atlases','params',atlases{:},'-v7.3')

	%% create symmetrized atlases
	% TODO
end 

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Helper Functions - now as standalone functions

%%% compression/expansion
%%simple compression method: only store non-zero values, as single precision. since ROIs are very sparse
%%for the most part, the space savings (thus load time improvement) is
%%large. E.g. aseg only 1.4% of voxels are non-zero
%function c = compressVol(vol)
	%s = size(vol);
	%if prod(s) >= 2^32, error('Volume too large to compress with uint32 indexes'), end
	%nz = cast(find(vol(:)~=0),'uint32'); %index of non-zero voxels
%
	%c = {s nz(:) single(vol(nz))}; %volume size, non-zero-voxel indices, non-zero-voxel values
%
%	fprintf('(compressed to %s%% of original)\n',num2str(100*2*length(nz)/prod(s),2)); %factor of 2 since we save index %and values
%
	%% companion expansion function when loading atlas. Not used here, but FYI.
	%%  e.g. 
	%% if iscell(aseg.prob)
	%%     aseg.prob = expandVol(aseg.prob);
	%% end
%end 
%
%function vol = expandVol(c)
	%vol = zeros(c{1},'single');
	%vol(c{2}) = c{3};
%end 