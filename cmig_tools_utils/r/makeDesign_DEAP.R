library(plyr)
library(dplyr)
library(tidyverse)
library(psych)
library(Matrix)
library(ordinal)
library(pracma)
library(data.table)
library(arrow)

makeDesign_DEAP <- function(jsonFile, outpath, outfile, study='ABCD6') 
{

	# jsonFile = input JSON spec from the front end
	# outpath  = location where the design matrix should be saved
	# outfile  = name of the design matrix
	# study    = either 'ABCD6' or 'HBCD'

	# To Dos:
	# Add data.table and jsonlite to yml
	# Handle transformations
	# Handle filtering, if required
	# Remember to do QC filtering for imaging variables
	# Add modality to JSON specification file so we know which QC filters to apply
	# Additionally include filtering by incidental findings both in the front end and back end
	# Remember to account for / handle non-FSE random effects
	# Remember to invoke MATLAB/R functions to create smooth basis functions

	# Load necessary libraries
	for (p in c('plyr','dplyr','tidyverse','psych','Matrix','ordinal','pracma', 'data.table', 'jsonlite')) {
	    if(!eval(parse(text=paste("require(",p,")")))) {
	            install.packages(p)
	            lapply(p,library,character.only=TRUE)
	    }
	}

	# Read the JSON job spec
	jobSpec <- read_json(jsonFile)

	# Work on fixed effects
	# ======================

	# Extract names of the variables that we need to read
	varMatrix <- sapply(jobSpec$params$fixed_effects$vars, "[")

	# Names of fixed effects
	varNames  <- as.character(varMatrix["name", ])

	# Which fixed effects are longitudinal
	isLong <- as.logical(varMatrix["longitudinal", ])
	delta  <- varNames[isLong]
	if (isempty(delta))
	{
		delta <- NULL
	}

	# Which fixed effects are to be smoothed; not currently handled
	isSmooth <- as.logical(varMatrix["smoothing", ])

	# Extract interaction terms
	interact <- paste(jobSpec$params$fixed_effects$interaction_terms)
	if (isempty(interact) | interact == "")
        {
                interact <- NULL
        }

	# Extract mediator information
	mediator <- jobSpec$params$fixed_effects$mediator
	if (isempty(mediator) | mediator == "")
        {
                mediator <- NULL
        }

	# Extract which variables need to be de-meaned
	demean <- as.logical(jobSpec$params$fixed_effects$mean_center)

	# Work on random effects
	# ======================

	# Extract random effect names excluding E
	RFXNames <- setdiff(rownames(as.matrix(jobSpec$params$random_effects)), "E")

	# Set family ID
	if ("F" %in% RFXNames)
	{
		familyID <- jobSpec$params$random_effects$F
	} else
	{
		familyID <- 'participant_id'
	}

	# Set participant ID
	if ("S" %in% RFXNames)
	{
		idvar <- jobSpec$params$random_effects$S
	} else
	{
		idvar <- 'participant_id'
	}

	# Create "NDA" data frame
	# =======================

	# Hard coding where the data should be read from
	# dirTabulated <- "/abcd/6_0/bids/rawdata/phenotype"
	dirTabulated <- "/space/ceph/1/ABCD/release/bids/rawdata/phenotype"

	# This part needs to be replaced with a call to the database / data dictionary and know which files to read
	# Currently hardcoded to read varNames
	tmp_dynamic_info <- read_parquet(file.path(dirTabulated, "ab_g_dyn.parquet"), as_data_frame=FALSE)
	tmp_static_info  <- read_parquet(file.path(dirTabulated, "ab_g_stc.parquet"), as_data_frame=FALSE)

	# Extract schema out and figure out which are non numeric columns
	# Define what counts as numeric types
    	numeric_types <- c("int8", "int16", "int32", "int64", "uint8", "uint16", "uint32", "uint64", "float", "double", "
    decimal")

    	# For dynamic info
	schema      <- tmp_dynamic_info$schema
	field_names <- names(schema)
	field_types <- sapply(schema$fields, function(f) f$type$ToString())
	
	# Find non-numeric columns
    	non_num_dynamic_cols <- field_names[!field_types %in% numeric_types]
	
	# Repeat for static cols
	schema_stc_info <- tmp_static_info$schema
	schema      <- tmp_static_info$schema
        field_names <- names(schema)
        field_types <- sapply(schema$fields, function(f) f$type$ToString())         
        
        # Find non-numeric columns
        non_num_static_cols <- field_names[!field_types %in% numeric_types]

	# Make a list of all non numeric categorical columns
	all_non_num_cols <- c(non_num_dynamic_cols, non_num_static_cols)

	# Join
	nda <- join(as.data.frame(tmp_dynamic_info), as.data.frame(tmp_static_info), by=c("participant_id"))
	# nda <- join(as.data.frame(tmp_dynamic_info), as.data.frame(tmp_static_info), by=c("participant_id", "session_id"))
	
	# This is where filtering needs to be done based on what imaging variables are selected and QC filters

	# create path to outfile if doesn't exist
	# outpath <- dirname(outfile)
	if ( ! dir.exists(outpath) ) {
		dir.create(outpath, recursive=TRUE)
	}

    if (study=='ABCD5') {
		agevar <- 'interview_age'
        # idvar <- 'src_subject_id'
        visitvar <- 'eventname'
        # familyID <- 'rel_family_id'
    } else if (study=='ABCD6') {
        agevar <- 'ab_g_dyn__visit_age'
        # idvar <- 'participant_id'
        visitvar <- 'session_id'
    } else if (study=='HBCD') {
            agevar <- 'candidate_age'
            # idvar <- 'participant_id'
            visitvar <- 'session_id'
            # familyID <- 'rel_family_id'
    }

    allvars <- varNames # c(contvar,catvar,delta)
    contvar <- setdiff(allvars, c(all_non_num_cols, delta))
    catvar  <- setdiff(allvars, c(contvar, delta))

    nda[,'age'] <- nda[, agevar]
    nda <- nda[,c(idvar, visitvar, familyID, 'age', allvars)]
    nda <- nda[complete.cases(nda),]
    time <- c("00A", "02A", "04A", "06A")
    idx_time <- grep(paste0(time, collapse='|'), nda[, visitvar])
    nda <- nda[idx_time,]	
    # get subject ids at that time point
    subjid <- as.character(nda[, idvar])

	# filter by specific variable and value
	# if ( !is.null(filtervar) ) {
	#	nda <- filter(nda, filtervar[1] == filtervar[2])
	# }
	
	# calculate deltas
	if ( !is.null(delta)) {
		uniq_subj <- unique(subjid)
		for (d in delta) {
			vec0 <- rep(NA,dim(nda)[1]);
			vecD <- rep(NA,dim(nda)[1]);
			for (i in 1:length(uniq_subj)) {
				idx_subj <- which(nda[, idvar]==uniq_subj[i])
				idx_tmp <- which(nda$age[idx_subj]==min(nda$age[idx_subj]))
                idx_subj0 <- idx_subj[idx_tmp]	
				vec0[idx_subj] <- nda[idx_subj0,d]
				vecD[idx_subj0] <- 0
				idx_fu <- setdiff(idx_subj, idx_subj0)
				vecD[idx_fu] <- nda[idx_fu,d] - nda[idx_subj0,d]
			}
			name_base <- paste0(d, '_base')
			name_delta <- paste0(d, '_delta')
			nda[, name_base] <- vec0
			nda[, name_delta] <- vecD
		}
	}
	
    # Transformations to be handled later
    quadratic <- NULL
	if (!is.null(quadratic)){
		quadvars=NULL
		for (i in 1:length(quadratic)){
			if (quadratic[i] %in% delta){
				nda[,paste0(quadratic[i],'_base:',quadratic[i],'_delta')]<-(nda[,paste0(quadratic[i],'_base')]*nda[,paste0(quadratic[i],'_base')])*2
				nda[,paste0(quadratic[i],'_sq')]<-nda[,paste0(quadratic[i],'_delta')]^2
				if (demean==TRUE) {
					nda[,paste0(quadratic[i],'_base:',quadratic[i],'_delta')]<-scale(nda[,paste0(quadratic[i],'_base:',quadratic[i],'_delta')], center=TRUE,scale=FALSE)
					nda[,paste0(quadratic[i],'_sq')]<-scale(nda[,paste0(quadratic[i],'_sq')], center=TRUE,scale=FALSE)
				}
				quadvars<-c(quadvars, paste0(quadratic[i],'_base:',quadratic[i],'_delta'), paste0(quadratic[i],'_sq'))
			} else {
				if ((quadratic[i] %in% contvar)==FALSE){
					warning (paste0(quadratic[i],' is not in contvar list. Linear effect NOT included in model. Please fix and rerun'))
				}
				nda[,paste0(quadratic[i],'_sq')]<-nda[,quadratic[i]]^2
				if (demean==TRUE){
					nda[,paste0(quadratic[i],'_sq')]<-scale(nda[,paste0(quadratic[i],'_sq')],center=TRUE,scale=FALSE)
				}
				quadvars<-c(quadvars, paste0(quadratic[i],'_sq'))
			}
		}
	}

	####################################
	# function to extract variables from nda	

	varextract = function(data, varname, index, dummy, int){
		if (dummy==1 && int==FALSE){
			varout <- dummy.code(data[,varname][index])
		} else if (dummy==1 && int==TRUE){
			varout <- dummy.code(data[,varname][index])
			if (ncol(varout)==2){
				dumname<-paste0(varname,'_',colnames(varout)[2])
				varout <- data.frame(varout[,2:ncol(varout)])
				colnames(varout)<-dumname
			} else {
				varout <- data.frame(varout[,2:ncol(varout)])
			}
		} else if (dummy==0){
			varout <- data[,varname][index]
		}
		return(varout)
	}

	####################################
	# extract variables 
	
	if ( !is.null(delta)) {
		name_base <- paste0(delta, '_base')
								name_delta <- paste0(delta, '_delta')
		d0 <- lapply(name_base, varextract, data=nda, index=idx_time, dummy=0, int=int)
		if (demean==TRUE){
			d0<-apply(data.frame(d0), 2, function(y) y - mean(y, na.rm=T))
		}
		dD <- lapply(name_delta, varextract, data=nda, index=idx_time, dummy=0)
		if (demean==TRUE){
			dD<-apply(data.frame(dD), 2, function(y) y - mean(y, na.rm=T))
		}
		nda[,name_base]<-data.frame(d0)
		nda[,name_delta]<-data.frame(dD)		
	}
	
	if (demean==TRUE){
		nda[,contvar]<-apply(data.frame(nda[,contvar]), 2, function(y) y - mean(y, na.rm=T))
	}
	
	if (!is.null(catvar)){
		for (i in 1:length(catvar)){
			if (!is.factor(nda[,catvar[i]])){
				warning (paste0(catvar[i],' is not a factor. Being transformed into a factor.'))
				nda[,catvar[i]]<-factor(nda[,catvar[i]])
			}
		}
	} 

	if (!is.null(quadratic)){
		quadvars=NULL
	  for (i in 1:length(quadratic)){
		if (demean==FALSE){ #Must demean variables before adding the quadratic
			warning ('Warning: quadratic in model, but demean=FALSE. Could lead to colinearity problems')
		}
		if (quadratic[i] %in% delta){
			nda[,paste0(quadratic[i],'_base:',quadratic[i],'_delta')]<-(nda[,paste0(quadratic[i],'_base')]*nda[,paste0(quadratic[i],'_base')])*2
			nda[,paste0(quadratic[i],'_sq')]<-nda[,paste0(quadratic[i],'_delta')]^2
		if (demean==TRUE) {
			nda[,paste0(quadratic[i],'_base:',quadratic[i],'_delta')]<-scale(nda[,paste0(quadratic[i],'_base:',quadratic[i],'_delta')], center=TRUE,scale=FALSE)
			nda[,paste0(quadratic[i],'_sq')]<-scale(nda[,paste0(quadratic[i],'_sq')], center=TRUE,scale=FALSE)
		}
		quadvars<-c(quadvars, paste0(quadratic[i],'_base:',quadratic[i],'_delta'), paste0(quadratic[i],'_sq'))
		} else {
			if ((quadratic[i] %in% contvar)==FALSE){
				warning (paste0(quadratic[i],' is not in contvar list. Linear effect NOT included in model. Please fix and rerun'))
			}
			if (demean==FALSE){ #Must demean variables before adding the quadratic
				warning ('Warning: quadratic in model, but demean=FALSE. Could lead to colinearity problems')
			}
				nda[,paste0(quadratic[i],'_sq')]<-nda[,quadratic[i]]^2
				if (demean==TRUE){
					nda[,paste0(quadratic[i],'_sq')]<-scale(nda[,paste0(quadratic[i],'_sq')],center=TRUE,scale=FALSE)
				}
				quadvars<-c(quadvars, paste0(quadratic[i],'_sq'))
			}
		}
	}
	
	
	## MAKING COVARIATE STRING FOR MODEL
	if (!is.null(delta)){
		covariate_str = paste(c(name_base,name_delta,contvar,catvar), collapse=' + ')
	} else {
		covariate_str = paste(c(contvar,catvar), collapse=' + ')
	}
	
	if (!is.null(quadratic)){
		covariate_str = paste0(covariate_str, ' + ', paste(c(quadvars), collapse=' + '))
	}
	
	if (!is.null(interact)){
		covariate_str = paste0(covariate_str, ' + ', paste(c(interact), collapse=' + '))
	}
	
	nda[,'dummy_var']<-rand(dim(nda)[1],1)
	modelmat<-model.matrix(as.formula(paste0('dummy_var ~ ',covariate_str)),nda)
	
	### Check for rank deficiency #####
	
	R <- rankMatrix(modelmat)
	
	if (R[1]<ncol(modelmat)){
		modelmat<-drop.coef(modelmat)
	}
	
	colnames(modelmat)[colnames(modelmat) == "(Intercept)"] = "intercept"
	
	modelmat_names = c(which(colnames(modelmat) != "intercept"), which(colnames(modelmat) == "intercept"))
	modelmat = modelmat[,modelmat_names]
	
	####################################

    # define the first four columns required by FEMA_wrapper
    outmat<-nda[,c(idvar, visitvar, familyID, 'age')]
	
    # add delta if supplied
	if (!is.null(delta)){
		delta_df<-data.frame(nda[,c(name_base, name_delta)])
		outmat<-data.frame(outmat, delta_df)

			modelmat<-modelmat[,-c(which(colnames(modelmat) %in% name_base))]
			modelmat<-modelmat[,-c(which(colnames(modelmat) %in% name_delta))]
	 }
	
	outmat<-cbind(outmat,modelmat)
	
    # add mediatator if supplied
	if (!is.null(mediator)){
		#Check mediator is included already in design matrix
		if (mediator %in% colnames(outmat)==FALSE){
			stop('ERROR! Mediator NOT in contvar list')
		}
		mediator_ind<-match(mediator,colnames(outmat))
		outmat_reduced<-outmat[,-c(mediator_ind)]
		outfile_reduced<-str_replace_all(outfile, '.txt', '_reduced.txt')
		write.table(outmat_reduced, outfile_reduced, col.names=TRUE, row.names=FALSE, sep='\t')
		
		ncols<-length(colnames(outmat))
		allind<-1:ncols
		allind<-c(allind[-c(mediator_ind,ncols)],mediator_ind,ncols)
		outmat_full<-outmat[,allind]
		
		outfile_full<-str_replace_all(outfile, '.txt', '_full.txt')
		write.table(outmat_full, outfile_full, col.names=TRUE, row.names=FALSE, sep='\t')
		
		outmat=list()
		outmat[[1]]<-outmat_reduced
		outmat[[2]]<-outmat_full
		
	}else{
		write.table(outmat, outfile, col.names=TRUE, row.names=FALSE, sep='\t') 
	}
	
return(outmat)
}	
