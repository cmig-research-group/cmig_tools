library(arrow)
library(furrr)
library(purrr)
library(plyr)
library(dplyr)
library(tidyverse)
library(psych)
library(Matrix)
library(ordinal)
library(pracma)

makeDesign <- function(data, outfile, time, contvar=NULL, catvar=NULL,	delta=NULL, interact=NULL, subjs=NULL, demean=FALSE, quadratic=NULL, mediator=NULL, familyID='ab_g_stc__design_id__fam') {

	#data = data frame with variables of interest or path to tabulated data or path to directory with complete tabulated data needed for design matrix
	#outfile = filepath and name to save design matrix to
	#time = 'eventname' the events that you want to include e.g. c('baseline_year_1_arm_1','2_year_follow_up_y_arm_1')
	#contvar = list of continuous variables e.g. c('interview_age','nihtbx_pattern_uncorrected')
	#catvar = list of categorical variables e.g. c('sex_at_birth','hisp','household.income','high.educ')
	#delta = list of variables to be divided into baseline and change scores e.g. c('bmi') - DO NOT ALSO INCLUDE IN OTHER LISTS
	#interact = list of pairwise interactions e.g. c('interview_age_delta*sex_at_birth','interview_age_base*sex_at_birth','interview_age_base*interview_age_delta')
	#subjs = path to text file with list of subjects to use if want to subsample NO HEADER
	#demean = if TRUE will demean all continuous variables in design matrix inc deltas
	#quadratic = list of variables that you want quadratic predictors for --> STILL BETA TESTING

	# FOR MEDIATION ANALYSIS:	 mediator = name of mediator
	#Include all variables as if creating for the full model.	Mediator MUST already be included in contvar, delta or interact
	#This function will make two design matrices:
	#	 '*_full.txt' = will include mediator as penultimate column (before intercept)
	#	 '*_red.txt' = will be nested model with same sample NOT including mediator
	#Currently only supported for continuous variables or interactions for mediation of moderation

	################################
	# The following R packages need to be loaded
	for (p in c('arrow','furrr','purrr','plyr','dplyr','tidyverse','psych','Matrix','ordinal','pracma')) {
	        if(!eval(parse(text=paste("require(",p,")")))) {
	                install.packages(p)
	                lapply(p,library,character.only=TRUE)
	        }
	}
	################################

	# check what variables have supplied 
	if (is.null(contvar) & is.null(catvar) & is.null(delta)) & is.null(allvar) {
		stop('ERROR! No variables supplied')
	}

	if (is.null(allvar)) { # TO DO: INLCLUDE DATA DICTIONARY TO CHECK WHETHER VARIABLES ARE CONTINUOUS OR CATEGORICAL 
		allvar <- c(contvar, catvar, delta)
	}

	# Define allowed values for damily ID
	valid_familyIDs <- c('ab_g_stc__design_id__fam', 'ab_g_stc__design_id__fam__gen', 'rel_family_id')
	if (!familyID %in% valid_familyIDs) {
		stop("Error: familyID must be either 'ab_g_stc__design_id__fam', 'ab_g_stc__design_id__fam__gen' or 'rel_family_id'.")
	  }

	# define the first four columns of the design matrix
	reqvar <- c(familyID,'ab_g_dyn__visit_age')

	# Function to join datasets while handling missing keys
	merge_dataframe <- function(df1, df2) {
		common_keys <- intersect(c('participant_id', 'session_id'), colnames(df1))
		if ('participant_id' %in% colnames(df1) & 'session_id' %in% colnames(df1) & 'participant_id' %in% colnames(df2) & 'session_id' %in% colnames(df2)) {
			full_join(df1, df2, by = c('participant_id', 'session_id'))
		} else {
			full_join(df1, df2, by = 'participant_id')
		}
	}

	if is.character(data) {
		print(paste0('Reading tabulated data from ', nda)) 
		# detect the file format used 
		data_files <- list.files(nda)
		ext <- unique(tools::file_ext(data_files))
		if 'parquet' %in% ext {
			# get only parquet files in case there are other file types in the directory
			parquet_files <- data_files[grep('parquet', data_files)]
			parquet_files <- paste0(nda, '/', parquet_files)
			# Function to extract variables from parquet files
			process_parquet_file <- function(file, extractvar) {
  				parquet_dataset <- open_dataset(file, format = "parquet")
				# Identify available columns in the file
  				parquetvar <- names(parquet_dataset)

				# Find which columns from `allvar` are present
				present_extractvar <- intersect(extractvar, parquetvar)
				# If no relevant columns, return NULL to avoid processing
  				if (length(present_extractvar) == 0) {
  					return(NULL)  # Or return an empty data frame: return(data.frame())
  				}
				# If at least one column from `allvar` is present, also extract `idvar`
				id_columns <- intersect(c('participant_id', 'session_id'), parquetvar)  # Only keep existing IDs
 				selected_columns <- unique(c(id_columns, present_extractvar))
  				# Select columns if available
  				extracted_data <- parquet_dataset %>% select(any_of(selected_columns))
  				# Convert to a data frame
  				return(as.data.frame(extracted_data))
			}
			# Number of parallel processes
			plan(multisession, workers = 4) 
			data_list <- furrr::future_map(parquet_files, process_parquet_file, allvar)
			data_list <- data_list[!sapply(data_list, is.null)]
			# Merge all data frames iteratively
			combined_df <- reduce(data_list, merge_dataframe)
		}

		# use csv and tsv files if present but not parquet
		if (('csv' %in% ext | 'tsv' %in% ext ) & !'parquet' %in% ext ) {
			# get only tsv files in case there are other file types in the directory
			tsv_files <- data_files[grep("\\.tsv$", data_files)]
			if (length(tsv_files) > 0) {  # Change condition to check if files exist
				tsv_files <- paste0(nda, '/', tsv_files)
			} else {
				tsv_files <- NULL
			}
			csv_files <- data_files[grep("\\.csv$", data_files)]
			if (length(csv_files) > 0) {  # Change condition to check if files exist
				csv_files <- paste0(nda, '/', csv_files)
			} else {
				csv_files <- NULL
			}
			all_files <- c(tsv_files, csv_files)
			# Function to read and select columns
			process_delim_file <- function(file, extractvar) {
				# Read the file
				if (grepl(".csv$", file)) {
				  df <- read.csv(file)
				} else {
				  df <- read.delim(file)
				}
				# Keep only relevant columns that exist in the file
 				# Identify available columns in the file
				available_columns <- colnames(df)
				# Find which columns from `extractvar` are present
				present_extractvar <- intersect(extractvar, available_columns)
				# If no relevant columns, return NULL to avoid processing
				if (length(present_extractvar) == 0) {
				  return(NULL)
				}
 				# If at least one column from `extractvar` is present, also extract `idvar`
				id_columns <- intersect(c('participant_id', 'session_id'), available_columns)  # Only keep existing IDs
				selected_columns <- unique(c(id_columns, present_extractvar))
				# Select columns if available
  				extracted_data <- df[, selected_columns, drop = FALSE]
  				return(extracted_data)
			}
			# Number of parallel processes
			plan(multisession, workers = 4) 
			# Process TSV and CSV files in parallel
			data_list <- furrr::future_map(all_files, process_delim_file, extractvar = allvar)
			data_list <- data_list[!sapply(data_list, is.null)]
			combined_df <- reduce(data_list, merge_dataframe)
		}
	nda <- combined_df
	}
	

	# create path to outfile if doesn't exist
	outpath <- dirname(outfile)
	if ( ! dir.exists(outpath) ) {
		dir.create(outpath, recursive=TRUE)
	}
	
	# if filtering by subjid 
	if ( !is.null(subjs)) {
		s_mat <- read.delim(subjs, header=FALSE)
		s <- data.frame(s_mat)
		names(s) <- 'idevent's
		inc_idevent <- lapply(FUN=grep, s$idevent, nda$idevent)
		inc_idevent <- unlist(inc_idevent)
		nda <- nda[inc_idevent,]
	}
	
	# SHOUDL YOU BE ALLOWED TO MIX AND MATH CSV/TDV/PARQUET FILES?
	# NEED TO ADD HOW TO ASSING ALLVAR TO CATEFORICAL OR CONTINUOUS
	
	# add particupant_id and session_id to reqvar 
	if ("participant_id" %in% names(nda) & "session_id" %in% names(nda)) {
		reqvar <- c('participant_id','session_id',reqvar)
	} else if ("src_subject_id" %in% names(nda) & "eventname" %in% names(nda)) {
		reqvar <- c('src_subject_id','eventname',reqvar)
	} else {
		stop('ERROR! No participant_id and session_id or src_subject_id and eventname in data')
	}

	if ("src_subject_id" %in% names(nda)) { # I DONT KNOW WHETHER TO KEEP THIS OR NOT 
		nda[,'age']<-nda$interview_age
		nda = nda[,c("src_subject_id","eventname","rel_family_id","age",allvar)]
		nda = nda[complete.cases(nda),]
		idx_time <-grep(paste0(time, collapse='|'), nda$eventname)
		# get subject ids at that time point
		subjid<-as.character(nda$src_subject_id[idx_time])
		#src_subject_id <-gsub('NDAR_','',subjid)
		src_subject_id <- subjid
		eventname <- nda$eventname[idx_time]
		rel_family_id<-nda$rel_family_id[idx_time]
		nda <- nda[idx_time,]
	} else {
		# make sure column order is 'participant_id','session_id','ab_g_stc__design_id__fam', 'ab_g_dyn__visit_age'
		# check that all required variables are present
		if ( !all(reqvar %in% colnames(nda))) {
			stop('ERROR! Not all required variables present')
		}
		# check order re-arrange if necessary
		if ( !identical(colnames(nda), c('participant_id','session_id','ab_g_stc__design_id__fam','ab_g_dyn__visit_age',allvar))) {
			nda <- nda[,c('participant_id','session_id','ab_g_stc__design_id__fam','ab_g_dyn__visit_age',allvar)]
		}
		nda <- nda[complete.cases(nda),]
		idx_time <-grep(paste0(time, collapse='|'), nda$session_id)
		# get subject ids at that time point
		subjid<-as.character(nda$participant_id[idx_time])
		#src_subject_id <-gsub('NDAR_','',subjid)
		participant_id <- subjid
		session_id <- nda$session_id[idx_time]
		ab_g_stc__design_id_fam<-nda$ab_g_stc__design_id_fam[idx_time]
		nda <- nda[idx_time,]
	}
	
	# calculate deltas
	if ( !is.null(delta)) {
		uniq_subj <- unique(subjid)
		for (d in delta) {
			vec0 <- rep(NA,dim(nda)[1]);
			vecD <- rep(NA,dim(nda)[1]);
			for (i in 1:length(uniq_subj)) {
				if ("eventname" %in% names(nda)) {
					idx_subj <- which(nda$src_subject_id==uniq_subj[i])
					idx_tmp <- which(nda$interview_age[idx_subj]==min(nda$interview_age[idx_subj]))	
				} else {
					idx_subj <- which(nda$participant_id==uniq_subj[i])
					idx_tmp <- which(nda$ab_g_dyn__visit_age[idx_subj]==min(nda$ab_g_dyn__visit_age[idx_subj]))
				}
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
	# function to extract variables from rds	

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
							 #d0 <- data.frame(src_subject_id, eventname, d0)
							 #names(d0) <- c('src_subject_id', 'eventname', name_base)
		dD <- lapply(name_delta, varextract, data=nda, index=idx_time, dummy=0)
		if (demean==TRUE){
			dD<-apply(data.frame(dD), 2, function(y) y - mean(y, na.rm=T))
		}
							 #dD <- data.frame(src_subject_id, eventname, dD)
							 #names(dD) <- c('src_subject_id', 'eventname', name_delta)
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
	
	R<-rankMatrix(modelmat)
	
	if (R[1]<ncol(modelmat)){
		modelmat<-drop.coef(modelmat)
	}
	
	colnames(modelmat)[colnames(modelmat) == "(Intercept)"] = "intercept"
	
	modelmat_names = c(which(colnames(modelmat) != "intercept"), which(colnames(modelmat) == "intercept"))
	modelmat = modelmat[,modelmat_names]
	
	####################################

	if ("src_subject_id" %in% names(nda)) {
		outmat<-nda[,c('src_subject_id','eventname','rel_family_id','age')]
	} else {
		outmat<-nda[,c('participant_id','session_id','ab_g_stc__design_id__fam','age')]
	}
	
	#if (fam==TRUE && incage==TRUE){
	#	outmat<-nda[,c('src_subject_id','eventname','rel_family_id','age')]
	#} else if (fam==TRUE && incage!=TRUE){
	#	outmat<-nda[,c('src_subject_id','eventname','rel_family_id')]
	#}else if (fam!=TRUE && incage==TRUE){
	#	outmat<-nda[,c('src_subject_id','eventname','age')]
	#}
	
	if (!is.null(delta)){
		delta_df<-data.frame(nda[,c(name_base, name_delta)])
		outmat<-data.frame(outmat, delta_df)

			modelmat<-modelmat[,-c(which(colnames(modelmat) %in% name_base))]
			modelmat<-modelmat[,-c(which(colnames(modelmat) %in% name_delta))]
	 }
	
	
	outmat<-cbind(outmat,modelmat)
	
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
