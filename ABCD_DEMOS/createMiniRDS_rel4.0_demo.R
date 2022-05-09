################################

# This example script describes how to combine variables from ABCD instruments 
# and save as an RDS file. It can be used as input to the script 
# run_makeDesign_example.R which describes how to create design
# matrices for the FEMA package. Users should save a local copy edit paths 
# appropriately. 

################################
# The following R packages need to be loaded

library(tidyverse)
library(psych)
library(plyr)
library(dplyr)
library(PerformanceAnalytics)
library(pracma)

################################
# This section defines input and output paths for files and functions called. 

# Define the path to the directory which contains the tabulated ABCD data 
inpath <- '/path/to/abcd-sync/4.0/tabulated/released'
# Define the path to the plink file which contains the genetic PCs 
# subject data. Depending on the research question(s), it may be favourable to
# covary for genetic ancestry, and this data is stored seperately from the rest
# of the ABCD data 
pcfile <- '/path/to/plink2.eigenvec'

# Define the full path to the output RDS file 
outpath <- '/path/to/your/fave/output/directory'
fname <- 'nda4.0_offrel.RDS'
outmatfile <- paste0(outpath, '/', fname)

# Define the path to tge cmig_utils/r directory, R needs to be able to 
# parse functions from this directory
funcpath <- '/path/to/cmig_tools_utils/r'
# The functionmakeDEAPdemos.R requires the path to the directory which 
# contains the tabulated ABCD data defined explicitly here
deappath <- inpath

# Define the file names for the instruments from which we need to pull
# variables. 
physfile <- 'abcd_ant01.txt'
img_tabfile1 <- 'abcd_smrip10201.txt'
MRIinfofile <- 'abcd_mri01.txt'
imgincfile <- 'abcd_imgincl01.txt';
# Define the full paths to these files 
physfile <- paste0(inpath, '/', physfile)
img_tabfile1 <- paste0(inpath, '/', img_tabfile1)
imgincfile <- paste0(inpath, '/', imgincfile)
MRIinfofile <- paste0(inpath, '/', MRIinfofile)

################################

# R needs to parse two functions from cmig_utils/r. One is load.txt
# to load the .txt files and remove unnecessary 2nd row, the other is 
# makeDEAPdemos.R which collates SES variables into the format used by
# DESP. 
source(paste0(funcpath, '/', 'loadtxt.R'))
source(paste0(funcpath, '/', 'makeDEAPdemos.R'))

################################
# Load the physical health instrument file  
phys <- loadtxt(physfile)
# Extract the variables of interes 
physvar <- c('src_subject_id', 'eventname', 'interview_age', 'sex')
phys <- phys[,physvar]
# Write to a dataframe 
outmat <- phys

################################
# Load the genetic PCs file (this does not require loadtxt.R!) 
pc_mat <- read.delim(pcfile)
# Get just the first 10 PCs and write to a dataframe  
pc <- data.frame(pc_mat[,c('IID','PC1','PC2','PC3','PC4','PC5', 'PC6', 'PC7', 'PC8', 'PC9', 'PC10')])
names(pc) <- c('src_subject_id','PC1','PC2','PC3','PC4','PC5', 'PC6', 'PC7', 'PC8', 'PC9', 'PC10')
# Combine with the physical health variables. 
outmat <- join(outmat, pc, by='src_subject_id', match = "all")

################################
# Load the MRI info instrument  and extract the device serial number and software version 
# variables which are always needed as covariates when using imaging data  
MRIinfo <- loadtxt(MRIinfofile)
MRIinfo <- MRIinfo[,c('src_subject_id','eventname', grep('mri_info', names(MRIinfo), value=TRUE))]
MRIinfo[,'idevent'] <- paste0(MRIinfo$src_subject_id, '_', MRIinfo$eventname)
MRIinfo <- MRIinfo[duplicated(MRIinfo$idevent)==FALSE,]
MRIinfo[which(MRIinfo$mri_info_deviceserialnumber==""),]<-NA
lvl <- unique(MRIinfo$mri_info_deviceserialnumber)
lvl <- lvl[is.na(lvl)==FALSE]
MRIinfo$mri_info_deviceserialnumber<-factor(MRIinfo$mri_info_deviceserialnumber, levels=lvl)
MRIinfo[which(MRIinfo$mri_info_softwareversion==""),]<-NA
lvl <- unique(MRIinfo$mri_info_softwareversion)
lvl <- lvl[is.na(lvl)==FALSE]
MRIinfo$mri_info_softwareversion<-factor(MRIinfo$mri_info_softwareversion, levels=lvl)
MRIinfo <- select(MRIinfo,  -c('idevent'))
# Combine with the previously extracted variables
outmat <- join(outmat, MRIinfo, by=c('src_subject_id', 'eventname'))

################################
# Create the SES variables as coded by DEAP
deap <- makeDEAPdemos(deappath)
deap <- deap[ , -which(names(deap) %in% c("sex", "interview_date", "interview_age"))]
# Combine with the previously extracted variables
outmat <- join(outmat, deap, by=c('src_subject_id', 'eventname'))

################################
# If you wish to work with the morphological imaging varibles (surface area, 
# cortical thickness, Jacobians etc) it is advisable to include the global 
# measures for these variable as covariates. 
# Load imaging data files from tabulated data 
img_tab1 <- loadtxt(img_tabfile1)
# Extract intracranial volume  mean thickness and surface area  
img_tabvar1 <-c('src_subject_id', 'eventname', 'smri_thick_cdk_mean', 'smri_area_cdk_total', 'smri_vol_scs_intracranialv')
img_tab1 <- img_tab1[, img_tabvar1]
# Combine with the previously extracted variables
outmat <- join(outmat, img_tab1, by=c('src_subject_id', 'eventname'))

################################
# Include the MRI QC include/exclude variable 
imginc <- loadtxt(imgincfile)
# Exctract the include/exclude variable for all imaging modalities 
imgincvar <- c('src_subject_id', 'eventname', grep('include', names(imginc), value=TRUE))
imginc <- imginc[, imgincvar]
# Combine with the previously extracted variables
outmat <- join(outmat, imginc, by=c('src_subject_id', 'eventname'))

################################
# Save the "outmat" as an RDS 

if ( ! dir.exists(outpath) ) {
        dir.create(outpath, recursive=TRUE)
}

saveRDS(outmat, file=outmatfile)


