################################

# This example script describes how to combine variables from ABCD instruments 
# and save as a .txt file. It can be used as input to the script 
# run_makeDesign_example.R which describes how to create design
# matrices for the FEMA package. Users should save a local copy and edit paths 
# appropriately. 

################################
# The following R packages need to be loaded

for (p in c("tidyverse", "psych", "plyr", "dplyr", "pracma", "PerformanceAnalytics", "splines")){
        if(!eval(parse(text=paste("require(",p,")")))) {
                install.packages(p)
                lapply(p,library,character.only=TRUE)
        }
}

################################
# This section defines input and output paths for files and functions called. 

# Define the path to the directory which contains the tabulated ABCD data 
inpath <- '/space/syn65/1/data/abcd-sync/5.0/tabulated/released/core'

# Define the path to the plink file which contains the genetic PCs 
# subject data. Depending on the research question(s), it may be favourable to
# covary for genetic ancestry, and this data is stored seperately from the rest
# of the ABCD data 
pcfile <- '/space/syn65/1/data/abcd-sync/5.0/genomics/abcd_gen_y_hat.tsv'

# Define the full path to the output RDS file 
outpath <- '/space/syn65/1/data/abcd-sync/5.0/support_files/designmat'
fname <- 'nda5.0_data'
outmatfile <- paste0(outpath, '/', fname)

# Define the path to tge cmig_utils/r directory, R needs to be able to 
# parse functions from this directory
funcpath <- '/space/amdale/1/tmp/ABCD_cache/Packages/FEMA/beta/cmig_tools_utils/r'

# The functionmakeDEAPdemos.R requires the path to the directory which 
# contains the tabulated ABCD data defined explicitly here
deappath <- '/space/syn65/1/data/abcd-sync/5.0/tabulated/released/core/abcd-general'

# Define the file names for the instruments from which we need to pull
# variables. 
img_thk_file <- 'imaging/mri_y_smr_thk_dsk.csv'
img_area_file <- 'imaging/mri_y_smr_area_dsk.csv'
img_vol_file <- 'imaging/mri_y_smr_vol_aseg.csv'
MRIinfofile <- 'imaging/mri_y_adm_info.csv'
imgincfile <- 'imaging/mri_y_qc_incl.csv'
motionfile <- 'imaging/mri_y_qc_motion.csv'
pdspfile <- 'physical-health/ph_p_pds.csv'
pdsyfile <- 'physical-health/ph_y_pds.csv'
physfile <- 'physical-health/ph_y_anthro.csv' 

# Define the full paths to these files 
img_thk_file <- paste0(inpath, '/', img_thk_file)
img_area_file <- paste0(inpath, '/', img_area_file)
img_vol_file <- paste0(inpath, '/', img_vol_file)
imgincfile <- paste0(inpath, '/', imgincfile)
MRIinfofile <- paste0(inpath, '/', MRIinfofile)
motionfile <- paste0(inpath, '/', motionfile)
pdspfile <- paste0(inpath, '/', pdspfile)
pdsyfile <- paste0(inpath, '/', pdsyfile) 
physfile <- paste0(inpath, '/', physfile)  

################################
# R needs to parse two functions: makeDEAPdemos.R which collates SES variables 
# into the format used by DEAP, and createBasis.R which is used to generate 
# the basis matrix used for modeling smooth functions
source(paste0(funcpath, '/', 'makeDEAPdemos.R'))
source(paste0(funcpath, '/', 'createBasis.R'))

################################
# Load the MRI info instrument  and extract the device serial number and software version 
# variables which are always needed as covariates when using imaging data  
MRIinfo <- read.csv(MRIinfofile)
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
outmat <- MRIinfo

################################
# Load the genetic PCs file 
pc_mat <- read.delim(pcfile)
# Get just the first 10 PCs and write to a dataframe
pc_names <- paste0('genetic_pc_',c(1:10))
pc <- data.frame(pc_mat[,c('src_subject_id',pc_names)])
names(pc) <- c('src_subject_id',paste0('PC',c(1:10)))
# Combine with the physical health variables. 
outmat <- join(outmat, pc, by='src_subject_id', match = "all")

################################
## Create the SES variables as coded by DEAP
deap <- makeDEAPdemos(deappath)
deap <- deap[ , -which(names(deap) %in% c("interview_date"))]
# Combine with the previously extracted variables
outmat <- join(outmat, deap, by=c('src_subject_id', 'eventname'))

################################
# If you wish to work with the morphological imaging varibles (surface area, 
# cortical thickness, etc) it is advisable to include the global 
# measures for these variable as covariates. 

# Load imaging data files from tabulated data 
img_thk <- read.csv(img_thk_file)
# Extract intracranial volume  mean thickness and surface area  
img_thk_vars <-c('src_subject_id', 'eventname', 'smri_thick_cdk_mean')
img_thk <- img_thk[, img_thk_vars]
# Combine with the previously extracted variables
outmat <- join(outmat, img_thk, by=c('src_subject_id', 'eventname'))

img_area <- read.csv(img_area_file)
# Extract total surface area 
img_area_vars <-c('src_subject_id', 'eventname', 'smri_area_cdk_total')
img_area <- img_area[, img_area_vars]
# Combine with the previously extracted variables
outmat <- join(outmat, img_area, by=c('src_subject_id', 'eventname'))

img_vol <- read.csv(img_vol_file)
# Extract total surface area 
img_vol_vars <-c('src_subject_id', 'eventname', 'smri_vol_scs_intracranialv')
img_vol <- img_vol[, img_vol_vars]
# Combine with the previously extracted variables
outmat <- join(outmat, img_vol, by=c('src_subject_id', 'eventname'))

################################
# Include the MRI QC include/exclude variable 
imginc <- read.csv(imgincfile)
# Exctract the include/exclude variable for all imaging modalities 
imgincvar <- c('src_subject_id', 'eventname', grep('include', names(imginc), value=TRUE))
imginc <- imginc[, imgincvar]

# Combine with the previously extracted variables
outmat <- join(outmat, imginc, by=c('src_subject_id', 'eventname'))

################################
# Include the MRI QC motion variable
motion <- read.csv(motionfile)
# Extract intracranial volume  mean thickness and surface area  
motion_vars <-c('src_subject_id', 'eventname', 'dmri_meanmotion', 'rsfmri_meanmotion')
motion <- motion[, motion_vars]
# Combine with the previously extracted variables
outmat <- join(outmat, motion, by=c('src_subject_id', 'eventname'))

################################
# Need to format design matrix so that first three columns are 
# src_subject_id, eventname, rel_family_id
colnames = c('src_subject_id', 'eventname', 'rel_family_id')
outmat = cbind(outmat[,colnames], outmat[,-which(names(outmat) %in% colnames)])

################################
# Pubertal Development, PDS
# average parent and youth report or use whichever report is available if only one informant
pds_y<-read.csv(pdsyfile)
pds_p<-read.csv(pdspfile)

#Tanner stage categories
# MALES Prepubertal = 3; early Pubertal = 4 or 5 (no 3-point responses); Midpubertal = 6,7, or 8 (no 4-point responses; Late pubertal = 9-11; Postpubertal = 12. 
# FEMALES Prepubertal = 3; Early Puberty = 3 and no menarche; Midpubertal = 4 and no menarche; Late Puberty = <=7 and menarche; Postpubertal = 8 and menarche

#youth
pds_y[,'pds_y_ss_category_all']<-coalesce(pds_y$pds_y_ss_female_category_2, pds_y$pds_y_ss_male_cat_2)
pds_y$pds_y_ss_category_all<-as.numeric(pds_y$pds_y_ss_category_all)

#parent
pds_p[,'pds_p_ss_category_all']<-coalesce(pds_p$pds_p_ss_female_category_2, pds_p$pds_p_ss_male_category_2)
pds_p$pds_p_ss_category_all<-as.numeric(pds_p$pds_p_ss_category_all)

pds<-join(pds_p,pds_y,c("src_subject_id","eventname"))

#take average from parent and youth reports; if one is missing, take the non-missing value
pds$pds_y_p_average <- ifelse(is.na(pds$pds_y_ss_category_all), pds$pds_p_ss_category_all, ifelse(is.na(pds$pds_p_ss_category_all), pds$pds_y_ss_category_all, (pds$pds_y_ss_category_all + pds$pds_p_ss_category_all) / 2))

pds_vars = c('src_subject_id', 'eventname', 'pds_y_p_average', 'pds_y_ss_category_all', 'pds_p_ss_category_all')
pds <- pds[, pds_vars] 

outmat <- join(outmat, pds, by=c('src_subject_id', 'eventname'))

################################
# Body Mass Index (BMI)
phys <- read.csv(physfile)
# calculate bmi and tmi
weightkg <- phys$anthroweightcalc*0.453592
heightm <- phys$anthroheightcalc*0.0254
bmi <- weightkg/(heightm^2)
tmi <- weightkg/(heightm^3)
phys$anthro_bmi_calc <- bmi
phys$anthro_tmi_calc <- tmi
# remove biologically implausible values
ulim <- 45
llim <- 11
rm_bmi <- which(phys$anthro_bmi_calc>ulim | phys$anthro_bmi_calc<llim)
anthro_bmi_corr <- phys$anthro_bmi_calc
anthro_tmi_corr <- phys$anthro_tmi_calc
phys$anthro_bmi_corr <- anthro_bmi_corr
phys$anthro_tmi_corr <- anthro_tmi_corr
phys[rm_bmi,'anthro_bmi_corr'] <- NA
phys[rm_bmi,'anthro_tmi_corr'] <- NA

phys_vars = c('src_subject_id', 'eventname', 'anthro_bmi_corr')
phys <- phys[, phys_vars] 

outmat <- join(outmat, phys, by=c('src_subject_id', 'eventname'))

################################
# compute basis functions for age and save corresponding file
agevec = seq(from=100,to=200,length=101) # create age vector (in months)

# save basis matrix
basis = createBasis(agevec, intercept = TRUE, demean = TRUE)
write.table(basis, file = paste0(outpath, '/basis_age.txt'), sep = "\t", row.names = FALSE)

# Apply function to each row of interview_age in outmat
basis_values_df = get_basis_values(outmat,basis,'interview_age')
colnames(basis_values_df) <- colnames(basis)

# add basis functions to outmat
outmat <- cbind(outmat, basis_values_df)

################################
# compute basis functions for puberty and save corresponding file
pdsvec = seq(from=0.5,to=5.5,length=11) # create PDS vector

# save basis matrix
basis_pds = createBasis(pdsvec, intercept = TRUE, demean = TRUE)
colnames(basis_pds) <- gsub("demean", "pdsavg", names(basis_pds))
write.table(basis_pds, file = paste0(outpath, '/basis_pds.txt'), sep = "\t", row.names = FALSE)

# Apply function to each row of pds in outmat
basis_pds_df = get_basis_values(outmat,basis_pds,'pds_y_p_average')
colnames(basis_pds_df) <- colnames(basis_pds)

# add basis functions to outmat
outmat <- cbind(outmat, basis_pds_df)

################################
# compute basis functions for BMI and save corresponding file
bmivec = seq(from=11,to=45,length=35) # create PDS vector

# save basis matrix
basis_bmi = createBasis(bmivec, intercept = TRUE, demean = TRUE)
colnames(basis_bmi) <- gsub("demean", "bmi", names(basis_bmi))
write.table(basis_bmi, file = paste0(outpath, '/basis_bmi.txt'), sep = "\t", row.names = FALSE)

# Apply function to each row of bmi in outmat
basis_bmi_df = get_basis_values(outmat,basis_bmi,'anthro_bmi_corr')
colnames(basis_bmi_df) <- colnames(basis_bmi)

# add basis functions to outmat
outmat <- cbind(outmat, basis_bmi_df)

################################
# Save the "outmat"

if ( ! dir.exists(outpath) ) {
        dir.create(outpath, recursive=TRUE)
}

write.table(outmat, file=paste0(outmatfile, '.txt'), sep = "\t", row.names = FALSE)
cat(paste0('File written to ', outmatfile, '.txt\n'))

