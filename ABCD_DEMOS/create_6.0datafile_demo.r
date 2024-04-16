################################
# WILL UPDATE AS MORE 6.0 DATA BECOMES AVAILABLE
# This example script describes how to combine variables from ABCD instruments for 6.0 data 
# and save as a .txt file. It can be used as input to the script 
# run_makeDesign_example.R which describes how to create design
# matrices for the FEMA package. Users should save a local copy and edit paths 
# appropriately. 
# ABCD 6.0 DATA: merging sociodemographic and behavioral (cognition both in and out of scanner, psychopathology) vars, physical vars, and other covariates of interest
# from here, select variables you want to create design matrix for FEMA
# Ali Rigby updating and combining 4.0 script from Dr. Carolina Makowski and script from Diana Smith
# April 2024

################################
# The following R packages need to be loaded
for (p in c("tidyverse", "psych", "plyr", "dplyr", "pracma", "argparse", "tableone", "data.table")){
        if(!eval(parse(text=paste("require(",p,")")))) {
                install.packages(p)
                lapply(p,library,character.only=TRUE)
        }
}

################################
# This section defines input and output paths for files and functions called. 

# Define the path to the directory which contains the tabulated ABCD data 
data5path <- '/space/syn65/1/data/abcd-sync/5.0/tabulated/released/core' 
data6path <- '/space/syn65/1/data/abcd-sync/6.0/tabulated/img' 

# Define pre-release data path
deap_6.0_file <- '/space/syn65/1/data/abcd-sync/6.0/data_deap_tabulated.csv'

# Define the path to the genetic PCs 
pcs_genesis_file <-'/space/syn65/1/data/abcd-sync/5.0/genomics/abcd_gen_y_hat.tsv'

# Define the full path to the output RDS file 
outpath<-'/space/syn65/1/data/abcd-sync/6.0/support_files' 
fname <- 'nda6.0_data'
outmatfile <- paste0(outpath, '/', fname)

# Define the path to tge cmig_utils/r directory, R needs to be able to 
# parse functions from this directory
funcpath <- '/home/arigby/cmig_github/cmig_tools_internal/cmig_tools_utils/r'

# The functionmakeDEAPdemos.R requires the path to the directory which 
# contains the tabulated ABCD data defined explicitly here
# deappath <- '/space/syn65/1/data/abcd-sync/5.0/tabulated/released/core/abcd-general'
# deappath <- data5path

# Define the file names for the instruments from which we need to pull
# variables. As more pre-release 6.0 becomes available, this file will be updated with the missing variables.
mri_info_file <- 'abcd_mri01.csv'
imgincfile <- 'abcd_imgincl01.csv'
motion_file_rsfmri <- 'abcd_betnet02.csv'
# motion_file_dmri <- 'abcd_dmdtifp202.csv'
img_thk_file <- 'abcd_smrip102.csv'
img_area_file <- 'abcd_smrip102.csv'
img_vol_file <- 'abcd_smrip102.csv'
# nih_tbx_file <- 'neurocognition/nc_y_nihtb.csv'
# ravlt_file <- 'neurocognition/nc_y_ravlt.csv'
# matrx_file <- 'neurocognition/nc_y_wisc.csv'
# nback_beh_file <- 'imaging/mri_y_tfmr_nback_beh.csv'
# dprime_file <- 'imaging/mri_y_tfmr_nback_rec_beh.csv'
# mid_beh_file <- 'imaging/mri_y_tfmr_mid_beh.csv'
# ssrt_file <- 'imaging/mri_y_tfmr_sst_beh.csv'
# cbcl_file <- 'mental-health/mh_p_cbcl.csv'
# upps_file <- 'mental-health/mh_y_upps.csv'
# bisbas_file <- 'mental-health/mh_y_bisbas.csv'


# Define the full paths to these files 
mri_info_file <- paste0(data6path, '/', mri_info_file)
imgincfile <- paste0(data6path, '/', imgincfile)
motion_file_rsfmri <- paste0(data6path, '/', motion_file_rsfmri)
# motion_file_dmri <- paste0(data6path, '/', motion_file_dmri)
img_thk_file <- paste0(data6path, '/', img_thk_file)
img_area_file <- paste0(data6path, '/', img_area_file)
img_vol_file <- paste0(data6path, '/', img_vol_file)
# nih_tbx_file <- paste0(data5path, '/', nih_tbx_file)
# ravlt_file <- paste0(data5path, '/', ravlt_file)
# matrx_file <- paste0(data5path, '/', matrx_file)
# nback_beh_file <- paste0(data5path, '/', nback_beh_file)
# dprime_file <- paste0(data5path, '/', dprime_file)
# mid_beh_file <- paste0(data5path, '/', mid_beh_file)
# ssrt_file <- paste0(data5path, '/', ssrt_file)
# cbcl_file <- paste0(data5path, '/', cbcl_file)
# upps_file <- paste0(data5path, '/', upps_file)
# bisbas_file <- paste0(data5path, '/', bisbas_file)

################################
# R needs to parse two functions from cmig_utils/r. One is load.txt
# to load the .txt files and remove unnecessary 2nd row, the other is 
# makeDEAPdemos.R which collates SES variables into the format used by
# DEAP. 
source(paste0(funcpath, '/', 'makeDesign.R'))
# source(paste0(funcpath, '/', 'makeDEAPdemos.R'))
# source(paste0(funcpath, '/', 'loadtxt.R'))

################################
# Create the SES variables as coded by DEAP - cannot use until we have complete 6.0 data 
# Will also need to update makeDEAPdemos function for 6.0
# deapdemos<-makeDEAPdemos(deappath)

################################
# Add imaging/MRI info data
# Load the MRI info instrument and extract the device serial number, manufacturer, and software version 
# variables which are always needed as covariates when using imaging data  
MRIinfo<-read.csv(mri_info_file)
MRIinfo<-MRIinfo[,c('src_subject_id','eventname','mri_info_deviceserialnumber','mri_info_manufacturer','mri_info_softwareversion')]
MRIinfo[,'idevent']<-paste0(MRIinfo$src_subject_id, '_', MRIinfo$eventname)
MRIinfo<-MRIinfo[duplicated(MRIinfo$idevent)==FALSE,]
MRIinfo[which(MRIinfo$mri_info_deviceserialnumber==""),]<-NA
lvl<-unique(MRIinfo$mri_info_deviceserialnumber)
lvl<-lvl[is.na(lvl)==FALSE]
MRIinfo$mri_info_deviceserialnumber<-factor(MRIinfo$mri_info_deviceserialnumber, levels=lvl)
MRIinfo[which(MRIinfo$mri_info_softwareversion==""),]<-NA
lvl<-unique(MRIinfo$mri_info_softwareversion)
lvl<-lvl[is.na(lvl)==FALSE]
MRIinfo$mri_info_softwareversion<-factor(MRIinfo$mri_info_softwareversion, levels=lvl)
MRIinfo <- select(MRIinfo,  -c('idevent'))
# Combine with the previously extracted variables
df <- MRIinfo

################################
# Read PC file, generated by genesis toolbox
pc_mat <- read.delim(pcs_genesis_file)
# Get just the first 20 PCs and write to a dataframe
pc_names <- paste0('genetic_pc_',c(1:10))
pc <- data.frame(pc_mat[,c('src_subject_id',pc_names)])
names(pc) <- c('src_subject_id',paste0('PC',c(1:10)))

# Joining DEAPdemos with PCs in a new dataframe
df = join(df, pc, by='src_subject_id', match = "all")

################################
# Include the MRI QC include/exclude variable 
imginc <- read.csv(imgincfile)

# Reformat "VisitID" variable to match src_subject_id and eventname
visitid = data.frame(str_split_fixed(imginc[,'VisitID'],"_",3))
visitid$src_subject_id = paste0('NDAR_',visitid[,'X2'])
visitid$eventname = case_match(visitid[,'X3'], 
'baseline' ~ 
'baseline_year_1_arm_1', 
'2year' ~ '2_year_follow_up_y_arm_1',
'4year' ~ '4_year_follow_up_y_arm_1',
'6year' ~ '6_year_follow_up_y_arm_1')

imginc = cbind(visitid[,c('src_subject_id', 'eventname')],imginc)

# Exctract the include/exclude variable for all imaging modalities 
imgincvar <- c('src_subject_id', 'eventname', grep('include', names(imginc), value=TRUE))
imginc <- imginc[, imgincvar]

# Combine with the previously extracted variables
df <- join(df, imginc, by=c('src_subject_id', 'eventname'))

# ################################
# # Include the MRI QC include/exclude variable 
# imgincl<-read.delim((imgincfile), header=TRUE, sep = ',')

# # Exctract the include/exclude variable for all imaging modalities 
# imgincl_vars<-c('imgincl_t1w_include','imgincl_t2w_include','imgincl_dmri_include','imgincl_mid_include','imgincl_nback_include','imgincl_sst_include')
# imgincl<-imgincl[,c('src_subject_id','eventname',imgincl_vars)]

# # Add MRI QC include/exclude variables to main dataframe
# df<-join(df,imgincl)

################################
# Load imaging data files from tabulated data 
# Add global measures (total surface area, mean cortical thickness, total intracranial volume for respective analyses)
img_thk <- read.csv(img_thk_file)
# Extract intracranial volume  mean thickness and surface area  
img_thk_vars <-c('src_subject_id', 'eventname', 'smri_thick_cdk_mean')
img_thk <- img_thk[, img_thk_vars]
# Combine with the previously extracted variables
df <- join(df, img_thk, by=c('src_subject_id', 'eventname'))

img_area <- read.csv(img_area_file)
# Extract total surface area 
img_area_vars <-c('src_subject_id', 'eventname', 'smri_area_cdk_total')
img_area <- img_area[, img_area_vars]
# Combine with the previously extracted variables
df <- join(df, img_area, by=c('src_subject_id', 'eventname'))

img_vol <- read.csv(img_vol_file)
# Extract total surface area 
img_vol_vars <-c('src_subject_id', 'eventname', 'smri_vol_scs_intracranialv')
img_vol <- img_vol[, img_vol_vars]
# Combine with the previously extracted variables
df <- join(df, img_vol, by=c('src_subject_id', 'eventname'))

################################
# Include the MRI QC motion variable
motion_rsfmri <- read.csv(motion_file_rsfmri)
# motion_dmri <- read.csv(motion_file_dmri)
# Extract intracranial volume  mean thickness and surface area  
motion_vars_rsfmri <-c('src_subject_id', 'eventname', 'rsfmri_c_ngd_meanmotion')
# motion_vars_dmri <-c('src_subject_id', 'eventname', 'dmri_meanmotion')
motion_rsfmri <- motion_rsfmri[, motion_vars_rsfmri]
# motion_dmri <- motion_dmri[, motion_vars_dmri]
# Combine with the previously extracted variables
df <- join(df, motion_rsfmri, by=c('src_subject_id', 'eventname'))
# df <- join(df, motion_dmri, by=c('src_subject_id', 'eventname'))

################################
# We don't yet have family ID for 6.0 data, so using a file 
# from 5.0 (since family ID doesn't change) 
# Will update once we have complete 6.0 data
famfile = '/space/syn65/1/data/abcd-sync/5.0/support_files/birth_id.txt'
fam = read.delim(famfile, sep = ' ')

# Get just the first 10 PCs and write to a dataframe
famnames <- c('pguid', 'update_family_id')
fam = fam[,famnames]
colnames(fam) = c('src_subject_id', 'rel_family_id')

# Combine with the physical health variables. 
df <- join(df, fam, by='src_subject_id', match = "all")

# Need to format design matrix so that first three columns are 
# src_subject_id, eventname, rel_family_id
colnames = c('src_subject_id', 'eventname', 'rel_family_id')
df = cbind(df[,colnames], df[,-which(names(df) %in% colnames)])

################################
# add variables from DEAP pre-release file
deap = read.csv(deap_6.0_file)

# rename subject id, eventname, age
deap$src_subject_id = deap$id_redcap
deap$eventname = deap$redcap_event_name
deap$interview_age = deap$age_visit

# sex
sextmp = data.frame(deap[deap$eventname=='baseline_year_1_arm_1',c('src_subject_id','demo_sex_v2b')])
sextmp$sex = recode(as.factor(sextmp$demo_sex_v2b), "1" = "M","2" = "F", "3" = "I")
deap<-join(deap,sextmp[,c('src_subject_id','sex')], by='src_subject_id', match = "all")

# parental education
deap[,'demo_prnt_ed_p']<-coalesce(deap$demo_prnt_ed_v2b,deap$demo_prnt_ed_v2_l)
deap[,'demo_prnt_ed_p']<-coalesce(deap$demo_prnt_ed_p,deap$demo_prnt_ed_v2_2yr_l)

# commented out because we don't have partner education for 6.0 yet
# deap$demo_prtnr_ed_v2_l = as.integer(deap$demo_prtnr_ed_v2_l)
# deap[,'demo_prtnr_ed_p']<-coalesce(deap$demo_prtnr_ed_v2,deap$demo_prtnr_ed_v2_l)
# deap[,'demo_prtnr_ed_p']<-coalesce(deap$demo_prtnr_ed_p,deap$demo_prtnr_ed_v2_2yr_l)

#highest education: 5 different levels. These levels correspond to the numbers published by the American Community Survey (ACS).
high.educ1 = deap$demo_prnt_ed_p
# high.educ2 = alldems$demo_prtnr_ed_p
high.educ1[which(high.educ1 == "999")] = NA
# high.educ2[which(high.educ2 == "999")] = NA
high.educ1[which(high.educ1 == "777")] = NA
# high.educ2[which(high.educ2 == "777")] = NA
high.educ1[which(high.educ1 == "22" | high.educ1=="23")] = 15 #22 and 23 = some college --> lower level than 18+
# high.educ2[which(high.educ2 == "22" | high.educ2=="23")] = 15
high.educ = pmax(as.numeric(as.character(high.educ1)), na.rm=T) # high.educ = pmax(as.numeric(as.character(high.educ1)), as.numeric(as.character(high.educ2)), na.rm=T)
idx <- which(high.educ %in% 0:12, arr.ind = TRUE)
high.educ[idx] = 1 # "< HS Diploma"
idx <- which(high.educ %in% 13:14, arr.ind = TRUE)
high.educ[idx] = 2 # "HS Diploma/GED"
idx <- which(high.educ %in% c(15:17,22:23), arr.ind = TRUE)
high.educ[idx] = 3 # "Some College"
idx <- which(high.educ == 18, arr.ind = TRUE)
high.educ[idx] = 4 # "Bachelor"
idx <- which(high.educ %in% 19:21, arr.ind = TRUE)
high.educ[idx] = 5 # "Post Graduate Degree"
high.educ[which(high.educ == "999")]=NA
high.educ[which(high.educ == "777")]=NA
deap$high.educ = factor( high.educ, levels= 1:5, labels = c("< HS Diploma","HS Diploma/GED","Some College","Bachelor","Post Graduate Degree") )

# household income
deap[,'demo_comb_income_p']<-coalesce(deap$demo_comb_income_v2b,deap$demo_comb_income_v2_l)

household.income = deap$demo_comb_income_p
household.income[deap$demo_comb_income_p == "1"] = 1 # "[<50K]"
household.income[deap$demo_comb_income_p == "2"] = 1 # "[<50K]"
household.income[deap$demo_comb_income_p == "3"] = 1 # "[<50K]"
household.income[deap$demo_comb_income_p == "4"] = 1 # "[<50K]"
household.income[deap$demo_comb_income_p == "5"] = 1 # "[<50K]"
household.income[deap$demo_comb_income_p == "6"] = 1 # "[<50K]"
household.income[deap$demo_comb_income_p == "7"] = 2 # "[>=50K & <100K]"
household.income[deap$demo_comb_income_p == "8"] = 2 # "[>=50K & <100K]"
household.income[deap$demo_comb_income_p == "9"] = 3 # "[>=100K]"
household.income[deap$demo_comb_income_p == "10"] = 3 # "[>=100K]"
household.income[deap$demo_comb_income_p == "777"] = NA
household.income[deap$demo_comb_income_p == "999"] = NA
household.income[household.income %in% c(NA, "999", "777")] = NA
deap$household.income = factor( household.income, levels= 1:3, labels = c("[<50K]", "[>=50K & <100K]", "[>=100K]") )

### Household income (continuous) - assign value based on middle of category
household.income_cont = deap$demo_comb_income_p
household.income_cont[deap$demo_comb_income_p == "1"] = 2500 # Less than $5,000
household.income_cont[deap$demo_comb_income_p == "2"] = 8500 # $5,000 through $11,999
household.income_cont[deap$demo_comb_income_p == "3"] = 14000 # $12,000 through $15,999
household.income_cont[deap$demo_comb_income_p == "4"] = 20500 # $16,000 through $24,999
household.income_cont[deap$demo_comb_income_p == "5"] = 30000 # $25,000 through $34,999;
household.income_cont[deap$demo_comb_income_p == "6"] = 42500 # $35,000 through $49,999
household.income_cont[deap$demo_comb_income_p == "7"] = 62500 # $50,000 through $74,999
household.income_cont[deap$demo_comb_income_p == "8"] = 87500 # $75,000 through $99,999
household.income_cont[deap$demo_comb_income_p == "9"] = 150000 # $100,000 through $199,999
household.income_cont[deap$demo_comb_income_p == "10"] = 250000 # $200,000 and greater
household.income_cont[deap$demo_comb_income_p == "777"] = NA # Refuse to answer
household.income_cont[deap$demo_comb_income_p == "999"] = NA # Don't know
household.income_cont[household.income_cont %in% c(NA, "999", "777")] = NA
deap$household.income_cont = household.income_cont

# Household income (10 level)
household.income_10level = deap$demo_comb_income_p
household.income_10level[deap$demo_comb_income_p == "777"] = NA # Refuse to answer
household.income_10level[deap$demo_comb_income_p == "999"] = NA # Don't know
household.income_10level[household.income_10level %in% c(NA, "999", "777")] = NA
deap$household.income_10level = household.income_10level

# calculate bmi and tmi
weightkg <- deap$anthro_weight_calc*0.453592
heightm <- deap$anthro_height_calc*0.0254
bmi <- weightkg/(heightm^2)
tmi <- weightkg/(heightm^3)
deap$anthro_bmi_calc <- bmi
deap$anthro_tmi_calc <- tmi
# remove biologically implausible values
ulim <- 45
llim <- 11
rm_bmi <- which(deap$anthro_bmi_calc>ulim | deap$anthro_bmi_calc<llim)
anthro_bmi_corr <- deap$anthro_bmi_calc
anthro_tmi_corr <- deap$anthro_tmi_calc
deap$anthro_bmi_corr <- anthro_bmi_corr
deap$anthro_tmi_corr <- anthro_tmi_corr
deap[rm_bmi,'anthro_bmi_corr'] <- NA
deap[rm_bmi,'anthro_tmi_corr'] <- NA

################################
# Pubertal Development, PDS
# average parent and youth report or use whichever report is available if only one informant

#Tanner stage categories
# MALES Prepubertal = 3; early Pubertal = 4 or 5 (no 3-point responses); Midpubertal = 6,7, or 8 (no 4-point responses; Late pubertal = 9-11; Postpubertal = 12. 

## FEMALES Prepubertal = 3; Early Puberty = 3 and no menarche; Midpubertal = 4 and no menarche; Late Puberty = <=7 and menarche; Postpubertal = 8 and menarche

#youth
deap[,'pds_y_ss_category_all']<-coalesce(deap$pds_y_ss_female_category_2, deap$pds_y_ss_male_category_2)
deap$pds_y_ss_category_all<-as.numeric(deap$pds_y_ss_category_all)

#parent
deap[,'pds_p_ss_category_all']<-coalesce(deap$pds_p_ss_female_category_2, deap$pds_p_ss_male_category_2)
deap$pds_p_ss_category_all<-as.numeric(deap$pds_p_ss_category_all)

#take average from parent and youth reports; if one is missing, take the non-missing value
deap$pds_y_p_average <- ifelse(is.na(deap$pds_y_ss_category_all), deap$pds_p_ss_category_all, ifelse(is.na(deap$pds_p_ss_category_all), deap$pds_y_ss_category_all, (deap$pds_y_ss_category_all + deap$pds_p_ss_category_all) / 2))

# merge deapvars with outmat
deap_vars = c('src_subject_id', 'eventname', 'sex', 'interview_age', 'high.educ', 'household.income', 
'household.income_cont', 'household.income_10level', 'anthro_bmi_corr', 'pds_y_p_average', 'pds_y_ss_category_all', 'pds_p_ss_category_all')
deap <- deap[, deap_vars] 

df <- join(df, deap, by=c('src_subject_id', 'eventname'))

################################
# Save the df

if ( ! dir.exists(outpath) ) {
        dir.create(outpath, recursive=TRUE)
}

write.table(df, file=paste0(outmatfile, '.txt'), sep = "\t", row.names = FALSE)
cat(paste0('File written to ', outmatfile, '.txt\n'))

# ################################
# # Load NIH toolbox data files from tabulated data 
# nihtbx<-read.csv(nih_tbx_file)

# # Note about longitudinal data: cardsort, list, and composite scores are not part of year2 data [fluid no longer calculated, and total comp would not be equivalent between baseline & 2yr])
# # could keep flanker for now but remote COVID assessments moved to another flanker task where scores aren't super comparable
# # using uncorrected scores because less seems to be missing? for instance pic vocab and reading; ask Frank difference between raw and uncorrected

# # Extract uncorrected scores from tests
# nihtbx_vars <- c('nihtbx_reading_uncorrected','nihtbx_flanker_uncorrected','nihtbx_cardsort_uncorrected','nihtbx_pattern_uncorrected',
#                'nihtbx_picture_uncorrected','nihtbx_picvocab_uncorrected','nihtbx_list_uncorrected','nihtbx_totalcomp_uncorrected',
#                'nihtbx_fluidcomp_uncorrected','nihtbx_cryst_uncorrected')

# nihtbx <- nihtbx[,c('src_subject_id','eventname',nihtbx_vars)]

# # Add variables to main dataframe
# df <- join(df,nihtbx)

# ################################
# # Load matrix reasoning and ravlt data files from tabulated data 
# matrx <- read.csv(matrx_file)
# ravlt <- read.csv(ravlt_file)
# # Extract wisc-v matrix reasoning score
# matrx <- matrx[,c('src_subject_id','eventname', 'pea_wiscv_trs')]
# # Create vector of columns to sum
# ind_pea_ravlt <- c('pea_ravlt_sd_trial_i_tc','pea_ravlt_sd_trial_ii_tc','pea_ravlt_sd_trial_iii_tc','pea_ravlt_sd_trial_iv_tc',
#                   'pea_ravlt_sd_trial_v_tc', 'pea_ravlt_sd_listb_tc','pea_ravlt_ld_trial_vii_tc');

# # Ensure variables are in a numeric format
# matrx[3] = apply(matrx[3],1,as.numeric)
# ravlt[,ind_pea_ravlt] = apply(ravlt[,ind_pea_ravlt],1,as.numeric)
# # Sum columns and put values into a new column
# ravlt$pea_ravlt_sd_trial_sum5trials = apply(ravlt[,ind_pea_ravlt],1,sum)

# # Join vars into new df
# pea <- join(matrx, ravlt)
# # Extract columns of interest
# pea_vars <- c('pea_wiscv_trs','pea_ravlt_sd_trial_i_tc','pea_ravlt_sd_trial_ii_tc','pea_ravlt_sd_trial_iii_tc','pea_ravlt_sd_trial_iv_tc',
#             'pea_ravlt_sd_trial_v_tc','pea_ravlt_sd_trial_sum5trials','pea_ravlt_sd_listb_tc','pea_ravlt_ld_trial_vii_tc')
# pea <- pea[,c('src_subject_id','eventname',pea_vars)]

# # Combine new variables into main dataframe
# df <- join(df, pea)

# ################################
# # Load nback in-scanner behavior data file from tabulated data
# nback_beh <- read.csv(nback_beh_file)

# # variables of interest
# # tfmri_nb_all_beh_c2b_rate - The rate of correct responses to 2 back stimuli during run 1 and run 2
# # tfmri_nb_all_beh_c2b_mrt - Average reaction time for all correct responses to 2 back stimuli during run 1 and run 2
# # tfmri_nb_r1_beh_c2b_rate; tfmri_nb_r2_beh_c2b_rate- rate of corr responses, run 1 or run 2 only
# # tfmri_nb_r1_beh_c2b_mrt; tfmri_nb_r2_beh_c2b_mrt - avg RT for correct responses, run 1 or run 2 only

# # Extract variables of interest
# nback_vars<-c('tfmri_nb_all_beh_c2b_rate','tfmri_nb_all_beh_c2b_mrt','tfmri_nb_r1_beh_c2b_rate',
#               'tfmri_nb_r2_beh_c2b_rate','tfmri_nb_r1_beh_c2b_mrt','tfmri_nb_r2_beh_c2b_mrt')

# nback_beh<-nback_beh[,c('src_subject_id','eventname',nback_vars)]
# # Ensure variables are in a numeric format
# nback_beh[3:8] <- sapply(nback_beh[3:8],as.numeric)
# # Combine variables with main dataframe
# df<-join(df,nback_beh)

# ################################
# # Load dprime behavior data file from tabulated data  ***CAN OMIT
# dprime<-read.csv(dprime_file)
# # Compute the average of scores and create new column
# dprime$tfmri_rec_all_beh_dp_avg<- mean(dprime$tfmri_rec_all_beh_posf_dpr + dprime$tfmri_rec_all_beh_neutf_dp + dprime$tfmri_rec_all_beh_negf_dp + dprime$tfmri_rec_all_beh_place_dp)
# dprime_vars<-c('tfmri_rec_all_beh_posf_dpr','tfmri_rec_all_beh_neutf_dp','tfmri_rec_all_beh_negf_dp','tfmri_rec_all_beh_place_dp','tfmri_rec_all_beh_dp_avg')
# dprime<-dprime[,c('src_subject_id','eventname',dprime_vars)]
# # Combine variables with main dataframe
# df<-join(df,dprime)

# ################################
# # Load mid behavior data file from tabulated data  ***CAN OMIT
# mid_beh<-read.csv(mid_beh_file)

# #variables of interest
# #tfmri_mid_all_beh_lrw_mrt - average reaction time for large reward trials, run 1 + run 2
# #tfmri_mid_r1_beh_lrw_mrt, tfmri_mid_r2_beh_lrw_mrt - avg RT for large reward trials, run1 + run2 sep
# #tfmri_mid_all_beh_ll_mrt - avg RT for large loss trials, run1 + run2
# #tfmri_mid_r1_beh_ll_mrt - avg RT for large loss trials, run1 + run2 sep

# # Extract variables of interest
# mid_vars<-c('tfmri_mid_all_beh_lrw_mrt','tfmri_mid_r1_beh_lrw_mrt','tfmri_mid_r2_beh_lrw_mrt',
#             'tfmri_mid_all_beh_ll_mrt','tfmri_mid_r1_beh_ll_mrt','tfmri_mid_r2_beh_ll_mrt','tfmri_mid_all_beh_t_earnings')
# mid_beh<-mid_beh[,c('src_subject_id','eventname',mid_vars)]
# # Combine variables with main dataframe
# df<-join(df,mid_beh)

# ################################
# # Load ssrt data file from tabulated data
# ssrt<-read.csv(ssrt_file)
# # mssrt is mean stop signal reaction time; issrt is integrated stop signal reaction time, they give very similar results with brain-behav associations..
# ssrt<-ssrt[,c('src_subject_id','eventname','tfmri_sst_all_beh_total_mssrt','tfmri_sst_all_beh_total_issrt')]

# # Ensure variables are in a numeric format
# ssrt$tfmri_sst_all_beh_total_mssrt<-as.numeric(ssrt$tfmri_sst_all_beh_total_mssrt)
# ssrt$tfmri_sst_all_beh_total_issrt<-as.numeric(ssrt$tfmri_sst_all_beh_total_issrt)
# ssrt_vars<-c('tfmri_sst_all_beh_total_mssrt','tfmri_sst_all_beh_total_issrt')

# #Combine new variables with main dataframe
# df<-join(df, ssrt)

################################
# Load cbcl data file from tabulated data
# cbcl<-read.csv(cbcl_file)

# # Extract variables of interest
# cbcl<-cbcl[,c("src_subject_id","eventname",'cbcl_scr_syn_totprob_r', 'cbcl_scr_syn_internal_r', 'cbcl_scr_syn_external_r',
#               'cbcl_scr_syn_aggressive_r','cbcl_scr_syn_anxdep_r','cbcl_scr_syn_rulebreak_r','cbcl_scr_syn_attention_r', 'cbcl_scr_syn_social_r','cbcl_scr_syn_thought_r','cbcl_scr_syn_somatic_r','cbcl_scr_syn_withdep_r')]
# cbcl_vars=grep("cbcl_", names(cbcl), value=TRUE)

# # Ensure variables are in a numeric format
# cbcl[,cbcl_vars] <- sapply(cbcl[,cbcl_vars],as.numeric)

# # Combine new variables with main dataframe
# df<-join(df,cbcl)

################################
# Might add p-factor once generated for 5.0 and other timepoints

#add in p-factor from CBCL (Wes' code) -ONLY AVAILABLE AT BASELINE
# pfactor<-read.delim(paste0(data4path,"ABCD_lavaan_pfactor.csv"),header=T,sep=",")
# names(pfactor)[1]<-c("src_subject_id")
# pfactor_vars<-c('PF10_lavaan','PF10_INT_lavaan','PF10_EXT_lavaan')
# pfactor[,pfactor_vars] <- sapply(pfactor[,pfactor_vars],as.numeric)

# df<-join(df,pfactor,by='src_subject_id') #because this is only baseline, the other timepoints will just be populated with 'baseline' data since it joins by subject

# ################################
# # Load other mental health data file from tabulated data
# upps<-read.csv(upps_file)
# bisbas<-read.csv(bisbas_file)

# # Extract variables of interest
# upps<-upps[,c('src_subject_id','eventname','upps_y_ss_lack_of_perseverance', 'upps_y_ss_lack_of_planning', 'upps_y_ss_positive_urgency', 
#                         'upps_y_ss_negative_urgency', 'upps_y_ss_sensation_seeking')]

# bisbas<-bisbas[,c('src_subject_id','eventname','bis_y_ss_bas_drive', 'bis_y_ss_bas_fs', 'bis_y_ss_bas_rr', 'bis_y_ss_bis_sum')]

# # Join vars into new df
# mh <- join(upps, bisbas)

# # Ensure variables are in numeric format
# mh_vars<-c('upps_y_ss_lack_of_perseverance', 'upps_y_ss_lack_of_planning', 'upps_y_ss_positive_urgency',
#            'upps_y_ss_negative_urgency', 'upps_y_ss_sensation_seeking', 'bis_y_ss_bas_drive', 'bis_y_ss_bas_fs', 
#            'bis_y_ss_bas_rr', 'bis_y_ss_bis_sum')
# mh[,mh_vars] <- sapply(mh[,mh_vars],as.numeric)

# # Combine new variables with main dataframe
# df<-join(df, mh)

################################
################################
# Save the "df"

# if ( ! dir.exists(outpath) ) {
#         dir.create(outpath, recursive=TRUE)
# }
# write.table(df,paste0(outpath, '/','ABCD_rel6.0_demos_PCs_covars_behvars_pds_phys.csv'),col.names=TRUE, row.names=FALSE, sep=',')

#List of vars you can select from that were defined above: 
#imgincl_vars (img_incl flags for each imaging modality)
#glob_vars (total surface area, mean cortical thickness, total intracranial volume)
#nihtbx_vars (note not all scales are collected longitudinally)
#pea_vars (ravlt + matrix reasoning, from Pearson assessments)
#nback_vars
#dprime_vars (can likely skip, related to n-back accuracy)
#mid_vars
#ssrt_vars
#cbcl_vars
#mh_vars (includes upps and bis bas)

#e.g.:
# df_filt<-df[,c('src_subject_id','eventname','rel_family_id','interview_age','household.income', 'high.educ', 
#               'race.4level', 'hisp', 'mri_info_deviceserialnumber', 'mri_info_softwareversion', paste0("genesis_PC",1:10),
     #          imgincl_vars, nihtbx_vars, cbcl_vars)]

##feed variables you want from df_filt you've created into makeDesign.R, filter for specific timepoints etc



