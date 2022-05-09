###############################################################

# This script describes how to create design matrices in the format 
# expected by the FEMA package. Users will need to save a local
# copy and modify paths and variables of interest appropriately.

###############################################################
# First need to make sure that R know where to find the makeDesign function 
# which makes (obv) the matrices 
source('/path/to/cmig_tools_utils/r/makeDesign.R')

###############################################################
# The function makeDesign reads variables from an R data frame. This 
# data frame can be loaded as the official ABCD RDS file or any other 
# data frame, for example if you have saved individual instruments as a 
# mini RDS as in the example script createMiniRDS.R

# Load the data 
ndafile <- '/path/to/nda4.0_offrel.RDS'
nda <- readRDS(ndafile)

# It is STRONGLY recommended that you only include subjects which pass QC 
# for the imaging variable you are interested in testing.  For example this step
# filters using the imgincl_dmri_include variable from the abcd_imgincl01.txt instrument. 
idx_dmri_inc <- which(nda$imgincl_dmri_include==1)
nda_dmri_inc <- nda[idx_dmri_inc,]

# Define the path to the directory where you would like to save out your design matrix 
outpath <- '/path/to/your/fave/output/directory'

###############################################################
# The following section describes how to set up a design matrix for CROSS-SECTIONAL 
# analysis, ie to model the inter-subject variability associated with your predictors 
# while accounting for repeated measures if you include multiple time points. 

# Define the name of your design matrix file 
fname <- 'designMat1.txt'
# and it's full path to your fave output directory
outfile <- paste0(outpath, '/', fname) 

# makeDesign encodes continuous and categorical variables differently, therefore they are 
# specified using different flags. Continuous variables are added using the "contvar" flag. 
# Here we include age and the genetic PCs as continuous variables.
contvar <- c('interview_age', 'PC1', 'PC2', 'PC3', 'PC4', 'PC5', 'PC6', 'PC7', 'PC8','PC9', 'PC10')

# Categorical variables are dummy coded and added using the "catvar" flag. makeDesign automatically
# includes an intercept. For each categorical variable one category is defined as the reference category 
# and that column is dropped. makeDesign also checks whether your matrix is rank deficient and if so
# will automatically drop further categories to avoid this. Here we include sex, scanner info and 
# SES demographics as categorical.
catvar <- c('sex', 'high.educ', 'hisp', 'household.income', 'mri_info_deviceserialnumber', 'mri_info_softwareversion')

# The time points for which we wish to extract data are specified. You do not have to specify multiple 
# time points, however if you do, specify them in chronoligical order ( this will be important for 
# longitudinal modelling)
time <- c('baseline_year_1_arm_1', '2_year_follow_up_y_arm_1') # order matters! start with baseline!

# You can specify whether you wish to demean your continuous variables or not. the default is set to 
# demean=TRUE. Other flags include "delta" for longitudinal modelling, "interact" to include interactions, 
# "subjs" if you wish to filter by subject and "quadratic" if you wish to include a quadratic term 
# (more details on these in the following sections). The defaults to these are set to null. 

# Now run makeDesign! 
makeDesign(nda_dmri_inc, outfile, time, contvar=contvar, catvar=catvar, delta=NULL, interact=NULL, subjs=NULL, demean=TRUE, quadratic=NULL)

############################################################
# The following describes how to set up a design matrix for LONGITUDINAL analysis, ie to model the both 
# the inter-subject and intra-subject variability associated with predictors of your choice across multiple 
# time points. The following example creates a design matrix to model longitudinal changes associated with age. 

# First specify a different name for your design matrix file 
fname <- 'designMat2.txt'
outfile <- paste0(outpath, '/', fname) 

# Seeing as we want to model variability associated with change in age, it is no longer included in 
# the "contvar" flag. Age is now specified using the "delta" flag, which separates the age variable 
# into its baseline value (saved as interview_age_base) and the change in age (saved as interview_age_delta)
delta <- c('interview_age'); 

# The other continuous variables are specified as before.
contvar <- c('PC1', 'PC2', 'PC3', 'PC4', 'PC5', 'PC6', 'PC7', 'PC8','PC9', 'PC10')

# Categorical variables are specified as before. 
catvar <- c('sex', 'high.educ', 'hisp', 'household.income', 'mri_info_deviceserialnumber', 'mri_info_softwareversion')

# The time points for which we wish to extract data have to be specified in chronoligical order, currently 
# it only works for two time points. 
time <- c('baseline_year_1_arm_1', '2_year_follow_up_y_arm_1') # order matters! start with baseline!

# Now run makeDesign! 
makeDesign(nda_dmri_inc, outfile, time, contvar=contvar, catvar=catvar, delta=delta, interact=NULL, subjs=NULL, demean=TRUE, quadratic=NULL)

############################################################
# The following section describes how to create a design matrix which includes INTERACTIONS. This e
# example uses a cross-sectional model and includes an interaction between sex and age

# Again we specify a new file name 
fname <- 'designMat3.txt'
outfile <- paste0(outpath, '/', fname) 

# We specify the continuous variables. Seeing as we are not modeling longitudinal change with age, 
# interview_age is back in with continuous variables. Both age and sex are included explicitly and 
# seperately from the interaction terms that includes them both   
contvar <- c('interview_age', 'PC1', 'PC2', 'PC3', 'PC4', 'PC5', 'PC6', 'PC7', 'PC8','PC9', 'PC10')

# The categorical variables are specified as before (sex is in there!)
catvar <- c('sex', 'high.educ', 'hisp', 'household.income', 'mri_info_deviceserialnumber', 'mri_info_softwareversion')

# The interaction term is specified using the convention "interview_age*sex" 
interact <- c('interview_age*sex'); 

# The time points which we wish to extract are specified as before 
time <- c('baseline_year_1_arm_1', '2_year_follow_up_y_arm_1') # order matters! start with baseline!

# Now run make Design! 
makeDesign(nda_dmri_inc, outfile, time, contvar=contvar, catvar=catvar, delta=NULL, interact=interact, subjs=NULL, demean=TRUE, quadratic=NULL)

############################################################

# The following section describes how to create design matrices for a MEDIATION analysis.
# This code will produce two design matrices one that is termed '_reduced' and the other which
# which is termed '_full'. These are nested models. Here the mediator is body mass index (BMI).
# The FULL model will include all of the independent variables specified in contvar and catvar.
# The REDUCED model will have all of these excluding the mediator (BMI).
# The mediator is always placed as the penultimate column (before the intercept) in the deisgn matrix.
# Using this code ensures the two design matrices have the same number of subjects and are nested.

# Again we specify a new file name 
fname <- 'designMat4.txt'
outfile <- paste0(outpath, '/', fname) 

# We specify the continuous variables. Here we have additionally included 
contvar <- c('interview_age', 'PC1', 'PC2', 'PC3', 'PC4', 'PC5', 'PC6', 'PC7', 'PC8','PC9', 'PC10','anthro_bmi_calc')

# The categorical variables are specified as before (sex is in there!)
catvar <- c('sex', 'high.educ', 'hisp', 'household.income', 'mri_info_deviceserialnumber', 'mri_info_softwareversion')

# The interaction term is specified using the convention "interview_age*sex" 
interact <- NULL; 

# The time points which we wish to extract are specified as before 
time <- c('baseline_year_1_arm_1', '2_year_follow_up_y_arm_1') # order matters! start with baseline!

# Now run make Design! 
makeDesign(nda_dmri_inc, outfile, time, contvar=contvar, catvar=catvar, delta=NULL, interact=interact, subjs=NULL, demean=TRUE, quadratic=NULL, mediator='anthro_bmi_calc')

