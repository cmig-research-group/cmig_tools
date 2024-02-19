makeDEAPdemos<-function(datapath){

  #Need to source('~/github/cmig_library/loadtxt.R')

library(dplyr)
library(data.table)

  files<-dir(datapath)

  txt_tmp<-grep('txt',files)
  csv_tmp<-grep('csv',files)
  if (!isempty(txt_tmp)){
    ext='txt';
  } else if (!isempty(csv_tmp)) {
    ext='csv'
  } else if (isempty(txt_tmp) & isempty(csv_tmp)) {
    warning('Files in datapath not csv or txt files. Script will not read files correctly.')
  }

  acsfile<-files[grep('acspsw03',files)]
  basedemsfile<-files[grep('pdem02',files)]
  longdemsfile<-files[grep('abcd_lpds01',files)]
  ltsitefile<-files[grep('abcd_lt01',files)]

  if (!isempty(grep('5.?',datapath))){
    longdems5.0file <- files[grep('abcd_p_demo',files)]
    ltsite5.0file <- files[grep('abcd_y_lt',files)]
    alldems <- read.table(paste0(datapath,'/',longdems5.0file), header=T, sep=",")
    lt_site <- read.table(paste0(datapath,'/',ltsite5.0file), header=T, sep=",")
  } else {
    if (ext=='csv'){
        acsdems<-read.table(paste0(datapath,'/',acsfile), header=T, sep=",")
        basedems<-read.table(paste0(datapath,'/',basedemsfile), header=T, sep=",")
        longdems<-read.table(paste0(datapath,'/',longdemsfile), header=T, sep=",")
        lt_site<-read.table(paste0(datapath,'/',ltsitefile), header=T, sep=",")
      } else if (ext=='txt'){
        acsdems<-loadtxt(paste0(datapath,'/',acsfile))
        acsdems<-acsdems[,-c(which(colnames(acsdems) %in% c('collection_id','dataset_id')))]
        basedems<-loadtxt(paste0(datapath,'/',basedemsfile))
        basedems<-basedems[,-c(which(colnames(basedems) %in% c('collection_id','dataset_id')))]
        longdems<-loadtxt(paste0(datapath,'/',longdemsfile))
        longdems<-longdems[,-c(which(colnames(longdems) %in% c('collection_id','dataset_id')))]
        lt_site<-loadtxt(paste0(datapath,'/',ltsitefile))
        lt_site<-lt_site[,-c(which(colnames(lt_site) %in% c('collection_id','dataset_id')))]
      }
    alldems<-join(longdems,basedems)
  }

alldems[,'demo_prnt_ed_p']<-coalesce(alldems$demo_prnt_ed_v2,alldems$demo_prnt_ed_v2_l)
alldems[,'demo_prnt_ed_p']<-coalesce(alldems$demo_prnt_ed_p,alldems$demo_prnt_ed_v2_2yr_l)
alldems$demo_prtnr_ed_v2_l = as.integer(alldems$demo_prtnr_ed_v2_l)
alldems[,'demo_prtnr_ed_p']<-coalesce(alldems$demo_prtnr_ed_v2,alldems$demo_prtnr_ed_v2_l)
alldems[,'demo_prtnr_ed_p']<-coalesce(alldems$demo_prtnr_ed_p,alldems$demo_prtnr_ed_v2_2yr_l)

alldems[,'demo_comb_income_p']<-coalesce(alldems$demo_comb_income_v2,alldems$demo_comb_income_v2_l)

alldems[,'demo_prnt_marital_p']<-coalesce(alldems$demo_prnt_marital_v2, alldems$demo_prnt_marital_v2_l)

alldems<-join(alldems,lt_site)

if (!isempty(grep('5.?',datapath))) {
  deap_demos<-data.frame(alldems[,c('src_subject_id','eventname','interview_date','interview_age')])
} else if (!isempty(grep('4.0',datapath))){
  deap_demos<-data.frame(alldems[,c('src_subject_id','eventname','interview_date','interview_age','sex','sched_delay','sched_hybrid')])
} else if (!isempty(grep('3.0',datapath))){
  deap_demos<-data.frame(alldems[,c('src_subject_id','eventname','interview_date','interview_age','sex')])
} else {
  deap_demos<-data.frame(alldems[,c('src_subject_id','eventname','interview_date','interview_age','sex')])
  warning('Datapath does not include release number. Hybrid variables only extracted for 4.0.')
}

deap_demos$abcd_site = alldems$site_id_l #site_id_l is in longitudianl tracking instrument

if (!isempty(grep('5.?',datapath))) {
  sextmp = data.frame(alldems[alldems$eventname=='baseline_year_1_arm_1',c('src_subject_id','demo_sex_v2')])
  sextmp$sex = recode(as.factor(sextmp$demo_sex_v2), "1" = "M","2" = "F", "3" = "I")
  deap_demos<-join(deap_demos,sextmp[,c('src_subject_id','sex')], by='src_subject_id', match = "all")
} else {
  deap_demos$sex[which(deap_demos$sex=="")]=NA
  deap_demos$sex=factor( deap_demos$sex, levels= c("F","M"))
}

### Household income
household.income = alldems$demo_comb_income_p
household.income[alldems$demo_comb_income_p == "1"] = 1 # "[<50K]"
household.income[alldems$demo_comb_income_p == "2"] = 1 # "[<50K]"
household.income[alldems$demo_comb_income_p == "3"] = 1 # "[<50K]"
household.income[alldems$demo_comb_income_p == "4"] = 1 # "[<50K]"
household.income[alldems$demo_comb_income_p == "5"] = 1 # "[<50K]"
household.income[alldems$demo_comb_income_p == "6"] = 1 # "[<50K]"
household.income[alldems$demo_comb_income_p == "7"] = 2 # "[>=50K & <100K]"
household.income[alldems$demo_comb_income_p == "8"] = 2 # "[>=50K & <100K]"
household.income[alldems$demo_comb_income_p == "9"] = 3 # "[>=100K]"
household.income[alldems$demo_comb_income_p == "10"] = 3 # "[>=100K]"
household.income[alldems$demo_comb_income_p == "777"] = NA
household.income[alldems$demo_comb_income_p == "999"] = NA
household.income[household.income %in% c(NA, "999", "777")] = NA
deap_demos$household.income = factor( household.income, levels= 1:3, labels = c("[<50K]", "[>=50K & <100K]", "[>=100K]") )

### Household income (continuous) - assign value based on middle of category
household.income_cont = alldems$demo_comb_income_p
household.income_cont[alldems$demo_comb_income_p == "1"] = 2500 # Less than $5,000
household.income_cont[alldems$demo_comb_income_p == "2"] = 8500 # $5,000 through $11,999
household.income_cont[alldems$demo_comb_income_p == "3"] = 14000 # $12,000 through $15,999
household.income_cont[alldems$demo_comb_income_p == "4"] = 20500 # $16,000 through $24,999
household.income_cont[alldems$demo_comb_income_p == "5"] = 30000 # $25,000 through $34,999;
household.income_cont[alldems$demo_comb_income_p == "6"] = 42500 # $35,000 through $49,999
household.income_cont[alldems$demo_comb_income_p == "7"] = 62500 # $50,000 through $74,999
household.income_cont[alldems$demo_comb_income_p == "8"] = 87500 # $75,000 through $99,999
household.income_cont[alldems$demo_comb_income_p == "9"] = 150000 # $100,000 through $199,999
household.income_cont[alldems$demo_comb_income_p == "10"] = 250000 # $200,000 and greater
household.income_cont[alldems$demo_comb_income_p == "777"] = NA # Refuse to answer
household.income_cont[alldems$demo_comb_income_p == "999"] = NA # Don't know
household.income_cont[household.income_cont %in% c(NA, "999", "777")] = NA
deap_demos$household.income_cont = household.income_cont

# Household income (10 level)
household.income_10level = alldems$demo_comb_income_p
household.income_10level[alldems$demo_comb_income_p == "777"] = NA # Refuse to answer
household.income_10level[alldems$demo_comb_income_p == "999"] = NA # Don't know
household.income_10level[household.income_10level %in% c(NA, "999", "777")] = NA
deap_demos$household.income_10level = household.income_10level

#highest education: 5 different levels. These levels correspond to the numbers published by the American Community Survey (ACS).
high.educ1 = alldems$demo_prnt_ed_p
high.educ2 = alldems$demo_prtnr_ed_p
high.educ1[which(high.educ1 == "999")] = NA
high.educ2[which(high.educ2 == "999")] = NA
high.educ1[which(high.educ1 == "777")] = NA
high.educ2[which(high.educ2 == "777")] = NA
high.educ1[which(high.educ1 == "22" | high.educ1=="23")] = 15 #22 and 23 = some college --> lower level than 18+
high.educ2[which(high.educ2 == "22" | high.educ2=="23")] = 15
high.educ = pmax(as.numeric(as.character(high.educ1)), as.numeric(as.character(high.educ2)), na.rm=T)
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
deap_demos$high.educ = factor( high.educ, levels= 1:5, labels = c("< HS Diploma","HS Diploma/GED","Some College","Bachelor","Post Graduate Degree") )
### Marital status
married = rep(NA, length(alldems$demo_prnt_marital_p))
married[alldems$demo_prnt_marital_p == 1] = 1
married[alldems$demo_prnt_marital_p %in% 2:6] = 0
deap_demos$married = factor( married, levels= 0:1, labels = c("No", "Yes") )
#Add another variable that also includes couples that just live together.
#married.livingtogether = rep(NA, length(alldems$demo_prnt_marital_p))
#married.livingtogether[alldems$demo_prnt_marital_p %in% c(1,6)] = 1
#married.livingtogether[alldems$demo_prnt_marital_p %in% 2:5] = 0
#deap_demos$married.or.livingtogether = factor( married.livingtogether, levels= 0:1, labels = c("no", "yes") )

#########################################

#Race of child not asked again longitudinally? Repeating this from baseline
race_cols = c("demo_ethn_v2", "demo_race_a_p___10", "demo_race_a_p___11","demo_race_a_p___12", "demo_race_a_p___13",
              "demo_race_a_p___14", "demo_race_a_p___15", "demo_race_a_p___16", "demo_race_a_p___17",
              "demo_race_a_p___18", "demo_race_a_p___19", "demo_race_a_p___20", "demo_race_a_p___21", "demo_race_a_p___22","demo_race_a_p___23",
              "demo_race_a_p___24", "demo_race_a_p___25",
              "demo_race_a_p___77", "demo_race_a_p___99")
if (isempty(grep('5.?',datapath))){
  racedf<-basedems[,match(c('src_subject_id','eventname',race_cols), names(basedems))]
  raceind<-which(colnames(alldems) %in% race_cols)
  alldems<-alldems[,-c(raceind)]

  alldems<-merge(alldems, racedf, by.x = 'src_subject_id', by.y = 'src_subject_id')
  racedf<-basedems[,match(c('src_subject_id','eventname',race_cols), names(basedems))] # OLD - from 4.0 and prior
} else {
  racedf<-alldems[,match(c('src_subject_id','eventname',race_cols), names(alldems))]
  raceind<-which(colnames(alldems) %in% race_cols)
  alldems<-alldems[,-c(raceind)]

  alldems<-join(alldems, racedf[racedf$eventname=='baseline_year_1_arm_1',], by='src_subject_id', match = "all")
  # racedf<-alldems[,match(c('src_subject_id','eventname',race_cols), names(alldems))]
}

raceind<-which(colnames(alldems) %in% race_cols)
dat = data.table(alldems[,raceind])

# White
dat[, white:= (demo_race_a_p___10 == 1)*1 ]

# Black
dat[, black:= (demo_race_a_p___11 == 1)*1 ]

# Asian
dat[, asian:= 0]
dat[ (demo_race_a_p___18 == 1 | demo_race_a_p___19 == 1 | demo_race_a_p___20 == 1 |
        demo_race_a_p___21 == 1 | demo_race_a_p___22 == 1 | demo_race_a_p___23 == 1 |
        demo_race_a_p___24==1), asian:= 1 ]

# AIAN: American Indian and Alaska Native
dat[, aian:= 0]
dat[ (demo_race_a_p___12 == 1 | demo_race_a_p___13 == 1), aian:=1 ]

#NHPI: Native Hawaiian and Other Pacific
dat[, nhpi:= 0]
dat[ demo_race_a_p___14 == 1 | demo_race_a_p___15 == 1 | demo_race_a_p___16 == 1 |
       demo_race_a_p___17 == 1, nhpi:= 1 ]

# Other
dat[, other:= 0 ]
dat[ demo_race_a_p___25 == 1, other:= 1 ]

# Mixed
dat[, mixed:= (white + black + asian + aian + nhpi + other)]
dat[, table(mixed, useNA = "if")]
dat[ mixed <= 1, mixed:= 0]
dat[ mixed > 1, mixed:= 1]
dat[, table(mixed, useNA = "if")]

# Race 4 level
dat[ white == 1, race.4level:= 1]
dat[ black == 1,race.4level:= 2]
dat[ asian == 1,race.4level:= 3]
dat[ aian == 1,race.4level:= 4]
dat[ nhpi == 1,race.4level:= 4]
dat[ other == 1,race.4level:= 4]
dat[ mixed == 1,race.4level:= 4]
dat[, table(race.4level, useNA = "if") ]

dat$race.4level<- factor(dat$race.4level,
                         levels = 1:4,
                         labels = c("White","Black","Asian","Other/Mixed"))

dat$race.eth[dat$race.4level==1] = "White"
dat$race.eth[dat$race.4level==2] = "Black"
dat$race.eth[dat$race.4level==3] ="Asian"
dat$race.eth[dat$race.4level==4] = "Other/Mixed"
dat[, table(race.4level, useNA = "if") ]

# Race 6 level
dat[white==1,race.6level:=1]
dat[black==1,race.6level:=2]
dat[asian==1,race.6level:=3]
dat[aian==1,race.6level:=4]
dat[nhpi==1,race.6level:=4]
dat[other==1,race.6level:=5]
dat[mixed==1,race.6level:=6]
dat[, table(race.6level,useNA="if") ]

dat$race.6level<- factor(dat$race.6level,
                         levels=1:6,
                         labels= c("White","Black","Asian","AIAN/NHPI","Other","Mixed"))
dat$race.eth[dat$race.6level==1] ="White"
dat$race.eth[dat$race.6level==2] ="Black"
dat$race.eth[dat$race.6level==3] ="Asian"
dat$race.eth[dat$race.6level==4] ="AIAN/NHPI"
dat$race.eth[dat$race.6level==5] ="Other"
dat$race.eth[dat$race.6level==6] ="Mixed"
dat[, table(race.6level,useNA="if") ]

# Hispanic
dat$hisp=NA;
indx.1=which(dat$demo_ethn_v2==1)
indx.0=which(dat$demo_ethn_v2==2)
dat$hisp[indx.1]=1;
dat$hisp[indx.0]=0;


dat$hisp<- factor(dat$hisp,
                  levels=0:1,
                  labels= c("No","Yes"))

dat$src_subject_id<-alldems$src_subject_id

if (isempty(grep('5.?',datapath))) {
  dat$eventname<-alldems$eventname.x
} else {
  dat$eventname<-alldems$eventname
}

racetmp<-data.frame(dat[,c('src_subject_id','eventname','race.4level','race.6level','hisp')])

deap_demos<-join(deap_demos,racetmp)

############################################

#Replicating family ID over all time points
if (!isempty(acsfile)){ # <5.0
  alldems<-merge(alldems, acsdems[acsdems$eventname=='baseline_year_1_arm_1',c('src_subject_id','rel_family_id','rel_group_id')])
  alldems$eventname<-alldems$eventname.x
  famtmp<-data.frame(alldems[,c('src_subject_id','eventname','rel_family_id','rel_group_id')])
  deap_demos<-join(deap_demos,famtmp)
  deap_demos$rel_group_id<-factor(deap_demos$rel_group_id)

} else{ #5.0
  famtmp<-data.frame(alldems[alldems$eventname=='baseline_year_1_arm_1',c('src_subject_id','rel_family_id','rel_birth_id')])
  deap_demos<-join(deap_demos,famtmp, by='src_subject_id', match = "all")
  deap_demos$rel_birth_id<-factor(deap_demos$rel_birth_id)
}

return(deap_demos)
}
