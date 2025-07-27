if (!require("glmmTMB")){
  install.packages("glmmTMB")
}
library(glmmTMB)
library(dplyr)
library(ggplot2)
library(reshape2)
library(stringr)
library(R.matlab)
source("Func.R")

file_path = "./Simulation/DATA/"

file_list = list.files(path = file_path, full.names = TRUE)


glmmTMB_results = matrix(NA, nrow=length(file_list), ncol=16)

colnames(glmmTMB_results) = c('Var_FFX','Var_FID','Var_IID','r',
                              'Est_BETA1','Est_BETA2','Est_BETA3','Est_BETA4','Est_BETA5',
                              'Est_SIG1','Est_SIG2',
                              'TRUE_BETA1','TRUE_BETA2','TRUE_BETA3','TRUE_BETA4','TRUE_BETA5')

start_time = Sys.time()

for (i in 1:length(file_list)){
  
  dat = readMat(file_list[i])[["test.dat"]]
  
  glmmTMB_results[i, "Var_FFX"] = dat[[8]]
  glmmTMB_results[i, "Var_FID"] = dat[[9]]
  glmmTMB_results[i, "Var_IID"] = dat[[10]]
  glmmTMB_results[i, "r"] = dat[[11]]
  glmmTMB_results[i, str_detect(colnames(glmmTMB_results),"TRUE")] = dat[[7]]
  
  dat_tmp = data.frame(y_binary=dat[[2]],X=dat[[3]],fid=unlist(dat[[4]]),iid=unlist(dat[[5]]))
  
  fit_model = glmmTMB(y_binary ~ -1+`X.1`+`X.2`+`X.3`+`X.4`+`X.5`+(1|fid)+(1|iid), data=dat_tmp, family=binomial)
  
  tmp = summary(fit_model)
  
  glmmTMB_results[i, str_detect(colnames(glmmTMB_results),"Est_BETA")] = tmp$coefficients$cond[,1]
  glmmTMB_results[i, "Est_SIG1"]   = tmp[["varcor"]][["cond"]][["fid"]][1,1]
  glmmTMB_results[i, "Est_SIG2"]   = tmp[["varcor"]][["cond"]][["iid"]][1,1]
  
}

end_time = Sys.time() 
time_taken = end_time - start_time
print(time_taken)

write.table(glmmTMB_results,"./Simulation/glmmTMB.txt",row.names=F,quote=F,sep="\t")


#################################################################################
# Comparison - FEMAlog VS true
#################################################################################
CpDat = read.table("./Simulation/FEMA_binary.txt",sep="\t",header=T)
########### FE
diff  = cbind(CpDat[, 5:9] - CpDat[, 13:17])
colnames(diff) = paste0("BETA",1:5)
set   = CpDat[, 1:3]
plot_dat = cbind(set, diff)

pdf(
  file = "./Simulation/FEMAb_true_fe_diff.pdf",
  # type = "cairo", # ??????
  # res = 1000, # 300ppi ?ֱ???
  width = 10, height = 10,
  bg = "white" # ͸??????
)

plot_fe_diff(plot_dat,title="FEMA logistic VS true value - Fixed effects")

dev.off()


plot_dat = plot_dat %>% group_by(Var_FFX, Var_FID, Var_IID) %>%
  summarize(across(starts_with("BETA"), ~ sqrt(mean(.x^2))))


pdf(
  file = "./Simulation/FEMAb_true_fe_rmse.pdf",
  # type = "cairo", # ??????
  # res = 1000, # 300ppi ?ֱ???
  width = 10, height = 10,
  bg = "white" # ͸??????
)

plot_fe_rmse(plot_dat,title="FEMA logistic VS true value - Fixed effects")

dev.off()

########### RE
diff  = cbind(CpDat[, 10:11] - CpDat[, 2:3])
colnames(diff) = c("FID","IID")
set   = CpDat[, 1:3]
plot_dat = cbind(set, diff)

pdf(
  file = "./Simulation/FEMAb_true_re_diff.pdf",
  # type = "cairo", # ??????
  # res = 1000, # 300ppi ?ֱ???
  width = 10, height = 10,
  bg = "white" # ͸??????
)

plot_re_diff(plot_dat,title="FEMA logistic VS true value - Random effects")

dev.off()


plot_dat = plot_dat %>% group_by(Var_FFX, Var_FID, Var_IID) %>%
  summarize(across(c("FID","IID"), ~ sqrt(mean(.x^2))))


pdf(
  file = "./Simulation/FEMAb_true_re_rmse.pdf",
  # type = "cairo", # ??????
  # res = 1000, # 300ppi ?ֱ???
  width = 10, height = 10,
  bg = "white" # ͸??????
)

plot_re_rmse(plot_dat,title="FEMA logistic VS true value - Random effects")

dev.off()



#################################################################################
# Comparison - glmmTMB VS true
#################################################################################
CpDat = read.table("./Simulation/glmmTMB.txt",sep="\t",header=T)
########### FE
diff  = cbind(CpDat[, 5:9] - CpDat[, 12:16])
colnames(diff) = paste0("BETA",1:5)
set   = CpDat[, 1:3]
plot_dat = cbind(set, diff)

pdf(
  file = "./Simulation/glmmTMB_true_fe_diff.pdf",
  # type = "cairo", # ??????
  # res = 1000, # 300ppi ?ֱ???
  width = 10, height = 10,
  bg = "white" # ͸??????
)

plot_fe_diff(plot_dat, title="glmmTMB VS true value - Fixed effects")

dev.off()


plot_dat = plot_dat %>% group_by(Var_FFX, Var_FID, Var_IID) %>%
  summarize(across(starts_with("BETA"), ~ sqrt(mean(.x^2))))


pdf(
  file = "./Simulation/glmmTMB_true_fe_rmse.pdf",
  # type = "cairo", # ??????
  # res = 1000, # 300ppi ?ֱ???
  width = 10, height = 10,
  bg = "white" # ͸??????
)

plot_fe_rmse(plot_dat,title="glmmTMB VS true value - Fixed effects")

dev.off()

########### RE
diff  = cbind(CpDat[, 10:11] - CpDat[, 2:3])
colnames(diff) = c("FID","IID")
set   = CpDat[, 1:3]
plot_dat = cbind(set, diff)

pdf(
  file = "./Simulation/glmmTMB_true_re_diff.pdf",
  # type = "cairo", # ??????
  # res = 1000, # 300ppi ?ֱ???
  width = 10, height = 10,
  bg = "white" # ͸??????
)

plot_re_diff(plot_dat,title="glmmTMB VS true value - Random effects")

dev.off()


plot_dat = plot_dat %>% group_by(Var_FFX, Var_FID, Var_IID) %>%
  summarize(across(c("FID","IID"), ~ sqrt(mean(.x^2))))


pdf(
  file = "./Simulation/glmmTMB_true_re_rmse.pdf",
  # type = "cairo", # ??????
  # res = 1000, # 300ppi ?ֱ???
  width = 10, height = 10,
  bg = "white" # ͸??????
)

plot_re_rmse(plot_dat,title="glmmTMB VS true value - Random effects")

dev.off()


#################################################################################
# Comparison - FEMAb VS glmmTMB
#################################################################################
FEMAb_results = read.table("./Simulation/FEMA_binary.txt",sep="\t",header=T)%>%arrange(VarFFX,VarFID,r)
glmmTMB_results = read.table("./Simulation/glmmTMB.txt",sep="\t",header=T)%>%arrange(VarFFX,VarFID,r)
########### FE
diff  = cbind(FEMAb_results[, 5:9] - glmmTMB_results[, 5:9])
# diff  = cbind(abs(CpDat[, 1:5]-FEMAb_results[,12:16]) - abs(FEMAb_results[,4:8]-FEMAb_results[, 12:16]))
colnames(diff) = paste0("BETA",1:5)
set   = FEMAb_results[, 1:3]
plot_dat = cbind(set, diff)

pdf(
  file = "./Simulation/glmmTMB_FEMAb_fe_diff.pdf",
  # type = "cairo", # ??????
  # res = 1000, # 300ppi ?ֱ???
  width = 10, height = 10,
  bg = "white" # ͸??????
)

plot_fe_diff(plot_dat,title="FEMA logistic VS glmmTMB - Fixed effects")

dev.off()


plot_dat = plot_dat %>% group_by(Var_FFX, Var_FID, Var_IID) %>%
  summarize(across(starts_with("BETA"), ~ sqrt(mean(.x^2))))


pdf(
  file = "./Simulation/glmmTMB_FEMAb_fe_rmse.pdf",
  # type = "cairo", # ??????
  # res = 1000, # 300ppi ?ֱ???
  width = 10, height = 10,
  bg = "white" # ͸??????
)

plot_fe_rmse(plot_dat,title="FEMA logistic VS glmmTMB - Fixed effects")

dev.off()

########### RE
# diff  = cbind(abs(FEMAb_results[, 10:11]-FEMAb_results[,2:3]) - abs(glmmTMB_results[, 6:7]-FEMAb_results[,2:3]))
diff  = cbind(FEMAb_results[, 10:11] - glmmTMB_results[, 10:11])
colnames(diff) = c("FID","IID")
set   = FEMAb_results[, 1:3]
plot_dat = cbind(set, diff)

pdf(
  file = "./Simulation/glmmTMB_FEMAb_re_diff.pdf",
  # type = "cairo", # ??????
  # res = 1000, # 300ppi ?ֱ???
  width = 10, height = 10,
  bg = "white" # ͸??????
)

plot_re_diff(plot_dat,title="FEMA logistic VS glmmTMB - Random effects")

dev.off()


plot_dat = plot_dat %>% group_by(Var_FFX, Var_FID, Var_IID) %>%
  summarize(across(c("FID","IID"), ~ sqrt(mean(.x^2))))


pdf(
  file = "./Simulation/glmmTMB_FEMAb_re_rmse.pdf",
  # type = "cairo", # ??????
  # res = 1000, # 300ppi ?ֱ???
  width = 10, height = 10,
  bg = "white" # ͸??????
)

plot_re_rmse(plot_dat,title="FEMA logistic VS glmmTMB - Random effects")

dev.off()
