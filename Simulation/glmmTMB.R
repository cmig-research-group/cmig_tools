if (!require("glmmTMB")){
  install.packages("glmmTMB")
}
library(glmmTMB)
library(dplyr)
library(ggplot2)
library(ggnewscale)
library(ggpubr)
library(grid)
library(reshape2)
library(stringr)
library(R.matlab)
library(doParallel)
library(parallel)
library(expression)


# file_path = "./FEMA/DATA/"
# 
# nXvars = 5
# nRepeats = 100
# nsamples = 10
# N = 10000
# 
# n_cores = 50
# 
# if (exists("cl")) {
#   stopCluster(cl)
#   rm(cl)
# }
# if (requireNamespace("doParallel", quietly = TRUE)) {
#   doParallel::stopImplicitCluster()
# }
# try({
#   parallel::stopCluster(parallel::makeCluster(0))
# }, silent = TRUE)
# 
# 
# cl = makeCluster(n_cores)
# registerDoParallel(cl)
# 
# start_time = Sys.time()
# 
# for (M in c(5,10,20)){
#   for (V in c(1,4,7)){
#     
#     cat(sprintf("Start processing M=%d, V=%d\n", M, V))
#     
#     glmmTMB_results = foreach(
#       r = 1:nRepeats,
#       .combine = "rbind",
#       .packages = c("R.matlab", "glmmTMB"),
#       .export = c("N", "M", "V", "nsamples", "nXvars", "file_path")
#     ) %dopar% {
#       
#       beta = matrix(NA, nsamples, nXvars)
#       beta_hat_sum = matrix(NA, nsamples, nXvars)
#       sig2mat_sum = matrix(NA, nsamples, 1)
#       
#       for (i in 1:nsamples) {
#         
#         fileName = paste0(file_path,sprintf("N%d_M%d_V%d/dat_%d-%d.mat",N,M,V,r,i))
#         
#         dat = readMat(fileName)[["data"]]
#         
#         dat_tmp = data.frame(
#           y_binary = dat[[2]],
#           X = dat[[3]],
#           iid = unlist(dat[[5]])
#         )
#         
#         fit_model = glmmTMB(
#           y_binary ~ -1 + `X.1` + `X.2` + `X.3` + `X.4` + `X.5` + (1 | iid),
#           data = dat_tmp,
#           family = binomial
#         )
#         
#         tmp = summary(fit_model)
#         
#         beta[i, ] = dat[[7]]
#         beta_hat_sum[i, ] = tmp[["coefficients"]][["cond"]][, 1]
#         sig2mat_sum[i, ] = tmp[["varcor"]][["cond"]][["iid"]][1, 1]
#       }
#       
#       c(colMeans(beta), colMeans(beta_hat_sum), mean(sig2mat_sum))
#       
#     }
#     
#     colnames(glmmTMB_results) = c(
#       'TRUE_BETA1', 'TRUE_BETA2', 'TRUE_BETA3', 'TRUE_BETA4', 'TRUE_BETA5',
#       'Est_BETA1', 'Est_BETA2', 'Est_BETA3', 'Est_BETA4', 'Est_BETA5', 'Est_Var_IID'
#     )
#     
#     save_path = sprintf("./FEMA/N%d_M%d_V%d.txt",N,M,V)
#     write.table(glmmTMB_results, save_path, row.names=F, quote=F, sep="\t")
#     
#     cat(sprintf("Completed M=%d, V=%d\n", M, V))
#     
#   }
# }
# 
# end_time = Sys.time() 
# time_taken = end_time - start_time
# print(time_taken)
# 
# stopCluster(cl)



#################################################################################
# Simulation 1 - summary results
#################################################################################
calculate_rmse = function(true_values, estimated_values) {
  sqrt(mean((estimated_values - true_values)^2))
}

nXvars = 5
N = 10000

beta_rmse = matrix(NA, nrow=9, ncol=5)
colnames(beta_rmse) = c("M","V","FEMA_true","glmmTMB_true","FEMA_glmmTMB")
beta_FEMA = list()
beta_TMB = list()
sigma_est = list()
sigma_rmse = matrix(NA, nrow=9, ncol=5)
colnames(sigma_rmse) = c("M","V","FEMA_true","glmmTMB_true","FEMA_glmmTMB")
idx = 1
for (M in c(5,10,20)){
  for (V in c(1,4,7)){
    
    # FEMA
    fileName = sprintf("./Simulation_1/FEMA_binary/N%d_M%d_V%d.txt",N,M,V)
    FEMA_dat = read.table(fileName,sep="\t",header=T)
    
    # glmmTMB
    fileName = sprintf("./Simulation_1/glmmTMB/N%d_M%d_V%d.txt",N,M,V)
    TMB_dat = read.table(fileName,sep="\t",header=T)
    
    #######################################################BETA
    #FEMA_true
    FEMA_true = numeric(5)
    
    for (i in 1:5) {
      FEMA_true[i] = calculate_rmse(
        FEMA_dat[[paste0("TRUE_BETA", i)]],
        FEMA_dat[[paste0("Est_BETA", i)]]
      )
    }
    
    #glmmTMB_true
    glmmTMB_true = numeric(5)
    
    for (i in 1:5) {
      glmmTMB_true[i] = calculate_rmse(
        TMB_dat[[paste0("TRUE_BETA", i)]],
        TMB_dat[[paste0("Est_BETA", i)]]
      )
    }
    
    #FEMA_glmmTMB
    FEMA_glmmTMB = numeric(5)
    
    for (i in 1:5) {
      FEMA_glmmTMB[i] = calculate_rmse(
        FEMA_dat[[paste0("Est_BETA", i)]],
        TMB_dat[[paste0("Est_BETA", i)]]
      )
    }
    
    # summary
    beta_rmse[idx, "M"] = M
    beta_rmse[idx, "V"] = V
    beta_rmse[idx, "FEMA_true"] = mean(FEMA_true)
    beta_rmse[idx, "glmmTMB_true"] = mean(glmmTMB_true)
    beta_rmse[idx, "FEMA_glmmTMB"] = mean(FEMA_glmmTMB)
    
    beta_FEMA[[idx]] = data.frame(
      M = M,
      V = V,
      unlist(FEMA_dat[,1:5]),
      unlist(FEMA_dat[,6:10])
    )
    colnames(beta_FEMA[[idx]])=c("M","V","True","FEMA_Est")
    
    beta_TMB[[idx]] = data.frame(
      M = M,
      V = V,
      unlist(TMB_dat[,1:5]),
      unlist(TMB_dat[,6:10])
    )
    colnames(beta_TMB[[idx]])=c("M","V","True","TMB_Est")
    
    #######################################################SIGMA
    #FEMA_true
    FEMA_true = numeric(5)
    
    for (i in 1:5) {
      FEMA_true[i] = calculate_rmse(
        FEMA_dat[["Est_Var_IID"]],
        V/10
      )
    }
    
    
    #glmmTMB_true
    glmmTMB_true = numeric(5)
    
    for (i in 1:5) {
      glmmTMB_true[i] = calculate_rmse(
        TMB_dat[["Est_Var_IID"]],
        V/10
      )
    }
    
    #FEMA_glmmTMB
    FEMA_glmmTMB = numeric(5)
    
    for (i in 1:5) {
      FEMA_glmmTMB[i] = calculate_rmse(
        FEMA_dat[["Est_Var_IID"]],
        TMB_dat[["Est_Var_IID"]]
      )
    }
    
    # summary
    sigma_rmse[idx, "M"] = M
    sigma_rmse[idx, "V"] = V
    sigma_rmse[idx, "FEMA_true"] = mean(FEMA_true)
    sigma_rmse[idx, "glmmTMB_true"] = mean(glmmTMB_true)
    sigma_rmse[idx, "FEMA_glmmTMB"] = mean(FEMA_glmmTMB)
    
    sigma_est[[idx]] = data.frame(
      M = M,
      V = V,
      FEMA = FEMA_dat$Est_Var_IID,
      TMB = TMB_dat$Est_Var_IID
    ) 
    
    
    idx = idx+1
    
  }
}


write.table(beta_rmse,"./Simulation_1/beta_rmse.txt",sep="\t",col.names=T,quote=F,row.names=F)
write.table(sigma_rmse,"./Simulation_1/sigma_rmse.txt",sep="\t",col.names=T,quote=F,row.names=F)

beta_FEMA = bind_rows(beta_FEMA)
beta_TMB = bind_rows(beta_TMB)
sigma_est = bind_rows(sigma_est)

############################### Plot
plot_dat = cbind(beta_FEMA,beta_TMB$TMB_Est)
colnames(plot_dat)[5]="TMB_Est"
plot_dat = melt(plot_dat,id.vars=c("M","V","True"))
plot_dat = plot_dat%>%arrange(desc(variable))

x_lim = c(-0.2,0.2)
y_lim = range(plot_dat$value, na.rm = TRUE)

plot_list=list()

for (M_val in c(5,10,20)){
  for (V_val in c(1,4,7)){
    
    plot_dat_tmp = plot_dat %>% filter(M==M_val, V==V_val)
    
    p = ggplot(data=plot_dat_tmp)+
      geom_point(aes(x=True,y=value,col=variable),alpha=0.5,pch=16,size=1.5)+
      scale_color_manual(values=c("TMB_Est"="#dba4a5","FEMA_Est"="#a3b9d6"),
                         breaks=c("TMB_Est","FEMA_Est"),
                         labels=c("glmmTMB","FEMA binary"))+
      guides(col="none")+
      new_scale_color()+
      geom_line(aes(x=True,y=value,col=variable),stat="smooth",se=F,lwd=1.5,alpha=0.8)+
      scale_color_manual(values=c("TMB_Est"="#A4514F","FEMA_Est"="#0367A6"),
                         breaks=c("TMB_Est","FEMA_Est"),
                         labels=c("glmmTMB","FEMA binary"))+
      theme_classic()+
      coord_cartesian(xlim = x_lim, ylim = y_lim, clip = "off") +
      labs(x="True value",y="Estimate")+
      theme(legend.position="bottom",
            strip.background = element_blank(),
            strip.text = element_blank())
    
    if (M_val == 5) {
      p = p + ggtitle(sprintf("Var(FE) = %.1f, Var(RE) = %.1f", 1-V_val/10,V_val/10)) +
        theme(plot.title = element_text(hjust = 0.5, face = "bold", size = 10))
    }
    
    
    if (M_val != 20) {
      p = p + labs(x = NULL)
    }
    
    if (V_val != 1) {
      p = p + labs(y = NULL)
    }
    
    if (V_val == 7) {
      p = p + 
        theme(plot.margin = margin(t = 10, r = 100, b = 10, l = 10))+
        annotation_custom(
        grob = textGrob(label = sprintf("Cluster Size = %d", M_val), 
                        hjust = 0, gp = gpar(fontface = "bold", fontsize = 10)),
        xmin = x_lim[2] + (x_lim[2] - x_lim[1]) * 0.1,
        xmax = x_lim[2] + (x_lim[2] - x_lim[1]) * 0.1,
        ymin = mean(y_lim), ymax = mean(y_lim)
      )
    }
    
    
    plot_list[[paste0("M", M_val, "_V", V_val)]] = p
    
  }
}



pdf(
  file = "./Simulation_1/Fig1_beta.pdf",
  # type = "cairo", # ??????
  # res = 300, # 300ppi ?ֱ???
  width = 10, height = 7,
  bg = "white" # ͸??????
)

combined_plot = wrap_plots(plot_list, ncol = 3, nrow = 3) + 
  plot_layout(guides = "collect") & 
  theme(legend.position = "bottom",
        legend.title = element_blank())

combined_plot

dev.off()


plot_dat = melt(sigma_est,id.vars=c("M","V"))
plot_dat$M = as.factor(plot_dat$M)
plot_dat$variable = factor(plot_dat$variable,levels=c("TMB","FEMA"))
y_lim = range(plot_dat$value, na.rm = TRUE)
plot_list = list()

for (V_val in c(1,4,7)){
  
  plot_dat_tmp = plot_dat%>%filter(V==V_val)
  p = ggplot(plot_dat_tmp)+
    geom_boxplot(aes(x=M,y=value,fill=variable),linewidth=0.25,
                 outlier.size = 0.7,outlier.color = "gray40")+
    geom_hline(yintercept = V_val/10, lty="dashed",col="gray70",lwd=0.7)+
    scale_fill_manual(values=c("TMB"="#A4514F","FEMA"="#0367A6"),
                      breaks=c("TMB","FEMA"),
                      labels=c("glmmTMB","FEMA binary"))+
    # guides(fill="none")+
    coord_cartesian(ylim = y_lim, clip = "off") +
    theme_classic()+
    labs(x="Cluster Size",y=expression(hat(sigma)^2))+
    ggtitle(sprintf("Var(FE) = %.1f, Var(RE) = %.1f", 1-V_val/10,V_val/10)) +
    theme(plot.title = element_text(hjust = 0.5, face = "bold", size = 10))+
    theme(legend.position="bottom",
          legend.title = element_blank())
  
  plot_list[[paste0("V_",V_val)]] = p
  
}


pdf(
  file = "./Simulation_1/Fig1_sigma.pdf",
  # type = "cairo", # ??????
  # res = 300, # 300ppi ?ֱ???
  width = 10, height = 3,
  bg = "white" # ͸??????
)

combined_plot = wrap_plots(plot_list, ncol = 3, nrow = 1) + 
  plot_layout(guides = "collect") & 
  theme(legend.position = "bottom",
        legend.title = element_blank())

combined_plot

dev.off()


#################################################################################
# Simulation 2 - Type I error control
#################################################################################
nXvars = 100
n_snps = 100000
p_lower_ci = qbeta(0.025, 1:n_snps, n_snps - (1:n_snps) + 1)
p_upper_ci = qbeta(0.975, 1:n_snps, n_snps - (1:n_snps) + 1)

neg_log10_ci_upper = sort(-log10(p_lower_ci))
neg_log10_ci_lower = sort(-log10(p_upper_ci))

plot_list = list()

for (V in c(1,4,7)){
  
  all_qq_data = data.frame()
  
  for (M in c(5,10,20)){
    
    p_values_matrix = read.table(sprintf("./Simulation_2/N10000_M%d_V%d.txt",M,V))
    
    p_values = melt(p_values_matrix[,2:(nXvars+1)])[,2]
    qq_values = p_values
    p_values = 10^(-p_values)
    
    observed_logp =sort(qq_values)
    expected_logp = -log10(ppoints(n_snps))
    qq_data = data.frame(
      V = V,
      M = M,
      observed = observed_logp,
      expected = rev(expected_logp),
      ci_upper = neg_log10_ci_upper,
      ci_lower = neg_log10_ci_lower)
    
    all_qq_data = rbind(all_qq_data, qq_data)
    
  }
  
  all_qq_data$M = as.factor(all_qq_data$M)
  
  p = ggplot(all_qq_data) +
    geom_point(aes(x = expected, y = observed, color = M),pch=16,
      alpha = 0.5, size = 2.5) +
    geom_abline(intercept = 0, slope = 1,
      color = "#be0101", linewidth = 1) +
    geom_line(aes(x = expected, y = ci_lower),lty="dashed",col="gray20",
      alpha = 0.8, linewidth = 0.8) +
    geom_line(aes(x = expected, y = ci_upper),lty="dashed",col="gray20",
              alpha = 0.8, linewidth = 0.8) +
    scale_color_manual(values = c("5"="#C9DBED","10"="#97B8DD","20"="#0367A6")) +
    guides(col=guide_legend(title="Cluster Size",override.aes = list(size = 3)))+
    theme_classic() +
    labs(
      x = expression("Expected" ~ -log[10](italic(P))),
      y = expression("Observed" ~ -log[10](italic(P)))) +
    ggtitle(sprintf("Var(FE) = %.1f, Var(RE) = %.1f", 1-V/10,V/10)) +
    theme(plot.title = element_text(hjust = 0.5, face = "bold", size = 10))+
    theme(legend.position="right")
  
  plot_list[[paste0("V_",V)]] = p
}


# chi_squared_values <- qchisq(1 - p_values, df = 1)
# median_observed_chi_sq <- median(chi_squared_values, na.rm = TRUE)
# median_expected_chi_sq <- qchisq(0.5, df = 1)
# lambda <- median_observed_chi_sq / median_expected_chi_sq

png(
  file = "./Simulation_2/test.png",
  type = "cairo", # ??????
  res = 300, # 300ppi ?ֱ???
  width = 3500, height = 1500,
  bg = "white" # ͸??????
)

combined_plot = wrap_plots(plot_list, ncol = 3, nrow = 1) + 
  plot_layout(guides = "collect") & 
  theme(legend.position = "right")

combined_plot

dev.off()

