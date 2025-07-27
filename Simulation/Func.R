library(dplyr)
library(reshape2)
library(ggplot2)

plot_fe_diff = function(plot_dat,title){

  plot_dat = melt(plot_dat,id.vars=c("Var_FFX","Var_FID","Var_IID"))

  plot_dat = plot_dat %>%mutate(Var_FFX = paste0("Var_FFX=",Var_FFX))

  ggplot(plot_dat)+
    geom_boxplot(aes(x=factor(Var_FID),y=value,fill=variable),alpha=0.8,outlier.size=0.7,lwd=0.3,width=0.7)+
    # facet_wrap(Var_FFX~., ncol=1)+
    facet_wrap(Var_FFX~., ncol=1)+
    scale_fill_manual(values=c("#A4514F","#6C91C2","#A7C0DE","#BFB1D0","#DCD7C1"))+
    scale_y_continuous(limits=c(-0.15,0.2))+
    labs(x="Var_FID",y="Difference",title=)+
    # scale_y_continuous(sec.axis = sec_axis(~ . , name = "Var_FFX", breaks = NULL, labels = NULL))+
    theme_bw()+
    theme(
      plot.title = element_text(hjust=0.5),
      legend.title = element_blank()
    )

}

plot_re_diff = function(plot_dat,title){

  plot_dat = melt(plot_dat,id.vars=c("Var_FFX","Var_FID","Var_IID"))

  plot_dat = plot_dat %>%mutate(Var_FFX = paste0("Var_FFX=",Var_FFX))

  ggplot(plot_dat)+
    geom_boxplot(aes(x=factor(Var_FID),y=value,fill=variable),alpha=0.8,outlier.size=0.7,lwd=0.3,width=0.7)+
    # facet_wrap(Var_FFX~., ncol=1)+
    facet_wrap(Var_FFX~., ncol=1)+
    scale_fill_manual(values=c("#fcdaba", "#e7e6d4"))+
    scale_y_continuous(limits=c(-0.4,0.7))+
    labs(x="Var_FID",y="Difference",title=title)+
    # scale_y_continuous(sec.axis = sec_axis(~ . , name = "Var_FFX", breaks = NULL, labels = NULL))+
    theme_bw()+
    theme(
      plot.title = element_text(hjust=0.5),
      legend.title = element_blank()
    )

}


plot_fe_rmse = function(plot_dat,title){

  plot_dat = melt(plot_dat,id.vars=c("Var_FFX","Var_FID","Var_IID"))

  plot_dat = plot_dat %>%mutate(Var_FFX = paste0("Var_FFX=",Var_FFX))

  ggplot(plot_dat)+
    geom_bar(aes(x=factor(Var_FID),y=value,fill=variable),col="gray30",alpha=0.8,stat="identity", position="dodge")+
    # facet_wrap(Var_FFX~., ncol=1)+
    facet_wrap(Var_FFX~., ncol=1)+
    scale_fill_manual(values=c("#A4514F","#6C91C2","#A7C0DE","#BFB1D0","#DCD7C1"))+
    scale_y_continuous(limits=c(0,0.07))+
    labs(x="Var_FID",y="RMSE",title=title)+
    # scale_y_continuous(sec.axis = sec_axis(~ . , name = "Var_FFX", breaks = NULL, labels = NULL))+
    theme_bw()+
    theme(
      plot.title = element_text(hjust=0.5),
      legend.title = element_blank()
    )

}



plot_re_rmse = function(plot_dat,title){
  
  plot_dat = melt(plot_dat,id.vars=c("Var_FFX","Var_FID","Var_IID"))
  
  plot_dat = plot_dat %>%mutate(Var_FFX = paste0("Var_FFX=",Var_FFX))
  
  ggplot(plot_dat)+
    geom_bar(aes(x=factor(Var_FID),y=value,fill=variable),col="gray30",alpha=0.8,
             stat="identity", position="dodge",width=0.7)+
    # facet_wrap(Var_FFX~., ncol=1)+
    facet_wrap(Var_FFX~., ncol=1)+
    scale_fill_manual(values=c("#fcdaba", "#e7e6d4"))+
    scale_y_continuous(limits=c(0,0.35))+
    labs(x="Var_FID",y="RMSE",title=title)+
    # scale_y_continuous(sec.axis = sec_axis(~ . , name = "Var_FFX", breaks = NULL, labels = NULL))+
    theme_bw()+
    theme(
      plot.title = element_text(hjust=0.5),
      legend.title = element_blank()
    )
  
}
