#################################################
############## Density Estimation ###############
##############   Chicago Crimes   ###############
#############   Post - Processing    ############
#################################################

library(ggplot2)
library(viridis)
library(colorspace)
library(grid)
library(gridExtra)
source("../Auxiliary/R_plot_graph.ggplot2.R")

boxplot_CV_error <-function(CV_errors,
                        title.size=20,
                        begin=0.95, #color
                        end=0.25,   #color
                        width =0.75,
                        title="CV error")
{
  
  model = rep(c("DE-PDE", "KDE-PDE", "KDE-2D", "VORONOI"), each=nrow(CV_errors))
  RMSE =  as.vector(CV_errors)
  dataFrame = data.frame(RMSE=RMSE, model = model)
  
  MyTheme <- theme(
    axis.text = element_text(size=title.size-5),
    axis.title = element_text(size=title.size),
    title = element_text(size=title.size),
    plot.title = element_text(hjust = 0.5),
    legend.text = element_text(size=title.size-5),
    legend.key.size = unit(1,"cm"),
    legend.key.height = unit(1,"cm"),
    legend.title = element_blank(),
    legend.background = element_rect(fill="white", color="black",
                                     size=c(1,0.5))
  )
  
  border_col = darken(viridis(ncol(CV_errors), begin=end,end=begin), amount=0.25)
  
  p<-ggplot(dataFrame)+
    geom_boxplot(aes(x=model,
                     y=RMSE, group=model,
                     fill=model,
                     color=model), width=width)+
    scale_x_discrete(limits=c("DE-PDE", "KDE-PDE", "KDE-2D", "VORONOI"))+
    labs(x="", y="",
         title=title)+
    scale_fill_viridis(begin = end,
                       end = begin,
                       option = "viridis", discrete=T) + #ok
    scale_color_manual(values=border_col) +
    scale_y_continuous(labels = function(x) format(x, scientific = TRUE))+
    MyTheme + 
    theme(#plot.title=element_blank(),
      axis.ticks.x = element_blank(),
      legend.position = "none")
  return(p)  
}

plot_estimates <-function(estimates, # list of estimates
                          palette = NULL, #palette = NULL (ggplot default colors), "viridis", "magma"
                          line.size=1,
                          names = c(rep("",times=length(estimates))) )
  {
  
  if(!is.null(palette)){
    if(palette=="viridis")
      p=viridis
    else if(palette=="magma")
      p=magma
  }else{
    p=jet.col
  }
  
  max.col = max(estimates[[2]]$coeff, estimates[[1]]$coeff)
  min.col = min(estimates[[2]]$coeff, estimates[[1]]$coeff)
  
  max.col = max(max.col, estimates[[3]]$coeff)
  min.col = min(min.col, estimates[[3]]$coeff)
  #max.col = max(max.col, estimates[[4]]$coeff)
  #min.col = min(min.col, estimates[[4]]$coeff)
  
  estimates.plot = list()
  
  for(i in 1:length(estimates)){
    
    estimates.plot[[i]] = R_plot_graph.ggplot2.2( estimates[[i]], # FEM object
                                                 line.size = line.size,
                                                 color.min = min.col,
                                                 color.max = max.col,
                                                 title = names[[i]],
                                                 return.ggplot.object = T,
                                                 palette=p,
                                                 legend.pos = "right")
    
  }
  
  # estimates.plot[[4]] = R_plot_graph.ggplot2.2( estimates[[4]], # FEM object
  #                                               line.size = line.size,
  #                                               #color.min = min.col,
  #                                               #color.max = max.col,
  #                                               title = names[[4]],
  #                                               return.ggplot.object = T,
  #                                               palette=viridis,
  #                                               legend.pos = "right")
   
  return(estimates.plot)
  
}

if(!dir.exists("img/")) {
  dir.create("img/")
}

pdf(paste("img/DE_chicago_img_",date_,".pdf",sep=""))
boxplot(CV_errors, names = c("DE-PDE", "KDE-PDE", "KDE-2D", "VORONOI"),
        col = c("blue4", "green4", "green", "yellow"),
        main = "CV error")
boxplot_CV_error(CV_errors = CV_errors)

estimates = list()
estimates[[1]] = DE_PDE.FEM
estimates[[2]] = KDE_PDE.FEM
estimates[[3]] = KDE_2D.FEM
estimates[[4]] = KDE_VORONOI.FEM

PLOTS <- plot_estimates(estimates, 
                        palette="viridis",
                        names = c("DE-PDE", "KDE-PDE", "KDE-2D", "VORONOI"),
                        line.size = 0.75)

for(i in 1:length(estimates)){
  print(PLOTS[[i]])
}

PLOTS <- plot_estimates(estimates, 
                        palette="magma",
                        names = c("DE-PDE", "KDE-PDE", "KDE-2D", "VORONOI"),
                        line.size = 0.75)

for(i in 1:length(estimates)){
  print(PLOTS[[i]])
}

PLOTS <- plot_estimates(estimates, 
                        names = c("DE-PDE", "KDE-PDE", "KDE-2D", "VORONOI"),
                        line.size = 0.75)

for(i in 1:length(estimates)){
  print(PLOTS[[i]])
}

dev.off()
