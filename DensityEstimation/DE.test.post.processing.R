#################################################
############## Density Estimation ###############
##############       Test 1       ###############
#############   Post - Processing    ############
#################################################

library(ggplot2)
library(viridis)
library(colorspace)
library(grid)
library(gridExtra)
source("../Auxiliary/R_plot_graph.ggplot2.R")

boxplot_RMSE <- function(RMSE,
                         methods,
                         methods.names,
                         nsim,
                         title.size=20,
                         begin=0.95, #color
                         end=0.25,   #color
                         width =0.75,
                         n = c(50, 100, 150, 250),
                         title="RMSE")
{
  RMSE = RMSE[, methods]
  
  METHODS = rep(methods.names[methods], each = nrow(RMSE)) 
  N = as.character(rep(rep(n, each = nsim), times = sum(methods)))
  
  RMSE =  as.vector(RMSE)
  dataFrame = data.frame(RMSE=RMSE, METHODS = METHODS, N = N)
  
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
  
  border_col = darken(viridis(sum(methods), begin=end,end=begin), amount=0.25)
  
  p<-ggplot(dataFrame)+
    geom_boxplot(aes(x=N,
                     y=RMSE, group=interaction(METHODS,N),
                     fill=METHODS))+
    scale_x_discrete(limits=as.character(n))+
    labs(x="", y="",
         title=title)+
    scale_fill_viridis(begin = end,
                       end = begin,
                       option = "viridis", discrete=T) + #ok
    scale_color_manual(values=border_col) +
    #scale_y_continuous(labels = function(x) format(x, scientific = TRUE))+
    MyTheme + 
    theme(#plot.title=element_blank(),
      axis.ticks.x = element_blank(),
      legend.position = c(0.85,0.85) )
  return(p)  
}

plot_estimates <-function(estimates, # list of estimates
                          methods,
                          methods.names =  c(rep("",times=length(estimates))), 
                          palette = NULL, #palette = NULL (ggplot default colors), "viridis", "magma"
                          line.size=1)
{
  
  if(!is.null(palette)){
    if(palette=="viridis")
      p=viridis
    else if(palette=="magma")
      p=magma
  }else{
    p=jet.col
  }
  
  max.col = -1e8
  min.col =  1e8
  
  for(i in 1:length(methods.names)){
    max.col = max(max.col, estimates[[i]]$coeff)
    min.col = min(min.col, estimates[[i]]$coeff)
  }
  
  estimates.plot = list()
  
  for(i in 1:length(estimates)){
    if( methods[i]){
    estimates.plot[[i]] = R_plot_graph.ggplot2.2( estimates[[i]], # FEM object
                                                  line.size = line.size,
                                                  color.min = min.col,
                                                  color.max = max.col,
                                                  title = methods.names[[i]],
                                                  return.ggplot.object = T,
                                                  palette=p,
                                                  legend.pos = "right")
    }
  }
  
  return(estimates.plot)
  
}



#load("data/test/DE_test1_RMSE_2022-08-25-15_45_45.RData")

folder.imgs = paste(folder.name,"img/",sep="")
if(!dir.exists(folder.imgs)) {
  dir.create(folder.imgs)
}

RMSE = matrix(nrow = nsim * length(n), ncol=length(methods.names)) 
RMSE[,1] = as.vector(rmse.DE_PDE)
RMSE[,2] = as.vector(rmse.KDE_PDE)
RMSE[,3] = as.vector(rmse.KDE_ES)
RMSE[,4] = as.vector(rmse.KDE_2D)
RMSE[,5] = as.vector(rmse.VORONOI)


pdf(paste(folder.imgs,"RMSE",".pdf",sep=""))

boxplot_RMSE(RMSE, 
             methods = methods,
             methods.names = methods.names, 
             nsim = nsim,
             title.size=20,
             begin=0.95, #color
             end=0.25,   #color
             width =0.75,
             n = n,
             title="RMSE")
dev.off()

pdf(paste(folder.imgs,"estimates",".pdf",sep=""))
estimates = list()
estimates[[1]] = FEM(DE_PDE.FEM, FEMbasis)
estimates[[2]] = FEM(KDE_PDE.FEM, FEMbasis)
estimates[[3]] = FEM(KDE_ES.FEM, FEMbasis)
estimates[[4]] = FEM(KDE_2D.FEM, FEMbasis)
estimates[[5]] = FEM(VORONOI.FEM, FEMbasis)

PLOTS <- plot_estimates(estimates,
                        methods = methods,
                        methods.names = methods.names,
                        palette="viridis",
                        line.size = 0.75)

for(i in 1:length(estimates)){
  print(PLOTS[[i]])
}

PLOTS <- plot_estimates(estimates, 
                        methods = methods,
                        methods.names = c("DE-PDE", "KDE-PDE", "KDE-2D", "VORONOI"),
                        palette="magma",
                        line.size = 0.75)

for(i in 1:length(estimates)){
  print(PLOTS[[i]])
}

PLOTS <- plot_estimates(estimates,
                        methods = methods,
                        methods.names = c("DE-PDE", "KDE-PDE", "KDE-2D", "VORONOI"),
                        line.size = 0.75)

for(i in 1:length(estimates)){
  print(PLOTS[[i]])
}

dev.off()
