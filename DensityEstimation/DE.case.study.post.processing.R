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
                            methods,
                            methods.names,
                        title.size=20,
                        begin=0.95, #color
                        end=0.25,   #color
                        width =0.75,
                        title="CV error")
{
  
  METHODS = rep(methods.names[methods], each=nrow(CV_errors))
  RMSE =  as.vector(CV_errors)
  dataFrame = data.frame(RMSE=RMSE, METHODS = METHODS)
  
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
  fill_col = viridis(ncol(CV_errors), begin=end, end=begin)
  
  BORDER = c()
  FILL = c()
  for(i in 1:length(methods)){
    if(methods[i]){ 
      FILL = append(FILL, fill_col[i])
      BORDER = append(BORDER, border_col[i])
    }
  }
  
  dataFrame$METHODS = factor(dataFrame$METHODS, 
                             levels=methods.names) 
  
  
  p<-ggplot(dataFrame)+
    geom_boxplot(aes(x=METHODS,
                     y=RMSE, group=METHODS,
                     fill=METHODS,
                     color=METHODS), width=width)+
    scale_x_discrete(limits=methods.names[methods])+
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

# estimates, lista di oggetti FEM ( fem = FEM(coeff, FEMbasis))
# palette, color palette (es. jet.col, viridis, magma, inferno, plasma)
# line.size, dimensione segmenti network
# titles, vettore di titoli dei plot 
plot_estimates <-function(estimates, # list of estimates
                          palette = jet.col, #palette = NULL (ggplot default colors), "viridis", "magma"
                          line.size=1,
                          titles = c(rep("",times=length(estimates))) )
  {
  
  # max.col = -1e8
  # min.col = 1e8
  # for(i in 1:length(estimates)){
  #   max.col = max(estimates[[i]]$coeff, max.col)
  #   min.col = min(estimates[[i]]$coeff, min.col)
  # }
  
  estimates.plot = list()
  
  for(i in 1:length(estimates)){
    
    estimates.plot[[i]] = R_plot_graph.ggplot2.2( estimates[[i]], # FEM object
                                                 line.size = line.size,
                                                #  color.min = min.col,
                                                #  color.max = max.col,
                                                 title = titles[[i]],
                                                 palette=palette,
                                                 legend.pos = "right")
    
  }
  
  return(estimates.plot)
  
}

folder.imgs = paste(folder.name,"img/",sep="")
if(!dir.exists(folder.imgs)) {
  dir.create(folder.imgs)
}

pdf(paste(folder.imgs,"CV_error.pdf",sep=""))
methods = c(T,T,F,T,T)
methods.names = c("DE-PDE", "KDE-PDE", "KDE-ES", "KDE-2D", "VORONOI")
boxplot_CV_error(CV_errors = CV_errors, 
                 methods = methods,
                 methods.names = methods.names)
dev.off()

pdf(paste(folder.imgs, "estimates.pdf",sep=""))
estimates = list()
estimates[[1]] = DE_PDE.FEM
estimates[[2]] = KDE_PDE.FEM
estimates[[3]] = KDE_2D.FEM
estimates[[4]] = KDE_VORONOI.FEM

PLOTS <- plot_estimates(estimates, 
                        palette=viridis,
                        titles = methods.names[methods],
                        line.size = 0.75)

for(i in 1:length(estimates)){
  print(PLOTS[[i]])
}

PLOTS <- plot_estimates(estimates, 
                        palette=magma,
                        titles = methods.names[methods],
                        line.size = 0.75)

for(i in 1:length(estimates)){
  print(PLOTS[[i]])
}

PLOTS <- plot_estimates(estimates, 
                        titles = methods.names[methods],
                        line.size = 0.75)

for(i in 1:length(estimates)){
  print(PLOTS[[i]])
}
dev.off()

# ref mesh 
estimates = list()
estimates[[1]] = DE_PDE.FEM
estimates[[2]] = KDE_PDE.FEM
estimates[[3]] = KDE_2D.FEM
estimates[[4]] = KDE_VORONOI.FEM

mesh.ref = refine.mesh.1.5D(DE_PDE.FEM$FEMbasis$mesh, delta = 10)
FEMbasis.ref = create.FEM.basis(mesh.ref)
locs = mesh.ref$nodes

for(i in 1:length(estimates))
  estimates[[i]] = FEM(eval.FEM(estimates[[i]], locations = locs),
                       FEMbasis.ref)


#pdf(paste(folder.imgs,"estimates_ref.pdf",sep=""))
pdf(paste(folder.imgs,"estimates_ref.pdf",sep=""))
PLOTS <- plot_estimates(estimates, 
                        palette = viridis,
                        titles = methods.names[methods],
                        line.size = 0.75)

for(i in 1:length(estimates)){
  print(PLOTS[[i]])
}

PLOTS <- plot_estimates(estimates, 
                        palette = inferno,
                        titles = methods.names[methods],
                        line.size = 0.75)

for(i in 1:length(estimates)){
  print(PLOTS[[i]])
}

PLOTS <- plot_estimates(estimates, 
                        palette=magma,
                        titles = methods.names[methods],
                        line.size = 0.75)

for(i in 1:length(estimates)){
  print(PLOTS[[i]])
}

PLOTS <- plot_estimates(estimates, 
                        palette=plasma,
                        titles = methods.names[methods],
                        line.size = 0.75)

for(i in 1:length(estimates)){
  print(PLOTS[[i]])
}
dev.off()

pdf(paste(folder.imgs, "point_pattern.pdf",sep=""))
plot(mesh, pch=".")
points(chicago$data$x, chicago$data$y, pch=16, col="red3")
dev.off()

pdf("estimates_plasma.pdf")
estimates = list()
estimates[[1]] = DE_PDE.FEM
estimates[[2]] = KDE_PDE.FEM
estimates[[3]] = KDE_2D.FEM
estimates[[4]] = KDE_VORONOI.FEM

PLOTS <- plot_estimates(estimates, 
                        palette=plasma,
                        titles = methods.names[methods],
                        line.size = 0.75)

for(i in 1:length(estimates)){
  print(PLOTS[[i]])
}
dev.off()



