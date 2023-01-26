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
  
  mesh = estimates[[1]]$FEMbasis$mesh
  
  num_edges= dim(mesh$edges)[1]
  coef=matrix(0, nrow= num_edges, ncol=length(estimates) )

  for(i in 1:length(estimates)){
  for(e in 1:num_edges){
    
    coef[e,i]= (estimates[[i]]$coeff[mesh$edges[e,1]] + estimates[[i]]$coeff[mesh$edges[e,2]])/2  
    
    }
  }
  
  max.col = max(coef)
  min.col = min(coef)
  
  estimates.plot = list()
  
  for(i in 1:length(estimates)){
    
    estimates.plot[[i]] = R_plot_graph.ggplot2.2( estimates[[i]], # FEM object
                                                 line.size = line.size,
                                                  color.min = min.col,
                                                  color.max = max.col,
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

date_ = "2023-01-25-16_20_21"
folder.name = paste("data/chicago/", date_,"/",sep="")
load(paste(folder.name, "CV_error.RData",sep=""))
load(paste(folder.name, "estimates.RData",sep=""))

pdf(paste(folder.imgs,"CV_error.pdf", sep=""))
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

PLOTS <- plot_estimates(estimates, 
                        palette=plasma,
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

mesh.ref = refine.by.splitting.mesh.1.5D(DE_PDE.FEM$FEMbasis$mesh)
mesh.ref = refine.by.splitting.mesh.1.5D(mesh.ref)
mesh.ref = refine.by.splitting.mesh.1.5D(mesh.ref)
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
points(spat.stat.linnet$data$x, spat.stat.linnet$data$y, pch=16, col="red3")
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



