library(viridis)
library(gridExtra)
library(ggplot2)
library(latex2exp)
source("C:/Users/Aldo/Desktop/simulationStudies/Auxiliary/R_plot_graph.ggplot2.R")

# length(n_data) == 4
# ncol(mean.field.fdaPDE) == 4
RegressionNoCovPlots<-function(imgfile,
                               true.field,
                               mean.field.fdaPDE,
                               observations,
                               FEMbasis = FEMbasis,
                               n_data = n_data,
                               RMSE,legend.pos.RMSE = "right",
                               palette = NULL, #palette = NULL (ggplot default colors), "viridis", "magma"
                               line.size=1.){
  
  if(!is.null(palette)){
    if(palette=="viridis")
      p=viridis
    else if(palette=="magma")
      p=magma
  }else{
    p=jet.col
  }
  max.col = max(true.field)
  min.col = min(true.field)
  
  true.spatial.field<- R_plot_graph.ggplot2.2(FEM(true.field, FEMbasis),
                         line.size = line.size,
                         color.min = min.col,
                         color.max = max.col,
                         title = bquote(f),
                         return.ggplot.object = T,
                         palette=p,
                         legend.pos = "right")
  
  rmse <- boxplot_RMSE(RMSE, n_data, model_ = c(T,T,T,T), 
                       names_ = c("SR-PDE","GWR","Lattice","RR-Krig"),
                       legend.pos = legend.pos.RMSE, palette=palette)

  mean.spatial.field.1 <- R_plot_graph.ggplot2.2(FEM(mean.field.fdaPDE[,1], FEMbasis),
                                           line.size = line.size,
                                           color.min = min.col,
                                           color.max = max.col,
                                           title = bquote(hat(f) ~ .(paste("(n=",n_data[1],")",sep=""))), #expression(hat(f) ~ paste("(n=",n_data[1],")",sep="")),
                                           
                                           return.ggplot.object = T, 
                                           palette=p,
                                           legend.pos = "right")
  mean.spatial.field.2 <- R_plot_graph.ggplot2.2(FEM(mean.field.fdaPDE[,2], FEMbasis),
                                               line.size = line.size,
                                               color.min = min.col,
                                               color.max = max.col,
                                               title = bquote(hat(f) ~ .(paste("(n=",n_data[2],")",sep=""))),
                                               return.ggplot.object = T, 
                                               palette=p,
                                               legend.pos = "right")
  mean.spatial.field.3 <- R_plot_graph.ggplot2.2(FEM(mean.field.fdaPDE[,3], FEMbasis),
                                               line.size = line.size,
                                               color.min = min.col,
                                               color.max = max.col,
                                               title = bquote(hat(f) ~ .(paste("(n=",n_data[3],")",sep=""))),
                                               return.ggplot.object = T,
                                               palette=p,
                                               legend.pos = "right")
  mean.spatial.field.4 <- R_plot_graph.ggplot2.2(FEM(mean.field.fdaPDE[,4], FEMbasis),
                                               line.size = line.size,
                                               color.min = min.col,
                                               color.max = max.col,
                                               title = bquote(hat(f) ~ .(paste("(n=",n_data[4],")",sep=""))),
                                               return.ggplot.object = T, 
                                               palette=p,
                                               legend.pos = "right")

  nnodes = nrow(FEMbasis$mesh$nodes)
  sample_ = sample(1:nnodes, n_data[1])
  points_ = FEMbasis$mesh$nodes[sample_,]
  mu_ =  observations[sample_]
  observations.example.1 <- R_plot_mesh.ggplot(mesh = FEMbasis$mesh,
                                               points_ = points_,
                                               mu = mu_, 
                                               color.min = min.col,
                                               color.max = max.col,
                                               line.size=line.size,
                                               palette=p,
                                               title = bquote(z[i] == f(bold(p)[i]) + epsilon[i] )) #bold(w)[i]^T * bold(beta)
                                                 
                                                 # paste("Sample ",
                                                      #      "(n=",n_data[1],")",sep=""))
  
  sample_ = sample(1:nnodes, n_data[2])
  points_ = FEMbasis$mesh$nodes[sample_,]
  mu_ =  observations[sample_]
  observations.example.2 <- R_plot_mesh.ggplot(mesh = FEMbasis$mesh,
                                               points_ = points_,
                                               mu = mu_, 
                                               color.min = min.col,
                                               color.max = max.col,
                                               line.size=line.size,
                                               palette=p,
                                               title = bquote(z[i] == f(bold(p)[i]) + epsilon[i]) ) 
                                                 #paste("Sample ",
                                                      #      "(n=",n_data[2],")",sep=""))
  
  sample_ = sample(1:nnodes, n_data[3])
  points_ = FEMbasis$mesh$nodes[sample_,]
  mu_ =  observations[sample_]
  observations.example.3 <- R_plot_mesh.ggplot(mesh = FEMbasis$mesh,
                                               points_ = points_,
                                               mu = mu_, 
                                               color.min = min.col,
                                               color.max = max.col,
                                               line.size=line.size,
                                               palette=p,
                                               title = bquote(z[i] == f(bold(p)[i]) + epsilon[i]) )
                                                 # paste("Sample ",
                                                      #      "(n=",n_data[3],")",sep=""))
  
  sample_ = sample(1:nnodes, n_data[4])
  points_ = FEMbasis$mesh$nodes[sample_,]
  mu_ =  observations[sample_]
  observations.example.4 <- R_plot_mesh.ggplot(mesh = FEMbasis$mesh,
                                               points_ = points_,
                                               mu = mu_, 
                                               color.min = min.col,
                                               color.max = max.col,
                                               line.size=line.size,
                                               palette=p,
                                               title = bquote(z[i] == f(bold(p)[i]) + epsilon[i]) ) 
                                                 #paste("Sample ",
                                                      #      "(n=",n_data[4],")",sep=""))
  
  
  pdf(imgfile)
  print(true.spatial.field)
  print(rmse)
  print(mean.spatial.field.1)
  print(mean.spatial.field.2)
  print(mean.spatial.field.3)
  print(mean.spatial.field.4)
  print(observations.example.1)
  print(observations.example.2)
  print(observations.example.3)
  print(observations.example.4)
  dev.off()
}
