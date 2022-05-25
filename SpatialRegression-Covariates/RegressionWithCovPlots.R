library(gridExtra)
library(ggplot2)
library(latex2exp)
source("../Auxiliary/R_plot_graph.ggplot2.R")

# length(n_data) == 4
# ncol(W) == 2
RegressionWithCovPlots<-function(imgfile,
                                 true.field,            # f 
                                 true.signal,           # f + W beta 
                                 mean.field.fdaPDE,
                                 observations,          # f + W beta + eps
                                 FEMbasis,
                                 n_data,
                                 W , betas,
                                 RMSE,legend.pos.RMSE = "right",
                                 palette = NULL,
                                 line.size=1.){
  
  if(!is.null(palette)){
    if(palette=="viridis")
      p=viridis
    else if(palette=="magma")
      p=magma
  }else{
    p=jet.col
  }
  
  max.col = max(true.field, mean.field.fdaPDE)
  min.col  = min(true.field, mean.field.fdaPDE)
  
  true.spatial.field<- R_plot_graph.ggplot2.2(FEM(true.field, FEMbasis),
                                              line.size = line.size,
                                              color.min = min.col,
                                              color.max = max.col,
                                              title = bquote(f),
                                              return.ggplot.object = T,
                                              palette=p,
                                              legend.pos = "right")
  
  true.spatial.signal<- R_plot_graph.ggplot2.2(FEM(true.signal, FEMbasis),
                                              line.size = line.size,
                                              title = TeX("$f + W\\beta$", italic = T),
                                              return.ggplot.object = T,
                                              palette=p,
                                              legend.pos = "right")
  
  rmse <- boxplot_RMSE(RMSE, n_data, model_ = c(T,T,F,F), 
                       names_ = c("SR-PDE","GWR","",""), palette = palette,
                       legend.pos = legend.pos.RMSE)
  
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
                                               line.size=line.size,
                                               palette=p,
                                               title = bquote(z[i] == bold(w)[i]^T * bold(beta) + f(bold(p)[i]) + epsilon[i] )) #bold(w)[i]^T * bold(beta)
  
  sample_ = sample(1:nnodes, n_data[2])
  points_ = FEMbasis$mesh$nodes[sample_,]
  mu_ =  observations[sample_]
  observations.example.2 <- R_plot_mesh.ggplot(mesh = FEMbasis$mesh,
                                               points_ = points_,
                                               mu = mu_, 
                                               line.size=line.size,
                                               palette=p,
                                               title = bquote(z[i] == bold(w)[i]^T * bold(beta) + f(bold(p)[i]) + epsilon[i] ))
  
  sample_ = sample(1:nnodes, n_data[3])
  points_ = FEMbasis$mesh$nodes[sample_,]
  mu_ =  observations[sample_]
  observations.example.3 <- R_plot_mesh.ggplot(mesh = FEMbasis$mesh,
                                               points_ = points_,
                                               mu = mu_, 
                                               line.size=line.size,
                                               palette=p,
                                               title = bquote(z[i] == bold(w)[i]^T * bold(beta) + f(bold(p)[i]) + epsilon[i] ))
  
  sample_ = sample(1:nnodes, n_data[4])
  points_ = FEMbasis$mesh$nodes[sample_,]
  mu_ =  observations[sample_]
  observations.example.4 <- R_plot_mesh.ggplot(mesh = FEMbasis$mesh,
                                               points_ = points_,
                                               mu = mu_, 
                                               line.size=line.size,
                                               palette=p,
                                               title = bquote(z[i] == bold(w)[i]^T * bold(beta) + f(bold(p)[i]) + epsilon[i] ))
  
  firstCov <- R_plot_graph.ggplot2.2(FEM(W[,1], FEMbasis),
                                     line.size = line.size,
                                     title = bquote(N(0.5, 0.25^2)), #"First Covariate", #expression(hat(f) ~ paste("(n=",n_data[1],")",sep="")),
                                     return.ggplot.object = T,
                                     palette=p,
                                     legend.pos = "right")
  secondCov <- R_plot_graph.ggplot2.2(FEM(W[,2], FEMbasis),
                                                 line.size = line.size,
                                                 title = "sinusoidal function", # second covariate
                                                 return.ggplot.object = T, 
                                                 palette=p,
                                                 legend.pos = "right")
  
  nnodes = nrow(FEMbasis$mesh$nodes)
  sample_ = sample(1:nnodes, n_data[1])
  points_ = FEMbasis$mesh$nodes[sample_,]
  firstCov.example.1 <- R_plot_mesh.ggplot(mesh = FEMbasis$mesh,
                                 points_ = points_, 
                                 mu = W[sample_,1],
                                 palette=p,
                                 title = paste("First Covariate",
                                               "(n=",n_data[1],")",sep=""))
  
  secondCov.example.1 <- R_plot_mesh.ggplot(mesh = FEMbasis$mesh,
                                  points_ = points_, 
                                  mu = W[sample_,2],
                                  palette = p,
                                  title = paste("Second Covariate",
                                                "(n=",n_data[1],")",sep=""))
  
  nnodes = nrow(FEMbasis$mesh$nodes)
  sample_ = sample(1:nnodes, n_data[2])
  points_ = FEMbasis$mesh$nodes[sample_,]
  firstCov.example.2 <- R_plot_mesh.ggplot(mesh = FEMbasis$mesh,
                                           points_ = points_, 
                                           mu = W[sample_,1], 
                                           palette=p,
                                           title = paste("First Covariate",
                                                         "(n=",n_data[2],")",sep=""))
  
  secondCov.example.2 <- R_plot_mesh.ggplot(mesh = FEMbasis$mesh,
                                            points_ = points_, 
                                            mu = W[sample_,2],
                                            palette=p,
                                            title = paste("Second Covariate",
                                                          "(n=",n_data[2],")",sep=""))
  
  nnodes = nrow(FEMbasis$mesh$nodes)
  sample_ = sample(1:nnodes, n_data[3])
  points_ = FEMbasis$mesh$nodes[sample_,]
  firstCov.example.3 <- R_plot_mesh.ggplot(mesh = FEMbasis$mesh,
                                           points_ = points_, 
                                           mu = W[sample_,1], 
                                           palette=p,
                                           title = paste("First Covariate",
                                                         "(n=",n_data[3],")",sep=""))
  
  secondCov.example.3 <- R_plot_mesh.ggplot(mesh = FEMbasis$mesh,
                                            points_ = points_, 
                                            mu = W[sample_,2],
                                            palette=p,
                                            title = paste("Second Covariate",
                                                          "(n=",n_data[3],")",sep=""))
  
  nnodes = nrow(FEMbasis$mesh$nodes)
  sample_ = sample(1:nnodes, n_data[4])
  points_ = FEMbasis$mesh$nodes[sample_,]
  firstCov.example.4 <- R_plot_mesh.ggplot(mesh = FEMbasis$mesh,
                                           points_ = points_, 
                                           mu = W[sample_,1], 
                                           palette=p,
                                           title = paste("First Covariate",
                                                         "(n=",n_data[4],")",sep=""))
  
  secondCov.example.4 <- R_plot_mesh.ggplot(mesh = FEMbasis$mesh,
                                            points_ = points_, 
                                            mu = W[sample_,2],
                                            palette=p,
                                            title = paste("Second Covariate",
                                                          "(n=",n_data[4],")",sep=""))
  
  
  
  
  pdf(imgfile)
  print(true.spatial.field)
  print(true.spatial.signal)
  print(rmse)
  print(mean.spatial.field.1)
  print(mean.spatial.field.2)
  print(mean.spatial.field.3)
  print(mean.spatial.field.4)
  print(observations.example.1)
  print(observations.example.2)
  print(observations.example.3)
  print(observations.example.4)
  print(firstCov)
  print(secondCov)
  print(firstCov.example.1)
  print(firstCov.example.2)
  print(firstCov.example.3)
  print(firstCov.example.4)
  print(secondCov.example.1)
  print(secondCov.example.2)
  print(secondCov.example.3)
  print(secondCov.example.4)
  dev.off()
}
