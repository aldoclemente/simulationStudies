
library(gridExtra)
library(ggplot2)
library(latex2exp)
library(viridis)
source("../Auxiliary/R_plot_graph.ggplot2.R")

GLRWithCovPlots<-function(imgfile,
                          true.field,            # f 
                          true.signal,           # f + W beta 
                          mean.field.fdaPDE,
                          param,          # f + W beta + eps
                          observations,
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
  
  mesh=FEMbasis$mesh
  num_edges= dim(mesh$edges)[1]
  estimates = mean.field.fdaPDE
  coef=matrix(0, nrow= num_edges, ncol=ncol(estimates) )

  for(i in 1:ncol(estimates)){
  for(e in 1:num_edges){
    
    coef[e,i]= (estimates[mesh$edges[e,1],i] + estimates[mesh$edges[e,2],i])/2  
    
    }
  }
  
  max.col = max(coef)
  min.col = min(coef)

  max.col.true = max(true.field)
  min.col.true = min(true.field)

  max.col = max(max.col, max.col.true)
  min.col = min(min.col, min.col.true)
  
  true.spatial.field<- R_plot_graph.ggplot2.2(FEM(true.field, FEMbasis),
                                              line.size = line.size,
                                              color.min = min.col,
                                              color.max = max.col,
                                              title = bquote(f),
                                             
                                              palette=p,
                                              legend.pos = "right")
  
  true.spatial.signal<- R_plot_graph.ggplot2.2(FEM(true.signal, FEMbasis),
                                               line.size = line.size,
                                               title = TeX("$f + W\\beta$", italic = T),
                                             
                                               palette=p,
                                               legend.pos = "right")
  
  rmse <- boxplot_RMSE(RMSE, n_data, model_ = c(T,T,F,F), 
                       names_ = c("GSR-PDE","GGWR","",""), palette = palette,
                       legend.pos = legend.pos.RMSE)
  
  mean.spatial.field.1 <- R_plot_graph.ggplot2.2(FEM(mean.field.fdaPDE[,1], FEMbasis),
                                                 line.size = line.size,
                                                 color.min = min.col,
                                                 color.max = max.col,
                                                 title = bquote(hat(f) ~ .(paste("(n=",n_data[1],")",sep=""))), #expression(hat(f) ~ paste("(n=",n_data[1],")",sep="")),
                                                 
                                             
                                                 palette=p,
                                                 legend.pos = "right")
  mean.spatial.field.2 <- R_plot_graph.ggplot2.2(FEM(mean.field.fdaPDE[,2], FEMbasis),
                                                 line.size = line.size,
                                                 color.min = min.col,
                                                 color.max = max.col,
                                                 title = bquote(hat(f) ~ .(paste("(n=",n_data[2],")",sep=""))),
                                             
                                                 palette=p,
                                                 legend.pos = "right")
  mean.spatial.field.3 <- R_plot_graph.ggplot2.2(FEM(mean.field.fdaPDE[,3], FEMbasis),
                                                 line.size = line.size,
                                                 color.min = min.col,
                                                 color.max = max.col,
                                                 title = bquote(hat(f) ~ .(paste("(n=",n_data[3],")",sep=""))),
                                               
                                                 palette=p,
                                                 legend.pos = "right")
  mean.spatial.field.4 <- R_plot_graph.ggplot2.2(FEM(mean.field.fdaPDE[,4], FEMbasis),
                                                 line.size = line.size,
                                                 color.min = min.col,
                                                 color.max = max.col,
                                                 title = bquote(hat(f) ~ .(paste("(n=",n_data[4],")",sep=""))),
                                              
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
                                               title = bquote(g(mu[i]) == bold(x)[i]^T * bold(beta) + f(bold(p)[i]))) #bold(x)[i]^T * bold(beta)
  
  sample_ = sample(1:nnodes, n_data[2])
  points_ = FEMbasis$mesh$nodes[sample_,]
  mu_ =  observations[sample_]
  observations.example.2 <- R_plot_mesh.ggplot(mesh = FEMbasis$mesh,
                                               points_ = points_,
                                               mu = mu_, 
                                               line.size=line.size,
                                               palette=p,
                                               title = bquote(g(mu[i]) == bold(x)[i]^T * bold(beta) + f(bold(p)[i])))
  
  sample_ = sample(1:nnodes, n_data[3])
  points_ = FEMbasis$mesh$nodes[sample_,]
  mu_ =  observations[sample_]
  observations.example.3 <- R_plot_mesh.ggplot(mesh = FEMbasis$mesh,
                                               points_ = points_,
                                               mu = mu_, 
                                               line.size=line.size,
                                               palette=p,
                                               title = bquote(g(mu[i]) == bold(x)[i]^T * bold(beta) + f(bold(p)[i])))
  
  sample_ = sample(1:nnodes, n_data[4])
  points_ = FEMbasis$mesh$nodes[sample_,]
  mu_ =  observations[sample_]
  observations.example.4 <- R_plot_mesh.ggplot(mesh = FEMbasis$mesh,
                                               points_ = points_,
                                               mu = mu_, 
                                               line.size=line.size,
                                               palette=p,
                                               title = bquote(g(mu[i]) == bold(x)[i]^T * bold(beta) + f(bold(p)[i])))
  
  firstCov <- R_plot_graph.ggplot2.2(FEM(W[,1], FEMbasis),
                                     line.size = line.size,
                                     title = bquote(X[{1}[i]]), #"First Covariate", #expression(hat(f) ~ paste("(n=",n_data[1],")",sep="")),
                                
                                     palette=p,
                                     legend.pos = "right")
  secondCov <- R_plot_graph.ggplot2.2(FEM(W[,2], FEMbasis),
                                      line.size = line.size,
                                      title = bquote(X[{2}]), # second covariate
                                    
                                      palette=p,
                                      legend.pos = "right")
  
  nnodes = nrow(FEMbasis$mesh$nodes)
  sample_ = sample(1:nnodes, n_data[1])
  points_ = FEMbasis$mesh$nodes[sample_,]
  firstCov.example.1 <- R_plot_mesh.ggplot(mesh = FEMbasis$mesh,
                                           points_ = points_, 
                                           mu = W[sample_,1],
                                           palette=p,
                                           line.size=line.size,
                                           title = bquote(X[{1}[i]]))
  
  secondCov.example.1 <- R_plot_mesh.ggplot(mesh = FEMbasis$mesh,
                                            points_ = points_, 
                                            mu = W[sample_,2],
                                            palette = p,
                                            line.size=line.size,
                                            title = bquote(X[{2}[i]]))
  
  nnodes = nrow(FEMbasis$mesh$nodes)
  sample_ = sample(1:nnodes, n_data[2])
  points_ = FEMbasis$mesh$nodes[sample_,]
  firstCov.example.2 <- R_plot_mesh.ggplot(mesh = FEMbasis$mesh,
                                           points_ = points_, 
                                           mu = W[sample_,1], 
                                           palette=p, line.size=line.size,
                                           title = bquote(X[{1}[i]]))
  
  secondCov.example.2 <- R_plot_mesh.ggplot(mesh = FEMbasis$mesh,
                                            points_ = points_, 
                                            mu = W[sample_,2],
                                            palette=p, line.size=line.size,
                                            title = bquote(X[{2}[i]]))
  
  nnodes = nrow(FEMbasis$mesh$nodes)
  sample_ = sample(1:nnodes, n_data[3])
  points_ = FEMbasis$mesh$nodes[sample_,]
  firstCov.example.3 <- R_plot_mesh.ggplot(mesh = FEMbasis$mesh,
                                           points_ = points_, 
                                           mu = W[sample_,1], 
                                           palette=p, line.size=line.size,
                                           title = bquote(X[{1}[i]]))
  
  secondCov.example.3 <- R_plot_mesh.ggplot(mesh = FEMbasis$mesh,
                                            points_ = points_, 
                                            mu = W[sample_,2],
                                            palette=p, line.size=line.size,
                                            title = bquote(X[{2}[i]]))
  
  nnodes = nrow(FEMbasis$mesh$nodes)
  sample_ = sample(1:nnodes, n_data[4])
  points_ = FEMbasis$mesh$nodes[sample_,]
  firstCov.example.4 <- R_plot_mesh.ggplot(mesh = FEMbasis$mesh,
                                           points_ = points_, 
                                           mu = W[sample_,1], 
                                           palette=p, line.size=line.size,
                                           title = bquote(X[{1}[i]]))
  
  secondCov.example.4 <- R_plot_mesh.ggplot(mesh = FEMbasis$mesh,
                                            points_ = points_, 
                                            mu = W[sample_,2],
                                            palette=p, line.size=line.size,
                                            title = bquote(X[{2}[i]]))
  
  
  
  
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

plotting.estimates <- function(imgfile_,
                               estimates,
                               FEMbasis,
                               W,
                               line.size=0.5){
  p = viridis
  points_ = estimates$locations
  observations = estimates$observations # 
  true.signal = estimates$true.signal
  
  max.col = max(true.signal, estimates$fdaPDE)
  min.col = min(true.signal, estimates$fdaPDE)
  
  max.col = max(max.col, estimates$GWR)
  min.col = min(min.col, estimates$GWR)
  
  
  mu_ = true.signal
  estimates.true <- R_plot_mesh.ggplot(mesh = FEMbasis$mesh,
                                       points_ = points_,
                                       mu = mu_, 
                                       color.min = min.col,
                                       color.max = max.col,
                                       line.size=line.size,
                                       palette=p,
                                       title = "TRUE")
  
  mu_ = observations
  estimates.obs <- R_plot_mesh.ggplot(mesh = FEMbasis$mesh,
                                      points_ = points_,
                                      mu = mu_, 
                                      #color.min = min.col,
                                      #color.max = max.col,
                                      line.size=line.size,
                                      palette=p,
                                      title = bquote(y[i]) ) 
  
  
  mu_ = estimates$fdaPDE
  estimates.fdaPDE <- R_plot_mesh.ggplot(mesh = FEMbasis$mesh,
                                         points_ = points_,
                                         mu = mu_, 
                                         color.min = min.col,
                                         color.max = max.col,
                                         line.size=line.size,
                                         palette=p,
                                         title = "SR-PDE")
  
  mu_ = estimates$GWR
  estimates.GWR <- R_plot_mesh.ggplot(mesh = FEMbasis$mesh,
                                      points_ = points_,
                                      mu = mu_, 
                                      color.min = min.col,
                                      color.max = max.col,
                                      line.size=line.size,
                                      palette=p,
                                      title = "GWR")
  
  mu_ = estimates$X1
  estimates.X1 <- R_plot_mesh.ggplot(mesh = FEMbasis$mesh,
                                     points_ = points_,
                                     mu = mu_, 
                                     #color.min = min.col,
                                     #color.max = max.col,
                                     line.size=line.size,
                                     palette=p,
                                     title = bquote(X[{1}[i]]))
  
  mu_ = estimates$X2
  estimates.X2 <- R_plot_mesh.ggplot(mesh = FEMbasis$mesh,
                                     points_ = points_,
                                     mu = mu_, 
                                     #color.min = min.col,
                                     #color.max = max.col,
                                     line.size=line.size,
                                     palette=p,
                                     title = bquote(X[{2}[i]]))
  
  firstCov <- R_plot_graph.ggplot2.2(FEM(W[,1], FEMbasis),
                                     line.size = line.size,
                                     title = bquote(X[1]), #"First Covariate", #expression(hat(f) ~ paste("(n=",n_data[1],")",sep="")),
                                     return.ggplot.object = T,
                                     palette=p,
                                     legend.pos = "right")
  
  secondCov <- R_plot_graph.ggplot2.2(FEM(W[,2], FEMbasis),
                                      line.size = line.size,
                                      title = bquote(X[2]), # second covariate
                                      return.ggplot.object = T, 
                                      palette=p,
                                      legend.pos = "right")
  
  secondCov.i <- R_plot_graph.ggplot2.2(FEM(W[,2], FEMbasis),
                                        line.size = line.size,
                                        title = bquote( X[{2}[i]] ), # second covariate
                                        return.ggplot.object = T, 
                                        palette=p,
                                        legend.pos = "right")
  
  
  pdf(imgfile_)
  estimates.true
  estimates.obs
  estimates.fdaPDE
  estimates.GWR
  estimates.X1
  estimates.X2
  firstCov
  secondCov
  secondCov.i
  dev.off()
  
}
