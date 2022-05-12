
library(gridExtra)
library(ggplot2)
library(latex2exp)
source("C:/Users/Aldo/Documents/SimulationStudies/Auxiliary/R_plot_graph.ggplot2.R")

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
                          line.size=1.){
  
max.col = max(true.field, mean.field.fdaPDE)
min.col = min(true.field, mean.field.fdaPDE)

true.spatial.field<- R_plot_graph.ggplot2.2(FEM(true.field, FEMbasis),
                                            line.size = line.size,
                                            color.min = min.col,
                                            color.max = max.col,
                                            title = TeX("$f$", italic = T),
                                            return.ggplot.object = T, 
                                            legend.pos = "right")

true.spatial.signal<- R_plot_graph.ggplot2.2(FEM(true.signal, FEMbasis),
                                             line.size = line.size,
                                             title = TeX("$f + W\\beta$", italic = T),
                                             return.ggplot.object = T, 
                                             legend.pos = "right")

rmse <- boxplot_RMSE(RMSE, n_data, model_ = c(T,T,F,F), 
                     names_ = c("fdaPDE","GWR","",""),
                     legend.pos = legend.pos.RMSE)

mean.spatial.field.1 <- R_plot_graph.ggplot2.2(FEM(mean.field.fdaPDE[,1], FEMbasis),
                                               line.size = line.size,
                                               color.min = min.col,
                                               color.max = max.col,
                                               title = bquote(hat(f) ~ .(paste("(n=",n_data[1],")",sep=""))), #expression(hat(f) ~ paste("(n=",n_data[1],")",sep="")),
                                               
                                               return.ggplot.object = T, 
                                               legend.pos = "right")
mean.spatial.field.2 <- R_plot_graph.ggplot2.2(FEM(mean.field.fdaPDE[,2], FEMbasis),
                                               line.size = line.size,
                                               color.min = min.col,
                                               color.max = max.col,
                                               title = bquote(hat(f) ~ .(paste("(n=",n_data[2],")",sep=""))),
                                               return.ggplot.object = T, 
                                               legend.pos = "right")
mean.spatial.field.3 <- R_plot_graph.ggplot2.2(FEM(mean.field.fdaPDE[,3], FEMbasis),
                                               line.size = line.size,
                                               color.min = min.col,
                                               color.max = max.col,
                                               title = bquote(hat(f) ~ .(paste("(n=",n_data[3],")",sep=""))),
                                               return.ggplot.object = T, 
                                               legend.pos = "right")
mean.spatial.field.4 <- R_plot_graph.ggplot2.2(FEM(mean.field.fdaPDE[,4], FEMbasis),
                                               line.size = line.size,
                                               color.min = min.col,
                                               color.max = max.col,
                                               title = bquote(hat(f) ~ .(paste("(n=",n_data[4],")",sep=""))),
                                               return.ggplot.object = T, 
                                               legend.pos = "right")

nnodes = nrow(FEMbasis$mesh$nodes)
sample_ = sample(1:nnodes, n_data[1])
points_ = FEMbasis$mesh$nodes[sample_,]
mu_ =  observations[sample_]
observations.example.1 <- R_plot_mesh.ggplot(mesh = FEMbasis$mesh,
                                             points_ = points_,
                                             mu = mu_, 
                                             line.size=line.size,
                                             title = paste("Sample ",
                                                           "(n=",n_data[1],")",sep=""))

sample_ = sample(1:nnodes, n_data[2])
points_ = FEMbasis$mesh$nodes[sample_,]
mu_ =  observations[sample_]
observations.example.2 <- R_plot_mesh.ggplot(mesh = FEMbasis$mesh,
                                             points_ = points_,
                                             mu = mu_, 
                                             line.size=line.size,
                                             title = paste("Sample ",
                                                           "(n=",n_data[2],")",sep=""))

sample_ = sample(1:nnodes, n_data[3])
points_ = FEMbasis$mesh$nodes[sample_,]
mu_ =  observations[sample_]
observations.example.3 <- R_plot_mesh.ggplot(mesh = FEMbasis$mesh,
                                             points_ = points_,
                                             mu = mu_, 
                                             line.size=line.size,
                                             title = paste("Sample ",
                                                           "(n=",n_data[3],")",sep=""))

sample_ = sample(1:nnodes, n_data[4])
points_ = FEMbasis$mesh$nodes[sample_,]
mu_ =  observations[sample_]
observations.example.4 <- R_plot_mesh.ggplot(mesh = FEMbasis$mesh,
                                             points_ = points_,
                                             mu = mu_, 
                                             line.size=line.size,
                                             title = paste("Sample ",
                                                           "(n=",n_data[4],")",sep=""))

firstCov <- R_plot_graph.ggplot2.2(FEM(W[,1], FEMbasis),
                                   line.size = line.size,
                                   title = "First Covariate", #expression(hat(f) ~ paste("(n=",n_data[1],")",sep="")),
                                   return.ggplot.object = T, 
                                   legend.pos = "right")
secondCov <- R_plot_graph.ggplot2.2(FEM(W[,2], FEMbasis),
                                    line.size = line.size,
                                    title = "Second Covariate",
                                    return.ggplot.object = T, 
                                    legend.pos = "right")

nnodes = nrow(FEMbasis$mesh$nodes)
sample_ = sample(1:nnodes, n_data[1])
points_ = mesh$nodes[sample_,]
firstCov.example.1 <- R_plot_mesh.ggplot(mesh = mesh,
                                         points_ = points_, 
                                         mu = W[sample_,1], 
                                         title = paste("First Covariate",
                                                       "(n=",n_data[1],")",sep=""))

secondCov.example.1 <- R_plot_mesh.ggplot(mesh = mesh,
                                          points_ = points_, 
                                          mu = W[sample_,2],
                                          title = paste("Second Covariate",
                                                        "(n=",n_data[1],")",sep=""))

nnodes = nrow(FEMbasis$mesh$nodes)
sample_ = sample(1:nnodes, n_data[2])
points_ = mesh$nodes[sample_,]
firstCov.example.2 <- R_plot_mesh.ggplot(mesh = mesh,
                                         points_ = points_, 
                                         mu = W[sample_,1], 
                                         title = paste("First Covariate",
                                                       "(n=",n_data[2],")",sep=""))

secondCov.example.2 <- R_plot_mesh.ggplot(mesh = mesh,
                                          points_ = points_, 
                                          mu = W[sample_,2],
                                          title = paste("Second Covariate",
                                                        "(n=",n_data[2],")",sep=""))

nnodes = nrow(FEMbasis$mesh$nodes)
sample_ = sample(1:nnodes, n_data[3])
points_ = mesh$nodes[sample_,]
firstCov.example.3 <- R_plot_mesh.ggplot(mesh = mesh,
                                         points_ = points_, 
                                         mu = W[sample_,1], 
                                         title = paste("First Covariate",
                                                       "(n=",n_data[3],")",sep=""))

secondCov.example.3 <- R_plot_mesh.ggplot(mesh = mesh,
                                          points_ = points_, 
                                          mu = W[sample_,2],
                                          title = paste("Second Covariate",
                                                        "(n=",n_data[3],")",sep=""))

nnodes = nrow(FEMbasis$mesh$nodes)
sample_ = sample(1:nnodes, n_data[4])
points_ = mesh$nodes[sample_,]
firstCov.example.4 <- R_plot_mesh.ggplot(mesh = mesh,
                                         points_ = points_, 
                                         mu = W[sample_,1], 
                                         title = paste("First Covariate",
                                                       "(n=",n_data[4],")",sep=""))

secondCov.example.4 <- R_plot_mesh.ggplot(mesh = mesh,
                                          points_ = points_, 
                                          mu = W[sample_,2],
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

