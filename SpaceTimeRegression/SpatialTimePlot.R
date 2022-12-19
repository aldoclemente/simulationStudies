library(grid)
library(gridExtra)
library(ggplot2)
library(latex2exp)
library(viridis)
source("../Auxiliary/R_plot_graph.ggplot2.R")
source("../SpatialRegression-NoCovariate/NoCovariatesCore.R")

SpaceTimePlots <- function(imgfile, 
                           time_locations,
                           true.field,            # f 
                           true.signal,           # f + W beta 
                           mean.field.fdaPDE,
                           observations,          # f + W beta + eps
                           FEMbasis,
                           n_data,
                           W , betas,
                           RMSE,legend.pos.RMSE = "right",
                           palette = NULL,
                           line.size=1){
  
  if(!is.null(palette)){
    if(palette=="viridis")
      p=viridis
    else if(palette=="magma")
      p=magma
  }else{
    p=jet.col
  }

  col.min = min(field, mean.field.fdaPDE)
  col.max = max(field, mean.field.fdaPDE)
  nnodes = nrow(FEMbasis$mesh$nodes)
  mesh=FEMbasis$mesh
  n.locs.time = length(time_locations)
  coef.ex = matrix(true.field, 
                   nrow=nnodes,
                   ncol=n.locs.time)
  
  coef.signal = matrix(true.signal, 
                       nrow=nnodes,
                       ncol=n.locs.time)
  
  coef.mean.1 = matrix(mean.field.fdaPDE[,1], 
                       nrow=nnodes,
                       ncol=n.locs.time)
  coef.mean.2 = matrix(mean.field.fdaPDE[,2], 
                       nrow=nnodes,
                       ncol=n.locs.time)
  coef.mean.3 = matrix(mean.field.fdaPDE[,3], 
                       nrow=nnodes,
                       ncol=n.locs.time)
  coef.mean.4 = matrix(mean.field.fdaPDE[,4], 
                       nrow=nnodes,
                       ncol=n.locs.time)
  
  
  true.spatial.signal  = list()
  true.spatial.field   = list()
  mean.field.fdaPDE.1  = list()
  mean.field.fdaPDE.2  = list()
  mean.field.fdaPDE.3  = list()
  mean.field.fdaPDE.4  = list()
  
  for( i in 1:n.locs.time){
    title_ = TeX(sprintf("$t = %f", round(time_locations[i],digits = 2)))
    
    true.spatial.field[[i]] <- R_plot_graph.ggplot2.2(FEM(coef.ex[,i], FEMbasis),
                                               line.size = line.size,
                                               color.max = col.max,
                                               color.min = col.min,
                                               title = title_,
                                               return.ggplot.object = T,
                                               legend.pos = "right",
                                               palette = p,
                                               title.size = 16) 
    
    true.spatial.signal[[i]] <- R_plot_graph.ggplot2.2(FEM(coef.signal[,i], FEMbasis),
                                                      line.size = line.size,
                                                      title = title_,
                                                      color.max = max(coef.signal),
                                                      color.min = min(coef.signal),
                                                      return.ggplot.object = T,
                                                      legend.pos = "right",
                                                      palette = p,
                                                      title.size = 16)
    
    mean.field.fdaPDE.1[[i]] <- R_plot_graph.ggplot2.2(FEM(coef.mean.1[,i], FEMbasis),
                                                       line.size = line.size,
                                                       color.max = col.max,
                                                       color.min = col.min,
                                                       title = title_,
                                                       return.ggplot.object = T,
                                                       legend.pos = "right",
                                                       palette = p,
                                                       title.size = 16)
    
    mean.field.fdaPDE.2[[i]] <- R_plot_graph.ggplot2.2(FEM(coef.mean.2[,i], FEMbasis),
                                                       line.size = line.size,
                                                       color.max = col.max,
                                                       color.min = col.min,
                                                       title = title_,
                                                       return.ggplot.object = T,
                                                       legend.pos = "right",
                                                       palette = p,
                                                       title.size = 16)
    
    mean.field.fdaPDE.3[[i]] <- R_plot_graph.ggplot2.2(FEM(coef.mean.3[,i], FEMbasis),
                                                       line.size = line.size,
                                                       color.max = col.max,
                                                       color.min = col.min,
                                                       title = title_,
                                                       return.ggplot.object = T,
                                                       legend.pos = "right",
                                                       palette = p,
                                                       title.size = 16)
    
    mean.field.fdaPDE.4[[i]] <- R_plot_graph.ggplot2.2(FEM(coef.mean.4[,i], FEMbasis),
                                                       line.size = line.size,
                                                       color.max = col.max,
                                                       color.min = col.min,
                                                       title = title_,
                                                       return.ggplot.object = T,
                                                       legend.pos = "right",
                                                       palette = p,
                                                       title.size = 16)
  }  
  
rmse <- boxplot_RMSE(RMSE, n_data, model_ = c(T,T,F,F), 
             names_ = c("ST-PDE","GTWR","",""),
             palette = palette,
             legend.pos = legend.pos.RMSE)

### obs ###

observations.example.1 = list()
observations.example.2 = list()
observations.example.3 = list()
observations.example.4 = list()

tmp.sample_ = sample(x=(1:nnodes), size=n_data[1])
sample_ = vector(mode="integer")
for(k in 1:n.locs.time){
    sample_ = cbind(sample_, tmp.sample_ + (k-1)*nnodes)
}
sample_ = as.vector(sample_)
  
points_ = FEMbasis$mesh$nodes[tmp.sample_,]
mu_ =  observations[sample_]
mu_.coeff = matrix(mu_, nrow=n_data[1], ncol=n.locs.time)
for( i in 1:n.locs.time){
observations.example.1[[i]] = R_plot_mesh.ggplot(mesh = FEMbasis$mesh,
                                                    points_ = points_,
                                                    mu = mu_, 
                                                    line.size=line.size,
                                                    title = TeX(sprintf("$t = %f", round(time_locations[i],digits = 2))),
                                                 title.size = 16,
                                                 palette = p,
                                                 points.size = 0.25)
}

tmp.sample_ = sample(x=(1:nnodes), size=n_data[2])
sample_ = vector(mode="integer")
for(k in 1:n.locs.time){
  sample_ = cbind(sample_, tmp.sample_ + (k-1)*nnodes)
}
sample_ = as.vector(sample_)

points_ = FEMbasis$mesh$nodes[tmp.sample_,]
mu_ =  observations[sample_]
mu_.coeff = matrix(mu_, nrow=n_data[2], ncol=n.locs.time)
for( i in 1:n.locs.time){
  observations.example.2[[i]] = R_plot_mesh.ggplot(mesh = FEMbasis$mesh,
                                                    points_ = points_,
                                                    mu = mu_, 
                                                    line.size=line.size,
                                                    title = TeX(sprintf("$t = %f", round(time_locations[i],digits = 2))),
                                                   title.size = 16,
                                                   palette = p,
                                                   points.size = 0.25)
}

tmp.sample_ = sample(x=(1:nnodes), size=n_data[3])
sample_ = vector(mode="integer")
for(k in 1:n.locs.time){
  sample_ = cbind(sample_, tmp.sample_ + (k-1)*nnodes)
}
sample_ = as.vector(sample_)

points_ = FEMbasis$mesh$nodes[tmp.sample_,]
mu_ =  observations[sample_]
mu_.coeff = matrix(mu_, nrow=n_data[3], ncol=n.locs.time)
for( i in 1:n.locs.time){
  observations.example.3[[i]] = R_plot_mesh.ggplot(mesh = FEMbasis$mesh,
                                                    points_ = points_,
                                                    mu = mu_, 
                                                    line.size=line.size,
                                                    title = TeX(sprintf("$t = %f", round(time_locations[i],digits = 2))),
                                                   title.size = 16,
                                                   palette = p,
                                                   points.size = 0.25)
}

tmp.sample_ = sample(x=(1:nnodes), size=n_data[4])
sample_ = vector(mode="integer")
for(k in 1:n.locs.time){
  sample_ = cbind(sample_, tmp.sample_ + (k-1)*nnodes)
}
sample_ = as.vector(sample_)

points_ = FEMbasis$mesh$nodes[tmp.sample_,]
mu_ =  observations[sample_]
mu_.coeff = matrix(mu_, nrow=n_data[4], ncol=n.locs.time)
for( i in 1:n.locs.time){
  observations.example.4[[i]] = R_plot_mesh.ggplot(mesh = FEMbasis$mesh,
                                                    points_ = points_,
                                                    mu = mu_, 
                                                    line.size=line.size,
                                                    title = TeX(sprintf("$t = %f", round(time_locations[i],digits = 2))),
                                                   title.size = 16,
                                                   palette = p,
                                                   points.size = 0.25)
}

###########
firstCov <- R_plot_graph.ggplot2.2(FEM(W[1:nnodes,1], FEMbasis),
                                   line.size = line.size,
                                   title = bquote(N(0,0.25^2)), #expression(hat(f) ~ paste("(n=",n_data[1],")",sep="")),
                                   return.ggplot.object = T, 
                                   legend.pos = "right")
secondCov <- R_plot_graph.ggplot2.2(FEM(W[1:nnodes,2], FEMbasis),
                                    line.size = line.size,
                                    title = "sinusoidal function",
                                    palette = p,
                                    return.ggplot.object = T, 
                                    legend.pos = "right")

sample_ = sample(1:nnodes, n_data[1])
points_ = mesh$nodes[sample_,]
firstCov.example.1 <- R_plot_mesh.ggplot(mesh = mesh,
                                         points_ = points_, 
                                         mu = W[sample_,1], 
                                         palette = p,
                                         title = paste("First Covariate",
                                                       "(n=",n_data[1],")",sep=""))

secondCov.example.1 <- R_plot_mesh.ggplot(mesh = mesh,
                                          points_ = points_, 
                                          mu = W[sample_,2],
                                          palette = p,
                                          title = paste("Second Covariate",
                                                        "(n=",n_data[1],")",sep=""))

sample_ = sample(1:nnodes, n_data[2])
points_ = mesh$nodes[sample_,]
firstCov.example.2 <- R_plot_mesh.ggplot(mesh = mesh,
                                         points_ = points_, 
                                         mu = W[sample_,1], 
                                         palette = p,
                                         title = paste("First Covariate",
                                                       "(n=",n_data[2],")",sep=""))

secondCov.example.2 <- R_plot_mesh.ggplot(mesh = mesh,
                                          points_ = points_, 
                                          mu = W[sample_,2],
                                          palette = p,
                                          title = paste("Second Covariate",
                                                        "(n=",n_data[2],")",sep=""))

sample_ = sample(1:nnodes, n_data[3])
points_ = mesh$nodes[sample_,]
firstCov.example.3 <- R_plot_mesh.ggplot(mesh = mesh,
                                         points_ = points_, 
                                         mu = W[sample_,1], 
                                         palette = p,
                                         title = paste("First Covariate",
                                                       "(n=",n_data[3],")",sep=""))

secondCov.example.3 <- R_plot_mesh.ggplot(mesh = mesh,
                                          points_ = points_, 
                                          mu = W[sample_,2],
                                          palette = p,
                                          title = paste("Second Covariate",
                                                        "(n=",n_data[3],")",sep=""))

sample_ = sample(1:nnodes, n_data[4])
points_ = mesh$nodes[sample_,]
firstCov.example.4 <- R_plot_mesh.ggplot(mesh = mesh,
                                         points_ = points_, 
                                         mu = W[sample_,1], 
                                         palette = p,
                                         title = paste("First Covariate",
                                                       "(n=",n_data[4],")",sep=""))

secondCov.example.4 <- R_plot_mesh.ggplot(mesh = mesh,
                                          points_ = points_, 
                                          mu = W[sample_,2],
                                          palette = p,
                                          title = paste("Second Covariate",
                                                        "(n=",n_data[4],")",sep=""))


pdf(imgfile) # h = 9 l = 9 ?
grid.arrange(grobs = true.spatial.field, 
             layout_matrix = rbind(c(1,2,3),
                                   c(4,5,6)),top=grid.text(expression(f)) )

grid.arrange(grobs = true.spatial.signal, 
             layout_matrix = rbind(c(1,2,3),
                                   c(4,5,6)),top=grid.text(TeX("$f + W\\beta$", italic = T)) )

print(rmse)

grid.arrange(grobs = mean.field.fdaPDE.1, 
                        layout_matrix = rbind(c(1,2,3),
                                              c(4,5,6)), 
             grid.text(bquote(hat(f) ~ .(paste("(n=",n_data[1],")",sep="")))) )

grid.arrange(grobs = mean.field.fdaPDE.2, 
             layout_matrix = rbind(c(1,2,3),
                                   c(4,5,6)), top=bquote(hat(f) ~ .(paste("(n=",n_data[3],")",sep=""))) )

grid.arrange(grobs = mean.field.fdaPDE.3, 
             layout_matrix = rbind(c(1,2,3),
                                   c(4,5,6)), bquote(hat(f) ~ .(paste("(n=",n_data[4],")",sep=""))) )

grid.arrange(grobs = mean.field.fdaPDE.4, 
             layout_matrix = rbind(c(1,2,3),
                                   c(4,5,6)), bquote(hat(f) ~ .(paste("(n=",n_data[4],")",sep=""))) )

grid.arrange(grobs = observations.example.1,
               layout_matrix = rbind(c(1,2,3),
                                     c(4,5,6)), top=paste("Sample ", "(n=",n_data[1],")",sep=""))
grid.arrange(grobs = observations.example.2,
             layout_matrix = rbind(c(1,2,3),
                                   c(4,5,6)), top=paste("Sample ", "(n=",n_data[2],")",sep=""))
grid.arrange(grobs = observations.example.3,
             layout_matrix = rbind(c(1,2,3),
                                   c(4,5,6)), top=paste("Sample ", "(n=",n_data[3],")",sep=""))
grid.arrange(grobs = observations.example.4,
             layout_matrix = rbind(c(1,2,3),
                                   c(4,5,6)), top=paste("Sample ", "(n=",n_data[4],")",sep=""))


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

for( i in 1:n.locs.time){
  print(true.spatial.field[[i]])
}

for( i in 1:n.locs.time){
  print(mean.field.fdaPDE.4[[i]])
}

dev.off()



}
