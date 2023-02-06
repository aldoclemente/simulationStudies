rm(list=ls())
library(fdaPDE)
library(spatstat)
setwd("SR-chicago-case-study/")

source("refine1D.R")
source("../utils.R")
source("../DensityEstimation/utils.R")
data("chicago")

vertices = cbind(chicago$domain$vertices$x, chicago$domain$vertices$y)
edges = cbind( chicago$domain$from, chicago$domain$to)

mesh = create.mesh.1.5D(nodes = vertices, edges = edges)

#mesh = normalize_mesh(mesh)
tmp = normalize_mesh_unit(mesh)
mesh = tmp$mesh

LN = as.spatstat.linnet.fdaPDE(mesh)

# x.m = mean(chicago$domain$vertices$x)
# y.m = mean(chicago$domain$vertices$y)
# 
# x.sd = sd(chicago$domain$vertices$x)
# y.sd = sd(chicago$domain$vertices$y)

# x.norm = (chicago$data$x - x.m)/x.sd
# y.norm = (chicago$data$y - y.m)/y.sd

x.min = tmp$x.min
x.max = tmp$x.max

y.min = tmp$y.min
y.max = tmp$y.max

x.norm = (chicago$data$x - x.min)/(x.max-x.min)
y.norm = (chicago$data$y - y.min)/(y.max-y.min)

chicago.norm = spatstat.linnet::lpp(X= ppp(x.norm, y=y.norm, 
                                        window= LN$window), 
                                 L= LN)
#chicago.norm = list(data= data.norm, domain= LN )

DATA = chicago.norm # chicago

{
x11()
plot(mesh, pch=".")
points(DATA$data$x, DATA$data$y, col = "red", pch=16, cex=1.5)
}
### FAMILY - poisson ###
delta = 0.045

mesh = refine.mesh.1.5D(mesh, delta=delta)
FEMbasis = create.FEM.basis(mesh)
#new_to_old  = refine1D(mesh$nodes, mesh$edges, delta)$new_to_old

#setting regions 
nnodes = nrow(mesh$nodes)
#nregion = 30 #8 #6 #10 
nedges = nrow(mesh$edges)
ndata = length(DATA$data$x)
data_ = cbind(DATA$data$x, DATA$data$y)

#PP = runiflpp(n=nregion, L = chicago.norm$domain) 
#centroids_ = cbind(PP$data$x, PP$data$y)
#centroids_ = projection.points.1.5D(mesh=mesh, location=centroids_)
# idxs = c(2, 5, 10, 20, 35, 45, 55, 60, 65, 70, 75, 80, 90, 
#          100, 110, 114, 120, 130, 140, 143, 150, 160, 170, 180,
#          200, 210, 215, 225, 235, 245, 250, 255, 270, 280, 285,
#          315, 320, 330)
set.seed(0)
idxs= sample(1:length(DATA$domain$vertices$x), size=ceiling(length(DATA$domain$vertices$x)/15))
centroids_ = mesh$nodes[idxs,]
centroids_ = rbind(centroids_, mesh$nodes[c(1, 23, 119:121, 209, 211), ])

result_ <- kmeans(x=data_, centers=30, iter.max = 100)
centroids_ = rbind( centroids_, projection.points.1.5D(mesh, locations= result_$centers) )

{
  x11()
  plot(mesh, pch=".")
  points(data_, col="red", pch=16, cex=2)
  points(centroids_, col="green4", pch=16, cex=2)
  legend("topright", legend=c("data", "centroids"),
          col=c("red", "green4"), pch = 16, cex=1.25 )
}


# # lines == true edges of the network / edges == edges of the discretized network
# result_ <- kmeans(x=data_, centers=nregion, iter.max = 100)
# centroids_ = rbind( centroids_, projection.points.1.5D(mesh, locations= result_$centers) )

# {
# x11()
# plot(mesh, pch=".")
# points(data_, col="red", pch=16, cex=2)
# points(result_$centers, col="blue", pch = 16, cex=2)
# points(centroids_, col="green4", pch=16, cex=2)
# legend("topright", legend=c("data", "2D", "1.5D"), 
#         col=c("red", "blue", "green4"), pch = 16, cex=1.25 )
# }

lines_to_region <- set_region(centroids_, mesh=mesh, LN = DATA)

nregion = nrow(centroids_)

incidence_matrix = matrix(0, nrow=nregion, ncol=nedges)

for(i in 1:nedges){
   incidence_matrix[ lines_to_region[i],i] = 1
}

mask_= rep(0, times=nregion)
for( i in 1:nregion){
  if( !sum(incidence_matrix[i,])){
    mask_[i] = 1
    print(mask_[i])
  }
}

if( sum(mask_)){
  incidence_matrix = incidence_matrix[-which(mask_==1),]
  nregion = nregion - sum(mask_)
  centroids_ = centroids_[-which(mask_==1),]
  lines_to_region <- set_region(centroids_, mesh=mesh, LN = DATA)
}


response = rep(0, times= nregion)

for( i in DATA$data$seg){
  response[lines_to_region[i]] = response[lines_to_region[i]] + 1
}
range(response)

{
x11()
plot_region(lines_to_region, response, LN=DATA,mesh=mesh, nregion = 1)
}

{ pdf("regions.pdf")
  for(i in 1:nregion)
    plot_region(lines_to_region, response, LN=DATA, mesh=mesh, nregion = i)
  dev.off()
}

{
  x11()
  plot_region_gradient_color(mesh, response, lines_to_region, line.size = 1.5)
}

lambda_GSR = 10^seq(from=-5,to=-1,length.out=40)

GSR_PDE <- smooth.FEM(observations = response,
                         covariates = NULL,
                         FEMbasis = FEMbasis,
                         incidence_matrix = incidence_matrix,
                         lambda = lambda_GSR,
                         lambda.selection.criterion = "grid",
                         lambda.selection.lossfunction = "GCV",
                         DOF.evaluation = "exact",
                         family="poisson")

lambda_opt <- GSR_PDE$optimization$lambda_position
{
x11()
plot(log10(lambda_GSR), GSR_PDE$optimization$GCV_vector, xlab="log10(lambda)", ylab="GCV")
}

Mass = fdaPDE::CPP_get.FEM.Mass.Matrix(FEMbasis)
coeff_ <- exp(GSR_PDE$solution$f[,lambda_opt])
coeff_ <- coeff_ / sum( Mass%*% coeff_ )
plot(FEM(coeff_, FEMbasis))
range(coeff_)

lambda_DE = 10^seq(from=-5, to=-2.5,length.out = 20)
DE_PDE = fdaPDE::DE.FEM(data = cbind(DATA$data$x, DATA$data$y), FEMbasis = FEMbasis,
                        lambda = lambda_DE,
                        preprocess_method ="RightCV",
                        nfolds = 10)
x11()
plot(log10(lambda_DE), DE_PDE$CV_err, xlab="log10(lambda)", ylab="CV")

sum(Mass%*% exp(DE_PDE$g))
plot(FEM(exp(DE_PDE$g), FEMbasis))
range(exp(DE_PDE$g))

sum( Mass%*% (exp(DE_PDE$g) - coeff_)^2 )

################################################################################
K = 10

dataList = set_Kfold_data(DATA$data, seed = 0)
n = nrow(chicago$data)

date_ = gsub(":","_",gsub(" ","-",Sys.time()))
if(!dir.exists("data/")) {
  dir.create("data/")
}

if( !dir.exists(paste("data/chicago/",sep=""))){
  dir.create(paste("data/chicago/",sep=""))
}

folder.name = paste("data/chicago/", date_,"/",sep="")

if(!dir.exists(folder.name)) {
  dir.create(folder.name)
}


CV_errors = matrix(0, nrow = K, ncol = 2)

for(i in 1:K){
  
  tmp = get_Kfold_data(dataList, iter = i)
  train_data = tmp$train_data
  test_data = tmp$test_data
  
  # DE-PDE 
  DE_PDE = fdaPDE::DE.FEM(data = cbind(train_data$x, train_data$y), FEMbasis = FEMbasis,
                        lambda = lambda_DE,
                        preprocess_method ="RightCV",
                        nfolds = 10)
  DE_PDE.FEM = FEM(coeff= exp(DE_PDE$g), FEMbasis= FEMbasis )
  
  CV_errors[i,1] = cv_error(FEM = DE_PDE.FEM, 
                            R0 = Mass, data.k = test_data)

  # GSR-PDE 
  response = rep(0, times= nregion)
  for(k in train_data$seg){
    response[lines_to_region[k]] = response[lines_to_region[k]] + 1
  }

  GSR_PDE <- smooth.FEM(observations = response,
                         covariates = NULL,
                         FEMbasis = FEMbasis,
                         incidence_matrix = incidence_matrix,
                         lambda = lambda_GSR,
                         lambda.selection.criterion = "grid",
                         lambda.selection.lossfunction = "GCV",
                         DOF.evaluation = "exact",
                         family="poisson")

  lambda_opt <- GSR_PDE$optimization$lambda_position
  
  coeff_ <- exp(GSR_PDE$solution$f[,lambda_opt])
  coeff_ <- coeff_ / sum( Mass%*% coeff_ )

  GSR_PDE.FEM = FEM(coeff_, FEMbasis) 
  CV_errors[i,2] = cv_error(FEM = GSR_PDE.FEM, 
                            R0 = Mass, data.k = test_data)


}

save(CV_errors, date_, folder.name,
     file = paste(folder.name, "CV_error.RData", sep=""))

boxplot(CV_errors, main = "CV error")


response = rep(0, times= nregion)
for( i in DATA$data$seg){
  response[lines_to_region[i]] = response[lines_to_region[i]] + 1
}
range(response)
GSR_PDE <- smooth.FEM(observations = response,
                         covariates = NULL,
                         FEMbasis = FEMbasis,
                         incidence_matrix = incidence_matrix,
                         lambda = lambda_GSR,
                         lambda.selection.criterion = "grid",
                         lambda.selection.lossfunction = "GCV",
                         DOF.evaluation = "exact",
                         family="poisson")

lambda_opt <- GSR_PDE$optimization$lambda_position

coeff_ <- exp(GSR_PDE$solution$f[,lambda_opt])
coeff_ <- coeff_ / sum( Mass%*% coeff_ )
GSR_PDE.FEM = FEM(coeff_, FEMbasis)

DE_PDE = fdaPDE::DE.FEM(data = cbind(DATA$data$x, DATA$data$y), FEMbasis = FEMbasis,
                        lambda = lambda_DE,
                        preprocess_method ="RightCV",
                        nfolds = 10)
DE_PDE.FEM = FEM(coeff= exp(DE_PDE$g), FEMbasis)

save(GSR_PDE.FEM, DE_PDE.FEM, 
      file = paste(folder.name, "estimates.RData", sep=""))

{
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

}

folder.imgs = paste(folder.name,"img/",sep="")
if(!dir.exists(folder.imgs)) {
  dir.create(folder.imgs)
}

{
pdf(paste(folder.imgs,"CV_error.pdf",sep=""))
methods = c(T,T,F,F)
methods.names = c("DE-PDE", "GSR-PDE")
print(boxplot_CV_error(CV_errors = CV_errors, 
                 methods = methods,
                 methods.names = methods.names))
dev.off()
}

{
line.size = 1.
pdf(paste(folder.imgs, "estimates.pdf",sep=""))

estimates = list()
estimates[[1]] = DE_PDE.FEM
estimates[[2]] = GSR_PDE.FEM

num_edges= dim(mesh$edges)[1]
coef=matrix(0, nrow= num_edges, ncol=length(estimates) )

for(i in 1:length(estimates)){
for(e in 1:num_edges){
    
  coef[e,i]= (estimates[[i]]$coeff[mesh$edges[e,1]] + estimates[[i]]$coeff[mesh$edges[e,2]])/2  
    
  }
}
  
max.col = max(coef)
min.col = min(coef)

PLOTS = list()
for(i in 1:length(estimates))
{
    PLOTS[[i]] = R_plot_graph.ggplot2.2( estimates[[i]], 
                                          title = methods.names[i],
                                          color.min = min.col,
                                          color.max = max.col,
                                           palette=viridis, line.size = line.size)
    
}

for(i in 1:length(estimates)){
  print(PLOTS[[i]])
}

PLOTS = list()
for(i in 1:length(estimates)){
    PLOTS[[i]] = R_plot_graph.ggplot2.2( estimates[[i]], # FEM object
                                          title = methods.names[[i]],
                                          color.min = min.col,
                                          color.max = max.col,
                                          palette=magma, line.size = line.size)
    
}
for(i in 1:length(estimates)){
  print(PLOTS[[i]])
}


PLOTS = list()
for(i in 1:length(estimates)){
    PLOTS[[i]] = R_plot_graph.ggplot2.2( estimates[[i]], # FEM object
                                          title = methods.names[[i]],
                                          color.min = min.col,
                                          color.max = max.col,
                                          palette=jet.col, line.size = line.size)
    
}
for(i in 1:length(estimates)){
  print(PLOTS[[i]])
}

dev.off()
}

# ref mesh 
estimates = list()
estimates[[1]] = DE_PDE.FEM
estimates[[2]] = GSR_PDE.FEM

mesh.ref = refine.mesh.1.5D(DE_PDE.FEM$FEMbasis$mesh, delta = 0.0125)
FEMbasis.ref = create.FEM.basis(mesh.ref)
locs = mesh.ref$nodes

num_edges= dim(mesh.ref$edges)[1]
coef=matrix(0, nrow= num_edges, ncol=length(estimates) )

for(i in 1:length(estimates)){
for(e in 1:num_edges){
    
  coef[e,i]= (estimates[[i]]$coeff[mesh.ref$edges[e,1]] + estimates[[i]]$coeff[mesh.ref$edges[e,2]])/2  
    
  }
}
  
max.col = max(coef)
min.col = min(coef)

for(i in 1:length(estimates))
  estimates[[i]] = FEM(eval.FEM(estimates[[i]], locations = locs),
                       FEMbasis.ref)

{
pdf(paste(folder.imgs, "estimates_ref.pdf",sep=""))
PLOTS = list()
for(i in 1:length(estimates))
{
    PLOTS[[i]] = R_plot_graph.ggplot2.2( estimates[[i]], 
                                          title = methods.names[i],
                                          color.min = min.col,
                                          color.max = max.col,
                                          palette=viridis, line.size = line.size)
    
}

for(i in 1:length(estimates)){
  print(PLOTS[[i]])
}

PLOTS = list()
for(i in 1:length(estimates)){
    PLOTS[[i]] = R_plot_graph.ggplot2.2( estimates[[i]], # FEM object
                                          title = methods.names[[i]],
                                          color.min = min.col,
                                          color.max = max.col,
                                          palette=magma, line.size = line.size)
    
}
for(i in 1:length(estimates)){
  print(PLOTS[[i]])
}


PLOTS = list()
for(i in 1:length(estimates)){
    PLOTS[[i]] = R_plot_graph.ggplot2.2( estimates[[i]], # FEM object
                                          title = methods.names[[i]],
                                          color.min = min.col,
                                          color.max = max.col,
                                          palette=jet.col, line.size = line.size)
    
}
for(i in 1:length(estimates)){
  print(PLOTS[[i]])
}

dev.off()
}
