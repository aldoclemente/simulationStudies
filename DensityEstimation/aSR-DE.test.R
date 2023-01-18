#################################################
##############  aSR-PDE / DE-PDE  ###############
##################   Test   #####################
#################################################

setwd("DensityEstimation/")

graphics.off()
rm(list=ls())

source("../utils.R")
source("../SR-chicago-case-study/refine1D.R")
source("setting.R")
source("utils.R")

nsim =  30#30 #1
ntest = 1
domains = c("simplenet", "ontario")

# methods[1] -> DE-PDE
# methods[2] -> GSR-PDE

methods = c(T,T) # mask
methods.names = c("DE-PDE", "GSR-PDE")
tests.names = c("test_1", "test_2")
 
mesh.tmp = create.mesh.1.5D(nodes = cbind(simplenet$vertices$x, simplenet$vertices$y),
                            edges = cbind(simplenet$from, simplenet$to))
mesh.tmp = refine.by.splitting.mesh.1.5D(mesh.tmp)
#mesh.tmp = refine.by.splitting.mesh.1.5D(mesh.tmp)
centroids_ = mesh.tmp$nodes
dim(centroids_)
PP = runiflpp(n=15, L = simplenet) 
centroids_ = rbind(centroids_, cbind(PP$data$x, PP$data$y))

sett = setting(domains[ntest]) 
mesh = sett$mesh
centroids_ = projection.points.1.5D(mesh=mesh, location=centroids_)

FEMbasis = sett$FEMbasis
nnodes = sett$nnodes

#nregion = 20 #8 #6 #10 
nregion = nrow(centroids_)
nedges = nrow(mesh$edges)

spat.stat.linnet = sett$spat.stat.linnet

locs.test = runiflpp(1000, spat.stat.linnet)
locs.test = cbind(locs.test$data$x, locs.test$data$y)
# point pattern #
if(ntest==1) n = c(250)#n = c(50, 100, 150, 250)
#if(ntest==2) n = c(100, 250, 500, 1000)
set.seed(1234)

if(ntest==1) sources = c(6,8) 
#if(ntest==2) sources = c(32, 185, 400) #63  7 150 250

auxiliary_test = aux_test[[ntest]]

DENSITY = linfun(auxiliary_test, spat.stat.linnet) # test 2

PP = rlpp(n=n, f = DENSITY)  
data = cbind(PP$data$x, PP$data$y)
result_ <- kmeans(x=data, centers=nregion, iter.max = 100)
#centroids_ = projection.points.1.5D(mesh, locations= result_$centers)

# mask_ = mesh$nodesmarkers[1:10]

# mid_points = matrix(0, nrow=0, ncol=2)
# for(i in 1:length(simplenet$from)){
#     n1 = simplenet$from[i]
#     n2 = simplenet$to[i]
#     x.m = mean(c(simplenet$vertices$x[n1], simplenet$vertices$x[n2]))
#     y.m = mean(c(simplenet$vertices$y[n1], simplenet$vertices$y[n2]))
#     mid_points = rbind(mid_points, c(x.m, y.m))
# }

# mid_points = projection.points.1.5D(mesh, mid_points)

# centroids_ = rbind(centroids_, mesh$nodes[mesh$nodesmarkers, ], mid_points)

# nregion = nregion + sum(mesh$nodesmarkers) + nrow(mid_points)
# {
# x11()
# plot(mesh, pch=".")
# points(data, col="red", pch=16, cex=1)
# points(result_$centers, col="blue", pch = 16, cex=1)
# points(centroids_, col="green4", pch=16, cex=2)
# legend("topright", legend=c("data", "2D", "1.5D"), 
#         col=c("red", "blue", "green4"), pch = 16, cex=1.25 )
#}

lines_to_region <- set_region(centroids_, mesh=mesh, LN = PP)
incidence_matrix = matrix(0, nrow=nregion, ncol=nedges)

for(i in 1:nedges){
  incidence_matrix[ lines_to_region[i],i] = 1
}

for(i in 1:nregion){
  print(sum(incidence_matrix[i,]))
}
sum(incidence_matrix)

mask_= matrix(0, nrow=nregion, ncol=1)
for( i in 1:nregion){
  if( !sum(incidence_matrix[i,])){
    mask_[i] = 1
    print(mask_[i])
  }
}

if( sum(mask_)){
  incidence_matrix = incidence_matrix[-mask_,]
  nregion = nregion - sum(mask_)
  centroids_ = centroids_[-mask_,]
  lines_to_region <- set_region(centroids_, mesh=mesh, LN = PP)
}

response = rep(0, times= nregion)
for( i in PP$data$seg){
  response[lines_to_region[i]] = response[lines_to_region[i]] + 1
}
range(response)

{
x11()
plot_region(lines_to_region, response, LN=DATA,mesh=mesh)
}

Mass = CPP_get.FEM.Mass.Matrix(FEMbasis)
true.density = DENSITY(x=mesh$nodes[,1], y=mesh$nodes[,2])
true.density.FEM = true.density / sum( Mass %*% true.density)
true.density = DENSITY(x=locs.test[,1], y=locs.test[,2]) / sum( Mass %*% true.density)

lambda = 10^seq(from=-4, to=1,length.out = 500)
{
start = Sys.time()
GSR_PDE   <- smooth.FEM(observations = response,
                         FEMbasis = FEMbasis,
                         incidence_matrix = incidence_matrix,
                         lambda = lambda,
                         lambda.selection.criterion = "grid",
                         lambda.selection.lossfunction = "GCV",
                         DOF.evaluation = "exact",
                         family="poisson")
cat(paste("- GSR-PDE DONE, time elapsed = ", difftime(Sys.time(),start, units = "mins")," mins \n", sep="") )
plot(log10(lambda), GSR_PDE$optimization$GCV_vector)
lambda_opt = GSR_PDE$optimization$lambda_position
coeff_ <- exp(GSR_PDE$solution$f[,lambda_opt])
coeff_ <- coeff_ / sum( Mass%*% coeff_ )
plot(FEM(coeff_, FEMbasis))
sqrt(mean( (true.density - eval.FEM(FEM(coeff=coeff_,FEMbasis),
                                              locations = locs.test))^2 ))
}

plot(FEM(true.density.FEM, FEMbasis))
range(true.density.FEM)
range(coeff_)

rmse.DE_PDE = matrix(nrow=nsim, ncol = length(n))
rmse.GSR_PDE = rmse.DE_PDE

DE_PDE.FEM = matrix(0,nrow=nnodes,ncol=1)
GSR_PDE.FEM = DE_PDE.FEM

date_ = gsub(":","_",gsub(" ","-",Sys.time()))
if(!dir.exists("data/")) {
  dir.create("data/")
}

if( !dir.exists(paste("data/", tests.names[ntest],"/",sep=""))){
  dir.create(paste("data/", tests.names[ntest],"/",sep=""))
}

folder.name = paste("data/", tests.names[ntest],"/",date_,"/",sep="")

if(!dir.exists(folder.name)) {
  dir.create(folder.name)
}


if(ntest==1)lambda = 10^seq(from=-4, to=-3,length.out = 20)
if(ntest==2)lambda = 10^seq(from=-8, to=-5,length.out = 20)


for(j in 1:length(n)){
  cat(paste("###############  n = ", n[j], "  ###############\n", sep="") ) 
for( i in 1:nsim){
  cat(paste("###############  ", i, " / ", nsim,"  ###############\n", sep="") )
  PP = rlpp(n=n[j], f = DENSITY)  
  data = cbind(PP$data$x, PP$data$y)
  
  # DE-PDE #
if(methods[1]){
  lambda = 10^seq(from=-4, to=-3,length.out = 20)
  start = Sys.time()
  invisible(capture.output( DE_PDE <-  fdaPDE::DE.FEM(data = data, 
                                                      FEMbasis = FEMbasis,
                                                      preprocess_method ="RightCV",
                                                      nfolds = 10,
                                                      lambda = lambda) ) )
  cat(paste("- DE-PDE DONE, time elapsed = ", difftime(Sys.time(),start, units = "mins")," mins \n", sep="") )
  #plot(DE_PDE$CV_err, main=paste("n = ", n[j], sep=""))
  rmse.DE_PDE[i,j] = sqrt(mean((true.density - eval.FEM(FEM(coeff=exp(DE_PDE$g),FEMbasis),
                                              locations = locs.test))^2 ))
  DE_PDE.FEM = DE_PDE.FEM + exp(DE_PDE$g) / nsim
}

# GSR-PDE
if(methods[2]){
  lines_to_region <- set_region(centroids_, mesh=mesh, LN = PP)
  # incidence_matrix = matrix(0, nrow=nregion, ncol=nedges)

  # for(k in 1:nrow(new_to_old)){
  #   incidence_matrix[ lines_to_region[new_to_old[k]],i] = 1
  # }

  response = rep(0, times= nregion)

  for( k in PP$data$seg){
    response[lines_to_region[k]] = response[lines_to_region[k]] + 1
  }

  lambda = 10^seq(from=-4, to=1,length.out = 500)
  start = Sys.time()
  GSR_PDE   <- smooth.FEM(observations = response,
                         FEMbasis = FEMbasis,
                         incidence_matrix = incidence_matrix,
                         lambda = lambda,
                         lambda.selection.criterion = "grid",
                         lambda.selection.lossfunction = "GCV",
                         DOF.evaluation = "exact",
                         family="poisson")
cat(paste("- GSR-PDE DONE, time elapsed = ", difftime(Sys.time(),start, units = "mins")," mins \n", sep="") )

plot(log10(lambda), GSR_PDE$optimization$GCV_vector)
lambda_opt <- GSR_PDE$optimization$lambda_position
#coeff_ <- GSR_PDE$solution$f[,lambda_opt] / sum( Mass%*% GSR_PDE$solution$f[,lambda_opt] )
coeff_ <- exp(GSR_PDE$solution$f[,lambda_opt])
coeff_ <- coeff_ / sum( Mass%*% coeff_ )
rmse.GSR_PDE[i,j] = sqrt(mean( (true.density - eval.FEM(FEM(coeff=coeff_,FEMbasis),
                                              locations = locs.test))^2 ))

GSR_PDE.FEM = GSR_PDE.FEM + coeff_ / nsim
}
  # spatstat returns the INTENSITY function. 
  # The estimation of DENSITY or INTENSITY is equivalent, if n is fixed,
  # INTENSITY(p) = n DENSITY(p) for all p. 

  # See 
  # McSwiggan, Greg, Adrian Baddeley, and Gopalan Nair. 
  # "Kernel density estimation on a linear network." 
  # Scandinavian Journal of Statistics 44.2 (2017): 324-345.

  cat(paste("###############  ############### ###############\n", sep="") )
}
}

boxplot( cbind(rmse.DE_PDE, rmse.GSR_PDE) )


save(rmse.DE_PDE, rmse.GSR_PDE,  
      methods.names, methods, n,
      folder.name, date_, file = paste(folder.name,"RMSE_",".RData", sep=""))

save(DE_PDE.FEM, GSR_PDE.FEM, true.density.FEM,
       FEMbasis,
       file = paste(folder.name,"estimates",".RData", sep=""))

folder.imgs = paste(folder.name,"img/",sep="")
if(!dir.exists(folder.imgs)) {
  dir.create(folder.imgs)
}
nsim = nrow(rmse.DE_PDE)
RMSE = matrix(nrow = nsim * length(n), ncol=length(methods.names)) 
RMSE[,1] = as.vector(rmse.DE_PDE)
RMSE[,2] = as.vector(rmse.GSR_PDE)

{
pdf(paste(folder.imgs,"RMSE",".pdf",sep=""), width=12)

print(boxplot_RMSE(RMSE, 
             methods = methods,
             methods.names = methods.names, 
             nsim = nsim,
             title.size=26,
             begin=1, #color
             end=0.25,   #color
             width =0.75,
             n = n,
             title="RMSE"))
dev.off()
}

{
pdf(paste(folder.imgs,"estimates",".pdf",sep=""))
estimates = list()
estimates[[1]] = FEM(DE_PDE.FEM, FEMbasis)
estimates[[2]] = FEM(GSR_PDE.FEM, FEMbasis)

PLOTS <- plot_estimates(estimates,
                        true.density = FEM(true.density.FEM,FEMbasis),
                        methods = methods,
                        methods.names = methods.names,
                        true.density.main = "Density",
                        palette="viridis",
                        line.size = 0.75)

print(PLOTS$density.plot)
for(i in 1:length(estimates)){
  if(methods[i]) print(PLOTS$estimates.plot[[i]])
}

PLOTS <- plot_estimates(estimates, 
                        true.density = FEM(true.density.FEM,FEMbasis),
                        methods = methods,
                        methods.names = methods.names,
                        true.density.main = "Density",
                        palette="magma",
                        line.size = 0.75)
print(PLOTS$density.plot)
for(i in 1:length(estimates)){
  if(methods[i])  print(PLOTS$estimates.plot[[i]])
}

PLOTS <- plot_estimates(estimates,
                        true.density = FEM(true.density.FEM,FEMbasis),
                        methods = methods,
                        methods.names = methods.names,
                        true.density.main = "Density",
                        line.size = 0.75)

print(PLOTS$density.plot)
for(i in 1:length(estimates)){
  if(methods[i])  print(PLOTS$estimates.plot[[i]])
}

dev.off()

pdf(paste(folder.imgs,"point_pattern",".pdf",sep=""))
for(i in n){
PP = rlpp(i, DENSITY)
plot(mesh, pch=".")
points(PP$data$x, PP$data$y, pch=16, col ="red3")
}
dev.off()

}