#################################################
############## Density Estimation ###############
##################   Test   #####################
#################################################

setwd("DensityEstimation/")

graphics.off()
rm(list=ls())
library(fdaPDE)
source("../utils.R")
source("../settings.R")
source("setting.R")
source("utils.R")
nsim =  30 #30
ntest = 1
domains = c("simplenet", "ontario")

# methods[1] -> DE-PDE
# methods[2] -> KDE-PDE
# methods[3] -> KDE-ES   (very slow !)
# methods[4] -> KDE-2D
# methods[5] -> VORONOI  (slow !)

methods = c(T,T,F,T,T) # mask
methods.names = c("DE-PDE", "KDE-PDE", "KDE-ES", "KDE-2D", "VORONOI")
tests.names = c("test_1", "test_2")
 
sett = setting(domains[ntest]) 
mesh = sett$mesh
FEMbasis = sett$FEMbasis
nnodes = sett$nnodes

spat.stat.linnet = sett$spat.stat.linnet

locs.test = runiflpp(1000, spat.stat.linnet)
locs.test = cbind(locs.test$data$x, locs.test$data$y)
# point pattern #
if(ntest==1) n = c(50, 100, 150, 250)
if(ntest==2) n = c(100, 250, 500, 1000)
set.seed(1234)

if(ntest==1) sources = c(6,8) 
if(ntest==2) sources = c(32, 185, 400) #63  7 150 250

auxiliary_test = aux_test[[ntest]]

DENSITY = linfun(auxiliary_test, spat.stat.linnet) # test 2

Mass = CPP_get.FEM.Mass.Matrix(FEMbasis)
true.density = DENSITY(x=mesh$nodes[,1], y=mesh$nodes[,2])
true.density.FEM = true.density / sum( Mass %*% true.density)

true.density = DENSITY(x=locs.test[,1], y=locs.test[,2]) / sum( Mass %*% true.density)

rmse.DE_PDE = matrix(nrow=nsim, ncol = length(n))
rmse.KDE_PDE = rmse.DE_PDE
rmse.KDE_ES = rmse.DE_PDE
rmse.KDE_2D = rmse.DE_PDE
rmse.VORONOI = rmse.DE_PDE

DE_PDE.FEM = matrix(0,nrow=nnodes,ncol=1)
KDE_PDE.FEM = DE_PDE.FEM
KDE_ES.FEM = DE_PDE.FEM
KDE_2D.FEM = DE_PDE.FEM
VORONOI.FEM = DE_PDE.FEM

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
if(ntest==2)lambda = 10^seq(from=-6, to=-3,length.out = 20)


for(j in 1:length(n)){
  cat(paste("###############  n = ", n[j], "  ###############\n", sep="") ) 
for( i in 1:nsim){
  cat(paste("###############  ", i, " / ", nsim,"  ###############\n", sep="") )
  PP = rlpp(n=n[j], f = DENSITY)  
  data = cbind(PP$data$x, PP$data$y)
  
  # DE-PDE #
  if(methods[1]){
  start = Sys.time()
  # invisible(capture.output( DE_PDE <-  fdaPDE::DE.FEM(data = data, FEMbasis = FEMbasis,
  #                                   lambda = lambda,
  #                                   preprocess_method ="RightCV",
  #                                   nfolds = 10) ))
  invisible(capture.output( DE_PDE <-  fdaPDE::DE.FEM(data = data, 
                                                      FEMbasis = FEMbasis,
                                                      lambda = lambda[1]) ) )
  cat(paste("- DE-PDE DONE, time elapsed = ", difftime(Sys.time(),start, units = "mins")," mins \n", sep="") )
  #plot(DE_PDE$CV_err, main=paste("n = ", n[j], sep=""))
  rmse.DE_PDE[i,j] = sqrt(mean((true.density - eval.FEM(FEM(coeff=exp(DE_PDE$g),FEMbasis),
                                              locations = locs.test))^2 ))
}
  # spatstat returns the INTENSITY function. 
  # The estimation of DENSITY or INTENSITY is equivalent, if n is fixed,
  # INTENSITY(p) = n DENSITY(p) for all p. 

  # See 
  # McSwiggan, Greg, Adrian Baddeley, and Gopalan Nair. 
  # "Kernel density estimation on a linear network." 
  # Scandinavian Journal of Statistics 44.2 (2017): 324-345.

  # KDE-PDE
  if(methods[2]){  
start = Sys.time()  
invisible(capture.output(bw <- bw.lppl(X = PP) ))
invisible(capture.output(KDE_PDE <- densityHeat(x = as.lpp(PP), sigma = as.numeric(bw), iterMax = 1e+9) )) 
cat(paste("- KDE-PDE DONE, time elapsed = ", difftime(Sys.time(),start, units = "mins")," mins \n", sep="") )
rmse.KDE_PDE[i,j] = sqrt(mean( (true.density - eval.FEM(FEM(coeff=as.linfun(KDE_PDE/n[j])(mesh$nodes[,1], mesh$nodes[,2]),
                                               FEMbasis),
                                               locations = locs.test))^2 ))
}

  # KDE-ES 
  if(methods[3]){
  start = Sys.time()
  invisible(capture.output(bw <- bw.lppl(X = PP) ))
  invisible(capture.output(KDE_ES <- densityEqualSplit(x = PP, sigma = bw) ))
  cat(paste("- KDE-ES DONE, time elapsed = ", difftime(Sys.time(),start, units = "mins")," mins \n", sep="") )

  rmse.KDE_ES[i,j] = sqrt(mean( (true.density - eval.FEM(FEM(coeff=as.linfun(KDE_ES/n[j])(mesh$nodes[,1], mesh$nodes[,2]),
                                        FEMbasis),
                                    locations = locs.test))^2 ))
 }
  
  # KDE-2D
  if(methods[4]){
  start = Sys.time()
  invisible(capture.output(bw <- bw.scott(X = PP) ))
  invisible(capture.output(KDE_2D <-  densityQuick.lpp(X = PP, sigma = bw) )) #, at = points)
  cat(paste("- KDE-2D DONE, time elapsed = ", difftime(Sys.time(),start, units = "mins")," mins \n", sep="") )

  rmse.KDE_2D[i,j] = sqrt(mean( (true.density - eval.FEM(FEM(coeff=as.linfun(KDE_2D/n[j])(mesh$nodes[,1], mesh$nodes[,2]),
                                                   FEMbasis),
                                               locations = locs.test))^2 ))
 }
  
  # KDE-VORONOI #
  if(methods[5]){
  start = Sys.time()
  invisible(capture.output(bw <- bw.voronoi(X = PP) ))
  invisible(capture.output(VORONOI <- densityVoronoi(X = PP, sigma = bw) ))
  cat(paste("- VORONOI DONE, time elapsed = ", difftime(Sys.time(),start, units = "mins")," mins \n", sep="") )

  rmse.VORONOI[i,j] = sqrt(mean( (true.density - eval.FEM(FEM(coeff=as.linfun(VORONOI/n[j])(mesh$nodes[,1], mesh$nodes[,2]),
                                                   FEMbasis),
                                               locations = locs.test))^2 ))
  }
  cat(paste("###############  ############### ###############\n", sep="") )

  if(j == 3){
    if(methods[1]) DE_PDE.FEM = DE_PDE.FEM + exp(DE_PDE$g) / nsim
    if(methods[2]) KDE_PDE.FEM = KDE_PDE.FEM + as.linfun(KDE_PDE/n[j])(mesh$nodes[,1], mesh$nodes[,2])/ nsim
    if(methods[3]) KDE_ES.FEM = KDE_ES.FEM + as.linfun(KDE_ES/n[j])(mesh$nodes[,1], mesh$nodes[,2])/ nsim
    if(methods[4]) KDE_2D.FEM = KDE_2D.FEM + as.linfun(KDE_2D/n[j])(mesh$nodes[,1], mesh$nodes[,2])/ nsim
    if(methods[5]) VORONOI.FEM = KDE_2D.FEM + as.linfun(VORONOI/n[j])(mesh$nodes[,1], mesh$nodes[,2])/ nsim
  
  save(DE_PDE.FEM, KDE_PDE.FEM, KDE_2D.FEM, VORONOI.FEM, KDE_ES.FEM, true.density.FEM,
       FEMbasis,
       file = paste(folder.name,"estimates",".RData", sep=""))
  
}

  save(rmse.DE_PDE, rmse.KDE_PDE, rmse.KDE_ES,
     rmse.KDE_2D, rmse.VORONOI, i ,j, methods.names, methods, n,
     folder.name, date_, file = paste(folder.name,"RMSE",".RData", sep=""))

}
}

source("DE.test.post.processing.R")

#########

pdf(paste(folder.imgs,"nodi.pdf",sep=""))
for(i in 1:nrow(mesh$nodes)){
  plot(mesh$nodes[,1],mesh$nodes[,2],main= paste("node ", i, sep=""))
  points(mesh$nodes[i,1],mesh$nodes[i,2], col="red3",pch=16)
  
}
dev.off()