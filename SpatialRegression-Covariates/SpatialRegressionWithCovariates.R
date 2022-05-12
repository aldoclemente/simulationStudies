############################
###  Spatial Regression  ###
###    With Covariates   ###
###  locations at nodes  ###
############################

setwd("C:/Users/Aldo/Documents/SimulationStudies/SpatialRegression-Covariates")
library(plotrix)
source("../utils.R")
source("WithCovariatesCore.R")
source("../SpatialRegression-NoCovariate/NoCovariatesCore.R")
source("settings.R")

func = "exp" #exp
domain = "estevan" #ontario

set = setting(network = domain)

mesh = set$mesh
FEMbasis = set$FEMbasis
nnodes = set$nnodes

mesh.2D = set$mesh.2D
FEMbasis.2D = set$FEMbasis.2D

# in ../utils.R 
ND_ = compute_NetworkDist_matrix(mesh, 1:nnodes)
ED_ = computed_EuclideanDist_matrix(mesh, 1:nnodes)
set.seed(1234)

#sample.idx = sample(1:nnodes,size=4)
field = fun(func)(ND_,source = 63, sigma=2) #f(ND_,source = 1760, sigma=2)  #63 350
plot(FEM(field, FEMbasis))

# Test 1 estevan exp (source=63, sigma = 2) # NB estevan è sempre test 1 per salvataggio dati 
sample.idx =c(2378, 1271, 1802, 1529, 2693)
field = fun(func)(ND_,source = 63, sigma=2)
field.2 = fun(func)(ND_,source = sample.idx[1], sigma=1.5)
field.3 = fun(func)(ND_,source = 1039, sigma=2)
field.4 = fun(func)(ND_,source = sample.idx[5], sigma=1.5)

# plot(FEM(field,FEMbasis))
# plot(FEM(field.2,FEMbasis))
# plot(FEM(field.3,FEMbasis))
# plot(FEM(field.4,FEMbasis))
field = field + field.2 + field.3 + field.4 
plot(FEM(field, FEMbasis))

W = matrix(nrow=nnodes, ncol=2)
W[,1] = rnorm(nnodes, mean=0.5,sd=0.25)
#W[,2] = rnorm(nnodes, mean=-0.5, sd=0.25)
W[,2] = 1/4*fun("sin")(ND_, source=63, sigma=7.5) #2.5
betas=c(1.,1.)
signal = field + W%*%betas
#observations = signal + rnorm(nnodes, mean=0, sd=0.05*diff(range(field)))
observations = signal + rnorm(nnodes, mean=0, sd=0.05*diff(range(signal)))

plot(FEM(signal,FEMbasis))
plot(FEM(W[,1],FEMbasis))
plot(FEM(W[,2],FEMbasis))
plot(FEM(field,FEMbasis))
lambda = 10^seq(from=-5,to=-1,length.out=25)
n_data = c(250, 500, 1000, 1500)
###

# Test 4 estevan sin (source=63, sigma = 0.75)
# W = matrix(nrow=nnodes, ncol=2)
# W[,1] = rnorm(nnodes, mean=0.5,sd=0.25)
# W[,2] = rbeta(nnodes, shape1=1, shape2=2)
# betas = c(1.,1.)
# signal = field + W%*%betas
# observations = signal + rnorm(nnodes, mean=0, sd=0.05*diff(range(field)))
# lambda = 10^seq(from=-7,to=-5,length.out=25) # -5 era molto buono
# ED_ = NULL 
# ###

plot(FEM((observations_), FEMbasis))

lambda = 10^seq(from=-3,to=-0.5,length.out=25)
lambda.2D = 10^seq(from=-3,to=-0.5,length.out=25)

n_sim=30
invisible(capture.output(
  results<- WithCovariatesCore(ND_, ED_, observations,
                              n_sim, n_data, lambda,lambda.2D, mesh,
                              FEMbasis, FEMbasis.2D,
                              true.signal = signal,
                              model_=c(T,T,T,F), betas, W) ))

results$tot.time
RMSE = results$RMSE
mean.field.fdaPDE = results$mean.field.fdaPDE

boxplot_RMSE(RMSE, n_data, model_ = c(T,T,T,F),
             names_ = c("fdaPDE","GWR","lattice","fdaPDE.2D"))


if(domain=="estevan"){
  head = "test-1"
}else{
  head = "test-2"
}
date_ = gsub(":","_",gsub(" ","-",Sys.time()))
tail_ = paste(domain, func, date_, sep="-")

filename_ = paste(paste(paste("data/",head,sep=""),tail_, sep="-"), ".RData", sep="")
imgfile_  = paste(paste(paste(paste("img/" ,head,sep=""),"plots",sep="-"),tail_,sep="-"),".pdf",sep="")

save(RMSE, 
     field,
     signal,
     observations, 
     n_data, 
     mean.field.fdaPDE, 
     imgfile_, 
     FEMbasis,
     W, betas,
     file = filename_)


####
setwd("C:/Users/Aldo/Documents/SimulationStudies/SpatialRegression-Covariates")
source("RegressionWithCovPlots.R")

RegressionWithCovPlots(imgfile=imgfile_,
                       true.field=field,            # f 
                       true.signal=signal,           # f + W beta 
                       mean.field.fdaPDE=mean.field.fdaPDE,
                       observations=observations,          # f + W beta + eps
                       FEMbasis=FEMbasis,
                       n_data=n_data,
                       W=W , betas=betas,
                       RMSE,legend.pos.RMSE = "right",
                       line.size=0.5)

