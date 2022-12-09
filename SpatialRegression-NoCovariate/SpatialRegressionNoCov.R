############################
###  Spatial Regression  ###
###    No Covariates     ###
###  locations at nodes  ###
############################

setwd("SpatialRegression-NoCovariate/")
library(KrigLinCaution)
library(plotrix)
source("../utils.R")
source("../create_knots.R")
source("settings.R")
source("../rank_reduced_kriging.R")
source("NoCovariatesCore.R")

  func = "exp" # 'exp' 'sin', see settings.R
  domain = "estevan" # 'estevan' 'ontario', see settings.R
  
  set = setting(network = domain)
  
  mesh = set$mesh
  FEMbasis = set$FEMbasis
  nnodes = set$nnodes
  
  nnodes = nrow(mesh$nodes)
  
  # in ../utils.R 
  ND_ = compute_NetworkDist_matrix(mesh, 1:nnodes)
  ED_ = computed_EuclideanDist_matrix(mesh, 1:nnodes)
  lambda=10^seq(from=-8,to=-4,length.out=20)
  plot(mesh, asp=1,pch=16)
  
  set.seed(1234)
  
  # Test 1: Estevan Exp 
  sample.idx =c(2378, 1271, 1802, 1529, 2693)
  lambda = 10^seq(from=-5,to=0.,length.out=20)
  field = fun(func)(ND_,source = 63, sigma=2)
  field.2 = fun(func)(ND_,source = sample.idx[1], sigma=1.5)
  field.3 = fun(func)(ND_,source = 1039, sigma=2)
  field.4 = fun(func)(ND_,source = sample.idx[5], sigma=1.5)
  field = field + field.2 + field.3 + field.4 
  observations_ = field + rnorm(nnodes, mean=0, sd = 0.05*diff(range(field)))
  n_data = c(50, 75, 150, 200)
  
plot(FEM(observations_, FEMbasis))

n_sim = 30
invisible(capture.output(
  results<-NoCovariatesCore(ND_, ED_, observations_,
                            n_sim=n_sim, n_data,lambda,mesh,
                            true.signal = field,
                            model_ = c(T,T,T,T))))

results$tot.time
RMSE = results$RMSE
mean.field.fdaPDE = results$mean.field.fdaPDE
estimates = results$estimates
boxplot_RMSE(RMSE, n_data, 
             model_ = c(T,T,T,T),
             names_ = c("fdaPDE","GWR","lattice","Krig"))


date_ = gsub(":","_",gsub(" ","-",Sys.time()))
tail_ = paste(domain, func, date_, sep="-")

if(domain=="estevan"){
  head = "test-1"
}else{
  head = "test-2"
}


filename_ = paste(paste(paste("data/",head,sep=""),tail_, sep="-"), ".RData", sep="")
imgfile_  = paste(paste(paste(paste("img/" ,head,sep=""),"plots",sep="-"),tail_,sep="-"),".pdf",sep="")

save(RMSE, field, observations_, n_data, mean.field.fdaPDE, imgfile_, FEMbasis, estimates,
     file = filename_)

#################################
# imgs #
source("../SpatialRegression-NoCovariate/RegressionNoCovPlots.R")

palette = "magma" # "viridis" "magma
imgfile_ = paste("img/SpaceRegression-NoCov-",palette,".pdf",sep="")

if(palette == "ggplot")
  palette=NULL

RegressionNoCovPlots(imgfile = imgfile_,
                     true.field = field, 
                     mean.field.fdaPDE = mean.field.fdaPDE,
                     observations = observations_,
                     RMSE = RMSE,
                     FEMbasis = FEMbasis,
                     n_data = n_data,
                     palette = palette,
                     legend.pos.RMSE = "right",
                     line.size=0.5)

colors = viridis(n=4, begin=0.95, end=0.25)

pdf("img/RMSE-NoCov-18.pdf", width=18)
boxplot_RMSE(RMSE, n_data, model_ = c(T,T,T,T), 
             names_ = c("SR-PDE","GWR","Lattice","RR-Krig"),
             legend.pos = c(0.95,0.8625), palette=palette,
             colors=colors)
dev.off()



