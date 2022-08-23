#################################################
############## Density Estimation ###############
#################  Case Study ###################
#################################################
graphics.off()
rm(list=ls())

setwd("DensityEstimation/")
source("../utils.R")
source("setting.R")

data("chicago")

spat.stat.linnet = chicago$domain
mesh = as.fdaPDE.spatstat.linnet(chicago$domain)
FEMbasis = create.FEM.basis(mesh)
Mass = CPP_get.FEM.Mass.Matrix(FEMbasis)

K = 10

dataList = set_Kfold_data(chicago$data)
CV_errors = matrix(0, nrow = K, ncol = 4)

for(i in 1:K){
  
  tmp = get_Kfold_data(dataList, iter = i)
  train_data = tmp$train_data
  test_data = tmp$test_data

  # DE-PDE 
  lambda = 10^seq(from=-2, to=1,length.out = 10)
  DE_PDE = fdaPDE::DE.FEM(data = cbind(train_data$x, train_data$y), FEMbasis = FEMbasis,
                          lambda = lambda,
                          preprocess_method ="SimplifiedCV",
                          nfolds = 10)
  
  CV_errors[i,1] = cv_error(FEM = FEM(coeff= exp(DE_PDE$g), FEMbasis), R0 = Mass, data.k = test_data)
    
  
}


