get_Kfold_data <- function(SpLinesDataSet_List, iter, K = 10){
  
  tmp_data = data.frame()
  for(i in 1:K){
    if( i == iter){
      test_data = SpLinesDataSet_List[[i]]
    }else{
      tmp_data = rbind(tmp_data, SpLinesDataSet_List[[i]]@data)
    }
  }
  
  train_data = SpatialPointsDataFrame(coords = cbind(tmp_data$X, tmp_data$Y),
                                      data = tmp_data)
  
  ret_list = list(train_data = train_data, test_data = test_data)
  return(ret_list)
}

setwd("LondonHousePricing/")
library(purrr)
source("../utils.R")

load("LHP_1.RData")
load("LHP_2.RData")
load("LHP_3.RData")
load("LHP_4.RData")

data.frame.igr = as_data_frame(igr, what="both")
mesh = create.mesh.1.5D(nodes = as.matrix(data.frame.igr$vertices)/10^3,
                        edges = as.matrix(data.frame.igr$edges[,1:2]))
FEMbasis = create.FEM.basis(mesh=mesh)

LN.data.adj = LN.data
LN.data.adj$NODES.IDX = nodes.idx
LN.data.adj$DATA.IDX = 1:nrow(LN.data.adj)

# shuffle data
set.seed(27182) # 31415 # 1234
LN.data.adj = LN.data.adj[sample(1:nrow(LN.data.adj)), ]

listData = list()
K = 10 # 10-folds Cross-Validation
num_data = round(nrow(LN.data.adj)/K)
for(i in 1:(K-1)){
  listData[[i]] = LN.data.adj[(1 + num_data*(i-1)):(num_data*i),]
  
}
listData[[K]] = LN.data.adj[(num_data*(K-1) + 1):nrow(LN.data.adj), ]

tmp = get_Kfold_data(listData,iter=9)

# prova 

RMSE.prediction = list( RMSE.fdaPDE = matrix(0, nrow=K, ncol=1),
                        RMSE.GWR.ND = matrix(0, nrow=K,ncol=1))


start=Sys.time()
for(i in 1:K){
  kIter = get_Kfold_data(listData, iter = i, K = K)
  train_data = kIter$train_data
  test_data  = kIter$test_data
  # GWR - ND #
  train_ND = ND[train_data$DATA.IDX, train_data$DATA.IDX] # dMat2
  cross_ND = ND[train_data$DATA.IDX,-train_data$DATA.IDX] # dMat1 ???
  bw.ND = bw.gwr(PURCHASE ~ FLOORSZ + PROF + BATH2, 
                 data = train_data, 
                 approach="AIC", 
                 kernel="gaussian",
                 dMat = train_ND)
  
  GWR.ND = gwr.predict(PURCHASE ~ FLOORSZ + PROF + BATH2, 
                       data = train_data, 
                       predictdata = test_data,
                       kernel = "gaussian",
                       bw = bw.ND,
                       dMat1 = cross_ND,
                       dMat2 = train_ND)
  
  RMSE.prediction$RMSE.GWR.ND[i] = sqrt( mean( (GWR.ND$SDF$prediction - test_data$PURCHASE)^2 ) )
  
  #  fdaPDE  #
  
  observations = train_data$PURCHASE
  locs = cbind(train_data$X, train_data$Y)/10^3 #/10^3
  
  W = cbind( train_data$FLOORSZ, 
             train_data$PROF,   #, 
             train_data$BATH2) #, 
  lambda = 10^seq(from=0.5,to=2.5,by=0.0725) #28
  output_CPP = smooth.FEM(observations = observations, 
                          locations = locs,
                          FEMbasis = FEMbasis,
                          covariates = W,
                          lambda = lambda,
                          lambda.selection.criterion = "grid",
                          lambda.selection.lossfunction = "GCV",
                          DOF.evaluation = "stochastic")
  
  plot(log10(lambda), output_CPP$optimization$GCV_vector)
  
  beta1 = output_CPP$solution$beta[1]
  beta2 = output_CPP$solution$beta[2]
  beta3 = output_CPP$solution$beta[3]
  
  locs_pred = cbind(test_data$X, test_data$Y)/10^3
  #locs_pred = projection.points.1.5D(mesh, cbind(test_data$X, test_data$Y)) #/10^3 )
  prediction = beta1*test_data$FLOORSZ + 
    beta2*test_data$PROF + 
    beta3*test_data$BATH2 +
    eval.FEM(output_CPP$fit.FEM, locs_pred)
  RMSE.prediction$RMSE.fdaPDE[i] = sqrt(mean( (prediction - test_data$PURCHASE)^2 ) )
  
}
end = difftime(Sys.time(), start, units="hours")
tot.time = end

filename = paste("fdaPDE_vs_GWR-", gsub(":","_",gsub(" ","-",Sys.time())), ".RData",sep="")
save(tot.time, RMSE.prediction, file=filename)
