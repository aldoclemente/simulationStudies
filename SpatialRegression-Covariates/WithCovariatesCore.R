#' @param ND_, Network distance matrix between ALL the nodes 
#' @param ED_, Euclidean distance matrix between ALL the nodes
#' @param observations_, observations at ALL the nodes
#' @param n_sim, number of simulation to be performed 
#' @param n_data, a vector which contains the number of obserservations 
#' to be considered
#' @param lambda, LN smoothing parameters vector (fdaPDE)
#' @param lambda.2D, 2D smoothing paramenters vector (fdaPDE)
#' @param mesh, fdaPDE mesh 
#' @param FEMbasis, Linear Networks FEMbasis
#' @param FEMbasis.2D, 2D FEMbasis
#' @param true.signal, true.signal = field + W beta
#' @param model, logical vec, model[1] -> fdaPDE
#'                            model[2] -> GWR.ND
#'                            model[3] -> GWR.ED
#'                            model[4] -> fdaPDE.2D
#' @param beta_, true 
#' @param W, covariates matrix                            
WithCovariatesCore.rmse.at.nodes <- function(ND_, ED_, observations_,
                             n_sim, n_data, 
                             lambda,lambda.2D,mesh,
                             FEMbasis,
                             FEMbasis.2D,
                             true.signal,
                             model_=c(T,T,T,F), beta_, W)
{
  
  RMSE.fdaPDE  = NULL
  RMSE.GWR.ND  = NULL
  RMSE.GWR.ED  = NULL
  RMSE.fdaPDE.2D = NULL  
  if(model_[1]) RMSE.fdaPDE  = matrix(0,nrow=n_sim, ncol=length(n_data))
  if(model_[2]) RMSE.GWR.ND  = matrix(0,nrow=n_sim, ncol=length(n_data))
  if(model_[3]) RMSE.GWR.ED  = matrix(0,nrow=n_sim, ncol=length(n_data))
  if(model_[4]) RMSE.fdaPDE.2D = matrix(0,nrow=n_sim, ncol=length(n_data))
    
  
  nnodes = nrow(ND_)
  mean.field.fdaPDE = matrix(0,nrow=nnodes, ncol=length(n_data))
  start.time = Sys.time()
  data_.test = data.frame(observations = observations_,
                          W1 = W[,1],
                          W2 = W[,2])
  Sp.data.test = SpatialPointsDataFrame(coords = mesh$nodes,
                                        data = data_.test)
  
  for(j in 1:length(n_data)){  
    for(i in 1:n_sim){
      
      sample_ = sample(1:nnodes, size=n_data[j])
      
      locations = mesh$nodes[sample_,]
    
      observations = observations_[sample_]
      
      #ND = ND_[sample_, sample_]
      #ED = ED_[sample_, sample_]
      
      ### fdaPDE ### 
      if(model_[1]){
        output_CPP = smooth.FEM(observations = observations,
                                covariates = W[sample_,],
                                locations = locations,
                                FEMbasis = FEMbasis,
                                lambda = lambda,
                                lambda.selection.criterion = "grid",
                                lambda.selection.lossfunction = "GCV",
                                DOF.evaluation = "stochastic") # "stochastic"
        
        plot(log10(lambda), output_CPP$optimization$GCV_vector,
             xlab="log10(lambda)", ylab="", main=paste("GCV:",i,j,sep=" "))
        prediction = eval.FEM(output_CPP$fit.FEM, locations = mesh$nodes) +
                       W%*%output_CPP$solution$beta
        RMSE.fdaPDE[i,j] = rmse(true.signal, prediction)
        mean.field.fdaPDE[,j] = mean.field.fdaPDE[,j] + output_CPP$fit.FEM$coeff / n_sim 
      }
      ### GWR ### 
      
      data_ = data.frame(observations = observations,
                         W1 = W[sample_,1],
                         W2 = W[sample_,2])
      Sp.data = SpatialPointsDataFrame(coords = locations,
                                       data = data_)
      # ND
      if(model_[2]){
        train_ND = ND_[sample_, sample_] #dMat2: a pre-specified sysmetric distance matrix between data points
        cross_ND = ND_[sample_, ]        #dMat1: a pre-specified distance matrix between data points and prediction locations (prediction at nodes!)
        
        bw.ND = bw.gwr(observations ~ W1 + W2, 
                       data = Sp.data, 
                       approach="AIC", 
                       kernel="gaussian",
                       dMat = train_ND)
        
        GWR.ND = gwr.predict(observations ~ W1 + W2, 
                             data = Sp.data, 
                             predictdata = Sp.data.test,
                             kernel = "gaussian",
                             bw = bw.ND,
                             dMat1 = cross_ND,
                             dMat2 = train_ND)
        # bw.ND = bw.gwr(observations ~ W1 + W2,
        #                data = Sp.data, 
        #                approach="AIC", 
        #                kernel="gaussian",
        #                dMat = ND)
        # 
         GWR.ND.basic = gwr.basic(observations ~ W1 + W2, 
                            data = Sp.data, 
                            kernel = "gaussian",
                            bw = bw.ND,
                            dMat = train_ND)
        
        RMSE.GWR.ND[i,j] = rmse(GWR.ND$SDF$prediction,true.signal)
      }
      # ED
      if(model_[3]){
        train_ED = ED_[sample_, sample_] #dMat2: a pre-specified sysmetric distance matrix between data points
        cross_ED = ED_[sample_, ]        #dMat1: a pre-specified distance matrix between data points and prediction locations (prediction at nodes!)
      
        bw.ED = bw.gwr(observations ~ W1 + W2, 
                       data = Sp.data, 
                       approach="AIC", 
                       kernel="gaussian",
                       dMat = train_ED)
        
        GWR.ED = gwr.predict(observations ~ W1 + W2, 
                             data = Sp.data, 
                             predictdata = Sp.data.test,
                             kernel = "gaussian",
                             bw = bw.ED,
                             dMat1 = cross_ED,
                             dMat2 = train_ND)
        # bw.ED = bw.gwr(observations ~ W1 + W2,
        #                data = Sp.data, 
        #                approach="AIC", 
        #                kernel="gaussian",
        #                dMat = ED)
        # 
        # GWR.ED = gwr.basic(observations ~ W1 + W2, 
        #                    data = Sp.data, 
        #                    kernel = "gaussian",
        #                    bw = bw.ED,
        #                    dMat = ED)
         RMSE.GWR.ED[i,j] = rmse(GWR.ED$SDF$prediction,true.signal)
      }
      ### fdaPDE 2D model ###
      if(model_[4]){
        output_CPP = smooth.FEM(observations = observations,
                                covariates = W[sample_,],
                                locations = locations,
                                FEMbasis = FEMbasis.2D,
                                lambda = lambda.2D,
                                lambda.selection.criterion = "grid",
                                lambda.selection.lossfunction = "GCV",
                                DOF.evaluation = "stochastic") # "stochastic"
        
        plot(log10(lambda.2D), output_CPP$optimization$GCV_vector, 
             xlab="log10(lambda)", ylab="", main=paste("GCV 2D:",i,j,sep=" "))
        prediction = eval.FEM(output_CPP$fit.FEM, locations = locations)+
                      W[sample_,]%*%output_CPP$solution$beta
        RMSE.fdaPDE.2D[i,j] = rmse(prediction, true.signal)
      }
      
    }
  }
  tot.time = difftime(Sys.time(), start.time)
  
  RMSE = list( RMSE.fdaPDE  = RMSE.fdaPDE,
               RMSE.GWR.ND  = RMSE.GWR.ND,
               RMSE.GWR.ED  = RMSE.GWR.ED,
               RMSE.fdaPDE.2D = RMSE.fdaPDE.2D)
  
  res_ = list(RMSE = RMSE, tot.time=tot.time, mean.field.fdaPDE=mean.field.fdaPDE)
  return(res_)
}

#' @param ND_, Network distance matrix between ALL the nodes 
#' @param ED_, Euclidean distance matrix between ALL the nodes
#' @param observations_, observations at ALL the nodes
#' @param n_sim, number of simulation to be performed 
#' @param n_data, a vector which contains the number of obserservations 
#' to be considered
#' @param lambda, LN smoothing parameters vector (fdaPDE)
#' @param lambda.2D, 2D smoothing paramenters vector (fdaPDE)
#' @param mesh, fdaPDE mesh 
#' @param FEMbasis, Linear Networks FEMbasis
#' @param FEMbasis.2D, 2D FEMbasis
#' @param true.signal, true.signal = field + W beta
#' @param model, logical vec, model[1] -> fdaPDE
#'                            model[2] -> GWR.ND
#'                            model[3] -> GWR.ED
#'                            model[4] -> fdaPDE.2D
#' @param beta_, true 
#' @param W, covariates matrix                            
WithCovariatesCore <- function(ND_, ED_, observations_,
                               n_sim, n_data, 
                               lambda,lambda.2D,mesh,
                               FEMbasis,
                               FEMbasis.2D,
                               true.signal,
                               model_=c(T,T,T,F), beta_, W)
{
  
  RMSE.fdaPDE  = NULL
  RMSE.GWR.ND  = NULL
  RMSE.lattice  = NULL
  RMSE.fdaPDE.2D = NULL  
  if(model_[1]) RMSE.fdaPDE  = matrix(0,nrow=n_sim, ncol=length(n_data))
  if(model_[2]) RMSE.GWR.ND  = matrix(0,nrow=n_sim, ncol=length(n_data))
  if(model_[3]) RMSE.lattice  = matrix(0,nrow=n_sim, ncol=length(n_data))
  if(model_[4]) RMSE.fdaPDE.2D = matrix(0,nrow=n_sim, ncol=length(n_data))
  
  
  nnodes = nrow(ND_)
  mean.field.fdaPDE = matrix(0,nrow=nnodes, ncol=length(n_data))
  start.time = Sys.time()
  data_.test = data.frame(observations = observations_,
                          W1 = W[,1],
                          W2 = W[,2])
  # lattice method
  tmp = as.lattice.fdaPDE(mesh)
  nodes.lattice = tmp$nodes.lattice
  adj_matrix = tmp$adj_matrix
  T_matrix = makeTranMatrix(adj_matrix, M = 0.5)
  ##############################################
  estimates = list()
  
  for(j in 1:length(n_data)){  
    for(i in 1:n_sim){
      
      sample_ = sample(1:nnodes, size=n_data[j])
      
      locations = mesh$nodes[sample_,]
      
      observations = observations_[sample_]
      
      if(i==1 && j == 1){ 
        estimates$locations = locations
        estimates$observations = observations
        estimates$true.signal = true.signal[sample_]
        estimates$X1 = W[sample_,1]
        estimates$X2 = W[sample_,2]
      }
      
      #ND = ND_[sample_, sample_]
      #ED = ED_[sample_, sample_]
      
      ### fdaPDE ### 
      if(model_[1]){
        output_CPP = smooth.FEM(observations = observations,
                                covariates = W[sample_,],
                                locations = locations,
                                FEMbasis = FEMbasis,
                                lambda = lambda,
                                lambda.selection.criterion = "grid",
                                lambda.selection.lossfunction = "GCV",
                                DOF.evaluation = "stochastic") # "stochastic"
        
        plot(log10(lambda), output_CPP$optimization$GCV_vector,
             xlab="log10(lambda)", ylab="", main=paste("GCV:",i,j,sep=" "))
        prediction = eval.FEM(output_CPP$fit.FEM, locations = locations) +
          W[sample_,]%*%output_CPP$solution$beta
        RMSE.fdaPDE[i,j] = rmse(true.signal[sample_], prediction)
        mean.field.fdaPDE[,j] = mean.field.fdaPDE[,j] + output_CPP$fit.FEM$coeff / n_sim
        
        if(i==1 && j == 1){ 
          estimates$fdaPDE = prediction 
        }
      }
      ### GWR ### 
      
      data_ = data.frame(observations = observations,
                         W1 = W[sample_,1],
                         W2 = W[sample_,2])
      Sp.data = SpatialPointsDataFrame(coords = locations,
                                       data = data_)
      # ND
      if(model_[2]){
        train_ND = ND_[sample_, sample_] #dMat: a pre-specified sysmetric distance matrix between data points
        
        bw.ND = bw.gwr(observations ~ W1 + W2, 
                       data = Sp.data, 
                       approach="AIC", 
                       kernel="gaussian",
                       dMat = train_ND)
        
        GWR.ND = gwr.basic(observations ~ W1 + W2, 
                                 data = Sp.data, 
                                 kernel = "gaussian",
                                 bw = bw.ND,
                                 dMat = train_ND)
        
        RMSE.GWR.ND[i,j] = rmse(GWR.ND$SDF$yhat,true.signal[sample_])
        
        if(i==1 && j == 1){ 
          estimates$GWR = GWR.ND$SDF$yhat 
        }
      }
      # lattice based method
      if(model_[3]){
        # dummy z-coordinates !
        locations.lattice = cbind(locations, rep(0, times=n_data[j]))
        output_press = crossvalNparSmoother(T=T_matrix,
                                            nodelocs = nodes.lattice,
                                            locs = locations.lattice,
                                            Z = observations,
                                            k_max=1500)
        output_lattice = nparSmoother(T=T_matrix,
                                      nodelocs = nodes.lattice,
                                      locs = locations.lattice,
                                      Z = observations,
                                      k = output_press$k)
        prediction.latt = eval.FEM(FEM(output_lattice[,4], FEMbasis), locations=locations)  
        RMSE.lattice[i,j] = rmse(true.signal[sample_], prediction.latt)
      }
      ### fdaPDE 2D model ###
      if(model_[4]){
        output_CPP = smooth.FEM(observations = observations,
                                covariates = W[sample_,],
                                locations = locations,
                                FEMbasis = FEMbasis.2D,
                                lambda = lambda.2D,
                                lambda.selection.criterion = "grid",
                                lambda.selection.lossfunction = "GCV",
                                DOF.evaluation = "stochastic") # "stochastic"
        
        plot(log10(lambda.2D), output_CPP$optimization$GCV_vector, 
             xlab="log10(lambda)", ylab="", main=paste("GCV 2D:",i,j,sep=" "))
        prediction = eval.FEM(output_CPP$fit.FEM, locations = locations)+
          W[sample_,]%*%output_CPP$solution$beta
        RMSE.fdaPDE.2D[i,j] = rmse(prediction, true.signal[sample_])
      }
      
    }
  }
  tot.time = difftime(Sys.time(), start.time)
  
  RMSE = list( RMSE.fdaPDE  = RMSE.fdaPDE,
               RMSE.GWR.ND  = RMSE.GWR.ND,
               RMSE.lattice  = RMSE.lattice,
               RMSE.fdaPDE.2D = RMSE.fdaPDE.2D)
  
  res_ = list(RMSE = RMSE, tot.time=tot.time, mean.field.fdaPDE=mean.field.fdaPDE, estimates = estimates)
  return(res_)
}


library(aod)
library(Matrix)
BetaInference <- function(FEMbasis, 
                          nobs, 
                          sample_, 
                          covariates, 
                          lambda, 
                          observations, 
                          fitted.values,
                          beta.hat){
  
  # Inference example #
  nnodes = nrow(FEMbasis$mesh$nodes)
  
  R1 = CPP_get.FEM.Stiff.Matrix(FEMbasis = FEMbasis)
  R0 = CPP_get.FEM.Mass.Matrix(FEMbasis = FEMbasis)
  
  P = t(R1)%*%solve(R0)%*%R1
  X = covariates # n x p
  Psi = Matrix(data=0, nrow = nobs , ncol = nnodes, sparse =TRUE) # n x N
  Psi <- as(Psi, "dgCMatrix")
  for(i in 1:nobs){
    Psi[i, sample_[i]] = 1
  }
  
  I = Matrix(data=0, nrow=nobs, ncol=nobs, sparse = TRUE)
  I <- as(I, "dgCMatrix")
  for(i in 1:nobs){
    I[i,i] = 1
  }
  
  invXtX = solve(t(X)%*%X) # (p x n) X (n x p) = p x p 

  Q = I - X %*% invXtX %*% t(X) # n x n
  Q <- as(Q, "dgCMatrix")
  S = Psi%*%solve(t(Psi)%*%Q%*%Psi + lambda*P)%*%t(Psi)%*%Q # (n x N) [(N x n) (n x n) (n x N) + N x N ]  (N x n) (n x n) = n x n
  z_hat =  fitted.values # n
  sigma2_hat = 1/(nnodes - ncol(X) - sum(diag(S)))*t((z_hat - observations))%*%(z_hat - observations)
  
  #SIGMA_hat = sigma2_hat * (I-Q-Q%*%S)%*%t((I-Q-Q%*%S))
  
  SIGMA_beta = sigma2_hat * ( invXtX + invXtX%*%t(X)%*%S%*%t(S)%*%X%*%invXtX) # p x p 
  print("Complete Wald Test:")
  print(wald.test(Sigma=SIGMA_beta, b = beta.hat , Terms = 1:ncol(X))) 
  
  for(i in 1:ncol(X)){
    print(paste("Wald Test for :", i , " coeff", sep=""))
    print(wald.test(Sigma=SIGMA_beta, b = beta.hat , Terms = i))
  }
  
}

