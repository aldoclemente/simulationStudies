#' @param ND_, Network distance matrix between ALL the nodes 
#' @param ED_, Euclidean distance matrix between ALL the nodes
#' @param observations_, observations at ALL the nodes
#' @param n_sim, number of simulation to be performed 
#' @param n_data, a vector which contains the number of obserservations 
#' to be considered
#' @param lambdaS, LN space smoothing parameters vector (fdaPDE)
#' @param lambdaT, LN time smoothing paramenters vector (fdaDPE)
#' @param lambdaS.2D, 2D space smoothing paramenters vector (fdaPDE)
#' @param lambdaT.2D, 2D time smoothing paramenters vector (fdaPDE)
#' @param mesh, fdaPDE mesh 
#' @param FEMbasis, Linear Networks FEMbasis
#' @param FEMbasis.2D, 2D FEMbasis
#' @param true.signal, signal = field + W beta
#' @param model, logical vec, model[1] -> fdaPDE
#'                            model[2] -> GWR.ND
#'                            model[3] -> GWR.ED
#'                            model[4] -> fdaPDE.2D

library(purrr)
SpaceTimeRegressionCore <- function(ND_, 
                                    ED_=NULL, 
                                    observations,
                                    n_sim, 
                                    n_data, 
                                    lambdaS,
                                    lambdaT,
                                    mesh,
                                    FEMbasis,
                                    FEMbasis.2D=NULL, 
                                    lambdaS.2D=NULL,
                                    lambdaT.2D=NULL,
                                    true.signal,
                                    model_=c(T,T,T,F), 
                                    W = NULL, betas = NULL, time_locations,
                                    legend.pos.RMSE = "right")
{
  
  RMSE.fdaPDE  = NULL
  RMSE.GWR.ND  = NULL
  RMSE.GWR.ED  = NULL
  RMSE.fdaPDE.2D = NULL  
  if(model_[1]) RMSE.fdaPDE  = matrix(0,nrow=n_sim, ncol=length(n_data))
  if(model_[2]) RMSE.GWR.ND  = matrix(0,nrow=n_sim, ncol=length(n_data))
  if(model_[3]) RMSE.GWR.ED  = matrix(0,nrow=n_sim, ncol=length(n_data))
  if(model_[4]) RMSE.fdaPDE.2D = matrix(0,nrow=n_sim, ncol=length(n_data))
  
  
  nnodes_ = nrow(mesh$nodes)
  
  start.time = Sys.time()
  for(j in 1:length(n_data)){  
    for(i in 1:n_sim){
      
      tmp.sample_ = sample(x=(1:nnodes_), size=n_data[j])
      sample_ = vector(mode="integer")
      for(k in 1:length(time_locations)){
        sample_ = cbind(sample_, tmp.sample_ + (k-1)*nnodes_)
      }
      sample_ = as.vector(sample_)
      
      locations = mesh$nodes[tmp.sample_,]
      
      observations_ = observations[sample_]
      
      ND = ND_[sample_, sample_]
      ED = ED_[sample_, sample_]
      
      if(is.null(W)){
        ### fdaPDE ### 
        
        if(model_[1] ){
          output_CPP = smooth.FEM.time(observations = matrix(observations_,
                                                             nrow = n_data[j],
                                                             ncol = length(time_locations)),
                                  time_locations = time_locations,                           
                                  locations = locations,
                                  FEMbasis = FEMbasis,
                                  lambdaS = lambdaS,
                                  lambdaT = lambdaT,
                                  lambda.selection.criterion = "grid",
                                  lambda.selection.lossfunction = "GCV",
                                  DOF.evaluation = "stochastic") 
          
          prediction = eval.FEM.time(FEM.time = output_CPP$fit.FEM.time, 
                                     locations=locations,
                                     time.instants=time_locations,
                                     lambdaS = output_CPP$optimization$lambda_solution[1],
                                     lambdaT = output_CPP$optimization$lambda_solution[2])
          RMSE.fdaPDE[i,j] = rmse(true.signal[sample_], prediction)
        }
        ### GWR ### 
        
        data_ = data.frame(observations = observations_)
        
        Sp.data = SpatialPointsDataFrame(coords = locations,
                                         data = data_)
        # ND
        if(model_[2]){
          bw.ND = bw.gtwr(observations ~ 1,
                          data = Sp.data, 
                          approach="AIC", 
                          kernel="gaussian",
                          st.dMat = ND)
          
          GWR.ND = gtwr(observations ~ 1,
                              data = Sp.data, 
                              kernel = "gaussian",
                              st.bw = bw.ND,
                              st.dMat = ND)
          RMSE.GWR.ND[i,j] = rmse(GWR.ND$SDF$yhat, true.signal[sample_])
        }
        # ED
        if(model_[3]){
          bw.ED = bw.gtwr(observations ~ 1,
                          obs.tv = rep(time_locations,each = n_data[j]),
                          data = Sp.data, 
                          approach="AIC", 
                          kernel="gaussian",
                          st.dMat = ED)
          
          GWR.ED = gtwr(observations ~ 1,
                              obs.tv = rep(time_locations,each = n_data[j]),
                              data = Sp.data, 
                              kernel = "gaussian",
                              st.bw = bw.ED,
                              stMat = ED)
          RMSE.GWR.ED[i,j] = rmse(GWR.ED$SDF$yhat, true.signal[sample_])
        }
        ### fdaPDE 2D model ###
        if(model_[4]){
          output_CPP = smooth.FEM(observations = matrix(observations_,
                                                         nrow = n_data[j],
                                                         ncol = length(time_locations)),
                                  time_locations = time_locations,
                                  locations = locations,
                                  FEMbasis = FEMbasis.2D,
                                  #lambdaS = lambdaS.2D,
                                  #lambdaT = lambdaT.2D,
                                  lambda.selection.criterion = "newton",
                                  lambda.selection.lossfunction = "GCV",
                                  DOF.evaluation = "stochastic") 
          
          prediction = eval.FEM.time(FEM.time = output_CPP$fit.FEM.time, 
                                     locations=locations,
                                     time.instants=time_locations,
                                     lambdaS = output_CPP$optimization$lambda_solution[1],
                                     lambdaT = output_CPP$optimization$lambda_solution[2])
          RMSE.fdaPDE.2D[i,j] = rmse(prediction, true.signal[sample_])
        }
      }
      else{
        if(model_[1] ){
          output_CPP = smooth.FEM.time(observations = matrix(observations_,
                                                        nrow = n_data[j],
                                                        ncol = length(time_locations)),
                                  time_locations = time_locations,
                                  covariates = W[sample_,],
                                  locations = locations,
                                  FEMbasis = FEMbasis,
                                  lambdaS = lambdaS,
                                  lambdaT = lambdaT,
                                  lambda.selection.criterion = "grid",
                                  lambda.selection.lossfunction = "GCV",
                                  DOF.evaluation = "stochastic")
          prediction = 
            eval.FEM.time(FEM.time = output_CPP$fit.FEM.time, 
                          locations=locations,
                          time.instants=time_locations,
                          lambdaS = output_CPP$optimization$lambda_solution[1],
                          lambdaT = output_CPP$optimization$lambda_solution[2]) +
             W[sample_,]%*%as.vector(output_CPP$beta)
          RMSE.fdaPDE[i,j] = rmse(true.signal[sample_], prediction)
        }
        ### GWR ### 
        
        data_ = data.frame(observations = observations_,
                           W1 = W[sample_,1],
                           W2 = W[sample_,2])
        Sp.data = SpatialPointsDataFrame(coords = matrix(rep(locations, 
                                                             times=length(time_locations)),
                                                             nrow = n_data[j]*length(time_locations),
                                                             ncol = 2, byrow = T),
                                         data = data_)
        # ND
        if(model_[2]){
        bw.ND = bw.gtwr(observations ~ W1 + W2,
                          obs.tv = rep(time_locations,each = n_data[j]),
                          data = Sp.data, 
                          approach="AIC", 
                          kernel="gaussian",
                          st.dMat = ND)
          
          GWR.ND = gtwr(observations ~ W1 + W2,
                              obs.tv = rep(time_locations,each = n_data[j]),
                              data = Sp.data, 
                              kernel = "gaussian",
                              st.bw = bw.ND,
                              st.dMat = ND)
          RMSE.GWR.ND[i,j] = rmse(GWR.ND$SDF$yhat, true.signal[sample_])
        }
        # ED
        if(model_[3]){
          bw.ED = bw.gtwr(observations ~ W1 + W2,
                          obs.tv = rep(time_locations,each = n_data[j]),
                          data = Sp.data, 
                          approach="AIC", 
                          kernel="gaussian",
                          st.dMat = ED)
          
          GWR.ED = gtwr(observations ~ W1 + W2,
                              obs.tv = rep(time_locations,each = n_data[j]),
                              data = Sp.data, 
                              kernel = "gaussian",
                              st.bw = bw.ED,
                              st.dMat = ED)
          RMSE.GWR.ED[i,j] = rmse(GWR.ED$SDF$yhat, true.signal[sample_])
        }
        ### fdaPDE 1.5D PARABOLIC ###
        if(model_[4]){
          output_CPP = smooth.FEM.time(observations = matrix(observations_,
                                                             nrow = n_data[j],
                                                             ncol = length(time_locations)),
                                       time_locations = time_locations,
                                       covariates = W[sample_,],
                                       locations = locations,
                                       FEMbasis = FEMbasis,
                                       lambdaS = lambdaS,
                                       lambdaT = lambdaT,
                                       lambda.selection.criterion = "newton",
                                       lambda.selection.lossfunction = "GCV",
                                       DOF.evaluation = "stochastic",
                                       FLAG_PARABOLIC = T) 
          plot(1:length(output_CPP$optimization$GCV_vector), output_CPP$optimization$GCV_vector,
               main="PARABOLIC")
          prediction = eval.FEM.time(FEM.time = output_CPP$fit.FEM.time, 
                                     locations=locations,
                                     time.instants=time_locations,
                                     lambdaS = output_CPP$optimization$lambda_solution[1],
                                     lambdaT = output_CPP$optimization$lambda_solution[2]) +
            W[sample_,]%*%as.vector(output_CPP$beta)
          RMSE.fdaPDE.2D[i,j] = rmse(prediction, true.signal[sample_])
        }
        
      }
      
    }
  }
  tot.time = difftime(Sys.time(), start.time)
  
  RMSE = list( RMSE.fdaPDE  = RMSE.fdaPDE,
               RMSE.GWR.ND  = RMSE.GWR.ND,
               RMSE.GWR.ED  = RMSE.GWR.ED,
               RMSE.fdaPDE.2D = RMSE.fdaPDE.2D)
  
  res_ = list(RMSE = RMSE, tot.time=tot.time)
  return(res_)
}
