#' @param ND_, Network distance matrix between ALL the nodes 
#' @param observations_, observations at ALL the nodes
#' @param n_sim, number of simulation to be performed 
#' @param n_data, a vector which contains the number of obserservations 
#' to be considered
#' @param lambdaS, LN space smoothing parameters vector (fdaPDE)
#' @param lambdaT, LN time smoothing paramenters vector (fdaDPE)
#' @param mesh, fdaPDE mesh 
#' @param FEMbasis, Linear Networks FEMbasis
#' @param true.signal, signal = field + W beta

library(purrr)
SpaceTimeGLRCore <- function(ND_, 
                         observations,
                         n_sim, 
                         n_data, 
                         lambdaS,
                         lambdaT,
                         mesh,
                         FEMbasis,
                         true.signal,
                         true.field,
                         W = NULL, betas = NULL, time_locations,
                         FAMILY,
                         legend.pos.RMSE = "right")
{
  
  
  RMSE.fdaPDE  = matrix(0,nrow=n_sim, ncol=length(n_data))
  
  nnodes_ = nrow(mesh$nodes)
  mean.field.fdaPDE = matrix(0,nrow=nnodes*length(time_locations), ncol=length(n_data))
  
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
      
      if(is.null(W)){
        ### fdaPDE ### 
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
                                       DOF.evaluation = "stochastic",
                                     family = FAMILY) 
          
          prediction = eval.FEM.time(FEM.time = output_CPP$fit.FEM.time, 
                                     locations=mesh$nodes,
                                     time.instants=time_locations,
                                     lambdaS = output_CPP$optimization$lambda_solution[1],
                                     lambdaT = output_CPP$optimization$lambda_solution[2])
          RMSE.fdaPDE[i,j] = rmse(true.field, prediction)
          
          mean.field.fdaPDE[,j] = mean.field.fdaPDE[,j] + eval.FEM.time(FEM.time = output_CPP$fit.FEM.time, 
                                                                        locations=mesh$nodes,
                                                                        time.instants=time_locations,
                                                                        lambdaS = output_CPP$optimization$lambda_solution[1],
                                                                        lambdaT = output_CPP$optimization$lambda_solution[2])/n_sim
      }
      else{
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
                                       DOF.evaluation = "stochastic",
                                       family = FAMILY)
          prediction = 
            eval.FEM.time(FEM.time = output_CPP$fit.FEM.time, 
                          locations=locations,
                          time.instants=time_locations,
                          lambdaS = output_CPP$optimization$lambda_solution[1],
                          lambdaT = output_CPP$optimization$lambda_solution[2]) +
            W[sample_,]%*%as.vector(output_CPP$beta)
          RMSE.fdaPDE[i,j] = rmse(true.signal[sample_], prediction)
          
          mean.field.fdaPDE[,j] = mean.field.fdaPDE[,j] + eval.FEM.time(FEM.time = output_CPP$fit.FEM.time, 
                                                                        locations=mesh$nodes,
                                                                        time.instants=time_locations,
                                                                        lambdaS = output_CPP$optimization$lambda_solution[1],
                                                                        lambdaT = output_CPP$optimization$lambda_solution[2])/n_sim
        
      }
      
    }
  }
  tot.time = difftime(Sys.time(), start.time)
  
  RMSE = list( RMSE.fdaPDE  = RMSE.fdaPDE)
  
  res_ = list(RMSE = RMSE, tot.time=tot.time, mean.field.fdaPDE=mean.field.fdaPDE)
  return(res_)
}
