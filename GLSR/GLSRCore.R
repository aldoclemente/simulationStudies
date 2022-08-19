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
#' @param true.signal, true.signal = field
#'                     true.signal = field + W beta
#' @param model, logical vec, model[1] -> fdaPDE
#'                            model[2] -> GWR.ND
#'                            model[3] -> GWR.ED
#'                            model[4] -> fdaPDE.2D
#' @param family, "binomail" or "poisson" NB) gwr can handles ONLY these cases

library(purrr)
gamCore <- function(ND_, ED_, observations_,
                               n_sim, n_data, 
                               lambda,mesh,
                               FEMbasis,
                               FEMbasis.2D=NULL, lambda.2D=NULL,
                               true.signal,
                               model_=c(T,T,T,F), 
                               W = NULL, betas = NULL,
                               FAMILY="poisson")
{
  if(FAMILY=="poisson"){
    l<-make.link("log")
    link<-l$linkfun
    inv.link<-l$linkinv
  }else if(FAMILY == "binomial"){
    logit <- function(x){qlogis(x)}
    inv.logit <- function(x){plogis(x)}
    link = logit
    inv.link = inv.logit
  }else{
    stop("FAMILY must be 'poisson' or 'binomial. Try again.")
  }
  RMSE.fdaPDE  = NULL
  RMSE.GWR.ND  = NULL
  RMSE.GWR.ED  = NULL
  RMSE.fdaPDE.2D = NULL  
  if(model_[1]) RMSE.fdaPDE  = matrix(0,nrow=n_sim, ncol=length(n_data))
  if(model_[2]) RMSE.GWR.ND  = matrix(0,nrow=n_sim, ncol=length(n_data))
  if(model_[3]) RMSE.GWR.ED  = matrix(0,nrow=n_sim, ncol=length(n_data))
  if(model_[4]) RMSE.fdaPDE.2D = matrix(0,nrow=n_sim, ncol=length(n_data))
  
  estimates = list()
  
  nnodes = nrow(ND_)
  mean.field.fdaPDE = matrix(0,nrow=nnodes, ncol=length(n_data))
  start.time = Sys.time()
  for(j in 1:length(n_data)){  
    for(i in 1:n_sim){
      
      sample_ = sample(1:nnodes, size=n_data[j])
      
      locations = mesh$nodes[sample_,]
      
      observations = observations_[sample_]  
      
      ND = ND_[sample_, sample_]
      ED = ED_[sample_, sample_]
      
      if(i==1 && j == 1){ 
        estimates$locations = locations
        estimates$observations = observations
        estimates$true.signal = inv.link(true.signal[sample_])
        estimates$X1 = W[sample_,1]
        estimates$X2 = W[sample_,2]
      }
      
      if(is.null(W)){
      ### fdaPDE ### 
      if(model_[1] ){
        output_CPP = smooth.FEM(observations = observations,
                                locations = locations,
                                FEMbasis = FEMbasis,
                                lambda = lambda,
                                lambda.selection.criterion = "grid",
                                lambda.selection.lossfunction = "GCV",
                                DOF.evaluation = "stochastic",
                                family = FAMILY) 
        
        plot(log10(lambda), output_CPP$optimization$GCV_vector,
             xlab="log10(lambda)", ylab="", main=paste("GCV:",i,j,sep=" "))
        prediction = inv.link(eval.FEM(FEM(output_CPP$solution$f, FEMbasis), 
                              locations = locations))
        RMSE.fdaPDE[i,j] = rmse(inv.link(true.signal[sample_]), prediction)
        mean.field.fdaPDE[,j] = mean.field.fdaPDE[,j] + output_CPP$solution$f / n_sim
        if(i==1 && j == 1){ 
          estimates$fdaPDE = prediction 
        }
        
        
        
        
      }
      ### GWR ### 
      
      data_ = data.frame(observations = observations)
      
      Sp.data = SpatialPointsDataFrame(coords = locations,
                                       data = data_)
      # ND
      if(model_[2]){
        bw.ND = bw.ggwr(observations ~ 1,
                       family=FAMILY,
                       data = Sp.data, 
                       approach="AIC", 
                       kernel="gaussian",
                       dMat = ND)
        
        GWR.ND = ggwr.basic(observations ~ 1,
                           family = FAMILY,
                           data = Sp.data, 
                           kernel = "gaussian",
                           bw = bw.ND,
                           dMat = ND)
        RMSE.GWR.ND[i,j] = rmse(GWR.ND$SDF$yhat, inv.link(true.signal[sample_]))
        
        
        if(i==1 && j == 1){ 
          estimates$GWR = GWR.ND$SDF$yhat 
        }
      }
      # ED
      if(model_[3]){
        bw.ED = bw.ggwr(observations ~ 1,
                       family = FAMILY,
                       data = Sp.data, 
                       approach="AIC", 
                       kernel="gaussian",
                       dMat = ED)
        
        GWR.ED = ggwr.basic(observations ~ 1,
                            family = FAMILY,
                           data = Sp.data, 
                           kernel = "gaussian",
                           bw = bw.ED,
                           dMat = ED)
        RMSE.GWR.ED[i,j] = rmse(GWR.ED$SDF$yhat, inv.link(true.signal[sample_]))
      }
      ### fdaPDE 2D model ###
      if(model_[4]){
        output_CPP = smooth.FEM(observations = observations,
                                locations = locations,
                                FEMbasis = FEMbasis.2D,
                                lambda = lambda.2D,
                                lambda.selection.criterion = "grid",
                                lambda.selection.lossfunction = "GCV",
                                DOF.evaluation = "stochastic",
                                family = FAMILY) 
        
        plot(log10(lambda.2D), output_CPP$optimization$GCV_vector, 
             xlab="log10(lambda)", ylab="", main=paste("GCV 2D:",i,j,sep=" "))
        prediction = inv.link(eval.FEM(FEM(output_CPP$solution$f, FEMbasis.2D), 
                              locations = locations))  
        RMSE.fdaPDE.2D[i,j] = rmse(prediction, inv.link(true.signal[sample_]))
      }
      }
      else{
        if(model_[1] ){
          output_CPP = smooth.FEM(observations = observations,
                                  covariates = W[sample_,],
                                  locations = locations,
                                  FEMbasis = FEMbasis,
                                  lambda = lambda,
                                  lambda.selection.criterion = "grid",
                                  lambda.selection.lossfunction = "GCV",
                                  DOF.evaluation = "stochastic",
                                  family = FAMILY,
                                  max.steps.FPIRLS = 50) 
          
          plot(log10(lambda), output_CPP$optimization$GCV_vector,
               xlab="log10(lambda)", ylab="", main=paste("GCV:",i,j,sep=" "))
          prediction = inv.link(
            eval.FEM(FEM(output_CPP$solution$f[,output_CPP$optimization$lambda_position], FEMbasis), 
                                locations = locations) 
                       + W[sample_,]%*%output_CPP$solution$beta[,output_CPP$optimization$lambda_position])
          RMSE.fdaPDE[i,j] = rmse(inv.link(true.signal[sample_]), prediction)
          mean.field.fdaPDE[,j] = mean.field.fdaPDE[,j] + output_CPP$solution$f[,output_CPP$optimization$lambda_position] / n_sim
        
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
          bw.ND = bw.ggwr(observations ~ W1 + W2,
                          family=FAMILY,
                          data = Sp.data, 
                          approach="AIC", 
                          kernel="gaussian",
                          dMat = ND)
          
          GWR.ND = ggwr.basic(observations ~ W1 + W2,
                              family = FAMILY,
                              data = Sp.data, 
                              kernel = "gaussian",
                              bw = bw.ND,
                              dMat = ND)
          RMSE.GWR.ND[i,j] = rmse(GWR.ND$SDF$yhat, inv.link(true.signal[sample_]))
          
          if(i==1 && j == 1){ 
            estimates$GWR = GWR.ND$SDF$yhat 
          }
        }
        # ED
        if(model_[3]){
          bw.ED = bw.ggwr(observations ~ W1 + W2,
                          family = FAMILY,
                          data = Sp.data, 
                          approach="AIC", 
                          kernel="gaussian",
                          dMat = ED)
          
          GWR.ED = ggwr.basic(observations ~ W1 + W2,
                              family = FAMILY,
                              data = Sp.data, 
                              kernel = "gaussian",
                              bw = bw.ED,
                              dMat = ED)
          RMSE.GWR.ED[i,j] = rmse(GWR.ED$SDF$yhat, inv.link(true.signal[sample_]))
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
                                  DOF.evaluation = "stochastic",
                                  family = FAMILY) 
          
          plot(log10(lambda.2D), output_CPP$optimization$GCV_vector, 
               xlab="log10(lambda)", ylab="", main=paste("GCV 2D:",i,j,sep=" "))
          prediction = inv.link(
                      eval.FEM(FEM(output_CPP$solution$f[,output_CPP$optimization$lambda_position], FEMbasis.2D), 
                               locations = locations) 
                      + W[sample_,]%*%output_CPP$solution$beta[,output_CPP$optimization$lambda_position])
          RMSE.fdaPDE.2D[i,j] = rmse(prediction, inv.link(true.signal[sample_]))
        }
        
      }
      
    }
  }
  tot.time = difftime(Sys.time(), start.time)
  
  RMSE = list( RMSE.fdaPDE  = RMSE.fdaPDE,
               RMSE.GWR.ND  = RMSE.GWR.ND,
               RMSE.GWR.ED  = RMSE.GWR.ED,
               RMSE.fdaPDE.2D = RMSE.fdaPDE.2D)
  
  res_ = list(RMSE = RMSE, tot.time=tot.time, mean.field.fdaPDE=mean.field.fdaPDE, estimates=estimates)
  return(res_)
}
