library(ggplot2)

rmse <- function(x,y){
  return( sqrt(mean((x - y)^2 )) )
}

#' @param ND_, Network distance matrix between ALL the nodes 
#' @param ED_, Euclidean distance matrix between ALL the nodes
#' @param observations_, observations at ALL the nodes
#' @param n_sim, number of simulation to be performed 
#' @param n_data, a vector which contains the number of obserservations 
#' to be considered
#' @param lambda, smoothing parameters vector (fdaPDE)
#' @param mesh, fdaPDE mesh 
#' @param true.signal, true.signal = field
#' @param model, logical vec, model[1] -> fdaPDE
#'                            model[2] -> GWR.ND
#'                            model[3] -> GWR.ED
#'                            model[4] -> RR.Krig
NoCovariatesCore <- function(ND_, ED_, observations_,
                             n_sim, n_data, lambda,mesh,
                             true.signal,
                             model_=c(T,T,T,F))
{

RMSE.fdaPDE  = NULL
RMSE.GWR.ND  = NULL
RMSE.lattice  = NULL
RMSE.rr.krig = NULL  
if(model_[1]) RMSE.fdaPDE  = matrix(0,nrow=n_sim, ncol=length(n_data))
if(model_[2]) RMSE.GWR.ND  = matrix(0,nrow=n_sim, ncol=length(n_data))
if(model_[3]) RMSE.lattice  = matrix(0,nrow=n_sim, ncol=length(n_data))
if(model_[4]) RMSE.rr.krig = matrix(0,nrow=n_sim, ncol=length(n_data))

nnodes = nrow(ND_)

mean.field.fdaPDE = matrix(0,nrow=nnodes, ncol=length(n_data))

# lattice method
tmp = as.lattice.fdaPDE(mesh)
nodes.lattice = tmp$nodes.lattice
adj_matrix = tmp$adj_matrix
T_matrix = makeTranMatrix(adj_matrix, M = 0.5)
##############################################

start.time = Sys.time()
for(j in 1:length(n_data)){  
  for(i in 1:n_sim){
    
    sample_ = sample(1:nnodes, size=n_data[j])
    
    locations = mesh$nodes[sample_,]
    observations = observations_[sample_]
    
    ND = ND_[sample_, sample_]
    ED = ED_[sample_, sample_]
    
    ### fdaPDE ### 
    if(model_[1]){
    output_CPP = smooth.FEM(observations = observations, 
                            locations = locations,
                            FEMbasis = FEMbasis,
                            lambda = lambda,
                            lambda.selection.criterion = "grid",
                            lambda.selection.lossfunction = "GCV",
                            DOF.evaluation = "stochastic") # "stochastic"
    
    plot(log10(lambda), output_CPP$optimization$GCV_vector)
    prediction = eval.FEM(output_CPP$fit.FEM, locations = locations)
    RMSE.fdaPDE[i,j] = rmse(true.signal[sample_], prediction)
    
    mean.field.fdaPDE[,j] = mean.field.fdaPDE[,j] + output_CPP$fit.FEM$coeff / n_sim 
    }
    ### GWR ### 
    
    data_ = data.frame(observations = observations)
    Sp.data = SpatialPointsDataFrame(coords = locations,
                                     data = data_)
    # ND
    if(model_[2]){
    bw.ND = bw.gwr(observations ~ 1,
                   data = Sp.data, 
                   approach="AIC", 
                   kernel="gaussian",
                   dMat = ND)
    
    GWR.ND = gwr.basic(observations ~ 1, 
                       data = Sp.data, 
                       kernel = "gaussian",
                       bw = bw.ND,
                       dMat = ND)
    RMSE.GWR.ND[i,j] = rmse(GWR.ND$SDF$yhat,true.signal[sample_])
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
    ### Rank Reduced Kriging (Ver Hoef) ###
    if(model_[4]){
    matNames = as.character(1:nrow(locations))
    rownames(ND) = matNames
    colnames(ND) = matNames
    rownames(ED) = matNames
    colnames(ED) = matNames
    
    knots = create_knots(locations=locations)
    
    RR.krig = rank_reduced_kriging(observations, 
                                   ND, 
                                   knots, 
                                   model=c(F,T,F,F))[[2]]
    
    RMSE.rr.krig[i,j] = rmse(RR.krig$cv.pred, true.signal[sample_])
    }
    
  }
}
  tot.time = difftime(Sys.time(), start.time)
  
  RMSE = list( RMSE.fdaPDE  = RMSE.fdaPDE,
               RMSE.GWR.ND  = RMSE.GWR.ND,
               RMSE.lattice  = RMSE.lattice,
               RMSE.rr.krig = RMSE.rr.krig)
  
  res_ = list(RMSE = RMSE, tot.time=tot.time, mean.field.fdaPDE = mean.field.fdaPDE)
  return(res_)
}

boxplot_RMSE <- function(RMSE, n_data, model_=c(T,T,T,F),
                                       names_=c("fdaPDE", "GWR.ND","GWR.ED", "RR.Krig"), legend.pos = "right"){
 
  M = nrow(RMSE[[1]]) #n_sim
  N = ncol(RMSE[[1]]) #length(n_data) 
  
  model.ticks = vector(mode="character")
  
  n_models = 0
  for(i in 1:length(model_)){
    if(model_[i]){
      n_models = n_models + 1
      model.ticks = append(model.ticks, names_[i])
    }
    
  }
  
  RMSE_ = matrix(nrow=N*M,ncol=n_models)
  col_idx = 0
  for(i in 1:length(names_)){
   if(model_[i]){
    col_idx = col_idx + 1 
    RMSE_[,col_idx] = cbind(as.vector(RMSE[[i]]))  
   }
  }
  
  RMSE_ = as.vector(RMSE_)
  model.ticks = rep(model.ticks, each = N*M)
  
  obs_ = as.character(n_data)
  obs_ = rep(obs_, each=M)
  obs_ = rep(obs_, times=n_models)
  
  data_ = data.frame(RMSE_=RMSE_, model_ = model.ticks, obs_=obs_)
  
  MyTheme <- theme(
    axis.text = element_text(size=26),
    axis.title = element_text(size=26),
    title = element_text(size=26),
    plot.title = element_text(hjust = 0.5),
    legend.text = element_text(size=24),
    legend.key.size = unit(1,"cm"),
    legend.key.height = unit(1,"cm"),
    legend.title = element_blank(),
    legend.background = element_rect(fill="white", color="black",
                                     size=c(1,0.5))
  )
  # legend.position = c(0.85,0.85) # in theme
  
  ggplot(data_) + 
    geom_boxplot(aes(x=obs_, y=RMSE_, 
                     #group=interaction(obs_,model_), 
                     fill=model_))+
    scale_x_discrete(limits=as.character(n_data))+
    labs(x="observations", y="",
         title="RMSE",)+
    MyTheme + 
    theme(plot.title=element_text(hjust=0.5),
          axis.ticks.x = element_blank(),
          legend.position = legend.pos)
  
}
