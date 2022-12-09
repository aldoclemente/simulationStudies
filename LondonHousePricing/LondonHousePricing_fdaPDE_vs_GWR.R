############################
### LONDON HOUSE PRICING ###
#####  fdaPDE vs GWR  ######
############################
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
load("LHP_preliminaries_1.RData")
load("LHP_preliminaries_2.RData")
load("LHP_preliminaries_3.RData")

mesh = create.mesh.1.5D(nodes = as.matrix(data.frame.igr$vertices)/10^3,
                        edges = as.matrix(data.frame.igr$edges[,1:2]))
FEMbasis = create.FEM.basis(mesh=mesh)

LN.data.adj = LN.data
LN.data.adj$NODES.IDX = nodes.idx
LN.data.adj$DATA.IDX = 1:nrow(LN.data.adj)

# shuffle data
set.seed(27182) 
LN.data.adj = LN.data.adj[sample(1:nrow(LN.data.adj)), ]

listData = list()
K = 10 # 10-folds Cross-Validation
num_data = round(nrow(LN.data.adj)/K)
for(i in 1:(K-1)){
  listData[[i]] = LN.data.adj[(1 + num_data*(i-1)):(num_data*i),]
  
}
listData[[K]] = LN.data.adj[(num_data*(K-1)+1):nrow(LN.data.adj), ]


# prova 
tmp = get_Kfold_data(listData, 1)

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
  
  W = cbind( train_data$FLOORSZ, train_data$PROF, train_data$BATH2) #, 
  lambda = 10^seq(from=1,to=2.5,by=0.0725) #28
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
  RMSE.prediction$RMSE.fdaPDE[i] = sqrt(mean( (prediction - LN.test_data$PURCHASE)^2 ) )
  }
end = difftime(Sys.time, start, units="hours")
tot.time=end

filename = paste("fdaPDE_vs_GWR-LocsProj", gsub(":","_",gsub(" ","-",Sys.time())), ".RData",sep="")
save(tot.time, RMSE.prediction, file=filename)

###

load("fdaPDE_vs_GWR.RData")
library(ggplot2)

RMSE_ = rbind(RMSE.prediction$RMSE.fdaPDE, 
              RMSE.prediction$RMSE.GWR.ND, 
              RMSE.prediction$RMSE.GWR.ED)

model_ = c("fdaPDE", "GWR.ND","GWR.ED")
model_ = rep(model_, each= length(RMSE.prediction$RMSE.fdaPDE))
data_ = data.frame(RMSE_=RMSE_, model_ = model_)

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
  geom_boxplot(aes( y=RMSE_, fill=model_))+
  labs(x="", y="",
       title="RMSE (predicted)",)+
  MyTheme + 
  theme(plot.title=element_text(hjust=0.5),
        axis.text.x = element_blank(),
        axis.ticks.x = element_blank())
  
data(LNHP)
class(LN.prop)
library(maptools)
library(igraph)
library(shp2graph)
library(GWmodel)
data(LondonBorough) # library(GWmodel)
library(RColorBrewer)
mypalette<-brewer.pal(6,"Spectral")
#x11()
nsa = list("SpatialPolygonsRescale", layout.north.arrow(), offset = c(561900,200900),
           scale = 500, col=1)
spplot(LN.prop, "PURCHASE", key.space = "right", cex=1.5, cuts=5,
       ylim=c(155840.8,200933.9), xlim=c(503568.2,561957.5),
       main="TRUE PURCHASE",
       col.regions=mypalette, sp.layout=list(nsa, londonborough))

# ok Facciamo sti grafici esempio #
{
  
  load("LondonHousePricing_3.RData")
  igr = data.list.3$igr
  nodes_idx = data.list.3$nodes_idx
  
  # Network Distance matrix
  ND = data.list.3$ND
  
  locations = data.list.3$locations
  LN.data = data.list.3$LN.data
  
  # fdaPDE mesh #
  data(LNNT)
  data(LNHP)
  London.Road.Networks = maptools::as.linnet.SpatialLines(LN.nt)
  London.Point.Pattern= maptools::as.ppp.SpatialPointsDataFrame(LN.data)
  
  London.Road.Networks.2 = rescale(X = London.Road.Networks, s = 1e3, 
                                   unitname=c("kilometer","kilometer") )
  London.Point.Pattern.2 = rescale(X=London.Point.Pattern, s=1e3,
                                   unitname = c("kilometer","kilometer"))
  
  nodes = cbind(London.Road.Networks.2$vertices$x, 
                London.Road.Networks.2$vertices$y)
  edges = cbind(London.Road.Networks.2$from, 
                London.Road.Networks.2$to)
  mesh = create.mesh.1.5D(nodes=nodes, edges=edges)
  FEMbasis = create.FEM.basis(mesh=mesh)
  ###############
  
  # tot data 
  n_data = nrow(ND)
  
  # tot training set
  n_train = round(0.9*n_data)
  # tot test set
  n_test = nrow(ND) - n_train
  
  RMSE.fdaPDE = 0.
  RMSE.GWR.ND = 0.
  RMSE.GWR.ED = 0.
  
  data_idx = 1:n_data
  
  sample_ = sample(data_idx,size=n_train)
  train_idx = nodes_idx[sample_]
  test_idx = nodes_idx[-sample_]
  
  # training data
  
  train_data = data.frame(PURCHASE   = LN.data$PURCHASE[sample_],
                          FLOORSZ    = LN.data$FLOORSZ[sample_],
                          TYPEDETCH  = LN.data$TYPEDETCH[sample_], 
                          TPSEMIDTCH = LN.data$TPSEMIDTCH[sample_], 
                          TYPETRRD   = LN.data$TYPETRRD[sample_],
                          TYPEBNGLW  = LN.data$TYPEBNGLW[sample_],
                          TYPEFLAT   = LN.data$TYPEFLAT[sample_],
                          BLDPWW1    = LN.data$BLDPWW1[sample_],
                          BLDPOSTW   = LN.data$BLDPOSTW[sample_],
                          BLD60S     = LN.data$BLD60S[sample_],
                          BLD70S     = LN.data$BLD70S[sample_],
                          BLD80S     = LN.data$BLD80S[sample_],
                          BLD90S     = LN.data$BLD90S[sample_],
                          BATH2      = LN.data$BATH2[sample_],
                          BED2       = LN.data$BED2[sample_],
                          GARAGE1    = LN.data$GARAGE1[sample_],
                          CENTHEAT   = LN.data$CENTHEAT[sample_],
                          UNEMPLOY   = LN.data$UNEMPLOY[sample_],
                          PROF       = LN.data$PROF[sample_],
                          BLDINTW    = LN.data$BLDINTW[sample_],
                          X          = LN.data$X[sample_],
                          Y          = LN.data$Y[sample_]
  )
  test_data = data.frame(PURCHASE   = LN.data$PURCHASE[-sample_],
                         FLOORSZ    = LN.data$FLOORSZ[-sample_],
                         TYPEDETCH  = LN.data$TYPEDETCH[-sample_], 
                         TPSEMIDTCH = LN.data$TPSEMIDTCH[-sample_], 
                         TYPETRRD   = LN.data$TYPETRRD[-sample_],
                         TYPEBNGLW  = LN.data$TYPEBNGLW[-sample_],
                         TYPEFLAT   = LN.data$TYPEFLAT[-sample_],
                         BLDPWW1    = LN.data$BLDPWW1[-sample_],
                         BLDPOSTW   = LN.data$BLDPOSTW[-sample_],
                         BLD60S     = LN.data$BLD60S[-sample_],
                         BLD70S     = LN.data$BLD70S[-sample_],
                         BLD80S     = LN.data$BLD80S[-sample_],
                         BLD90S     = LN.data$BLD90S[-sample_],
                         BATH2      = LN.data$BATH2[-sample_],
                         BED2       = LN.data$BED2[-sample_],
                         GARAGE1    = LN.data$GARAGE1[-sample_],
                         CENTHEAT   = LN.data$CENTHEAT[-sample_],
                         UNEMPLOY   = LN.data$UNEMPLOY[-sample_],
                         PROF       = LN.data$PROF[-sample_],
                         BLDINTW    = LN.data$BLDINTW[-sample_],
                         X          = LN.data$X[-sample_],
                         Y          = LN.data$Y[-sample_]
  )
  train_coords = cbind(train_data$X, train_data$Y)
  test_coords = cbind(test_data$X, test_data$Y)
  LN.train_data = SpatialPointsDataFrame(coords = train_coords, data = train_data)
  LN.test_data =  SpatialPointsDataFrame(coords = test_coords,  data = test_data)
  
  # GWR - ND #
  train_ND = ND[sample_, sample_] # dMat2
  cross_ND = ND[sample_, -sample_]  # dMat1 ???
  bw.ND = bw.gwr(PURCHASE ~ FLOORSZ + PROF + BATH2, 
                 data = LN.train_data, 
                 approach="AIC", 
                 kernel="gaussian",
                 dMat = train_ND)
  
  GWR.ND = gwr.predict(PURCHASE ~ FLOORSZ + PROF + BATH2, 
                       data = LN.train_data, 
                       predictdata = LN.test_data,
                       kernel = "gaussian",
                       bw = bw.ND,
                       dMat1 = cross_ND,
                       dMat2 = train_ND)
  
  RMSE.GWR.ND = sqrt( mean( (GWR.ND$SDF$prediction - LN.test_data$PURCHASE)^2 ) )
  # GWR - ED #
  
  bw.ED = bw.gwr(PURCHASE ~ FLOORSZ + PROF + BATH2, 
                 data = LN.train_data, 
                 approach="AIC", 
                 kernel="gaussian")
  
  GWR.ED = gwr.predict(PURCHASE ~ FLOORSZ + PROF + BATH2, 
                       data = LN.train_data,
                       predictdata = LN.test_data,
                       kernel = "gaussian",
                       bw = bw.ED)
  
  RMSE.GWR.ED = sqrt( mean( (GWR.ED$SDF$prediction - LN.test_data$PURCHASE)^2 ) )
  
  #  fdaPDE  (with intercept) #
  
  observations = train_data$PURCHASE
  locs = cbind(train_data$X, train_data$Y)/10^3
  locs = projection.points.1.5D(mesh,locs)
  
  W = cbind( train_data$FLOORSZ, train_data$PROF, train_data$BATH2)
  lambda = 10^seq(from=0,to=2,by=0.0725) #28
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
  
  locs_pred = projection.points.1.5D(mesh, cbind(test_data$X, test_data$Y)/10^3 )
  prediction = beta1*test_data$FLOORSZ + 
    beta2*test_data$PROF + 
    beta3*test_data$BATH2 +
    eval.FEM(output_CPP$fit.FEM, locs_pred)
  RMSE.fdaPDE = sqrt(mean( (prediction - LN.test_data$PURCHASE)^2 ) )
  sp.fdaPDE = SpatialPointsDataFrame(coords=cbind(test_data$X, test_data$Y),
                                     data=data.frame(prediction=prediction))
}

  pdf(file="fdaPDE_vs_GWR_imgs.pdf")
  
  plot(mesh,type="n")
  points(LN.data$X/10^3, LN.data$Y/10^3, pch=16, col="red",cex=0.75)
  
  ggplot(data_) + 
    geom_boxplot(aes( y=RMSE_, fill=model_))+
    labs(x="", y="",
         title="RMSE (predicted)",)+
    MyTheme + 
    theme(plot.title=element_text(hjust=0.5),
          axis.text.x = element_blank(),
          axis.ticks.x = element_blank())
  
  spplot(LN.test_data, "PURCHASE", key.space = "right", cex=1.5, cuts=6,
         ylim=c(155840.8,200933.9), xlim=c(503568.2,561957.5),
         main="TRUE PURCHASE",
         col.regions=mypalette, sp.layout=list(nsa, londonborough))
  
  spplot(GWR.ND$SDF, "prediction", key.space = "right", cex=1.5, cuts=6,
         ylim=c(155840.8,200933.9), xlim=c(503568.2,561957.5),
         main="GWR.ND PURCHASE predictions",
         col.regions=mypalette, sp.layout=list(nsa, londonborough))
  
  spplot(GWR.ED$SDF, "prediction", key.space = "right", cex=1.5, cuts=6,
         ylim=c(155840.8,200933.9), xlim=c(503568.2,561957.5),
         main="GWR.ED PURCHASE predictions",
         col.regions=mypalette, sp.layout=list(nsa, londonborough))
  
  spplot(sp.fdaPDE, "prediction", key.space = "right", cex=1.5, cuts=6,
         ylim=c(155840.8,200933.9), xlim=c(503568.2,561957.5),
         main="fdaPDE PURCHASE predictions",
         col.regions=mypalette, sp.layout=list(nsa, londonborough))

 dev.off()
 save(GWR.ED, GWR.ND, output_CPP, RMSE.fdaPDE, RMSE.GWR.ED, RMSE.GWR.ND,
      file = "fdaPDE_vs_GWR_example.RData")
 