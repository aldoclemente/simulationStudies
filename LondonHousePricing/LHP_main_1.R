source("../utils.R")
source("LNH_utils.R")
data(LNNT)
data(LNHP)

#mesh = as.fdaPDE.SpatialLinesDataFrame.shp2graph(LN.nt)
#save(mesh, file="data/mesh_LN_shp2graph.RData")

spat.stat.linnet = maptools::as.linnet.SpatialLines(LN.nt)

mu_x = mean(spat.stat.linnet$vertices$x)
mu_y = mean(spat.stat.linnet$vertices$y)
coords_ = cbind(spat.stat.linnet$vertices$x-mu_x, 
                spat.stat.linnet$vertices$y-mu_y)/10^3
Windows_ = owin(xrange=c((min(spat.stat.linnet$vertices$x)-mu_x)/10^3,
                         (max(spat.stat.linnet$vertices$x)-mu_x)/10^3), 
                yrange=c((min(spat.stat.linnet$vertices$y)-mu_y)/10^3,
                         (max(spat.stat.linnet$vertices$y)-mu_y)/10^3)
)

spat.stat.linnet = linnet(vertices=as.ppp(coords_, W = Windows_), 
                          edges = cbind(spat.stat.linnet$from, spat.stat.linnet$to),
                          sparse = T)

locs = LN.prop@coords
#which.duplicated = which(duplicated(locs))
which.duplicated = which(duplicated(LN.prop@data))

locs = locs[-which.duplicated, ]
locs = cbind(locs[,1]-mu_x, locs[,2]-mu_y)/10^3

LPP = lpp(locs, spat.stat.linnet)

# Network distance matrix
ND <- pairdist.lpp(LPP) 

# fdaPDE mesh 
nodes = cbind(spat.stat.linnet$vertices$x, spat.stat.linnet$vertices$y)
edges = cbind(spat.stat.linnet$from, spat.stat.linnet$to)
mesh = create.mesh.1.5D(nodes=nodes,
                        edges=edges)
FEMbasis = create.FEM.basis(mesh)

# data
dataFrame = LN.prop[-which.duplicated,]
dataFrame@coords = cbind(LPP$data$x, LPP$data$y)
dataFrame$DATA.IDX = 1:nrow(dataFrame)
dataFrame$PURCHASE = dataFrame$PURCHASE/10^3 # k pounds

# 10-folds Cross Validation
K = 10
listData = set_Kfold_data(dataFrame, K = K)
RMSE.prediction = list( RMSE.fdaPDE = matrix(0, nrow=K, ncol=1),
                        RMSE.GWR.ND = matrix(0, nrow=K,ncol=1))
betas =list( beta.fdaPDE = matrix(0,nrow=K, ncol=3),
             beta.GWR = matrix(0,nrow=K, ncol=3))

mean.field.estimate = matrix(0,nrow=nrow(FEMbasis$mesh$nodes), ncol=1)
field.estimate  = matrix(0,nrow=nrow(FEMbasis$mesh$nodes), ncol=K)


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
  
  X = cbind( train_data$FLOORSZ, 
             train_data$PROF,   #, 
             train_data$BATH2) #, 
  lambda = seq(from=3,to=10,length.out=20) 
  output_CPP = smooth.FEM(observations = train_data$PURCHASE, 
                          locations = train_data@coords,
                          FEMbasis = FEMbasis,
                          covariates = X,
                          lambda = lambda,
                          lambda.selection.criterion = "grid",
                          lambda.selection.lossfunction = "GCV",
                          DOF.evaluation = "stochastic")
  
  plot(log10(lambda), output_CPP$optimization$GCV_vector)
  
  beta1 = output_CPP$solution$beta[1]
  beta2 = output_CPP$solution$beta[2]
  beta3 = output_CPP$solution$beta[3]
  
  prediction = beta1*test_data$FLOORSZ + 
    beta2*test_data$PROF + 
    beta3*test_data$BATH2 +
    eval.FEM(output_CPP$fit.FEM, test_data@coords)
  RMSE.prediction$RMSE.fdaPDE[i] = sqrt(mean( (prediction - test_data$PURCHASE)^2 ) )
  
  mean.field.estimate = mean.field.estimate + output_CPP$fit.FEM$coeff / K
  field.estimate[,i] = output_CPP$fit.FEM$coeff
  betas$beta.fdaPDE[i,] = c(beta1, beta2, beta3)
}
end = difftime(Sys.time(), start, units="hours")
tot.time = end

if(!dir.exists("data/")) {
  dir.create("data/")
}

filename = paste("data/LondonHousePrice-", gsub(":","_",gsub(" ","-",Sys.time())), ".RData",sep="")
save(tot.time, RMSE.prediction, betas,
     mean.field.estimate, field.estimate,
     FEMbasis, file=filename)


source("../Auxiliary/R_plot_graph.ggplot2.R")
library(viridis)
p = viridis
line.size=0.125
Field.plot <- R_plot_graph.ggplot2.2(FEM(mean.field.estimate,FEMbasis),
                                     line.size = line.size,
                                     #          color.min = min.col,
                                     #         color.max = max.col,
                                     title = bquote(hat(f)), #expression(hat(f) ~ paste("(n=",n_data[1],")",sep="")),
                                     return.ggplot.object = T,
                                     palette=p,
                                     legend.pos = "right")

rmse.plot <- boxplot_RMSE(rmse.SR.PDE = RMSE.prediction$RMSE.fdaPDE,
                          rmse.GWR =  RMSE.prediction$RMSE.GWR.ND,
                          title.size = 26, title="CV-RMSE")
pdf("img/LHP_1e3_images.pdf")
rmse.plot
Field.plot
for(i in 1:K){
  print( R_plot_graph.ggplot2.2(FEM(field.estimate[,i],FEMbasis),
                                line.size = line.size,
                                #          color.min = min.col,
                                #         color.max = max.col,
                                title = bquote(hat(f)), #expression(hat(f) ~ paste("(n=",n_data[1],")",sep="")),
                                return.ggplot.object = T,
                                palette=p,
                                legend.pos = "right")
  )
}
dev.off()


